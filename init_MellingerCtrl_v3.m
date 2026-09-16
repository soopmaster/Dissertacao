function [ctrlSt] = init_MellingerCtrl_v3(N,Ts,vehi)
%init_mellingerCtrl_v2 creates mellinger ctrlSt ( for parametope)
% dynamically generates trajectories
%
% ctrlSt : struct (struct with components of prmtp Mellinger)
% MAIN (interfaces with main.m)
%   .auxSt : contains auxiliary variables unique to this controller that
%       may need to be modified online
%   .contrlFunc
% AUXILIARY (does not interface with main.m)

% This version has a varying arrival time

import casadi.*

X = SX.sym('X',3*5+1,1); % States over optimization variable [m2v(W);tf]
%tf = SX.sym('tf');
%Y = MX.sym('om_y',3,8);
P = SX.sym('P',3*4+3+1,1); % parameters which include x0 and reference [pD(t0);vD(t0);aD(t0);jD(t0);t0;pA(tf)]

W = v2m(X(1:3*5),3);
tf = X(end);
td = X(end); % delta t
p0 = P(1:3);
v0 = P(4:6);
a0 = P(7:9);
j0 = P(10:12);
t0 = P(13);
pf = P(4*3+1+1:end); % Final point
amax = 10; % max acceleration
td_min = (1/amax)*(-norm(v0)+sqrt(norm(v0)^2-2*amax*(norm(p0-pf)))); % Minimum td possible for maximum acceleration

%C = [zeros(1,7);eye(7)];
%TN = dpolyN((N-1)*Ts,8,0);
%T0 = dpolyN(t0,8,0);
%Tf = polyN(tf*Ts,7);

% Equality Constraints vector
g = [W*dpolyN(0,5,0)-p0;
     W*dpolyN(td,5,0)-pf;
     W*dpolyN(0,5,1)-v0;
     W*dpolyN(0,5,2)-a0;
     W*dpolyN(0,5,3)-j0; % Initial condition constraints
    td-td_min]; % Inequality Constraint (dt > td_min)

% % Equality Constraints vector
% g = [W*dpolyN(t0,5,0)-p0;
%      W*dpolyN(tf,5,0)-pf;
%      W*dpolyN(t0,5,1)-v0;
%      W*dpolyN(t0,5,2)-a0;
%      W*dpolyN(t0,5,3)-j0; % Initial condition constraints
%     tf-t0]; % Inequality Constraint
%g = [g; W*N-Y*T0]; % Final condition constraints
%g = tf-t0; % Inequality Constraint
Q = eye(3);
R = 0.00001;
obj = jk(p0,pf,v0,a0,j0,t0,W,td,Q,R);
% make the decision vars one collumn vector
OPT_variables = [m2v(W);td];

nlp_prob = struct('f',obj,'x',OPT_variables, 'g',g,'p',P);
opts = struct;
opts.ipopt.max_iter = 100;
opts.ipopt.print_level = 5;
%opts.ipopt.print = 1;
opts.print_out = 0;
opts.print_in = 0;
opts.print_time = 0;
%opts.ipopt.acceptable_iter = 5;
% opts.ipopt.acceptable_tol = 1e-8;
% opts.ipopt.acceptable_obj_change_tol = 1e-6;
%opts.ipopt.acceptable_tol = 1e-1;
%opts.ipopt.acceptable_obj_change_tol = 1e-1;
solver = nlpsol('solver', 'ipopt', nlp_prob, opts);

args = struct;

% Equality Constraints
args.lbg(1:length(g)-1,1) = 0;
args.ubg(1:length(g)-1,1) = 0;
% InEquality Constraints (td-td_min(amax)>0)
args.lbg(length(g),1) = 0;
args.ubg(length(g),1) = Inf;

% States' bounds (-inf <= W <= +inf; 0 <= tf <= N*Ts)
% -inf <= W <= +inf
% args.lbx(1:3*5,1) = -inf;
% args.ubx(1:3*5,1) = inf;
% % % 0 <= tf <= N*Ts
% args.lbx(3*5+1,1) = 0;
% args.ubx(3*5+1,1) = N*Ts;

% === Initialize Auxiliary variables for mellinger controller === %
K_val = 1;
KK.Kp = 1*K_val.*eye(3); % dddp
KK.Kv = 1*K_val.*eye(3); % ddp
KK.KR = 1*K_val.*eye(3);   % dp
KK.Kom = 1*K_val.*eye(3);        % p
auxSt.KK = KK;

% Initial Conditions
auxSt.t0 = 0*Ts;
auxSt.tf = 1*Ts;
auxSt.W0 = zeros(3,5); % Init states decision variable
auxSt.dr = vehi.x0(4:6);
auxSt.ddr = zeros(3,1);
auxSt.dddr = zeros(3,1);

ctrlSt.auxSt = auxSt;

% === Specific Function Handles === %
ctrlSt.ctrlFunc = @(x0,xs,auxSt) mellinger_v3(x0,xs,auxSt,N,vehi,solver,args,Ts);

end

function [u, auxSt] = mellinger_v3(x0,xs,auxSt,pred_hor,vehi,solver,args,Ts)
% xs : set of A predictions
%nx = vehi.nx;

pred_hor = size(xs,2)-1;

% Path Planning - Generate 5th order trajectory between p0 and pf
t0 = auxSt.t0;
%tf0 = auxSt.tf0;
W0 = auxSt.W0;
v0 = auxSt.dr;
a0 = auxSt.ddr;
j0 = auxSt.dddr;
pf = xs(1:3,end);
p0 = x0(1:3);

args.p = [p0; v0; a0; j0; t0; pf]; % set values of parameter vector
args.x0 = [m2v(W0);100]; % init value of optimization variables
[Waux,td] = solverPolyTraj(solver,args);
tf = t0 + td;

% Sample 3 points from trajectory
t1 = t0 + (tf-t0)/4; t2 = t0 + 2*(tf-t0)/4; t3 = t0 + 3*(tf-t0)/4;
p1 = Waux*dpolyN(t1,5,0);
p2 = Waux*dpolyN(t2,5,0);
p3 = Waux*dpolyN(t3,5,0);

% Solve QP for 7th order trajectory
% P = [p0; p1; p2; p3; pf; v0; a0; j0];
% T = [tdim(t0,1,8,0); tdim(tf,1,8,0); tdim(t1,1,8,0); tdim(t2,1,8,0); tdim(t3,1,8,0); tdim(t0,1,8,1); tdim(t0,1,8,2); tdim(t0,1,8,3)];
%W = (T*P)/(T'*T);
vnan = [NaN;NaN;NaN];
W = trajectory3Dv4([t0,t1,t2,t3,tf],[p0,p1,p2,p3,pf],[v0,vnan,vnan,vnan,vnan],[a0,vnan,vnan,vnan,vnan],[j0,vnan,vnan,vnan,vnan]);

% Compute pnxt, vnxt, anxt, jnxt
% Compute Mellinger input
xpred = zeros(length(x0),pred_hor);
xpred(:,1) = x0;
tnxt = t0;
for ind=1:pred_hor
    tnxt = tnxt+Ts;
    [u(:,ind),rNxt] = trajNxtK(xpred(:,ind),W,tnxt,vehi,auxSt);
    xpred(:,ind+1) = runge_kutta4(@(y) vehi.partDynFunc(y,u(:,ind)), Ts, xpred(:,ind));
    if ind==1
        auxSt.dr = rNxt{2}; auxSt.ddr = rNxt{3}; auxSt.dddr = rNxt{4};
        auxSt.t0 = tnxt;
    end
end
 
auxSt.tf = tf;
auxSt.W0 = Waux;
end

function tvec = tdim(t,dim,polyOrder,dorder)
% Get t dpolyN vec of dim
tvec = repmat(dpolyN(t,polyOrder,dorder),dim,1);
end

function [W,td] = solverPolyTraj(solver,args)
% sol = solver('x0',args.x0,'lbx',args.lbx,'ubx',args.ubx,'lbg',args.lbg,'ubg',args.ubg,'p',args.p);
sol = solver('x0',args.x0,'lbg',args.lbg,'ubg',args.ubg,'p',args.p);

W = v2m(full(sol.x(1:3*5)),3);
td = full(sol.x(end));
end

function Y = minPolyTraj(xs,t0,Ts)
% generate A trajectory from 8 samples of xs
k_step = floor(size(xs,2)/8);
t = (0:k_step:size(xs,2))';
if length(t)>8
    t = t(1:end-1);
end
Y = trajectory3Dv3(t0+Ts.*t,Ts,xs(1:3,(t+1)));
Y = Y{1};
end

function cost = jk(p0,pf,v0,a0,j0,t0,W,td,Q,R)
% J(W,tf) = ||P-W*T||^2*Q + tf^2*R
% P = [p0; pf; v0; a0; j0]
% T = [dpolyN(t0,5,0); dpolyN(tf,5,0); dpolyN(t0,5,1); dpolyN(t0,5,2); dpolyN(t0,5,3)] 

% P = [p0; pf; v0; a0; j0];
P = [p0];
% T = [tdim(t0,1,5,0); tdim(tf,1,5,0); tdim(t0,1,5,1); tdim(t0,1,5,2); tdim(t0,1,5,3)];
T = [tdim(0,1,5,0)];
%Wk = kron(eye(5),W);

cost = (P-W*T)'*Q*(P-W*T) + td^2*R;
end

function [u,rNxt] = trajNxtK(x0,W,tnxt,vehi,auxSt)
nx = vehi.nx;

for k=1:4
    rNxt{k} = W*[dpolyN(tnxt,8,k-1)];
    % rNxt{k} = W*tdim(tnxt,3,8,k-1)<
end

yaw_ref= 0;
dyaw_ref= 0;
u = mellinger_controller(auxSt.KK,rNxt{1},rNxt{2},rNxt{3},rNxt{4},yaw_ref,dyaw_ref,x0(1:nx),vehi.I,vehi.M);

% Clamp u to bounds
u(1) = min(max(u(1),vehi.Fmin),vehi.Fmax);
tauBnorm = norm(u(2:3));
if tauBnorm == 0
    u(2:4) = [0;0;0];
elseif tauBnorm > vehi.tauBmax
    scale = vehi.tauBmax/tauBnorm;
    u(2:4) = scale.*u(2:4);
end
end

function [u, auxSt] = quadController(x0, xs, auxSt, tf, vehi)
    % x0 : current state [x0; alph0]
    % xs : desired state [ps; vs; ~]
    % auxSt : contains acceleration and jerk

    Ts = vehi.Ts;
    nx = vehi.nx;

    r = xs(1:3);
    dr = xs(4:6);
    % In the future, let these be unknowns
    ddr = zeros(3,1);
    dddr = zeros(3,1);
    yaw_ref = 0;
    dyaw_ref = 0;

    
    x0_cell = {x0(1:3), x0(4:6), auxSt.ddr, auxSt.dddr};
    xs_cell = {r,dr,ddr,dddr};
    t0 = auxSt.t0;
    [rNxt,drNxt,ddrNxt,dddrNxt] = trajNxtK(x0_cell, xs_cell, t0, tf, Ts);
    auxSt.ddr = ddrNxt; auxSt.dddr = dddrNxt; auxSt.t0 = t0+Ts;

    u = mellinger_controller(auxSt.KK,rNxt,drNxt,ddrNxt,dddrNxt,yaw_ref,dyaw_ref,x0(1:nx),vehi.I,vehi.M);

    % Clamp u to bounds
    u(1) = min(max(u(1),vehi.Fmin),vehi.Fmax);
    tauBnorm = norm(u(2:3));
    if tauBnorm == 0
        u(2:4) = [0;0;0];
    elseif tauBnorm > vehi.tauBmax
        scale = vehi.tauBmax/tauBnorm;
        u(2:4) = scale.*u(2:4);
    end

end