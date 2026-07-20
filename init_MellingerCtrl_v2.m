function [ctrlSt] = init_MellingerCtrl_v2(N,Ts,vehi)
%init_mellingerCtrl_v2 creates mellinger ctrlSt ( for parametope)
% dynamically generates trajectories
%
% ctrlSt : struct (struct with components of prmtp Mellinger)
% MAIN (interfaces with main.m)
%   .auxSt : contains auxiliary variables unique to this controller that
%       may need to be modified online
%   .contrlFunc
% AUXILIARY (does not interface with main.m)

import casadi.*

X = SX.sym('X',3*8,1); % States over optimization variable [m2v(W)]
%tf = SX.sym('tf');
%Y = MX.sym('om_y',3,8);
P = SX.sym('P',3*4+3*8+1); % parameters which include x0 and reference [pD(t0);vD(t0);aD(t0);jD(t0);t0;m2v(Y)]

W = v2m(X(1:3*8),3);
%tf = X(end);
pD0 = P(1:3);
vD0 = P(4:6);
aD0 = P(7:9);
jD0 = P(10:12);
t0 = P(13);
Y = v2m(P(4*3+1+1:end),3); % Omega matrix for A

%C = [zeros(1,7);eye(7)];
%TN = dpolyN((N-1)*Ts,8,0);
%T0 = dpolyN(t0,8,0);
%Tf = polyN(tf*Ts,7);

% Equality Constraints vector
g = [W*dpolyN(t0,8,0)-pD0;
     W*dpolyN(t0,8,1)-vD0;
     W*dpolyN(t0,8,2)-aD0;
     W*dpolyN(t0,8,3)-jD0]; % Initial condition constraints
%g = [g; W*N-Y*T0]; % Final condition constraints
obj = 0; % objective function
for k=0:N
    Q1 = 0.1.*eye(3);
    if k==N
        Q1=1.*eye(3);
    end
    Q2 = 1.*eye(3);
    if k==N
        Q2=10.*eye(3);
    end
    obj = jk(W,Y,t0+(k)*Ts,N,Q1,Q2,Ts);
end
% make the decision vars one collumn vector
OPT_variables = [m2v(W)];

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

% Equality Constraints (Initial and final conditions)
args.lbg(1:length(g),1) = 0;
args.ubg(1:length(g),1) = 0;

% States' bounds (-inf <= W <= +inf; 0 <= tf <= N*Ts)
% -inf <= W <= +inf
args.lbx(1:3*8,1) = -inf;
args.ubx(1:3*8,1) = inf;
% % 0 <= tf <= N*Ts
% args.lbx(3*8+1,1) = 0;
% args.ubx(3*8+1,1) = N*Ts;

% === Initialize Auxiliary variables for mellinger controller === %
K_val = 1;
KK.Kp = 1*K_val.*eye(3); % dddp
KK.Kv = 1*K_val.*eye(3); % ddp
KK.KR = 1*K_val.*eye(3);   % dp
KK.Kom = 1*K_val.*eye(3);        % p
auxSt.KK = KK;

% Initial Conditions
auxSt.t0 = 0*Ts;
auxSt.W0 = zeros(3,8); % Init states decision variable
auxSt.dr = vehi.x0(4:6);
auxSt.ddr = zeros(3,1);
auxSt.dddr = zeros(3,1);

ctrlSt.auxSt = auxSt;

% === Specific Function Handles === %
ctrlSt.ctrlFunc = @(x0,xs,auxSt) mellinger_v2(x0,xs,auxSt,N,vehi,solver,args,Ts);

end

function [u, auxSt] = mellinger_v2(x0,xs,auxSt,pred_hor,vehi,solver,args,Ts)
% xs : set of A predictions
%nx = vehi.nx;
t0 = auxSt.t0;
%tf0 = auxSt.tf0;
W0 = auxSt.W0;
v0 = auxSt.dr;
a0 = auxSt.ddr;
j0 = auxSt.dddr;

Y = minPolyTraj(xs,t0,Ts);

args.p = [x0(1:3); v0; a0; j0; t0; m2v(Y)]; % set values of parameter vector
args.x0 = [m2v(W0)]; % init value of optimization variables
W = solverPolyTraj(solver,args);

xpred = zeros(length(x0),pred_hor);
xpred(:,1) = x0;
tnxt = t0;
for ind=1:pred_hor
    tnxt = tnxt+Ts;
    [u(:,ind),rNxt] = trajNxtK(xpred(:,ind),W,tnxt,vehi,auxSt);
    xpred(:,ind+1) = runge_kutta4(@(y) vehi.dynFunc(y,u(:,ind)), Ts, xpred(:,ind));
    if ind==1
        auxSt.dr = rNxt{2}; auxSt.ddr = rNxt{3}; auxSt.dddr = rNxt{4};
        auxSt.t0 = tnxt;
    end
end
 
%auxSt.tf0 = tf;
auxSt.W0 = W;
end

function W = solverPolyTraj(solver,args)
sol = solver('x0',args.x0,'lbx',args.lbx,'ubx',args.ubx,'lbg',args.lbg,'ubg',args.ubg,'p',args.p);

W = v2m(full(sol.x(1:3*8)),3);
%tf = full(sol.x(end));
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

function cost = jk(W,Y,t,N,Q1,Q2,Ts)
% J(W,t) = 1/N*||sum(||r_D(t)-r_A(j*Ts)||^2, j=1:N)||^2*Q1 
%   + ||r_D(t)-r_A(t)||^2*Q2
% r_D(t) = W*polyN(t,7);    r_A(t) = Y*polyN(t,7);

rD = W*dpolyN(t,8,0);
drD = W*dpolyN(t,8,1); % Velocity
ddrD = W*dpolyN(t,8,2); % Acceleration
dddrD = W*dpolyN(t,8,3); % Jerk
%drD = W*dpolyN(t,8,1);
rA = Y*dpolyN(t,8,0);
RA = [];
for j=0:N
    rA_j = Y*dpolyN(j*Ts,8,0);
    %if drD*(rA_j-rD)>=0
        RA = [RA, rA_j];
    %end
end

Q3 = 1.*eye(3); % Penalty on velocity

e1 = 1/N*sum(rD-RA,2); % Average distance between rD(k) and and rA
e2 = rD-rA;
% cost = e1'*Q1*e1 + e2'*Q2*e2 + drD'*Q3*drD;
cost = e1'*Q1*e1 + e2'*Q2*e2 + dddrD'*Q3*dddrD;
end

function [u,rNxt] = trajNxtK(x0,W,tnxt,vehi,auxSt)
nx = vehi.nx;

for k=1:4
    rNxt{k} = W*[dpolyN(tnxt,8,k-1)];
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