function [ctrlSt] = init_MellingerCtrl_v2(N,Ts)
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

X = SX.sym('X',3*8+1,1); % States over optimization variable [m2v(W);tf]
%tf = SX.sym('tf');
%Y = MX.sym('om_y',3,8);
P = SX.sym('P',3+3*8); % parameters which include x0 and reference [p_D(0);m2v(Y)]

W = v2m(X(1:3*8),3);
tf = X(end);
Y = v2m(P(3+1:end),3); % Omega matrix for A

C = [zeros(1,7);eye(7)];
TN = polyN((N-1)*Ts,7);
T0 = polyN(0*Ts,7);
Tf = polyN(tf*Ts,7);

g = []; % Equality Constraints vector
g = [g; W(:,1)-P(1:3)]; % Initial condition constraints
g = [g; W*Tf-Y*T0]; % Final condition constraints
obj = 0; % objective function
for k=1:N
    Q = 1.*eye(3);
    if k==N
        Q=15.*eye(3);
    end
    R = 1.*eye(1);
    obj = jk(W,tf,Y,(k-1)*Ts,Q,R,C);
end
% make the decision vars one collumn vector
OPT_variables = [m2v(W);tf];

nlp_prob = struct('f',obj,'x',OPT_variables, 'g',g,'p',P);
opts = struct;
opts.ipopt.max_iter = 100;
opts.ipopt.print_level = 5;
%opts.ipopt.print = 1;
opts.print_out = 0;
opts.print_in = 0;
opts.print_time = 0;
opts.ipopt.acceptable_iter = 5;
% opts.ipopt.acceptable_tol = 1e-8;
% opts.ipopt.acceptable_obj_change_tol = 1e-6;
opts.ipopt.acceptable_tol = 1e-1;
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
% 0 <= tf <= N*Ts
args.lbx(3*8+1,1) = 0;
args.ubx(3*8+1,1) = N*Ts;

% Initial Conditions
auxSt.t0 = 0*Ts;
auxSt.tf0 = N*Ts;
auxSt.W0 = ones(3,8); % Init states decision variable

ctrlSt.auxSt = auxSt;

% === Specific Function Handles === %
ctrlSt.ctrlFunc = @(x0,xs,auxSt) NMPC_prmtp_v2(x0,xs,auxSt,solver,args,Ts);

end

function [u_out, auxSt] = NMPC_prmtp_v2(x0,xs,auxSt,solver,args,Ts)
% xs : set of A predictions
% generate A trajectory from 8 samples of xs
k_step = floor(size(xs,2)/8);
t = (0:k_step:size(xs,2))';
if length(t)>8
    t = t(1:end-1);
end
Y = trajectory3Dv3(Ts.*t,Ts,xs(1:3,(t+1)));
Y = Y{1};

t0 = auxSt.t0;
tf0 = auxSt.tf0;
W0 = auxSt.W0;
args.p = [x0(1:3); m2v(Y)]; % set values of parameter vector
args.x0 = [m2v(W0); tf0]; % init value of optimization variables
sol = solver('x0',args.x0,'lbx',args.lbx,'ubx',args.ubx,'lbg',args.lbg,'ubg',args.ubg,'p',args.p);

W = v2m(full(sol.x(1:3*8)),3);
tf = full(sol.x(end));


auxSt.t0 = t0+Ts;
auxSt.tf0 = tf;
auxSt.W0 = W;
end

function res = polyN(t,N) 
    res = (1./factorial(0:N).*t.^(0:N))';
end

function cost = jk(X,tf,Y,t,Q,R,C)
dT = polyN(t,6);

ev = crossprodmat(Y*C*dT)*X*C*dT; % Derivative error

cost = ev'*Q*ev + tf'*R*tf;
end

function varargout = trajNxtK(x0_cell, xs_cell, t0, tf, Ts)
    % calculate trajectory and return next point

    r0 = x0_cell{1}; rd = xs_cell{1};
    dr0 = x0_cell{2}; drd = xs_cell{2};
    ddr0 = x0_cell{3}; ddrd = xs_cell{3};
    dddr0 = x0_cell{4}; dddrd = xs_cell{4};

    % tf is some function of the distance between points
    %tf = floor(norm(rd-r0)/2/Ts)*Ts; % Assume average speed
    t1 = t0+Ts;
    tf = max(tf,t1);
    w_aux = trajectory3Dv2([t0;tf], Ts, [r0,rd], [dr0,drd], [ddr0,ddrd], [dddr0,dddrd]);
    w_aux = w_aux{1};

    traj_order = 4;
    prod_vec = ones(1,2*traj_order);
    for k=1:nargout
        pow_vec = 2*traj_order-k:-1:1;
        rx = polyval(prod_vec.*w_aux(1,end:-1:k),t1);
        ry = polyval(prod_vec.*w_aux(2,end:-1:k),t1);
        rz = polyval(prod_vec.*w_aux(3,end:-1:k),t1);
        
        prod_vec = prod_vec(1:length(pow_vec)).*pow_vec;
        r = [rx; ry; rz]; % Coordinate of trajectory position
        varargout{k} = r;
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