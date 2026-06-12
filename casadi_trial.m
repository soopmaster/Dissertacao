close all;
clearvars;

import casadi.*

simT = 5; % Sim time seconds
Ts = 0.01;
knum = simT/Ts;

N = 5; % MPC Time horizon
nx = 18;
nx_full = nx+nx^2;
ny = nx;
ny_full = nx_full;
nu = 4;

vehi = init_QuadDroneMellinger(Ts);

% Using github workshop
% (https://github.com/MMehrez/MPC-and-MHE-implementation-in-MATLAB-using-Casadi/blob/master/workshop_github/MPC_MHE_slides.pdf)

% Multiple shooting OCP to transform to NLP

x = SX.sym('x',nx_full);
u = SX.sym('u',nu);
rk4 = runge_kutta4(@(xfull) fullDynFunc(@(x) vehi.dynFunc(x,u), @(alph,x) prmtpU(alph,-vehi.jacFunc(x,u),nx),xfull,nx),Ts,x); % rk4 system r.h.s

f = Function('f',{x,u},{rk4},{'xk','uk'},{'xk_1'});
U = SX.sym('U',nu,N); % decision variables (control)
P = SX.sym('P',nx_full+ny_full); % parameters which include x0 and reference [x0;ref]
X = SX.sym('X',nx_full,N+1); % States over optimization problem

obj = 0; % objective function
g = []; % Constraints vector
Qk = 10.*eye(3);
Qmk = 10.*eye(3);
Rk = 5.*eye(4);

xk = X(:,1); % initial state
g = [g; xk-P(1:nx_full)]; % Initial condition constraints
ref = P(nx_full+1:end); % reference
for k=1:N
    xk = X(:,k); uk = U(:,k);
    obj = obj + jk(xk,uk,ref,Qk,Qmk,Rk,nx,ny); % Calc objective
    xk_1 = X(:,k+1); % xk+1
    fxk = f(xk,uk); % f(xk,uk)
    g = [g; xk_1-fxk]; % [x0 = xinit; ...; xk+1 = f(xk,uk)]
end

% make the decision vars one collumn vector
OPT_variables = [m2v(X); m2v(U)]; % = [x0; x1; ... ; xN+1; u0; u1; ... ; uN]

nlp_prob = struct('f',obj,'x',OPT_variables, 'g',g,'p',P);
opts = struct;
opts.ipopt.max_iter = 100;
opts.ipopt.print_level = 0;
opts.print_time = 0;
opts.ipopt.acceptable_tol = 1e-8;
opts.ipopt.acceptable_obj_change_tol = 1e-6;
solver = nlpsol('solver', 'ipopt', nlp_prob, opts);

args = struct;

% Equality Constraints
args.lbg(1:nx_full*(N+1),1) = 0;
args.ubg(1:nx_full*(N+1),1) = 0;

% States' bounds
bx = inf.*ones(nx_full*(N+1),1);
args.lbx(1:length(bx),1) = -bx;
args.ubx(1:length(bx),1) = bx;

% Control bounds
lbu = repmat(vehi.vbds_init(end-nu+1:end,1),N,1);
ubu = repmat(vehi.vbds_init(end-nu+1:end,2),N,1);
args.lbx(length(bx)+1:length(bx)+length(lbu),1) = lbu;
args.ubx(length(bx)+1:length(bx)+length(lbu),1) = ubu;

% SIMULATION LOOP
t0 = 0;
x0 = [vehi.x0; m2v(vehi.alph0)]; % initial condition
ps = [5;5;5];
xs = [[ps;vehi.x0(length(ps)+1:end)]; m2v(vehi.alph0)]; % reference

xx(:,1) = x0; % history of states
t(1) = t0;
u0 = zeros(nu,N);
X0 = repmat(x0,1,N+1); % Init states decision variable

% Start MPC
k = 1; xx1 = []; u_cl = [];
while (norm(proj3(x0)-proj3(xs),2)>1e-2 && k<knum)
    args.p = [x0;xs]; % set values of parameter vector
    args.x0 = [m2v(X0); m2v(u0)]; % init value of optimization variables
    sol = solver('x0',args.x0,'lbx',args.lbx,'ubx',args.ubx,'lbg',args.lbg,'ubg',args.ubg,'p',args.p);
    u = v2m(full(sol.x(nx_full*(N+1)+1:end)),nu); % control input from solution
    X0 = v2m(full(sol.x(1:nx_full*(N+1))),nx_full); % history of solution TRAJECTORY accross horizon
    xx1(1:nx_full,:,k) = X0;
    u_cl = [u_cl; u(1,:)];
    t(k+1) = t0;
    [t0,x0,u0] = mpc_shift(Ts,t0,x0,u,f);
    xx(:,k+1) = x0;
    X0 = [X0(:,2:end), X0(:,end)];
    k
    k = k+1;
end


return

% Using docs
% Formulate DAE
x = SX.sym('x',nx_full);
u = SX.sym('u',nu);
y = SX.sym('y',nx_full);
fsym = runge_kutta4(@(xfull) fullDynFunc(@(x) vehi.dynFunc(x,u), @(alph,x) prmtpU(alph,-vehi.jacFunc(x,u),nx),xfull,nx),Ts,x);
%f = fullDynFunc(@(x) vehi.dynFunc(x,u), @(alph,x) prmtpU(alph,-vehi.jacFunc(x,u),nx),x,nx);
%h = 

% Create Solver instance
Qk = eye(3);
Qmk = eye(3);
Rk = eye(4);
%jksym = jk(x,u,y,Qk,Qmk,Rk,nx,nx);
F = Function('F',{x,u},{fsym},{'x','u'},{'x_next'});
Jk = @(Qk,Qmk,Rk) Function('Jk',{x,u,y},{jk(x,u,y,Qk,Qmk,Rk,nx,nx)},{'x','u','y'}, {'Jk'});

% Formulate as
% min J(w)
% s.t. G(w) = 0, lbw <= w <= ubw

% Empty NLP
w = {}; lbw = []; ubw = [];
G = {}; J = 0;

% Initial Conditions
Xk = MX.sym('X0',nx_full);
w{end+1} = Xk;
lbx = -inf.*ones(nx_full,1);
ubx = inf.*ones(nx_full,1);
lbu = vehi.vbds_init(end-nu+1:end,1);
ubu = vehi.vbds_init(end-nu+1:end,2);
lbw = [lbw; lbx];
ubw = [ubw; ubx];

for k=1:N
    % Local control
    Uname = ['U' num2str(k-1)];
    Uk = MX.sym(Uname, nu);
    w{end+1} = Uk;
    lbw = [lbw; lbu];
    ubw = [ubw; ubu];

    Yname = ['Y' num2str(k-1)];
    Yk = MX.sym(Yname, nx_full);
    w{end+1} = Yk;
    lbw = [lbw; lbx];
    ubw = [ubw; ubx];

    % New Local State
    Xname = ['X' num2str(k)];
    Xk = MX.sym(Xname, nx_full);
    w{end+1} = Xk;
    lbw = [lbw; lbx];
    ubw = [ubw; ubx];

    % Call Integrator
    Fk = F('x',Xk,'u',Uk);
    Jk_aux = Jk(Qk,Qmk,Rk);
    J = J + Jk_aux(Fk.x_next,Uk,Yk);

    % Continuity Constraint
    G{end+1} = Fk.x_next - Xk;
end

% Create NLP Solver
nlp.f = J;
nlp.g = vertcat(G{:});
nlp.x = vertcat(w{:});
S = nlpsol('S','blocksqp',nlp);

% Solve NLP
% x0_val [x0; xalph0; u0; y0; yalph0]
x0_val = [vehi.x0; m2v(vehi.alph0); vehi.u0; vehi.x0; m2v(vehi.alph0)];
r = S('lbx',lbw,'ubx',ubw,'x0',0,'lbg',0,'ubg',0);




% Using tutorial
x = MX.sym('x',nx); % States
u = MX.sym('u',nu); % Free Parameter

xfull = MX.sym('x',nx_full); % States
ufull = MX.sym('u',nu); % Free Parameter

%ode = vehi.dynFunc(x,u);
fsym = runge_kutta4(@(xrk4) vehi.dynFunc(xrk4,u),Ts,x);
fsym2 = runge_kutta4(@(xfull) fullDynFunc(@(x) vehi.dynFunc(x,ufull), @(alph,x) prmtpU(alph,-vehi.jacFunc(x,ufull),nx),xfull,nx),Ts,xfull);
fk = Function('fk',{x,u},{fsym},{'x','u'},{'x_next'});
fk2 = Function('fk_full',{xfull,ufull},{fsym2},{'x','u'},{'x_next'});

sim = fk.mapaccum(knum);

opti = casadi.Opti();

x = opti.variable(nx,knum+1);
u = opti.variable(nu,knum);
x0 = opti.parameter(nx,1);

xfull = opti.variable(nx_full,N+1);
yfull = opti.variable(nx_full,N+1);
ufull = opti.variable(nu,N);
xfull0 = opti.parameter(nx_full,1);

% Make into cells
Qk = eye(3);
Qmk = eye(3);
Rk = eye(4);
%Jk(xfull,ufull,yfull,Qk,Qmk,Rk,nx,nx)
opti.minimize(Jk(xfull,ufull,yfull,Qk,Qmk,Rk,nx,nx));

for k=1:N
    opti.subject_to(xfull(:,k+1)==fk2(xfull(:,k),ufull(:,k)));
end

% Choose a concrete solver
opti.solver('sqpmethod',struct('qpsol','qrqp'));

% And choose a concrete value for p
xfull0_val = [vehi.x0; m2v(vehi.alph0)];
opti.set_value(xfull0,xfull0_val);
sol = opti.solve();

function [t0,x0,u0] = mpc_shift(Ts,t0,x0,u,f)
% Evolve system one timestep forward
xk = x0;
uk = u(:,1);
xk_1 = f(xk,uk);
x0 = full(xk_1);

t0 = t0+Ts;
u0 = [u(:,2:size(u,2)), u(:,size(u,2))];
end

function cost = jk(X,u,Y,Q,Qm,R,nx,ny)
    xind_step = nx^2+nx;
    yind_step = ny^2+ny;
    %N = size(u,2);
    N = 1;

    for k=1:N

        xnx = X(1:nx,k);
        x = proj3(xnx);
    
        yny = Y(1:ny,k);
        y = proj3(yny);
    
        Dnx = v2m(X(nx+1:xind_step,k),nx);
        D = proj3(Dnx);
    
        Any = v2m(Y(ny+1:yind_step,k),ny);
        A = proj3(Any);
    
        if k==1
            cost = -m2v(A)'*m2v((Qm*D)') + (y-x)'*Q*(y-x) + u(:,k)'*R*u(:,k);
        else
            cost = cost + -m2v(A)'*m2v((Qm*D)') + (y-x)'*Q*(y-x) + u(:,k)'*R*u(:,k);
        end
    end
    % xind_step = nx^2+nx;
    % yind_step = ny^2+ny;
    % Cx = [eye(3), zeros(3,nx-3)];
    % Cy = [eye(3), zeros(3,ny-3)];
    % 
    % Dbar = [];
    % xbar = [];
    % Abar = [];
    % ybar = [];
    % for k=0:N
    %     Dnxk = reshape(X(k*xind_step+(nx+1:xind_step)),nx,[]);
    %     Dk = proj3(Dnxk);
    %     Dbar = [Dbar; Dk];
    % 
    %     Anyk = reshape(Y(k*yind_step+(ny+1:yind_step)),ny,[]);
    %     Ak = proj3(Anyk);
    %     Abar = [Abar, Ak];
    % 
    %     xnxk = X(k*xind_step+(1:nx));
    %     xk = Cx*xnxk;
    %     xbar = [xbar; xk];
    % 
    %     ynyk = Y(k*yind_step+(1:ny));
    %     yk = Cy*ynyk;
    %     ybar = [ybar; yk];
    % end
    % 
    % 
    % cost = -m2v(Abar)'*m2v((Qbarm*Dbar)') + (ybar-xbar)'*Qbar*(ybar-xbar) + U'*Rbar*U;
end