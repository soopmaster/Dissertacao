function [ctrlSt] = init_NMPCctrl(ny,N,Ts,vehi)
%INIT_NMPC outputs ctrlStruct of NMPC (for parametope)
%   
% ctrlSt : struct (struct with components of prmtp NMPC)
% MAIN (interfaces with main.m)
%   .auxSt : contains auxiliary variables unique to this controller that
%       may need to be modified online
%   .contrlFunc
% AUXILIARY (does not interface with main.m)

% Using github workshop
% (https://github.com/MMehrez/MPC-and-MHE-implementation-in-MATLAB-using-Casadi/blob/master/workshop_github/MPC_MHE_slides.pdf)

% Multiple shooting OCP to transform to NLP

% Notes: I should be tracking error
% Maybe have several control modes; first to approach target, then to align
% Note: This is failing! It may be because problem is degenerate, i.e.
% information is lost

import casadi.*

nx = vehi.nx;
nx_full = nx+nx^2;
ny_full = ny+ny^2;
nu = vehi.nu;

x = SX.sym('x',nx_full);
u = SX.sym('u',nu);
rk4 = runge_kutta4(@(xfull) vehi.dynFunc(xfull,u),Ts,x); % rk4 system r.h.s

f = Function('f',{x,u},{rk4},{'xk','uk'},{'xk_1'});
U = SX.sym('U',nu,N); % decision variables (control)
P = SX.sym('P',nx_full+ny_full); % parameters which include x0 and reference [x0;ref]
X = SX.sym('X',nx_full,N+1); % States over optimization problem

obj = 0; % objective function
g = []; % Constraints vector

xk = X(:,1); % initial state
g = [g; xk-P(1:nx_full)]; % Initial condition constraints
ref = P(nx_full+1:end); % reference
for k=1:N
    %Qk = 10*exp((k-1)/2).*eye(3);
    Qk = 0.1.*eye(3);
    if k==N
        Qk=15.*eye(3);
    end
    %Qmk = 1/exp(k).*eye(3);
    %Qmk = 0.1*exp((k-1)/3).*eye(1);
    Qmk = 0.00000001;
    %Qmk = 1;
    Rk = 1.*eye(4);
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

% Initial Conditions
auxSt.u0 = repmat(vehi.u0,1,N);
auxSt.X0 = repmat(vehi.x0,1,N+1); % Init states decision variable

ctrlSt.auxSt = auxSt;

% === Specific Function Handles === %
ctrlSt.ctrlFunc = @(x0,xs,auxSt) NMPC_prmtp(x0,xs,auxSt,solver,args,N,nx_full,nu);

end

function [u_out, auxSt] = NMPC_prmtp(x0,xs,auxSt,solver,args,N,nx_full,nu)
u0 = auxSt.u0;
X0 = auxSt.X0;
args.p = [x0;xs]; % set values of parameter vector
args.x0 = [m2v(X0); m2v(u0)]; % init value of optimization variables
sol = solver('x0',args.x0,'lbx',args.lbx,'ubx',args.ubx,'lbg',args.lbg,'ubg',args.ubg,'p',args.p);
u = v2m(full(sol.x(nx_full*(N+1)+1:end)),nu); % control input from solution
u_out = u(:,1);
auxSt.X0 = v2m(full(sol.x(1:nx_full*(N+1))),nx_full); % history of solution TRAJECTORY accross horizon
auxSt.u0 = [u(:,2:end), u(:,end)];
end

function cost = jk(X,u,Y,Q,Qm,R,nx,ny)
xind_step = nx^2+nx;
yind_step = ny^2+ny;

xnx = X(1:nx);
x = proj3(xnx);

yny = Y(1:ny);
y = proj3(yny);

% Maximize inner product of the principal axis of ellipsoids (eigen_vecs)
Dnx = v2m(X(nx+1:xind_step),nx);
%Dnx = inv(Dnx);
D = proj3(Dnx'*Dnx);

Any = v2m(Y(ny+1:yind_step),ny);
%Any = inv(Any);
A = proj3(Any'*Any);

%eR = (m2v(A)-m2v(D)); % Ellipsoid Matrix error
eR = -m2v(A)'*m2v(D);

%cost = -m2v(A)'*m2v((Qm*D)') + (y-x)'*Q*(y-x) + u(:,k)'*R*u(:,k);
cost = Qm*eR + (y-x)'*Q*(y-x) + u'*R*u;
%cost =  + (y-x)'*Q*(y-x) + u'*R*u;
end