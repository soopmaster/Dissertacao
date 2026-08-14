function [ctrlSt] = init_NMPCctrl_v3(ny,N,Ts,vehi)
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

% Using documentation:
%   (https://systemanalysisdpt-cmc-msu.github.io/ellipsoids/doc/chap_ellcalc.html#affine-transformation)

% Generates control inputs based on ellipsoid error
% This version only propagates main system!! Not full ellipsoid

% Multiple shooting OCP to transform to NLP

import casadi.*

nx = vehi.nx;
nx_full = nx+nx^2;
ny_full = ny+ny^2;
nu = vehi.nu;

x = SX.sym('x',nx);
u = SX.sym('u',nu);
rk4 = runge_kutta4(@(xfull) vehi.partDynFunc(xfull,u),Ts,x); % rk4 system r.h.s

f = Function('f',{x,u},{rk4},{'xk','uk'},{'xk_1'});
U = SX.sym('U',nu,N); % decision variables (control)
% P = SX.sym('P',nx+3*N+ny_full*N); % parameters which include x0 and reference [x0;alphD0;alphD1...;ref1;ref2;...;refN]
P = SX.sym('P',nx+3*N+(ny+3^2)*N); % parameters which include x0 and reference [x0;alphD0;alphD1...;ref1;ref2;...;refN]
X = SX.sym('X',nx,N+1); % States over optimization problem
D_KL = SX.sym('DKL',N);

g = []; % Constraints vector

xk = X(:,1); % initial state
g = [g; xk-P(1:nx)]; % Initial condition constraints
% Y = v2m(P(nx+3*N+1:end),ny_full); % reference
Y = v2m(P(nx+3*N+1:end),(ny+3^2)); % reference
scale_vec = v2m(P(nx+1:nx+3*N),3); % Parametope alph

Obj = Jfunc(X,Y,U,scale_vec,D_KL,ny); % Get Cost Function
obj = Obj.obj_d + Obj.obj_u + Obj.obj_DKL;

for k=1:N
    xk = X(:,k);
    uk = U(:,k);

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
%opts.ipopt.acceptable_iter = 5;
% opts.ipopt.acceptable_tol = 1e-8;
% opts.ipopt.acceptable_obj_change_tol = 1e-6;
%opts.ipopt.acceptable_tol = 1e-1;
%opts.ipopt.acceptable_obj_change_tol = 1e-1;
solver = nlpsol('solver', 'ipopt', nlp_prob, opts);

args = struct;

% Equality Constraints
args.lbg(1:nx*(N+1),1) = 0;
args.ubg(1:nx*(N+1),1) = 0;

% States' bounds
bx = inf.*ones(nx*(N+1),1);
args.lbx(1:length(bx),1) = -bx;
args.ubx(1:length(bx),1) = bx;
args.lbx(3:nx:end) = 0; % Height must be above 0

% Control bounds
lbu = repmat(vehi.vbds_init(end-nu+1:end,1),N,1);
ubu = repmat(vehi.vbds_init(end-nu+1:end,2),N,1);
args.lbx(length(bx)+1:length(bx)+length(lbu),1) = lbu;
args.ubx(length(bx)+1:length(bx)+length(lbu),1) = ubu;

% Initial Conditions
auxSt.error_track = zeros(6,1);
auxSt.u0 = repmat(vehi.u0,1,N);
auxSt.X0 = repmat(vehi.x0,1,N+1); % Init states decision variable
auxSt.scale_factor = ones(3,N); % Scaling factors for axis of reduced parametopes

ctrlSt.auxSt = auxSt;

% === Specific Function Handles === %
ctrlSt.ctrlFunc = @(x0,xs,auxSt) NMPC_prmtp(x0,xs,auxSt,solver,args,N,nx,nu,ny);

end

function [u, auxSt] = NMPC_prmtp(x0,xs,auxSt,solver,args,N,nx,nu,ny)
u0 = auxSt.u0;
X0 = auxSt.X0;
scale = auxSt.scale_factor;
paramWA = xs(1:(ny+3^2),:);

Cy = eye(ny);
Cy = Cy(:,4:end);
By = zeros(size(Cy,2),1);

% Condition input data! and 
% scale entries should be positive. 
% all xs alphs should be full rank.
tol = 1.5e-15;
Tikh_coef = 1e-5;
for k=1:size(xs,2)
    scale(:,k) = 0.5.*[0.1;0.1;1];
    if any(scale(:,k)<=0)
        scale(:,k) = scale(:,k) + (-min(scale(:,k)) + Tikh_coef).*ones(size(scale(:,k)));
    end

    alphA = v2m(xs(ny+1:end,k),ny);

    QA = alphA'*alphA;
    [~,WA] = ellipIntersectPoly(zeros(ny,1),QA,Cy,By);

    % [evecA,evalA] = eig(QA);
    [evecA,evalA] = eig(WA);
    evalA1 = diag(evalA);
    if 1
        % any(evalA1<=0)
        %evalAnew = evalA1 + (-min(evalA1) + Tikh_coef).*ones(size(evalA1));
        evalAnew = 0.5.*[0.1;0.1;1];
        % QA = evecA*diag(evalAnew)*evecA';
        WA = evecA*diag(evalAnew)*evecA';
        %sc = (-min(evalA) + Tikh_coef);
        %QA = QA+sc.*eye(ny);
    end
    
    % %closeToSing = (rank(alphA) < size(alphA,1));
    % alphA_aux = alphA;
    % closeToSing = rcond(alphA_aux'*alphA_aux) < tol;
    % Tikh_coef = 1e-15;
    % while closeToSing
    %     % Use Ridge Regression (Tikhonov Regression) if QA is close to singular
    %     Tikh_coef = Tikh_coef*1e1;
    %     alphA_aux = alphA + Tikh_coef.*eye(ny);
    %     closeToSing = rcond(alphA_aux'*alphA_aux) < tol;
    % end
    % paramQA(ny+1:end,k) = m2v(QA);
    paramWA(ny+1:end,k) = m2v(WA);
end

args.p = [x0(1:nx); m2v(scale); m2v(paramWA)]; % set values of parameter vector
args.x0 = [m2v(X0(1:nx,:)); m2v(u0)]; % init value of optimization variables
sol = solver('x0',args.x0,'lbx',args.lbx,'ubx',args.ubx,'lbg',args.lbg,'ubg',args.ubg,'p',args.p);
u = v2m(full(sol.x(nx*(N+1)+1:end)),nu); % control input from solution
u_out = u(:,1);

% Record Error/Cost!
D_KL = [];
E = Jfunc(X0,paramWA,u0,scale,D_KL,ny);
auxSt.error_track = [E.obj_d; E.obj_u; E.obj_DKL; E.unsObj_d; E.unsObj_u; E.unsObj_DKL];

auxSt.X0 = v2m(full(sol.x(1:nx*(N+1))),nx); % history of solution TRAJECTORY accross horizon
auxSt.u0 = [u(:,2:end), u(:,end)];
end

function KK = ctrl_gains(k,N)
% Return control gains at k prediction horizon step
K0 = 1; % Start gain
KN = 150; % End Gain
steep = 2; % Steepness of exponential

%Qd = exp(steep*(-k+log(K0-KN)/steep)) + KN; % Starts big, goes down
Qd = exp(steep*(k-N+log(KN-K0)/steep)) + K0; % starts small, goes up
R = 0.01;
Qm = 2;

KK.Qd = Qd; KK.R = R; KK.Qm = Qm;
end

function [Obj] = Jfunc(X,Y,U,scale,D_KL,ny)
% Given X [nx, N+1], Y [ny_full, N] and U [nu, N] compute cost function for
% NMPC along prediction horizon
% Obj : Struct that contains scaled and unscaled objective function

% J = sum(||qD-qA||^2*Qk + ||u||^2*Rk, k=1:N) + Qm*LSE(-1, [D_KL(QD,QA)^2, ...])

N = size(U,2);
nx = size(X,1);
nu = size(U,1);

obj_d = 0; % distance objective function
obj_u = 0; % u objective function
unsObj_d = 0; % distance objective function
unsObj_u = 0; % u objective function
for k=1:N
    xk = X(:,k);
    ref = Y(:,k);
    uk = U(:,k);

    % 1) Intersect Rn ellipsoids with R3 hyperplane to get R3 ellipsoids
    Cx = eye(nx);
    Cx = Cx(:,4:end);
    Bx = zeros(size(Cx,2),1);
    % Cy = eye(ny);
    % Cy = Cy(:,4:end);
    % By = zeros(size(Cy,2),1);

    %alphA = v2m(ref(ny+1:end),ny);
    %QA = alphA'*alphA;
    % QA = v2m(ref(ny+1:end),ny);
    WA = v2m(ref(ny+1:end),3);

    % % Use Ridge Regression (Tikhonov Regression) if QA is close to singular
    % Tikh_coef = 1e-2;
    % QA = QA + Tikh_coef.*eye(ny);

    % % Sometimes eigenvalue is negative because previous matrix is close to singular. Just switch sign for now
    % [e_vecA,e_valA] = eig(QA);
    % QA = e_vecA*(sign(e_valA).*e_valA)*e_vecA';

    % [~,WA] = ellipIntersectPoly(zeros(ny,1),QA,Cy,By);

    rotD = v2m(xk(7:3*3+7-1),3); % Extract Rotation Matrix of sys D
    S = diag(scale(:,k));
    WD = rotD*S*rotD'; % Apply scaling

    wD = xk(1:3);
    wA = ref(1:3);

    KK = ctrl_gains(k,N);
    Qk = KK.Qd;
    Rk = KK.R.*eye(nu);

    % OR compute Kullback-Leibler Divergence
    D_KL(k) = kullback_leibler_div(wD,WD,wA,WA);
    %D_KL(k) = correlation_matrix_distance(wA,WA,wD,WD)^2;
    
    obj_u = obj_u + uk'*Rk*uk;
    unsObj_u = unsObj_u + uk'*uk; 
    obj_d = obj_d + (wA-wD)'*Qk*(wA-wD);
    unsObj_d = unsObj_d + (wA-wD)'*(wA-wD);
end
Qm = KK.Qm;
LSE_coef = -0.1;
%LSE_coef = -1;
%LSE_coef = -20;
unsObj_DKL = LSE(LSE_coef,D_KL);
obj_DKL = Qm*unsObj_DKL;

Obj.obj_d = obj_d;
Obj.obj_u = obj_u;
Obj.obj_DKL = obj_DKL;

Obj.unsObj_d = unsObj_d;
Obj.unsObj_u = unsObj_u;
Obj.unsObj_DKL = unsObj_DKL;

end