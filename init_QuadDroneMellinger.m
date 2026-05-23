function vehiStruct = init_QuadDroneMellinger(Ts)
%INIT_QUADDRONE outputs vehiStruct of quadDrone
%   
% vehiStruct : struct (struct with components of quadDrone)
% MAIN (interfaces with main.m)
%   .nx = nx;
%   .nu = nu;
%   .ny = ny;
%   .x0 = x0;
%   .u0 = u0;
%   .ysim0 = rad_sqr;
%   .alph0 = alph0;
%   .rStrt = range_struct;
%   .rConstDfMx = constDfMx;
%   .rCode = code;
%   .rHasVarMx = hasVarMx;
%   .vbds_init
%   .symvars
%   .dynFunc
%   .jacFunc
%   .contrlFunc
%   .bodyFunc
% AUXILIARY (does not interface with main.m)
%   .Ts
%   .M : num                       (mass)
%   .L : num                        (arm length)
%   .I : []                    (angular momentum matrix)
%   .niu_v : num                (drag coefficient for translational velocity)
%   .niu_w : num                (drag coefficient for angular velocity)
%   .F
%   .dF
%   .Fmax
%   .Fmin
%   .tauBmax
%   .K1..6
%   .trajFunc

% == Variables == %

% Extended quadDroneODE with F and dF to accomodate for Snap controller
% Variables for ODE: x = [p_cm; v; omegaB; A_vec; F; dF]
% p_cm: Position of center mass, relative to inertial
% v: Velocity of center mass, relative to inertial
% omegaB: Angular velocity, relative to body
% A_vec: DCM matrix of attitude in vector form

% y = [px; py; pz; omegaB3]

nx = 18;
nu = 4;
ny = 6;
vehiStruct.nx = nx;
vehiStruct.nu = nu;
vehiStruct.ny = ny;

vehiStruct.Ts = Ts;

% === Parameters of Dynamics === %
M = 0.5;
vehiStruct.M = M;
g = 9.81;
L = 0.25;
vehiStruct.L = L;
I = [0.0196; 0.0196; 0.0264];
vehiStruct.I = I;
niu_v = 0;
niu_w = 0;
vehiStruct.niu_v = niu_v;
vehiStruct.niu_w = niu_w;

% == Initial Conditions == %
p_cm0 = [0;0;0];
v0 = [0;0;0];
omega0B = [0;0;0];
A0 = eye(3);
A0_vec = reshape(A0,[],1);

x0 = [p_cm0; v0; A0_vec; omega0B];
vehiStruct.x0 = x0;

u0 = [0; zeros(3,1)];
vehiStruct.u0 = u0;

% === Symbolic ===%
syms psym vsym  omegaBsym tauBsym Isym [3 1]
syms Fsym dFsym ddFsym niu_vsym niu_wsym Msym
syms Asym  [3 3]
Avecsym = reshape(Asym,[],1);
usym = [Fsym; tauBsym];
xsym=[psym; vsym; Avecsym; omegaBsym];
symvars_ext = [xsym;usym];
vehiStruct.symvars = symvars_ext;

dxsym = quadDroneODE(xsym,usym,I,M,niu_v,niu_w);
%dxsym = quadDroneSnapODE(xsym,usym,Isym,Msym,niu_vsym,niu_wsym);

Df = jacobian(dxsym, xsym);

% === Bounds of state and input variables == %
Fmax = 20;
vehiStruct.Fmax = Fmax;
Fmin = 0.0001;
vehiStruct.Fmin = Fmin;
varmin = [-0.*ones(3,1); -0.*ones(3,1); -0.*ones(9,1); -0.*ones(3,1)];
varmax = [0.*ones(3,1); 0.*ones(3,1); 0.*ones(9,1); 0.*ones(3,1)];
varbounds = [varmin, varmax];
%tauBmax = 0.0007*pi;
tauBmax = pi;
vehiStruct.tauBmax = tauBmax;
ubounds_max = [0; tauBmax.*ones(3,1)];
ubounds_min = [0; -tauBmax.*ones(3,1)];
varbounds_ext = [varbounds; [ubounds_min, ubounds_max]];
varbounds_ext_init = varbounds_ext;
vehiStruct.vbds_init = varbounds_ext_init;

% === Initialize rangeStruct for parametope === %
[range_struct,constDfMx,code,hasVarMx] = initRangeStruct(Df,nx,symvars_ext,varbounds_ext);
vehiStruct.rStrt = range_struct;
vehiStruct.rConstDfMx = constDfMx;
vehiStruct.rCode = code;
vehiStruct.rHasVarMx = hasVarMx;

% Initialize parametope parameters (alph0, ysim0)
delta_p = [0.0808; 0.0800; 0.1069];
del_pinv = [1/delta_p(1); 1/delta_p(2); 1/delta_p(2)];

delta_v = 10^-3.*[4.6354e-04; 4.6195e-04; 0.2676*2];
del_vinv = [1/delta_v(1); 1/delta_v(2); 1/delta_v(2)];

delta_A = [0.0000	0.0000	0.0092;
           0.0000	0.0000	0.0091;
           0.0092	0.0091	0.0000]+1.0e-5;
del_Ainv = [1/delta_A(1,1) 1/delta_A(1,2) 1/delta_A(1,2);
            1/delta_A(2,1) 1/delta_A(2,2) 1/delta_A(2,2);
            1/delta_A(3,1) 1/delta_A(3,2) 1/delta_A(3,2);];

delta_w = [1.8262; 1.8325; 0.0000]+1.0e-5.*ones(3,1);
del_winv = [1/delta_w(1); 1/delta_w(2); 1/delta_w(2)];

alph0_pre = diag(2.*[del_pinv; del_vinv; reshape(del_Ainv,[],1); del_winv])^2;

rad_sqr = min(eig(alph0_pre))^-1;
vehiStruct.ysim0 = rad_sqr;

alph0=rad_sqr.*alph0_pre;
vehiStruct.alph0 = alph0;

% === Initialize Auxiliary variables for mellinger controller === %
K_val = 1;
KK.Kp = 1*K_val.*eye(3); % dddp
KK.Kv = 1*K_val.*eye(3); % ddp
KK.KR = 1*K_val.*eye(3);   % dp
KK.Kom = 1*K_val.*eye(3);        % p
vehiStruct.KK = KK;

% === Specific Function Handles === %
vehiStruct.dynFunc = @(y,u) quadDroneODE(y,u,I,M,niu_v,niu_w);
vehiStruct.jacFunc = @(x,u) quadDrone_jac([x;u],I,M,niu_v,niu_w);
vehiStruct.trajFunc = @(t,pts) quadDrone_traj(t,pts,vehiStruct);
vehiStruct.contrlFunc = @(x,ref) quadController(x,ref,vehiStruct);
vehiStruct.bodyFunc = @(x,p_vec) kron(eye(4),reshape(x(7:15),3,[]))*p_vec + [x(1:3); x(1:3); x(1:3); x(1:3)];
end

function ref = quadDrone_traj(t,pts,vehi)
    order = size(pts,1)/3;
    switch order
        case 2
            [r,dr,ddr,dddr] = trajectory3Dv2(t, vehi.Ts, pts(1:3,:), pts(4:6,:));
        case 4
            [r,dr,ddr,dddr] = trajectory3Dv2(t, vehi.Ts, pts(1:3,:), pts(4:6,:), pts(7:9,:), pts(10:12,:));
    end
    %[r,dr,ddr,dddr] = trajectory3Dv2(t, vehi.Ts, pts(1:3,:), pts(4:6,:), pts(7:9,:), pts(10:12,:));
    ref = [r; dr; ddr; dddr];
end

function u = quadController(x, ref, vehi)
    r = ref(1:3);
    dr = ref(4:6);
    ddr = ref(7:9);
    dddr = ref(10:12);
    yaw_ref = ref(13);
    dyaw_ref = ref(14);

    u = mellinger_controller(vehi.KK,r,dr,ddr,dddr,yaw_ref,dyaw_ref,x,vehi.I,vehi.M);

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