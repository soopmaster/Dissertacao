function vehiStruct = init_QuadDroneNMPC(p0,Ts)
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

% Variables for ODE: x = [p_cm; v; A_vec; omegaB]
% p_cm: Position of center mass, relative to inertial
% v: Velocity of center mass, relative to inertial
% omegaB: Angular velocity, relative to body
% A_vec: DCM matrix of attitude in vector form

% y = [px; py; pz; omegaB3]

nx = 18;
nx_full = nx+nx^2;
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
p_cm0 = p0;
v0 = [0;0;0];
omega0B = [0;0;0];
A0 = eye(3);
A0_vec = reshape(A0,[],1);

x0 = [p_cm0; v0; A0_vec; omega0B];

% Initialize parametope parameters (alph0, ysim0)
delta_p = 0.5e-1.*[0.0808; 0.0800; 0.1069];

delta_v = 10^-3.*[4.6354e-04; 4.6195e-04; 0.2676*2];

delta_A = [0.0000	0.0000	0.0092;
           0.0000	0.0000	0.0091;
           0.0092	0.0091	0.0000]+1.0e-5;
% delta_A = [0.0000	0.0000	0.0092;
%            0.0000	0.0000	0.0091;
%            0.0092	0.0091	0.0000];

delta_w = [1.8262; 1.8325; 0.0000]+1.0e-5.*ones(3,1);
% delta_w = [1.8262; 1.8325; 0.0000];

% alph0_pre = diag(2.*[del_pinv; del_vinv; reshape(del_Ainv,[],1); del_winv])^2;
alph0inv = diag([delta_p; delta_v; m2v(delta_A); delta_w]);
alph0_pre = inv(alph0inv);
% alph'*alph = Q^-1 => (alph'*alph)^-1 = Q, where Q's eigenvals are the
% squared axis scaling coefs
% therefore alph's eigenvals are the inverse of the axis scaling coefs
% x'*(alph0_pre'*alph0_pre)*x <= 1

%eigAlph0 = eig(alph0_pre);
%rad = max(eigAlph0(1:3))^-1;
rad = max(delta_p);
vehiStruct.ysim0 = rad^2;

alph0=rad.*alph0_pre;
vehiStruct.alph0 = alph0;
% x'*(alph0'*alph0)*x <= ysim0

% vehiStruct.x0 = [x0;m2v(alph0)];
vehiStruct.x0 = x0;

u0 = [M*g; zeros(3,1)];
vehiStruct.u0 = u0;

% === Plot Variables === %
L = 0.5;
p_NB = [L;0;0];
p_SB = [-L;0;0];
p_EB = [0;L;0];
p_OB = [0;-L;0];
p_vec = [p_NB; p_EB; p_SB; p_OB];
Xvec = [p_vec(1:3:12)];
Yvec = [p_vec(2:3:12)];
Zvec = [p_vec(3:3:12)];

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
% Fmax = 20;
Fmax = 30;
vehiStruct.Fmax = Fmax;
Fmin = 0.0001;
vehiStruct.Fmin = Fmin;
% varmin = [-0.*ones(3,1); -0.*ones(3,1); -0.*ones(9,1); -0.*ones(3,1)];
% varmax = [0.*ones(3,1); 0.*ones(3,1); 0.*ones(9,1); 0.*ones(3,1)];
varmin = [-delta_p; -delta_v; -0.*ones(9,1); -0.*ones(3,1)];
varmax = [delta_p; delta_v; 0.*ones(9,1); 0.*ones(3,1)];
varbounds = [varmin, varmax];
%tauBmax = 0.0007*pi;
tauBmax = 2*pi;
% tauBmax = 50;
vehiStruct.tauBmax = tauBmax;
ubounds_max = [Fmax; tauBmax.*ones(3,1)];
ubounds_min = [Fmin; -tauBmax.*ones(3,1)];
varbounds_ext = [varbounds; [ubounds_min, ubounds_max]];
varbounds_ext_init = varbounds_ext;
vehiStruct.vbds_init = varbounds_ext_init;

% === Initialize rangeStruct for parametope === %
[range_struct,constDfMx,code,hasVarMx] = initRangeStruct(Df,nx,symvars_ext,varbounds_ext);
vehiStruct.rStrt = range_struct;
vehiStruct.rConstDfMx = constDfMx;
vehiStruct.rCode = code;
vehiStruct.rHasVarMx = hasVarMx;

% === Specific Function Handles === %
ODEdynFunc = @(y,u) quadDroneODE(y,u,I,M,niu_v,niu_w);
localJacFunc = @(x,u) quadDrone_jac([x;u],I,M,niu_v,niu_w);
localDynFunc = @(y,u) fullDynFunc(@(x) ODEdynFunc(x,u), @(alph,x) prmtpU(alph,-localJacFunc(x,u)),y,nx); % rk4 system r.h.s
vehiStruct.dynFunc = @(y,u) localDynFunc(y,u);
vehiStruct.partDynFunc = @(y,u) quadDroneODE(y,u,I,M,niu_v,niu_w);
vehiStruct.jacFunc = @(x,u) quadDrone_jac([x(1:nx);u],I,M,niu_v,niu_w);
vehiStruct.trajFunc = @(t,pts) quadDrone_traj(t,pts,vehiStruct);
vehiStruct.contrlFunc = @(x,ref) quadController(x,ref,vehiStruct);
vehiStruct.bodyPlot = @() initBodyPlot(Xvec,Yvec,Zvec,p0);
vehiStruct.plotBodyFunc = @(x,bPlotSt) quadDrone_plotBody(x,bPlotSt,p_vec);
end

function bodyPlotSt = initBodyPlot(Xvec,Yvec,Zvec,p0)
bodyPlotSt.patch = patch(Xvec+p0(1), Yvec+p0(2), Zvec+p0(3), [0 0 1]);
bodyPlotSt.quiver = quiver3(p0(1),p0(2),p0(3), 0,0,1, 'LineWidth',1,'Color','k');
end

function bodyPlotSt = quadDrone_plotBody(x,bodyPlotSt,p_vec)
R = v2m(x(7:6+3^2),3); % Rotation Matrix
p_vecnew = kron(eye(4),R)*p_vec + [x(1:3); x(1:3); x(1:3); x(1:3)];
Xvec = [p_vecnew(1:3:12)];
Yvec = [p_vecnew(2:3:12)];
Zvec = [p_vecnew(3:3:12)];

% Update patch
bodyPlotSt.patch.XData = Xvec;
bodyPlotSt.patch.YData = Yvec;
bodyPlotSt.patch.ZData = Zvec;

% Update quiver
bodyPlotSt.quiver.XData = x(1);
bodyPlotSt.quiver.YData = x(2);
bodyPlotSt.quiver.ZData = x(3);

uvw = R*[0;0;1];
bodyPlotSt.quiver.UData = uvw(1);
bodyPlotSt.quiver.VData = uvw(2);
bodyPlotSt.quiver.WData = uvw(3);

end

function ref = quadDrone_traj(t,pts,vehi)
    %pts = pts(1:end-2,:); % Ignore Yaw for now...
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