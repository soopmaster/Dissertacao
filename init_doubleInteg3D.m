function vehiStruct = init_doubleInteg3D(p0in,Ts)
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

% Variables for ODE: x = [p; v; ang]
% p: Position of center mass, relative to inertial
% v: Velocity of center mass, relative to inertial
% ang: Angular position, relative to body

% y = [px; py; pz; omegaB3]

nx = 3*3;
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
p0 = p0in;
v0 = [0;0;0];
ang0 = [0;0;0];

x0 = [p0; v0; ang0];

% Initialize parametope parameters (alph0, ysim0)
% delta_p = [0.0808; 0.0800; 0.1069];
% del_pinv = [1/delta_p(1); 1/delta_p(2); 1/delta_p(2)];
% 
% delta_v = 10^-3.*[4.6354e-04; 4.6195e-04; 0.2676*2];
% del_vinv = [1/delta_v(1); 1/delta_v(2); 1/delta_v(2)];
% 
% delta_A = [0.0000	0.0000	0.0092;
%            0.0000	0.0000	0.0091;
%            0.0092	0.0091	0.0000]+1.0e-5;
% del_Ainv = [1/delta_A(1,1) 1/delta_A(1,2) 1/delta_A(1,2);
%             1/delta_A(2,1) 1/delta_A(2,2) 1/delta_A(2,2);
%             1/delta_A(3,1) 1/delta_A(3,2) 1/delta_A(3,2);];
% 
% delta_w = [1.8262; 1.8325; 0.0000]+1.0e-5.*ones(3,1);
% del_winv = [1/delta_w(1); 1/delta_w(2); 1/delta_w(2)];
% 
% % alph0_pre = diag(2.*[del_pinv; del_vinv; reshape(del_Ainv,[],1); del_winv])^2;
% alph0_pre = diag(2.*[del_pinv; del_vinv; reshape(del_Ainv,[],1); del_winv]);
alph0_pre = inv(0.5.*diag([1,1,1, 0.5,0.5,0.5, pi/6,pi/6,pi/6]));

Q = alph0_pre'*alph0_pre; % Ellipsoid shape matrix (x'*Q^-1*x <= 1)
rad_sqr = max(eig(inv(Q)))^-1;
vehiStruct.ysim0 = rad_sqr;

alph0=sqrt(rad_sqr).*alph0_pre;
vehiStruct.alph0 = alph0;

vehiStruct.x0 = [x0;m2v(alph0)];

u0 = [g; zeros(3,1)];
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
syms psym vsym  angsym u_angsym [3 1]
syms u_asym
usym = [u_asym; u_angsym];
xsym=[psym; vsym; angsym];
symvars_ext = [xsym;usym];
vehiStruct.symvars = symvars_ext;

dxsym = dIntegrator3DODE(xsym,usym);
%dxsym = quadDroneSnapODE(xsym,usym,Isym,Msym,niu_vsym,niu_wsym);

Df = jacobian(dxsym, xsym);

% === Bounds of state and input variables == %
Fmax = 20;
vehiStruct.Fmax = Fmax;
Fmin = 0.0001;
vehiStruct.Fmin = Fmin;
varmin = [-0.*ones(3,1); -0.*ones(3,1); -0.*ones(3,1)];
varmax = [0.*ones(3,1); 0.*ones(3,1); 0.*ones(3,1)];
varbounds = [varmin, varmax];
%tauBmax = 0.0007*pi;
tauBmax = pi;
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
ODEdynFunc = @(y,u) dIntegrator3DODE(y,u);
localJacFunc = @(x,u) dInteg3D_jac([x;u]);
localDynFunc = @(y,u) fullDynFunc(@(x) ODEdynFunc(x,u), @(alph,x) prmtpU(alph,-localJacFunc(x,u),nx),y,nx); % rk4 system r.h.s
vehiStruct.dynFunc = @(y,u) localDynFunc(y,u);
vehiStruct.jacFunc = @(x,u) dInteg3D_jac([x(1:nx);u]);
vehiStruct.trajFunc = @(t,pts) quadDrone_traj(t,pts,vehiStruct);
vehiStruct.contrlFunc = @(x,ref) quadController(x,ref,vehiStruct);
vehiStruct.bodyPlot = @() patch(Xvec, Yvec, Zvec, [0 0 1]);
vehiStruct.plotBodyFunc = @(x,patch) quadDrone_plotBody(x,patch,p_vec);
end

function p = quadDrone_plotBody(x,p,p_vec)
p_vecnew = kron(eye(4),reshape(x(7:15),3,[]))*p_vec + [x(1:3); x(1:3); x(1:3); x(1:3)];
Xvec = [p_vecnew(1:3:12)];
Yvec = [p_vecnew(2:3:12)];
Zvec = [p_vecnew(3:3:12)];

p.XData = Xvec;
p.YData = Yvec;
p.ZData = Zvec;

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