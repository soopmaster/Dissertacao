% Simple Rocket Simulation
% Definition of Parametope:
%   ((xo,alph,y))g := {x \in R^n : g(alph,x-xo)<=y}
% Controlled Embeddings:
%   dxo = f(xo,wo),   dalph = U(xo,alph,y), dy = E(xo,alph,y)
% Elipsoidal Reachable Sets:
%   g(alph,x-xo) = norm(alph*(x-xo))
%   U(xo,alph,y) = -alph*Df(xo)
%   E(xo,alph,y) = c*y
% where
%   c = max_i(max(eig((alph*(Mi-Df(xo))*alph^-1)'+ alph*(Mi-Df(xo))*alph^-1))
% and
%   Mi \in M (set of interval Mixed Jacobian Matrices)

close all;
clearvars;

% == Euler Axis-Angle == %
euler_ang = @(R) acos(1/2*(trace(R)-1)); % Euler Angle
euler_axis = @(R,Theta) 1/(2*sin(Theta))*(R'-R); % Euler Axis cross

% == DCM variation in time == %
% Adot = -crossprod(omegaB)*A

% == Parameters == %
    %            x
    %           ^
    %           |
    %           O
    %           x
    %           x
    % ---Oxxxxx<o>xxxxxO---> y
    %           O
    %           x
    %           x
    %           |

    % System is quadRotor

% === Sim Settings === %
ops = sdpsettings('verbose',0);
ops.plot.shade = 0.15;
ops.plot.linewidth=0.05;
ops.plot.wirestyle=':';
ops.plot.edgecolor="#EB4B2D";
simT = 5; % Sim time seconds
Ts = 0.01;
knum = simT/Ts;

% == Variables == %

% Variables for ODE: x = [p_vec; v; A_vec; omegaB]
% p_vec: Position of all points, relative to inertial
% v: Velocity of center mass, relative to inertial
% omegaB: Angular velocity, relative to body
% A_vec: DCM matrix of attitude in vector form 

% y = [px; py; pz; omegaB3]

% === Symbolic ===%
syms psym vsym  omegaBsym tauBsym [3 1]
syms Fsym 
syms Asym  [3 3]
Avecsym = reshape(Asym,[],1);
usym = [Fsym; tauBsym];
symvars=[psym; vsym; Avecsym; omegaBsym];
symvars_ext = [symvars;usym];

varmin = [-0.*ones(3,1); -0.*ones(3,1); -0.*ones(9,1); -0.*ones(3,1)];
varmax = -varmin;
varbounds = [varmin, varmax];
varbounds_init = varbounds;
Fmax = 20;
Fmin = 0.1;
tauBmax = 0.0007*pi;
ubounds_max = [Fmax; tauBmax.*ones(3,1)];
ubounds_min = [Fmin; -tauBmax.*ones(3,1)];
varbounds_ext = [varbounds; [ubounds_min, ubounds_max]];
varbounds_ext_init = varbounds_ext;

% % Info from "Design of Autonomous Quadcopter Using Orientation Sensor with
% %   Variations in Load Fulcrum Point"
% %   https://www.researchgate.net/publication/335719947_Design_of_Autonomous_Quadcopter_Using_Orientation_Sensor_with_Variations_in_Load_Fulcrum_Point
% L = 0.127;
% I = [0.008513; 0.008513; 0.015579]; % moment of inertia
% M = 0.990; % total mass of drone
% J = diag(I); % Inertia Matrix
% g = 9.8;
% 
% niu_v = 0.01;
% niu_w = 0.7;

M = 0.5;
g = 9.81;
L = 0.25;
I = [0.0196; 0.0196; 0.0264];
J = diag(I);
niu_v = 0;
niu_w = 0;

dxsym = quadDroneODE(symvars,usym,I,M,niu_v,niu_w);

nx = 18;
nu = 4;
ny = 6;

Df = jacobian(dxsym, symvars);
Df_vec = reshape(Df,1,[]);

% == Initial Conditions == %

p_cmB = [0;0;0]; % Position of centermass relative to body

p_cm0 = [0;0;0];
v0 = [0;0;0];
omega0B = [0;0;0];
A0 = eye(3);
A = A0;
A0_vec = reshape(A0,[],1);

x0 = [p_cm0; v0; A0_vec; omega0B];

p_NB = [L;0;0];
p_SB = [-L;0;0];
p_EB = [0;L;0];
p_OB = [0;-L;0];

p_vecB = [p_NB; p_EB; p_SB; p_OB];

x = zeros(size(x0,1),knum);
x(:,1) = x0;

p_cm = x(1:3,:);
A_vecs = x(7:15,:);

% == Control Action == %
% u = [FB; tauB]

F = M*g; % total thrust force [N]
%F = Fmin;
dF = 0;
tauB = [0.1;0;0]; % input torque [N/m]

u_val = [F; tauB];

u = u_val.*ones(nu,knum);

%===v PD CONTROLLER PARAMS v===%
K_val = 50;
K1 = 1*K_val.*eye(3); % dddp
K2 = 1*K_val.*eye(3); % ddp
K3 = 1*K_val.*eye(3);   % dp
K4 = 1*K_val.*eye(3);        % p
K5 = K_val.*eye(1); % dyaw
K6 = K_val.*eye(1); % yaw

pinit_ref = x0(1:3);
p1_ref = [0;0;1];
p2_ref = [5;-5;2];
p_ref = [p1_ref, p2_ref];
r_vec = p1_ref-pinit_ref;
r_norm = norm(r_vec);

% === Trajectory ===
%[r,w_star] = trajectory3D([0 0 0; 2 0 15; 30 0 30]',[0.1 0 1; 0.2 0 1; 0.1 0 0.1]',[0 0.5 1],100);
v_aux = p2_ref-pinit_ref;
v_aux = [v_aux(1:2)/norm(v_aux(1:2)); 1];
v_aux = 2.*v_aux/norm(v_aux);
v_aux = [0;0;0];
a_aux = p2_ref-pinit_ref;
a_aux = [a_aux(1:2)/norm(a_aux(1:2)); 0.5];
a_aux = 0.5.*a_aux/norm(a_aux);
a_aux = [0;0;0];
j_aux = p2_ref-pinit_ref;
j_aux = [j_aux(1:2)/norm(j_aux(1:2)); 0.1];
j_aux = 0.1.*j_aux/norm(j_aux);
j_aux = [0;0;0];
[r,dr,ddr,dddr] = trajectory3Dv2([0 1 4],Ts,[pinit_ref, p_ref(:,1), p_ref(:,2)],[[0;0;0.1], [0;0;0], v_aux], [[0;0;0.1], [0;0;0], a_aux], [[0;0;0.1], [0;0;0], j_aux]);
r = [r, r(:,end).*ones(size(r,1), knum-size(r,2))];
dr = [dr, dr(:,end).*ones(size(dr,1), knum-size(dr,2))];
ddr = [ddr, ddr(:,end).*ones(size(ddr,1), knum-size(ddr,2))];
dddr = [dddr, dddr(:,end).*ones(size(dddr,1), knum-size(dddr,2))];

yaw_ref = 0;
dyaw_ref = 0;

p1_ref = zeros(3,knum);
dp_ref = zeros(3,knum);
ddp_ref = zeros(3,knum);
dddp_ref = zeros(3,knum);

ddddp_des = ones(3,knum);
dddp = ones(3,knum);
ddp = ones(3,knum);

% === Parametope === %
enable_prmtp = 1;
inner_step = 20;
% Initialize range_struct
[range_struct,constDfMx,code,hasVarMx] = initRangeStruct(Df,nx,symvars_ext,varbounds_ext);

% Initialize alph
%delta_p = 1.0e-05.*([0.1161; 0.1157; 0]+1.0e-5.*ones(3,1));
delta_p = [0.0808; 0.0800; 0.1069];
del_pinv = [1/delta_p(1); 1/delta_p(2); 1/delta_p(2)];

delta_v = 10^-3.*[4.6354e-04; 4.6195e-04; 0.2676*2];
%delta_v = [1.4307; 1.4161; 1.0666];
del_vinv = [1/delta_v(1); 1/delta_v(2); 1/delta_v(2)];

delta_A = [0.0000	0.0000	0.0092;
           0.0000	0.0000	0.0091;
           0.0092	0.0091	0.0000]+1.0e-5;
% delta_A = [0.5910	0.5930	0.9966;
%            0.5930	0.5950	0.9932;
%            0.9966	0.9932	0.5929];
del_Ainv = [1/delta_A(1,1) 1/delta_A(1,2) 1/delta_A(1,2);
            1/delta_A(2,1) 1/delta_A(2,2) 1/delta_A(2,2);
            1/delta_A(3,1) 1/delta_A(3,2) 1/delta_A(3,2);];

delta_w = [1.8262; 1.8325; 0.0000]+1.0e-5.*ones(3,1);
%delta_w = [17.2183; 17.1595; 0.00001];
del_winv = [1/delta_w(1); 1/delta_w(2); 1/delta_w(2)];

alph0_pre = diag(2.*[del_pinv; del_vinv; reshape(del_Ainv,[],1); del_winv])^2;
%alph = diag([1/2.*delta_p; 1/2.*delta_v; 1/2.*reshape(delta_A,[],1); 1/2.*delta_w]);
%alph = diag([1/0.1.*ones(3,1); 1/0.01.*ones(3,1); reshape(1/0.1.*eye(3)+1/(1.0e-5).*ones(3,3),[],1); 1/2.*del_winv]);

rad_sqr = min(eig(alph0_pre))^-1;

alph0=rad_sqr.*alph0_pre;

alph_vec = zeros(nx^2,knum/inner_step);
alph_vec(:,1) = reshape(alph0,[],1);

ysim = zeros(1,knum/inner_step);
ysim(1) = rad_sqr;

% === Start Figure === %
enable_disp = 0;
figure;
hold on;
grid on
axis equal
view(3)
xlabel('$x$','Interpreter','latex');
ylabel('$y$','Interpreter','latex');
zlabel('$z$','Interpreter','latex');
%plot3(NaN, NaN, NaN, 'DisplayName','Rocket Body','Color','b'); % Rocket Body
%plot3(NaN, NaN, NaN, 'DisplayName','Thruster Plume','Color','m'); % Thruster Plume
plot3(p_ref(1,:),p_ref(2,:),p_ref(3,:),'xr','DisplayName','target');
%axis([-20 20 -20 20 0 20]);
%grid on;
plot3(r(1,:), r(2,:), r(3,:) ,'r','DisplayName','r_T');
quiver3(r(1,:), r(2,:), r(3,:), dr(1,:), dr(2,:), dr(3,:),1,'m','DisplayName','dr_T');
quiver3(r(1,:), r(2,:), r(3,:), ddr(1,:), ddr(2,:), ddr(3,:),1,'g','DisplayName','ddr_T');
%legend;
uvw = [0;0;0];

Xvec = [p_vecB(1:3:12)];
Yvec = [p_vecB(2:3:12)];
Zvec = [p_vecB(3:3:12)];
bodyPlot = patch(Xvec, Yvec, Zvec, [0 0 1],'DisplayName','Drone Body');
thrustPlot = quiver3(p_cm(1,1), p_cm(1,2), p_cm(1,3), 0,0,0,1/tauBmax,'m','DisplayName','F');
dthrustPlot = quiver3(p_cm(1,1), p_cm(1,2), p_cm(1,3), 0,0,0,1,'m','DisplayName','tau');
%vPlot = quiver3(p_cm(1,1), p_cm(1,2), p_cm(1,3), 0,0,0,1,'r','DisplayName','v');
%vdesPlot = quiver3(p_cm(1,1), p_cm(1,2), p_cm(1,3), 0,0,0,1,'g','DisplayName','v_{des}');
pathPlot = animatedline(p_cm(1,1),p_cm(2,1), p_cm(3,1),'LineWidth',2.5, 'Color','b','LineStyle','-','DisplayName','Path');
legend;

% 18-Sphere Coords
%sph_coords = nsphere(nx,10);

% === Simulation === %
telapse_total = tic;
for k=1:knum-1 
    disp(['=== Computing ', num2str(k*Ts), '/', num2str(simT), 's ===']);
    % === Calculate Control action
    
    [ddddp_des(:,k),dddp(:,k),ddp(:,k),ddF_des,domega_des,tauB] = snap_controller(K1,K2,K3,K4,K5,K6,r(:,k),dr(:,k),ddr(:,k),dddr(:,k),yaw_ref,dyaw_ref,x(:,k),F,dF,I,M);

    % Integrate ddF_des to obtain dF_des and F_des
    Fout = runge_kutta4(@(y) [0 1; 0 0]*y+[0;1]*ddF_des, Ts, [F;dF]);
    F = Fout(1);
    dF = Fout(2);

    % Clamp u to bounds
    F=min(max(F,Fmin),Fmax);
    tauBnorm = norm(tauB);
    if tauBnorm == 0
        tau = [0;0;0];
    elseif tauBnorm > tauBmax
        scale = tauBmax/tauBnorm;
        tauB = scale.*tauB;
    end

    u(1,k) = F;
    u(2:end,k) = tauB;

    % if k>=1
    %    u(1,k) = Fmax;
    %    u(2:end,k) = zeros(3,1);
    % end

    % === Propagate Parametope === %
    if enable_prmtp && ~mod(k+inner_step-1,inner_step)
        inner_k = (k+inner_step-1)/inner_step;
        telapse_k =tic;
        if enable_disp
            disp('----------------');
            disp(['Computing ', num2str(k), ' time interval ']);
        end

        %Dfxo = reshape(quadDrone_jac([x(:,k);u(:,k)]),1,[]);
        Dfxo = quadDrone_jac([x(:,k);u(:,k)],I,M,niu_v,niu_w);

        alph = reshape(alph_vec(:,inner_k),nx,[]);
        alphk_1 = runge_kutta4(@(y) prmtpU(y,-Dfxo,nx),Ts*inner_step,alph);
        alph_vec(:,inner_k+1) = reshape(alphk_1,[],1);

        varbounds_ext(:,1) = varbounds_ext_init(:,1)+[x(:,k);u(:,k)];
        varbounds_ext(:,2) = varbounds_ext_init(:,2)+[x(:,k);u(:,k)];

        % Calculate c just once outside of runge-kutta
        [alph_err,parametope_c] = prmtpEllips_c(range_struct,symvars_ext,varbounds_ext,-Dfxo,alph_vec(:,inner_k),nx,code,hasVarMx,constDfMx);
        %[a_eV,a_eE] = eig(alph_err);
        ysim(inner_k+1) = runge_kutta4(@(y) prmtpE(y,parametope_c^-2), Ts*inner_step, ysim(inner_k));
        %ysim(inner_k+1) = ysim(inner_k);

        disp(['Computing time for interval: ', num2str(toc(telapse_k)), 's']);
        if enable_disp
            %X.G = rand(2,m);
            rad_sqr = ysim(inner_k);
            alph_pre = rad_sqr^-1.*alph;
            X.G = inv(sqrtm(alph_pre));
            %X.c = zeros(2,1);
            X.c = x(:,k);
            X.A = zeros(0,nx);
            %X.A = ones(1,m);
            X.b = zeros(0,1);
            %X.b = ones(1,1);
            X.idx = [nx];
            X.type = [2];
            [Fc,p] = compileCCG(X);
            legend('off');
            plot(Fc,p(1:3),ops.plot.edgecolor,300,ops);
            plt_cell=plot(Fc,p(1:3));
            drawnow;
        end
    end

    %   === Calculate System ===
    x(:,k+1) = runge_kutta4(@(y) quadDroneODE(y,u(:,k),I,M,niu_v,niu_w), Ts, x(:,k));

    p_cm = x(1:3,k);
    v = x(4:6,k);
    A_vec = x(7:15,k);
    omegaB = x(16:end,k);

    A = reshape(A_vec,3,[]);

    p_vec = kron(eye(4),A)*p_vecB + [p_cm; p_cm; p_cm; p_cm];
    p_N = p_vec(1:3);
    p_E = p_vec(4:6);
    p_S = p_vec(7:9);
    p_O = p_vec(10:12);

    uvw_thrust = A*u(2:end,k);
    uvw_dthrust = A*u(1,k)*[0;0;1/(M*g*2)];

    Xvec = [p_vec(1:3:12)];
    Yvec = [p_vec(2:3:12)];
    Zvec = [p_vec(3:3:12)];
    bodyPlot.XData = Xvec;
    bodyPlot.YData = Yvec;
    bodyPlot.ZData = Zvec;

    set(thrustPlot,'XData',p_cm(1),'YData',p_cm(2),'ZData',p_cm(3),'UData',uvw_thrust(1),'VData',uvw_thrust(2),'WData',uvw_thrust(3));
    set(dthrustPlot,'XData',p_cm(1),'YData',p_cm(2),'ZData',p_cm(3),'UData',uvw_dthrust(1),'VData',uvw_dthrust(2),'WData',uvw_dthrust(3));
    %set(vPlot,'XData',p_cm(1),'YData',p_cm(2),'ZData',p_cm(3),'UData',v(1),'VData',v(2),'WData',v(3));
    %set(vdesPlot,'XData',p_cm(1),'YData',p_cm(2),'ZData',p_cm(3),'UData',dp_ref(1,k),'VData',dp_ref(2,k),'WData',dp_ref(3,k));
    addpoints(pathPlot,p_cm(1), p_cm(2), p_cm(3));
    drawnow;
    %pause(Ts);
end

%plot3(x(1,1:knum-1),x(2,1:knum-1), x(3,1:knum-1),'-og',p(1,1:knum-1),p(2,1:knum-1), p(3,1:knum-1), 'xr',tip(1,1:knum-1), tip(2,1:knum-1), tip(3,1:knum-1),'^k');
%plot3(p_ref(1,:),p_ref(2,:),p_ref(3,:),'-r');
%plot3(p2_ref(1),p2_ref(2),p2_ref(3),'-r');
%legend('Location','best','Interpreter','latex');
hold off;
disp(['Total sim computing time: ', num2str(toc(telapse_total)), 's']);

% Plot snap_controller data for tuning (p, v, p_ref, dp_ref)
tplot = 0:Ts:simT-Ts;
figure;
hold on;
plot(tplot,x(3,:),"DisplayName","x(3)");
%plot(tplot,p1_ref(3).*ones(1,knum),"DisplayName","ref");
plot(tplot,r(3,:),"DisplayName","p_{ref}");
legend
hold off;

