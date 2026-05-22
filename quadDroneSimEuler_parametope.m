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

simT = 2; % Sim time seconds
Ts = 0.01;
knum = simT/Ts;

% === Symbolic ===%
syms psym vsym ssym omegaBsym tauBsym [3 1]
syms Fsym 
usym = [Fsym; tauBsym];
symvars=[psym; vsym; ssym; omegaBsym];
symvars_ext = [symvars;usym];

varmin = [-0.*ones(3,1); -0.*ones(3,1); -0.*ones(3,1); -0.*ones(3,1)];
varmax = -varmin;
varbounds = [varmin, varmax];
varbounds_init = varbounds;
Fmax = 15;
tauBmax = pi/4;
ubounds = [Fmax; tauBmax.*ones(3,1)];
varbounds_ext = [varbounds; [-ubounds, ubounds]];
varbounds_ext_init = varbounds_ext;

% L = 1; % arm length
% R = 0.25;
% m = 0.1; % area mass of drone
% I = [2*m*R^2*((L+R/2)^3-(R/2)^3)/3; 2*m*R^2*((L+R/2)^3-(R/2)^3)/3; m*R*(2*L+R)*2*(R/2)^3/3]; % moment of inertia
% M = m*R^2*(2*L+R); % total mass of drone
% J = diag(I); % Inertia Matrix
% g = 9.8;

% Info from "Design of Autonomous Quadcopter Using Orientation Sensor with
%   Variations in Load Fulcrum Point"
%   https://www.researchgate.net/publication/335719947_Design_of_Autonomous_Quadcopter_Using_Orientation_Sensor_with_Variations_in_Load_Fulcrum_Point
L = 0.127;
I = [0.008513; 0.008513; 0.015579]; % moment of inertia
M = 0.990; % total mass of drone
J = diag(I); % Inertia Matrix
g = 9.8;

niu_v = 0.01;
niu_w = 0.7;

% aBsym = 1/M*Fsym.*[0;0;1]; %  Linear acceleration in relation to body
%
% dpsym = vsym;
% dvsym = Asym*aBsym+[0;0;-g]-niu_v.*vsym;
% dAsym = crossprodmat(omegaBsym)*Asym;
% domegaBsym = J\(tauBsym-crossprodmat(omegaBsym)*J*omegaBsym)-niu_w.*omegaBsym;
% 
% dAsym_vec = reshape(dAsym,[],1);
% 
% dxsym = [dpsym; dvsym; domegaBsym; dAsym_vec];


SOv = rotM(ssym);
SOome = [1, sin(ssym(1))*tan(ssym(2)), cos(ssym(1))*tan(ssym(2));
             0, cos(ssym(1)), -sin(ssym(1));
             0, sin(ssym(1))/cos(ssym(2)), cos(ssym(1))/cos(ssym(2))]; % Rotation Matrix for omega
fxsym = [SOv*vsym;
          crossprodmat(vsym)*omegaBsym+SOv'*[0;0;-g]-niu_v.*vsym;
          SOome*omegaBsym;
          -J\(crossprodmat(omegaBsym)*J*omegaBsym)-niu_w.*omegaBsym];

gFsym = 1/M.*SOv*[0;0;1];

gxsym = [zeros(3,4);
    gFsym, zeros(3,3);
    zeros(3,4);
    zeros(3,1), J^-1];

dxsym = fxsym + gxsym*usym;

nx = 12;
nu = 4;
ny = 6;

Df = jacobian(dxsym, symvars);
Df_vec = reshape(Df,1,[]);

dim = nx;

% == Variables == %

% == Variables == %
% Variables for ODE: x = [p_vec; v; omegaB; A_vec]
% p_vec: Position of all points, relative to inertial
% v: Velocity of center mass, relative to inertial
% s: euler angles of drone (roll, pitch, yaw)
% omegaB: Angular velocity, relative to body

% y = [px; py; pz; omegaB3]

Cvec = [1;1;1; 0;0;0; 0;0;1; zeros(9,1)];
Cmat = zeros(ny,nx);
Cmat(1,1) = 1;
Cmat(2,2) = 1;
Cmat(3,3) = 1;
Cmat(4,4) = 1;
Cmat(5,5) = 1;
Cmat(6,6) = 1;

% == Initial Conditions == %

p_cmB = [0;0;0]; % Position of centermass relative to body

p_cm0 = [0;0;0];
v0 = [0;0;0];
s0 = [0;0;0];
omega0B = [0;0;0];

x0 = [p_cm0; v0; s0; omega0B];

p_NB = [L;0;0];
p_SB = [-L;0;0];
p_EB = [0;L;0];
p_OB = [0;-L;0];

p_vecB = [p_NB; p_EB; p_SB; p_OB];

ref = [5;-5;20];

% == Control Action == %
% u = [FB; tauB]

FB = 15; % total thrust force [N]
tauB = [0.1;0;0]; % input torque [N/m]

u_val = [FB; tauB];

x = zeros(size(x0,1),knum);
x(:,1) = x0;

u = u_val.*ones(nu,knum);

p_cm = x(1:3,:);

% % === Trajectory ===
% r = trajectory3D([0 0 0; 2 0 15; 30 0 30]',[0.1 0 1; 0.2 0 1; 0.1 0 0.1]',[0 50 100],100);
% 
% % === Error ===
% ref_theta = zeros(1,knum);
% ref_phi = zeros(1,knum);
% ref_psi = zeros(1,knum);
% error = zeros(4,knum);
% 
% Kp  = [0.5;1];
% Kd = [.05;.05];

%===v PD CONTROLLER PARAMS v===%

Kp = 1*ones(nu,ny/2);
%Kp(3,:) = 0.001;
Kp(4,:) = 0.01;
Kd = 1*ones(nu,ny/2);
%Kd(3,:) = 0.001;
Kd(4,:) = 0.01;

KpF = 0.9;
KdF = 0.1;
KpA = 0.01;
KdA = 0.5;

% === Parametope === %
enable_prmtp = 1;
inner_step = 20;
% Initialize range_struct
[range_struct,constDfMx,code,hasVarMx] = initRangeStruct(Df,dim,symvars_ext,varbounds_ext);

% Initialize alph
delta_p = 1.0e-05.*([0.1161; 0.1157; 0]+1.0e-5.*ones(3,1));
del_pinv = [1/delta_p(1); 1/delta_p(2); 1/delta_p(2)];

delta_v = 10^-3.*([4.6354e-04; 4.6195e-04; 2.6699e-07]+1.0e-8.*ones(3,1));
del_vinv = [1/delta_v(1); 1/delta_v(2); 1/delta_v(2)];

delta_s = 1.0e-05.*([0.1161; 0.1157; 0]+1.0e-5.*ones(3,1));
del_sinv = [1/delta_p(1); 1/delta_p(2); 1/delta_p(2)];

delta_w = [1.8262; 1.8325; 0.0000]+1.0e-5.*ones(3,1);
del_winv = [1/delta_w(1); 1/delta_w(2); 1/delta_w(2)];

alph = diag([1/2.*del_pinv; 1/2.*del_vinv; 1/2.*del_sinv; 1/2.*del_winv]);
%alph = diag(1.0e-9.*ones(dim,1)+[1/2.*delta_p; 1/2.*delta_v; 1/2.*delta_s; 1/2.*delta_w]);
%alph = diag([1/0.1.*ones(3,1); 1/0.01.*ones(3,1); reshape(1/0.1.*eye(3)+1/(1.0e-5).*ones(3,3),[],1); 1/2.*del_winv]);

initial_offset=1^(1/2);
initial_alph=(alph)^(1/2);

alph_vec = zeros(dim^2,knum/inner_step);
alph_vec(:,1) = reshape(initial_alph,[],1);

ysim = zeros(1,knum/inner_step);
ysim(1) = initial_offset;

% === Start Figure === %
enable_disp = 1;
figure;
hold on;
grid on
axis equal
view(3)
xlabel('$x$','Interpreter','latex');
ylabel('$y$','Interpreter','latex');
zlabel('$z$','Interpreter','latex');
plot3(NaN, NaN, NaN, 'DisplayName','Rocket Body','Color','b'); % Rocket Body
plot3(NaN, NaN, NaN, 'DisplayName','Thruster Plume','Color','m'); % Thruster Plume
%axis([-20 20 -20 20 0 20]);
%grid on;
%plot3(r(1,:), r(2,:), r(3,:) ,'r');
uvw = [0;0;0];

Xvec = [p_vecB(1:3:12)];
Yvec = [p_vecB(2:3:12)];
Zvec = [p_vecB(3:3:12)];
bodyPlot = patch(Xvec, Yvec, Zvec, [0 0 1]);
thrustPlot = quiver3(p_cm(1,1), p_cm(1,2), p_cm(1,3), 0,0,0,1,'m');
dthrustPlot = quiver3(p_cm(1,1), p_cm(1,2), p_cm(1,3), 0,0,0,1,'m');
pathPlot = animatedline(p_cm(1,1),p_cm(2,1), p_cm(3,1), 'Color','r','LineStyle','-');

% === Simulation === %
telapse_total = tic;
for k=1:knum-1 
    % % === Calculate Control action
    % n_vec = A*[0;0;1]; % Normal vector
    % 
    % e_pos = ref-x(1:3,k); % positional error
    % de_pos = -x(4:6,k);
    % 
    % e_att = e_pos./norm(e_pos)-n_vec; % attitude error
    % de_att = -x(16:end,k);
    % 
    % F = KpF*e_pos(3) + KdF*de_pos(3);
    % tauB = KpA*norm(e_att)*crossprodmat([0;0;e_pos(3)])*[e_pos(1:2);0]+KdA.*de_att;
    % 
    % tauBnorm = norm(tauB);
    % 
    % % Clamp u to bounds
    % F=min(max(F,-Fmax),Fmax);
    % if tauBnorm > tauBmax
    %     scale = tauBmax/tauBnorm;
    %     tauB = scale.*tauB;
    % end

    u(1,k) = FB;
    u(2:end,k) = tauB;
    
    % === Calculate System ===

    % [t,y] = ode45(@(t,y) quadDroneODE(y,u(:,k),I,M),[0 Ts],x(:,k));
    % x(:,k+1) = y(end,:)';
    x(:,k+1) = runge_kutta4(@(y) quadDroneEulerODE(y,u(:,k),I,M,niu_v,niu_w), Ts, x(:,k));

    p_cm = x(1:3,k);
    v = x(4:6,k);
    s = x(7:9,k);
    omegaB = x(10:end,k);
    %A_vec = x(10:end,k);

    %A = reshape(A_vec,3,[]);
    A = rotM(s);

    p_vec = kron(eye(4),A)*p_vecB + [p_cm; p_cm; p_cm; p_cm];
    p_N = p_vec(1:3);
    p_E = p_vec(4:6);
    p_S = p_vec(7:9);
    p_O = p_vec(10:12);

    uvw_thrust = A*u(2:end,k);
    uvw_dthrust = A*u(1,k)*[0;0;1/(M*g)];

    Xvec = [p_vec(1:3:12)];
    Yvec = [p_vec(2:3:12)];
    Zvec = [p_vec(3:3:12)];
    bodyPlot.XData = Xvec;
    bodyPlot.YData = Yvec;
    bodyPlot.ZData = Zvec;

    % === Propagate Parametope === %
    if enable_prmtp && ~mod(k+inner_step-1,inner_step)
        inner_k = (k+inner_step-1)/inner_step;
        telapse_k =tic;
        if enable_disp
            disp('----------------');
            disp(['Computing ', num2str(k), ' time interval ']);
        end

        %Dfxo = reshape(quadDrone_jac([x(:,k);u(:,k)]),1,[]);
        Dfxo = quadDroneEuler_jac([x(:,k);u(:,k)],I,M,niu_v,niu_w);

        alph = reshape(alph_vec(:,inner_k),dim,[]);
        alphk_1 = runge_kutta4(@(y) prmtpU(y,-Dfxo,dim),Ts*inner_step,alph);
        alph_vec(:,inner_k+1) = reshape(alphk_1,[],1);

        varbounds_ext(:,1) = varbounds_ext_init(:,1)+[x(:,k);u(:,k)];
        varbounds_ext(:,2) = varbounds_ext_init(:,2)+[x(:,k);u(:,k)];

        % Calculate c just once outside of runge-kutta
        parametope_c = prmtpEllips_c(range_struct,symvars_ext,varbounds_ext,Dfxo,alph_vec(:,inner_k),dim,code,hasVarMx,constDfMx);
        ysim(inner_k+1) = runge_kutta4(@(y) prmtpE(y,parametope_c), Ts*inner_step, ysim(inner_k));
        %ysim(inner_k+1) = ysim(inner_k);

        disp(['Computing time for interval: ', num2str(toc(telapse_k)), 's']);
        if enable_disp
            m = dim;
            %X.G = rand(2,m);
            rad = ysim(inner_k);
            rad= 1;
            X.G = inv(sqrtm(1/(rad^2).*alph));
            %X.c = zeros(2,1);
            X.c = x(:,k);
            X.A = zeros(0,m);
            X.b = zeros(0,1);
            X.idx = [m];
            X.type = [2];
            [Fc,p] = compileCCG(X);
            plot(Fc,p(1:3));
            % plt_cell = plot(Fc,p(1:3));
            % plt_pt = plt_cell{1};
            % plot3(plt_pt(1,:),plt_pt(2,:),plt_pt(3,:),':g');
            %plot_ellipsoid_offset(x(1:3,k),alph(1:3,1:3),ysim(inner_k),'r');
            drawnow;
        end
    end

    set(thrustPlot,'XData',p_cm(1),'YData',p_cm(2),'ZData',p_cm(3),'UData',uvw_thrust(1),'VData',uvw_thrust(2),'WData',uvw_thrust(3));
    set(dthrustPlot,'XData',p_cm(1),'YData',p_cm(2),'ZData',p_cm(3),'UData',uvw_dthrust(1),'VData',uvw_dthrust(2),'WData',uvw_dthrust(3));
    addpoints(pathPlot,p_cm(1), p_cm(2), p_cm(3));
    drawnow;
    pause(Ts);
end

%plot3(x(1,1:knum-1),x(2,1:knum-1), x(3,1:knum-1),'-og',p(1,1:knum-1),p(2,1:knum-1), p(3,1:knum-1), 'xr',tip(1,1:knum-1), tip(2,1:knum-1), tip(3,1:knum-1),'^k');
%legend('Location','best','Interpreter','latex');
hold off;
disp(['Total sim computing time: ', num2str(toc(telapse_total)), 's']);
