% Simple Rocket Simulation

close all;
clearvars;

% == Euler Angles == %
Rx = @(x)[1 0 0;
          0 cos(x) -sin(x);
          0 sin(x) cos(x)]; % Rotation Matrix
Ry = @(x)[cos(x) 0 sin(x);
          0 1 0;
          -sin(x) 0 cos(x)]; % Rotation Matrix
Rz = @(x)[cos(x) -sin(x) 0;
          sin(x) cos(x) 0;
          0 0 1]; % Rotation Matrix

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


L = 1; % arm length
R = 0.25;
m = 0.1; % area mass of drone
I = [-(2*R^2*m*(R^3/8 - (L + R/2)^3))/3; -(2*R^2*m*(R^3/8 - (L + R/2)^3))/3; (R^4*m*(2*L + R))/12]; % moment of inertia
M = m*R^2*(2*L+R); % total mass of drone
J = diag(I); % Inertia Matrix
g = 9.8;

% == Variables == %

% Variables for ODE: x = [p_vec; v; omegaB; A_vec]
% p_vec: Position of all points, relative to inertial
% v: Velocity of center mass, relative to inertial
% s: euler angles of drone (roll, pitch, yaw)
% omegaB: Angular velocity, relative to body

% y = [px; py; pz; omegaB3]

nx = 18;
nu = 4;
ny = 6;

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
%A0 = [0 0 -1; 0 1 0; 1 0 0];
A0 = eye(3);
A = A0;
A0_vec = reshape(A0,[],1);

x0 = [p_cm0; v0; s0; omega0B];

p_NB = [L+R/2;0;0];
p_SB = [-L-R/2;0;0];
p_EB = [0;L+R/2;0];
p_OB = [0;-L-R/2;0];

p_vecB = [p_NB; p_EB; p_SB; p_OB];

ref = [5;-5;20];

% == Control Action == %
% u = [FB; tauB]

FB = 1; % total thrust force [N]
tauB = [0.1;0.1;0]; % input torque [N/m]

Fmax = 2;
tauBmax = 5;

u_val = [FB; tauB];

% === Sim Settings === %

simT = 5; % Sim time seconds
Ts = 0.01;
knum = simT/Ts;

x = zeros(size(x0,1),knum);
x(:,1) = x0;

u = u_val.*ones(nu,knum);

p_cm = x(1:3,:);
%A_vecs = x(10:end,:);

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

KpF = 0.2;
KdF = 0.1;
KpA = 0.1;
KdA = 0.5;

%===v NMPC BATCH APPROACH v===%

%Q = kron(eye(2),diag([10;0.1]));
Q = diag([10;10;10; 10]);
R = 1*eye(nu);
P = 1*Q;
N = 5;
Qbar = blkdiag(kron(Q,eye(N)),P);
Rbar = kron(R,eye(N));

CUmat = blkdiag(Cmat, Cmat, Cmat, Cmat, Cmat, Cmat, eye(nu*N));


Fmpc = @(X,U,X0) [X(1:nx)-X0;               % x0 = init
            quadDroneODE(X(1:nx),U(1:nu),I,M)-X(1*nx+(1:nx));   % f(x0,u0) = x1
            quadDroneODE(X(1*nx+(1:nx)),U(1*nu+(1:nu)),I,M)-X(2*nx+(1:nx));   % f(x1,u1) = x2
            quadDroneODE(X(2*nx+(1:nx)),U(2*nu+(1:nu)),I,M)-X(3*nx+(1:nx));   % f(x2,u2) = x3
            quadDroneODE(X(3*nx+(1:nx)),U(3*nu+(1:nu)),I,M)-X(4*nx+(1:nx));   % f(x3,u3) = x4
            quadDroneODE(X(4*nx+(1:nx)),U(4*nu+(1:nu)),I,M)-X(5*nx+(1:nx));];  % f(x4,u4) = x5

REF = [kron(ones(N+1,1),ref);zeros(N*nu,1)];

cost = @(XU) (CUmat*XU-REF)'*blkdiag(Qbar,Rbar)*(CUmat*XU-REF);
Ceq = @(XU,X0) Fmpc(XU(1:(N+1)*nx,1),XU((N+1)*nx+(1:(N*nu)),1),X0);

X0 = zeros((N+1)*nx+N*nu,1);

function [c,ceq] = stateconstraints(XU,Ceq,X0)
    c = []; 
    ceq = Ceq(XU,X0);
end

% === Start Figure === %
figure;
xlabel('$x$','Interpreter','latex');
ylabel('$y$','Interpreter','latex');
zlabel('$z$','Interpreter','latex');
plot3(NaN, NaN, NaN, 'DisplayName','Rocket Body','Color','b'); % Rocket Body
plot3(NaN, NaN, NaN, 'DisplayName','Thruster Plume','Color','m'); % Thruster Plume
axis([-10 10 -10 10 -10 10]);
grid on;
hold on;
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
for k=1:knum-1 
    % === Calculate Control action
    % 
    % % Compute the control action with NMPC (batch notation)
    % nonlcon = @(XU) stateconstraints(XU,Ceq,x(:,k));
    % 
    % XU0 = [x(:,k); zeros(N*nx,1);u(:,k); zeros((N-1)*nu,1)];
    % X = fmincon(cost,XU0,[],[],[],[],[],[],nonlcon);
    % 
    % u(:,k) = X((N+1)*nx+(1:nu),1);
    % u(:,k)

    % % Compute the control action with PD

    % Cxk = Cmat*x(:,k);
    % x(1:9,k)
    % uPD = PDfunc(Cxk,ref,Kp,Kd)
    % u(:,k) = uPD;

    % n_vec = A*[0;0;1] % Normal vector
    % 
    % e_pos = ref-x(1:3,k) % positional error
    % de_pos = -x(4:6,k)
    % 
    % e_att = e_pos./norm(e_pos)-n_vec % attitude error
    % de_att = -x(7:9,k)
    % 
    % F = KpF*e_pos(3) + KdF*de_pos(3);
    % tauB = KpA*norm(e_att)*crossprodmat([0;0;e_pos(3)])*[e_pos(1:2);0]+KdA.*de_att;
    % 
    % tauBnorm = norm(tauB);
    % 
    % % Clamp u to bounds
    % 
    % if F > Fmax
    %     F = Fmax
    % elseif F < -Fmax
    %     F = -Fmax
    % end
    % 
    % if tauBnorm > tauBmax
    %     scale = tauBmax/tauBnorm;
    %     tauB = scale.*tauB
    % end

    u(1,k) = FB;
    u(2:end,k) = tauB;
    
    % === Calculate System ===

    % [t,y] = ode45(@(t,y) quadDroneODE(y,u(:,k),I,M),[0 Ts],x(:,k));
    % x(:,k+1) = y(end,:)';
    x(:,k+1) = runge_kutta4(@(y) quadDroneEulerODE(y,u(:,k),I,M),Ts,x(:,k));

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
    uvw_dthrust = A*[0;0;u(1,k)];

    Xvec = [p_vec(1:3:12)];
    Yvec = [p_vec(2:3:12)];
    Zvec = [p_vec(3:3:12)];
    bodyPlot.XData = Xvec;
    bodyPlot.YData = Yvec;
    bodyPlot.ZData = Zvec;
    set(thrustPlot,'XData',p_cm(1),'YData',p_cm(2),'ZData',p_cm(3),'UData',uvw_thrust(1),'VData',uvw_thrust(2),'WData',uvw_thrust(3));
    set(dthrustPlot,'XData',p_cm(1),'YData',p_cm(2),'ZData',p_cm(3),'UData',uvw_dthrust(1),'VData',uvw_dthrust(2),'WData',uvw_dthrust(3));
    addpoints(pathPlot,p_cm(1), p_cm(2), p_cm(3));
    drawnow;
    pause(Ts);
end

%plot3(x(1,1:knum-1),x(2,1:knum-1), x(3,1:knum-1),'-og',p(1,1:knum-1),p(2,1:knum-1), p(3,1:knum-1), 'xr',tip(1,1:knum-1), tip(2,1:knum-1), tip(3,1:knum-1),'^k');
%legend('Location','best','Interpreter','latex');
hold off;
