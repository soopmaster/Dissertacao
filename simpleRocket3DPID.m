% Simple Rocket Simulation

close all;
clearvars;

Rx = @(x)[1 0 0;
          0 cos(x) -sin(x);
          0 sin(x) cos(x)]; % Rotation Matrix
Ry = @(x)[cos(x) 0 sin(x);
          0 1 0;
          -sin(x) 0 cos(x)]; % Rotation Matrix
Rz = @(x)[cos(x) -sin(x) 0;
          sin(x) cos(x) 0;
          0 0 1]; % Rotation Matrix

Ts = 0.01;
Tsim = 5; % Sim time seconds
knum = Tsim/Ts; % Number of samples

% u1: thruster angle actuation around y axis (theta) [-1:1]
% u2: thruster angle actuation around z axis (psi) [-1:1]
% u3: thruster power actuation (between 0 and 1)

u1_lock = pi/4; % angle max [rad]
u2_lock = pi;
u3_lock = 100; % power max actuation [m^2/s]

u1 = zeros(1,knum);
u2 = zeros(1,knum);
u3 = zeros(1,knum);

%u1(30:40) = 1;
%u2(30:40) = -0.5;
%u1(30:35)=-1;
%u2(40:50)=-0.25;
% u1(70:75)=-1;
u3(1:knum) = 0.5;
%u3(30:35) = 1;
%u3(50:100) = 1;
u = [u1_lock; u2_lock; u3_lock].*[u1; u2; u3];

udim = size(u,1);

% ap = [apx; apy; apz]: acceleration (a) given by thruster (p) to rocket [m^2/s]
% omega_dotdot = [theta_dotdot; phi_dotdot; psi_dotdot]: angular acceleration of given by thruster to rocket attitude relative to referential[rad^2/s]

ap = zeros(3,knum);
omega_dotdot = zeros(3,knum);

L = 2; % distance from thruster to center of gravity
R = 0.25;
m = 10; % mass of rocket
I = [1/4*m*R^2+1/12*m*(2*L)^2; 1/4*m*R^2+1/12*m*(2*L)^2; 1/2*m*R^2]; % angular inertia
g = 9.8;

% ar = [arx; ary; arz]: thruster acceleration relative to referencial
% at = [atx; aty; atz]: total acceleration relative to referencial
% x = [x1 x2 x3 vx1 vx2 vx3 theta phi psi theta_dot phi_dot psi_dot]'
% u_cm = [at; omega_dot] (actuation on center of mass)

ar = zeros(3,knum);

u_cm = zeros(6,knum);
u_cmdim = size(u_cm,1);
x = zeros(12,knum);
xdim = size(x,1);

x(:,1) = [0 0 L 0 0 0 0 0 0 0 0 0]'; % center mass initial value
z(:,1) = [0 0 L 0 0 0 0 0 0 0 0 0]';

% x_dot = F*x + B*u + w
% z = H*x + v

z = zeros(xdim,knum);

p = zeros(3,knum); % thruster positions

tip = zeros(3,knum); % tip positions

% === State-Space System === %
Fcnt = kron([0 1 0 0; 0 0 0 0; 0 0 0 1; 0 0 0 0], eye(3));
Bcnt = kron([0 0; 1 0; 0 0; 0 1],eye(3));
Hcnt = eye(xdim);
Dcnt = zeros(xdim,u_cmdim);

sysCnt = ss(Fcnt,Bcnt,Hcnt,Dcnt);

sysk = c2d(sysCnt,Ts);

Fk = sysk.A; Bk = sysk.B; Hk = sysk.C; Dk = sysk.D;

% === Trajectory ===
r = trajectory3D([0 0 0; 2 0 15; 30 0 30]',[0.1 0 1; 0.2 0 1; 0.1 0 0.1]',[0 50 100],100);

% === Error ===
ref_theta = zeros(1,knum);
ref_phi = zeros(1,knum);
ref_psi = zeros(1,knum);
error = zeros(4,knum);

Kp  = [0.5;1];
Kd = [.05;.05];

% === Start Figure === %
figure;
xlabel('$x$','Interpreter','latex');
ylabel('$y$','Interpreter','latex');
zlabel('$z$','Interpreter','latex');
%plot3(NaN, NaN, NaN, 'DisplayName','Rocket Body','Color','b'); % Rocket Body
%plot3(NaN, NaN, NaN, 'DisplayName','Thruster Plume','Color','m'); % Thruster Plume
axis([-10 10 -10 10 0 10]);
grid on;
hold on;
plot3(r(1,:), r(2,:), r(3,:) ,'r');
uvw = [0;0;0];
bodyPlot = quiver3(p(1,1), p(1,2), p(1,3), uvw(1), uvw(2), uvw(3),2*L,'b');
thrustPlot = quiver3(p(1,1), p(1,2), p(1,3), 0,0,0,1,'m');

% === Simulation === %
for k=1:knum-1
    % Clamp Attitude Angles
    if x(7,k) > pi
        x(7,k) = -pi;
    elseif x(7,k) < -pi
        x(7,k) = pi;
    end
    if x(8,k) > pi
        x(8,k) = -pi;
    elseif x(8,k) < -pi
        x(8,k) = pi;
    end
    theta = x(7,k);
    phi = x(8,k);
    
    % === Calculate Control action
    tip(:,k) = x(1:3,k) + Rx(phi)*Ry(theta)*[0; 0; L];
    ind = dsearchn(r',tip(:,k)'); % Search nearest point of r to tip
    r_star = r(:,ind); 
    ref_vec = r_star-z(1:3,k);

    ref_theta(k) = asin(ref_vec(1)/norm([ref_vec(1); ref_vec(3)]));
    ref_phi(k) = asin(ref_vec(2)/norm(ref_vec));
    ref_psi(k) = atan2(ref_vec(2),ref_vec(1));

    if k>1
        ref_thetadot = (ref_theta(k)-ref_theta(k-1))/Ts;
        ref_phidot = (ref_phi(k)-ref_phi(k-1))/Ts;
        ref_psidot = (ref_psi(k)-ref_psi(k-1))/Ts;
        error_thetadot = ref_thetadot-z(10,k);
        error_psidot = ref_psidot-z(12,k);
    else
        ref_thetadot = 0;
        ref_phidot= 0;
        ref_psidot = 0;
        error_thetadot = 0;
        error_psidot = 0;
    end

    error_theta = ref_theta(k)-z(7,k);
    error_psi = ref_psi(k)-z(9,k)-pi;

    error(:,k) = 1/pi*[error_theta; error_psi; error_thetadot; error_psidot];
    %error(1,k) = clamp(error(1,k),-1,1);
    %error(2,k) = clamp(error(2,k),-1,1);

    uP = Kp.*error(1:2,k);
    uD = Kd.*error(3:4,k);

    u_ang = uP + uD;
    u_ang(1) = clamp(u_ang(1),-1,1);
    u_ang(2) = clamp(u_ang(2),-1,1);

    u(:,k) = [u1_lock; u2_lock; u3_lock].*[u_ang; u3(k)];

    % === Calculate System ===

    ap(:,k) = Rz(u(2,k))*Ry(u(1,k))*[0; 0; u(3,k)]; % Calculate thrust vector relative to rocket

    omega_dotdot(:,k) = L*m./I.*[-1;1;0].*ap(:,k);  % Calculate Angular acceleration induced by thrust

    ar(:,k) = Rx(phi)*Ry(theta)*ap(:,k); % Calculate thrust relative to ground

    at = ar(:,k) + [0; 0; -g]; % Calculate total thrust

    u_cm(:,k) = [at; omega_dotdot(:,k)];

    x(:,k+1) = Fk*x(:,k) + Bk*u_cm(:,k);
    z(:,k+1) = Hk*x(:,k);

    p(:,k) = x(1:3,k) + Rx(phi)*Ry(theta)*[0; 0; -L]; 
    p_vec = (1/u3_lock).*ar(:,k);

    uvw = Rx(phi)*Ry(theta-pi/2);
    set(bodyPlot,'XData',p(1,k),'YData',p(2,k),'ZData',p(3,k),'UData',uvw(1),'VData',uvw(2),'WData',uvw(3));
    set(thrustPlot,'XData',p(1,k),'YData',p(2,k),'ZData',p(3,k),'UData',p_vec(1),'VData',p_vec(2),'WData',p_vec(3));
    drawnow;
    pause(Ts);
end

plot3(x(1,1:knum-1),x(2,1:knum-1), x(3,1:knum-1),'-og',p(1,1:knum-1),p(2,1:knum-1), p(3,1:knum-1), 'xr',tip(1,1:knum-1), tip(2,1:knum-1), tip(3,1:knum-1),'^k');
%legend('Location','best','Interpreter','latex');
hold off;
