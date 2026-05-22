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
R = @(x)[cos(x(2))*cos(x(3)) cos(x(1))*sin(x(3))+sin(x(1))*sin(x(2))*cos(x(3)) sin(x(1))*sin(x(3))-cos(x(1))*sin(x(2))*cos(x(3));
        -cos(x(2))*sin(x(3)) cos(x(1))*cos(x(3))-sin(x(1))*sin(x(2))*sin(x(3)) cos(x(1))*sin(x(2))*sin(x(3))+sin(x(1))*cos(x(3));
        sin(x(2)) -sin(x(1))*cos(x(2)) cos(x(1))*cos(x(2))]; % Rotation Matrix (1-2-3) or (x-y-z)

% == Euler Axis-Angle == %
euler_ang = @(R) acos(1/2*(trace(R)-1)); % Euler Angle
euler_axis = @(R,Theta) 1/(2*sin(Theta))*(R'-R); % Euler Axis cross

% == DCM variation in time == %
% Adot = -crossprod(omegaB)*A

% == Parameters == %

L = 2; % distance from thruster to center of gravity
R = 0.25;
m = 10; % mass of rocket
I = [1/4*m*R^2+1/12*m*(2*L)^2; 1/4*m*R^2+1/12*m*(2*L)^2; 1/2*m*R^2]; % angular inertia
J = diag(I); % Inertia Matrix
g = 9.8;

% == Variables == %

% Variables for ODE: x = [p_vec; v; omegaB; A_vec]
% p_vec: Position of all points, relative to inertial
% v: Velocity of center mass, relative to inertial
% omegaB: Angular velocity, relative to body
% A_vec: DCM matrix of attitude in vector form 

% ux: thruster actuation x axis
% uy: thruster actuation y axis

% == Initial Conditions == %

p_bB = [0;0;-L]; % Position of bottom relative to body
p_cmB = [0;0;0]; % Position of centermass relative to body
p_tB = [0;0;L]; % Position of top relative to body
p_vecB = [p_bB; p_cmB; p_tB];

%A0 = [0 0 -1; 0 1 0; 1 0 0];
A0 = eye(3);
A0_vec = reshape(A0,[],1);

p_b0 = [0;0;0];
p_cm0 = [0;0;L];
p_t0 = [0;0;2*L];
v0 = [0;0;0];
omega0B = [0;0;0];

x0 = [p_cm0; v0; omega0B; A0_vec];

% == Control Action == %

u1_lock = pi/4; % angle max [rad]
u2_lock = pi;
u3_lock = 100; % power max actuation [m^2/s]

u_max = 10; % power max actuation [m^2/s]
aB_z = 100; % thruster acceleration [m^2/s]

%u = [u1_lock; u2_lock; u3_lock].*[1;0;1];
u_dir = [0;1];
u = u_max*u_dir/norm(u_dir);

% === Dynamics === %

simT = 5; % Sim time seconds
Ts = 0.01;
knum = simT/Ts;

x = zeros(size(x0,1),knum);
x(:,1) = x0;

p_cm = x(1:3,:);
p_b = zeros(3,knum);
p_t = zeros(3,knum);
A_vecs = x(10:end,:);

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

% === Start Figure === %
figure;
xlabel('$x$','Interpreter','latex');
ylabel('$y$','Interpreter','latex');
zlabel('$z$','Interpreter','latex');
plot3(NaN, NaN, NaN, 'DisplayName','Rocket Body','Color','b'); % Rocket Body
plot3(NaN, NaN, NaN, 'DisplayName','Thruster Plume','Color','m'); % Thruster Plume
axis([-10 10 -10 10 0 10]);
grid on;
hold on;
%plot3(r(1,:), r(2,:), r(3,:) ,'r');
uvw = [0;0;0];
bodyPlot = quiver3(p_cm(1,1), p_cm(1,2), p_cm(1,3), uvw(1), uvw(2), uvw(3),2*L,'b');
thrustPlot = quiver3(p_cm(1,1), p_cm(1,2), p_cm(1,3), 0,0,0,1,'m');
dthrustPlot = quiver3(p_cm(1,1), p_cm(1,2), p_cm(1,3), 0,0,0,1,'m');
pathPlot = animatedline(p_cm(1,1),p_cm(2,1), p_cm(3,1), 'Color','r','LineStyle','-');

% === Simulation === %
for k=1:knum-1
    % % Clamp Attitude Angles
    % if x(7,k) > pi
    %     x(7,k) = -pi;
    % elseif x(7,k) < -pi
    %     x(7,k) = pi;
    % end
    % if x(8,k) > pi
    %     x(8,k) = -pi;
    % elseif x(8,k) < -pi
    %     x(8,k) = pi;
    % end
    % theta = x(7,k);
    % phi = x(8,k);
    
    % % === Calculate Control action
    % tip(:,k) = x(1:3,k) + Rx(phi)*Ry(theta)*[0; 0; L];
    % ind = dsearchn(r',tip(:,k)'); % Search nearest point of r to tip
    % r_star = r(:,ind); 
    % ref_vec = r_star-z(1:3,k);
    % 
    % ref_theta(k) = asin(ref_vec(1)/norm([ref_vec(1); ref_vec(3)]));
    % ref_phi(k) = asin(ref_vec(2)/norm(ref_vec));
    % ref_psi(k) = atan2(ref_vec(2),ref_vec(1));
    % 
    % if k>1
    %     ref_thetadot = (ref_theta(k)-ref_theta(k-1))/Ts;
    %     ref_phidot = (ref_phi(k)-ref_phi(k-1))/Ts;
    %     ref_psidot = (ref_psi(k)-ref_psi(k-1))/Ts;
    %     error_thetadot = ref_thetadot-z(10,k);
    %     error_psidot = ref_psidot-z(12,k);
    % else
    %     ref_thetadot = 0;
    %     ref_phidot= 0;
    %     ref_psidot = 0;
    %     error_thetadot = 0;
    %     error_psidot = 0;
    % end
    % 
    % error_theta = ref_theta(k)-z(7,k);
    % error_psi = ref_psi(k)-z(9,k)-pi;
    % 
    % error(:,k) = 1/pi*[error_theta; error_psi; error_thetadot; error_psidot];
    % %error(1,k) = clamp(error(1,k),-1,1);
    % %error(2,k) = clamp(error(2,k),-1,1);
    % 
    % uP = Kp.*error(1:2,k);
    % uD = Kd.*error(3:4,k);
    % 
    % u_ang = uP + uD;
    % u_ang(1) = clamp(u_ang(1),-1,1);
    % u_ang(2) = clamp(u_ang(2),-1,1);
    % 
    % u(:,k) = [u1_lock; u2_lock; u3_lock].*[u_ang; u3(k)];

    % === Calculate System ===

    [t,y] = ode45(@(t,y) simpleRocketODE(t,y,[u;aB_z],I,m,p_tB),[0 Ts],x(:,k));
    x(:,k+1) = y(end,:)';

    p_cm = x(1:3,k);
    v = x(4:6,k);
    omegaB = x(7:9,k);
    A_vec = x(10:end,k);

    A = reshape(A_vec,3,[]);

    p_vec = kron(eye(3),A)*p_vecB + [p_cm; p_cm; p_cm];
    p_b = p_vec(1:3);
    p_t = p_vec(7:9);

    %aI = A*R1(u(2))*R3(u(1))*[u(3)/u3_lock;0;0] + p_b;

    uvw_body = A*[0;0;1]
    %uvw_thrust = A*R1(u(2))*R3(u(1))*[u(3)/u3_lock;0;0]
    uvw_thrust = A*[u/u_max;0];
    uvw_dthrust = A*[0;0;1];

    set(bodyPlot,'XData',p_b(1),'YData',p_b(2),'ZData',p_b(3),'UData',uvw_body(1),'VData',uvw_body(2),'WData',uvw_body(3));
    set(thrustPlot,'XData',p_t(1),'YData',p_t(2),'ZData',p_t(3),'UData',uvw_thrust(1),'VData',uvw_thrust(2),'WData',uvw_thrust(3));
    set(dthrustPlot,'XData',p_b(1),'YData',p_b(2),'ZData',p_b(3),'UData',uvw_dthrust(1),'VData',uvw_dthrust(2),'WData',uvw_dthrust(3));
    addpoints(pathPlot,p_cm(1), p_cm(2), p_cm(3));
    drawnow;
    pause(Ts*10);
end

%plot3(x(1,1:knum-1),x(2,1:knum-1), x(3,1:knum-1),'-og',p(1,1:knum-1),p(2,1:knum-1), p(3,1:knum-1), 'xr',tip(1,1:knum-1), tip(2,1:knum-1), tip(3,1:knum-1),'^k');
%legend('Location','best','Interpreter','latex');
hold off;
