% Simple Rocket Simulation

close all;
clearvars;

Rot = @(x)[cos(x) -sin(x);
           sin(x) cos(x)]; % Rotation Matrix

Ts = 0.01;
Tsim = 5; % Sim time seconds
knum = Tsim/Ts; % Number of samples

% u1: thruster angle actuation (between -1 and 1)
% u3: thruster power actuation (between 0 and 1)

u1_lock = pi/4; % angle max [rad]
u3_lock = 100; % power max actuation [m^2/s]

u1 = zeros(1,knum);
u3 = zeros(1,knum);

u1(1:30) = 1;
u1(30:35) = -1;
%u1(50:55)=1;
%u1(60:65)=-1;
%u1(70:75)=-1;
u3(1:300) = 1;
u = [u1_lock; u3_lock].*[u1; u3];

udim = size(u,1);

% ap_spi: acceleration (a) given by thruster (p) along spine (spi) [m^2/s]
% ap_lat: acceleration (a) given by thruster (p) along lateral (lat) [m^2/s]
% omega1_dot: angular acceleration of given by thruster to rocket attitude [rad^2/s]

L = 1; % distance from thruster to center of gravity
m = 5; % mass of rocket
I = 10; % angular inertia
g = 9.8;

ap_spi = u(2,:).*cos(u(1,:));
ap_lat = -u(2,:).*sin(u(1,:));
ap = [ap_lat; ap_spi];
omega1_dot = L*m/I*u(2,:).*sin(u(1,:));

% atx: total acceleration along x axis
% atz: total acceleration along z axis
% theta: rocket attitude [rad/s]
% x = [x1 x3 vx1 vx3 theta omega]'
% u_cm = [atx atz omega1_dot] (actuation on center of mass)

u_cm = zeros(3,knum);
u_cmdim = size(u_cm,1);
x = zeros(6,knum);
xdim = size(x,1);

x(2,1) = L; % center mass initial value

% x_dot = F*x + B*u + w
% z = H*x + v

z = zeros(xdim,knum);

p = zeros(2,knum); % thruster positions

l = zeros(1,knum); % L size

Fcnt = [0 0 1 0 0 0;
        0 0 0 1 0 0;
        0 0 0 0 0 0;
        0 0 0 0 0 0;
        0 0 0 0 0 1;
        0 0 0 0 0 0];
Bcnt = [0 0 0;
        0 0 0;
        1 0 0;
        0 1 0;
        0 0 0;
        0 0 1];
Hcnt = eye(xdim);
Dcnt = zeros(xdim,u_cmdim);

sysCnt = ss(Fcnt,Bcnt,Hcnt,Dcnt);

sysk = c2d(sysCnt,Ts);

Fk = sysk.A; Bk = sysk.B; Hk = sysk.C; Dk = sysk.D;

figure;
xlabel('$x$','Interpreter','latex');
ylabel('$z$','Interpreter','latex');
%legend('Center Mass','Thruster','Interpreter','latex');
axis([-10 10 0 30]);
grid on;
hold on;
for k=1:knum-1
    if x(5,k) > pi
        x(5,k) = -pi;
    elseif x(5,k) < -pi
        x(5,k) = pi;
    end
    theta = x(5,k);

    at = [-cos(theta) sin(theta); sin(theta) cos(theta)]*ap(:,k) + [0; -g];

    u_cm(:,k) = [at; omega1_dot(k)];

    x(:,k+1) = Fk*x(:,k) + Bk*u_cm(:,k);
    z(:,k+1) = Hk*x(:,k);

    p(:,k) = z(1:2,k) - L.*[cos(pi/2-theta); sin(pi/2-theta)]; 
    l(k) = norm(p(:,k)-z(1:2,k));
    plot([p(1,k) z(1,k)], [p(2,k) z(2,k)], 'r');
end

plot(z(1,1:knum-1),z(2,1:knum-1),'og',p(1,1:knum-1),p(2,1:knum-1),'xr');
hold off;
