function dx = simpleDroneODE(x,u)
    % Variables for ODE: x = [p_cm; v; s; omega]
    % p_cm = [x,y,z]: Position of center mass, relative to inertial
    % v = [u, v, w]: Velocity of center mass, relative to inertial
    % s = [phi(roll), theta(pitch), psi(yaw)]: Rotation relative to inertial in Euler angle
    % omega = [p, q, r]: Angular velocity relative to body in Euler

    % Input variables:
    % x: state vector
    % u = [FB; tauB] (input vector, FB = total thrust force value, tauB = input
    %   torques vector)
    % I: Inertial momentum matrix diagonal
    % m: mass of drone
    % p_bB: Point of application of lateral jets (u)

    % dx = f(x) + g(x)*u
    
    p = x(1:3);
    v = x(4:6);
    s = x(7:9);
    
    % Clamp s
    % s(1) = mod(s(1) + pi/2, pi) - pi/2;
    % s(2) = mod(s(1) + pi/2, pi) - pi/2;
    % s(3) = wrapToPi(s(1));
    s(1) = wrapTo2Pi(s(1));
    s(2) = wrapTo2Pi(s(2));
    s(3) = wrapTo2Pi(s(3));

    m = 0.5; % mass of drone

    g = [0;0;-9.8];
    niu_v = 0.01;

    fx = [v;
          g-niu_v*v;
          zeros(3,1)];

    gx = [zeros(3,4);
          1/m.*rotM(s)*[0;0;1], zeros(3,3);
          zeros(3,1), eye(3)];

    dx = fx + gx*u;

end