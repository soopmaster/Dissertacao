function dx = quadDroneEulerODE(x,u,I,m,niu_v,niu_w)
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
    omega = x(10:end);

    % Clamp s
    % s(1) = mod(s(1) + pi/2, pi) - pi/2;
    % s(2) = mod(s(1) + pi/2, pi) - pi/2;
    % s(3) = wrapToPi(s(1));
    s(1) = wrapTo2Pi(s(1));
    s(2) = wrapTo2Pi(s(2));
    s(3) = wrapTo2Pi(s(3));

    F = u(1);
    tau = u(2:end);
 
    c_roll = cos(s(1));
    c_pitch = cos(s(2));
    c_yaw = cos(s(3));
    s_roll = sin(s(1));
    s_pitch = sin(s(2));
    s_yaw = sin(s(3));
    t_roll = tan(s(1));
    t_pitch = tan(s(2));
    t_yaw = tan(s(3));
    
    SOv = [c_pitch*c_yaw, s_roll*s_pitch*c_yaw-c_roll*s_yaw, c_roll*s_pitch*c_yaw+s_roll*s_yaw;
          c_pitch*s_yaw, s_roll*s_pitch*s_yaw+c_roll*c_yaw, c_roll*s_pitch*s_yaw-s_pitch*c_yaw;
          s_pitch, -s_roll*c_pitch, -c_roll*c_pitch]; % Rotation Matrix for v
    SOv = rotM(s);

    SOome = [1, s_roll*t_pitch, c_roll*t_pitch;
             0, c_roll, -s_roll;
             0, s_roll/c_pitch, c_roll/c_pitch]; % Rotation Matrix for omega

    g = [0;0;-9.8];

    J = diag(I);

    fx = [SOv*v;
          crossprodmat(v)*omega+SOv'*g-niu_v.*v;
          SOome*omega;
          -J\(crossprodmat(omega)*J*omega)-niu_w.*omega];

    gF = 1/m.*SOv*[0;0;1];

    gx = [zeros(3,4);
          gF, zeros(3,3);
          zeros(3,4);
          zeros(3,1), J^-1];

    % gx = [zeros(3,4);
    %       (-c_roll*s_pitch*c_yaw-s_roll*s_yaw)/m 0 0 0;
    %       (-c_roll*s_pitch*s_yaw+s_roll*s_yaw)/m 0 0 0;
    %       (-c_roll*c_pitch)/m 0 0 0;
    %       zeros(3,4);
    %       zeros(3,1), J^-1];

    dx = fx + gx*u;

end