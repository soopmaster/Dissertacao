function dxdt = dIntegrator3DODE(x,u)
    % Variables for ODE: x = [p; v; ang]
    % p: Position of center mass, relative to inertial
    % v: Velocity of center mass, relative to inertial
    % ang: Angular position, relative to body

    % Input variables:
    % x: state vector
    % u = [a; tau] (input vector, a = acceleration, tauB = input
    %   angular velocity vector)
    
    g = 9.81;

    p = x(1:3);
    v = x(4:6);
    ang = x(7:end);

    a = u(1);
    tau = u(2:end); % Torque vector in relation to body

    aB = a.*[0;0;1]; %  Linear acceleration in relation to body

    p_dot = v;
    vdot = rotM(ang)*aB+[0;0;-g];
    angdot = tau;

    dxdt = [p_dot; vdot; angdot];

end