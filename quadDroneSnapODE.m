function dxdt = quadDroneSnapODE(x,u,I,m,niu_v,niu_w)
    % Extended quadDroneODE with F and dF to accomodate for Snap controller
    % Variables for ODE: x = [p_cm; v; omegaB; A_vec; F; dF]
    % p_cm: Position of center mass, relative to inertial
    % v: Velocity of center mass, relative to inertial
    % omegaB: Angular velocity, relative to body
    % A_vec: DCM matrix of attitude in vector form 
    % F: Thrust force
    % dF: Derivative of thrust force

    % Input variables:
    % x: state vector
    % u = [ddF; tauB] (input vector, ddF = second derivative of total thrust force value, tauB = input
    %   torques vector)
    % I: Inertial momentum matrix diagonal
    % m: mass of drone
    % p_bB: Point of application of lateral jets (u)
    
    g = 9.81;

    p_cm = x(1:3);
    v = x(4:6);
    A_vec = x(7:15);
    omegaB = x(16:18);
    F = x(19);
    dF = x(20);

    ddF = u(1);
    tauB = u(2:end); % Torque vector in relation to body

    A = reshape(A_vec,3,[]);
    %det(A)
    J = diag(I);
    omegaB_cr = crossprodmat(omegaB);
 
    aB = 1/m*F.*[0;0;1]; %  Linear acceleration in relation to body

    p_cmdot = v;
    vdot = A*aB+[0;0;-g]-niu_v.*v;
    Adot = omegaB_cr*A;
    omegaBdot = J\(tauB-omegaB_cr*(J*omegaB))-niu_w.*omegaB;
    Fdot = dF;
    Fddot = ddF;

    Adot_vec = reshape(Adot,[],1);

    dxdt = [p_cmdot; vdot; Adot_vec; omegaBdot; Fdot; Fddot];

end