function dxdt = robotArmODE(t,x)
    % Variables for ODE:

    % Input variables:
    
    u1 = 2;
    u2 = 1;
    m = 1;
    M = 1;
    L = sqrt(3);
    kp1 = 2;
    kp2 = 1;
    kd1 = 2;
    kd2 = 1;
    
    dxdt = [x(3);
         x(4); 
         1/(m*x(2)^2+M*L^2/3)*(-2*m*x(2)*x(3)*x(4)-kd1*x(3)+kp1*(u1-x(1)));
         x(2)*x(3)^2+1/m*(-kd2*x(4)+kp2*(u2-x(2)))];

end