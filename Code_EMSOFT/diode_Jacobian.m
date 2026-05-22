function J = diode_Jacobian(x)

J = zeros(length(x),length(x));
C = 1;
L = 1;
R = 0.2;
V_in = 0.3;


dI_dV_d = 5*803.712*x(1)^4-4*1086.288*x(1)^3+3*551.088*x(1)^2-2*124.548*x(1)+10.656;

% dx(1) = (-I_d + x(2))/C;
% dx(2) = (-x(1) - R*x(2) + V_in)/L;

J = [-dI_dV_d/C, 1/C; -1/L, -R/L];