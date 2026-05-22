function J = bruss_jac(y)
J = [-2.5+2*y(1)*y(2),y(1)^2; 1.5-2*y(1)*y(2), -y(1)^2];