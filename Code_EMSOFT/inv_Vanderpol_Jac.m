function  J = inv_Vanderpol_Jac(x)

J = [0 -1; 1+2*x(1)*x(2) x(1)^2-1];