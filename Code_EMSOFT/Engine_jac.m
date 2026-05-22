function J = Engine_jac(y)
%J = zeros(2,2);
J = [-3*y(1)-3.0*(y(1)^2)/2, -1;3,-1];
