function Df = robot_jac(x)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
q1=x(1);
q2=x(2);
z1=x(3);
z2=x(4);
Df=[1, 0, 0, 0; 
 0, 1, 0, 0;
 -(2*q2*z2 + 2)/(q2^2 + 1), -(2*q2*z1)/(q2^2 + 1), -2/(q2^2 + 1), (2*q2*(2*q1 + 2*z1 + 2*q2*z1*z2 - 4))/(q2^2 + 1)^2 - (2*z1*z2)/(q2^2 + 1); 
 2*q2*z1, -1, 0, z1^2 - 1];
end