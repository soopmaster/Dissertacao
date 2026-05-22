function Df = simpleDrone_jac(x)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
p1=x(1);
p2=x(2);
p3=x(3);
v1=x(4);
v2=x(5);
v3=x(6);
s1=x(7);
s2=x(8);
s3=x(9);
F = x(10);
w1 = x(11);
w2 = x(12);
w3 = x(13);

Df=[(0), (0), (0), (1), (0), (0), (0), (0), (0);
    (0), (0), (0), (0), (1), (0), (0), (0), (0); 
    (0), (0), (0), (0), (0), (1), (0), (0), (0); 
    (0), (0), (0), -(7/10), (0), (0), F*(2*cos(s1)*sin(s3) - 2*cos(s3)*sin(s1)*sin(s2)), 2*F*cos(s1)*cos(s2)*cos(s3), F*(2*cos(s3)*sin(s1) - 2*cos(s1)*sin(s2)*sin(s3)); 
    (0), (0), (0), (0), -(7/10), (0), -F*(2*cos(s1)*cos(s3) + 2*sin(s1)*sin(s2)*sin(s3)), 2*F*cos(s1)*cos(s2)*sin(s3), F*(2*sin(s1)*sin(s3) + 2*cos(s1)*cos(s3)*sin(s2)); 
    (0), (0), (0), (0), (0), -(7/10), -2*F*cos(s2)*sin(s1), -2*F*cos(s1)*sin(s2), (0); 
    (0), (0), (0), (0), (0), (0), (0), (0), (0); 
    (0), (0), (0), (0), (0), (0), (0), (0), (0); 
    (0), (0), (0), (0), (0), (0), (0), (0), (0)];
end