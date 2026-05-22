function Df = quadDroneEuler_jac(x,I,M, niu_v, niu_w)
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
omegaB1 = x(10);
omegaB2 = x(11);
omegaB3 = x(12);

F = x(13);
w1 = x(14);
w2 = x(15);
w3 = x(16);

I1 = I(1);
I2 = I(2);
I3 = I(3);

Df = [(0), (0), (0), cos(s2)*cos(s3), cos(s3)*sin(s1)*sin(s2) - cos(s1)*sin(s3), sin(s1)*sin(s3) + cos(s1)*cos(s3)*sin(s2), v2*(sin(s1)*sin(s3) + cos(s1)*cos(s3)*sin(s2)) + v3*(cos(s1)*sin(s3) - cos(s3)*sin(s1)*sin(s2)), v3*cos(s1)*cos(s2)*cos(s3) - v1*cos(s3)*sin(s2) + v2*cos(s2)*cos(s3)*sin(s1), v3*(cos(s3)*sin(s1) - cos(s1)*sin(s2)*sin(s3)) - v2*(cos(s1)*cos(s3) + sin(s1)*sin(s2)*sin(s3)) - v1*cos(s2)*sin(s3), (0), (0), (0); 
    (0), (0), (0), cos(s2)*sin(s3), (cos(s1)*cos(s3) + sin(s1)*sin(s2)*sin(s3)), cos(s1)*sin(s2)*sin(s3) - cos(s3)*sin(s1), - v2*(cos(s3)*sin(s1) - cos(s1)*sin(s2)*sin(s3)) - v3*(cos(s1)*cos(s3) + sin(s1)*sin(s2)*sin(s3)), v3*cos(s1)*cos(s2)*sin(s3) - v1*sin(s2)*sin(s3) + v2*cos(s2)*sin(s1)*sin(s3), v3*(sin(s1)*sin(s3) + cos(s1)*cos(s3)*sin(s2)) - v2*(cos(s1)*sin(s3) - cos(s3)*sin(s1)*sin(s2)) + v1*cos(s2)*cos(s3), (0), (0), (0); 
    (0), (0), (0), -sin(s2), cos(s2)*sin(s1), cos(s1)*cos(s2), v2*cos(s1)*cos(s2) - v3*cos(s2)*sin(s1), - v1*cos(s2) - v3*cos(s1)*sin(s2) - v2*sin(s1)*sin(s2), (0), (0), (0), (0); 
    (0), (0), (0), -niu_v, omegaB3, -omegaB2, (F*(cos(s1)*sin(s3) - cos(s3)*sin(s1)*sin(s2)))/M, (49*cos(conj(s2)))/5 + (F*cos(s1)*cos(s2)*cos(s3))/M, (F*(cos(s3)*sin(s1) - cos(s1)*sin(s2)*sin(s3)))/M, (0), -v3, v2; 
    (0), (0), (0), -omegaB3, -niu_v, omegaB1, - (49*cos(conj(s1))*cos(conj(s2)))/5 - (F*(cos(s1)*cos(s3) + sin(s1)*sin(s2)*sin(s3)))/M, (49*sin(conj(s1))*sin(conj(s2)))/5 + (F*cos(s1)*cos(s2)*sin(s3))/M, (F*(sin(s1)*sin(s3) + cos(s1)*cos(s3)*sin(s2)))/M, v3, (0), -v1; 
    (0), (0), (0), omegaB2, -omegaB1, -niu_v, (49*cos(conj(s2))*sin(conj(s1)))/5 - (F*cos(s2)*sin(s1))/M, (49*cos(conj(s1))*sin(conj(s2)))/5 - (F*cos(s1)*sin(s2))/M, (0), -v2, v1, (0); 
    (0), (0), (0), (0), (0), (0), omegaB2*cos(s1)*tan(s2) - omegaB3*sin(s1)*tan(s2), omegaB3*cos(s1)*(tan(s2)^2 + 1) + omegaB2*sin(s1)*(tan(s2)^2 + 1), (0), (1), sin(s1)*tan(s2), cos(s1)*tan(s2); 
    (0), (0), (0), (0), (0), (0), - omegaB3*cos(s1) - omegaB2*sin(s1), (0), (0), (0), cos(s1), -sin(s1); 
    (0), (0), (0), (0), (0), (0), (omegaB2*cos(s1))/cos(s2) - (omegaB3*sin(s1))/cos(s2), (omegaB3*cos(s1)*sin(s2))/cos(s2)^2 + (omegaB2*sin(s1)*sin(s2))/cos(s2)^2, (0), (0), sin(s1)/cos(s2), cos(s1)/cos(s2); 
    (0), (0), (0), (0), (0), (0), (0), (0), (0), -niu_w, (I2*omegaB3 - I3*omegaB3)/I1, (I2*omegaB2 - I3*omegaB2)/I1; 
    (0), (0), (0), (0), (0), (0), (0), (0), (0), -(I1*omegaB3 - I3*omegaB3)/I2, -niu_w, -(I1*omegaB1 - I3*omegaB1)/I2; 
    (0), (0), (0), (0), (0), (0), (0), (0), (0), (I1*omegaB2 - I2*omegaB2)/I3, (I1*omegaB1 - I2*omegaB1)/I3, -niu_w];

% Df=[(0), (0), (0), (1), (0), (0), (0), (0), (0), (0), (0), (0), (0), (0), (0), (0), (0), (0); 
%     (0), (0), (0), (0), (1), (0), (0), (0), (0), (0), (0), (0), (0), (0), (0), (0), (0), (0); 
%     (0), (0), (0), (0), (0), (1), (0), (0), (0), (0), (0), (0), (0), (0), (0), (0), (0), (0); 
%     (0), (0), (0), -(1/100), (0), (0), (0), (0), (0), (0), (0), (0), (640*F)/9, (0), (0), (0), (0), (0); 
%     (0), (0), (0), (0), -(1/100), (0), (0), (0), (0), (0), (0), (0), (0), (640*F)/9, (0), (0), (0), (0); 
%     (0), (0), (0), (0), (0), -(1/100), (0), (0), (0), (0), (0), (0), (0), (0), (640*F)/9, (0), (0), (0); 
%     (0), (0), (0), (0), (0), (0), (0), -omegaB3, omegaB2, (0), (0), (0), (0), (0), (0), (0), A(3,1), -A(2,1); 
%     (0), (0), (0), (0), (0), (0), omegaB3, (0), -omegaB1, (0), (0), (0), (0), (0), (0), -A(3,1), (0), A(1,1); 
%     (0), (0), (0), (0), (0), (0), -omegaB2, omegaB1, (0), (0), (0), (0), (0), (0), (0), A(2,1), -A(1,1), (0); 
%     (0), (0), (0), (0), (0), (0), (0), (0), (0), (0), -omegaB3, omegaB2, (0), (0), (0), (0), A(3,2), -A(2,2); 
%     (0), (0), (0), (0), (0), (0), (0), (0), (0), omegaB3, (0), -omegaB1, (0), (0), (0), -A(3,2), (0), A(1,2); 
%     (0), (0), (0), (0), (0), (0), (0), (0), (0), -omegaB2, omegaB1, (0), (0), (0), (0), A(2,2), -A(1,2), (0); 
%     (0), (0), (0), (0), (0), (0), (0), (0), (0), (0), (0), (0), (0), -omegaB3, omegaB2, (0), A(3,3), -A(2,3); 
%     (0), (0), (0), (0), (0), (0), (0), (0), (0), (0), (0), (0), omegaB3, (0), -omegaB1, -A(3,3), (0), A(1,3); 
%     (0), (0), (0), (0), (0), (0), (0), (0), (0), (0), (0), (0), -omegaB2, omegaB1, (0), A(2,3), -A(1,3), (0); 
%     (0), (0), (0), (0), (0), (0), (0), (0), (0), (0), (0), (0), (0), (0), (0), -(7/10), (719*omegaB3)/728, (719*omegaB2)/728; 
%     (0), (0), (0), (0), (0), (0), (0), (0), (0), (0), (0), (0), (0), (0), (0), -(719*omegaB3)/728, -(7/10), -(719*omegaB1)/728; 
%     (0), (0), (0), (0), (0), (0), (0), (0), (0), (0), (0), (0), (0), (0), (0), (0), (0), -(7/10)];
end