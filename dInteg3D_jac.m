function Df = dInteg3D_jac(x)
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
s3 = x(9);

u1 = x(10);
u2 = x(11);
u3 = x(12);
u4 = x(13);

Df = [0, 0, 0, 1, 0, 0,                                               0,                          0,                                              0;
0, 0, 0, 0, 1, 0,                                               0,                          0,                                              0;
0, 0, 0, 0, 0, 1,                                               0,                          0,                                              0;
0, 0, 0, 0, 0, 0,  u1*(cos(s1)*sin(s3) - cos(s3)*sin(s1)*sin(s2)), u1*cos(s1)*cos(s2)*cos(s3), u1*(cos(s3)*sin(s1) - cos(s1)*sin(s2)*sin(s3));
0, 0, 0, 0, 0, 0, -u1*(cos(s1)*cos(s3) + sin(s1)*sin(s2)*sin(s3)), u1*cos(s1)*cos(s2)*sin(s3), u1*(sin(s1)*sin(s3) + cos(s1)*cos(s3)*sin(s2));
0, 0, 0, 0, 0, 0,                             -u1*cos(s2)*sin(s1),        -u1*cos(s1)*sin(s2),                                              0;
0, 0, 0, 0, 0, 0,                                               0,                          0,                                              0;
0, 0, 0, 0, 0, 0,                                               0,                          0,                                              0;
0, 0, 0, 0, 0, 0,                                               0,                          0,                                              0];

end