function T = orthoAlignmentMx(v1,v2)
%orthoAlignmentMx Transformation matrix T that is orthonormal basis of
%vectors v1 and v2
%   Detailed explanation goes here
n = length(v1);

v1n = (1/norm(v1)).*v1;
v2n = (1/norm(v2)).*v2;

d = v1n'*v2n;
s = (1-d^2)^(1/2);

q1 = v1n;
q2 = (1/s).*(v2n-d*v1n);
Q1 = [q1,q2];

S = [d s; -s d];

T = eye(n) + Q1*(S-eye(2))*Q1';

end