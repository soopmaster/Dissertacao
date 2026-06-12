function At = proj3(A)
%proj18 project NxN ellipsoid matrix or Nx1 vector A into R3
%   Detailed explanation goes here

nx = size(A);

T = [eye(3), zeros(3,nx(1)-3)];

Tinv = pinv(T);
%Tinv = conj(T)/(T*conj(T));

%At = inv((T*sqrt(inv(A))*Tinv)^2);
if nx(2) == 1
    At = T*A;
else
    At = T*A*T';
end
end