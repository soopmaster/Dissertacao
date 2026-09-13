function V = ellipVol(Q)
%ellipVolcalculate volume of ellipsoid
%   INPUT:
% Q : ellipsoid shape matrix
%   OUTPUT:
% V : scalar volume

% Volume of unit ball
n = size(Q,1);

if ~mod(n,2)
    Svol = pi^(n/2)/factorial(n/2);
else
    Svol = (2^n*pi^((n-1)/2)*factorial((n-1)/2))/factorial(n);
end

V = Svol*sqrt(det(Q));
end