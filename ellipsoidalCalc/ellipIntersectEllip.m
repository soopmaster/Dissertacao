function [qI,QI] = ellipIntersectEllip(q1,Q1,q2,Q2)
%ellipIntersectEllip return ellipsoid resulting from the intersect of
%ellipsoid with ellipsoid
%   INPUT:
% q1,2 : ellipsoid center
% Q1,2 : ellipsoid shape matrix
%   OUTPUT:
% qI : resulting ellipsoid center
% QI : resulting ellipsoid shape matrix

% E(q,Q) = {x in R^n : (x-q)'*Q-1*(x-q)<=1}
% E_I(qI,QI) = E(q1,Q1) intersect E(q1,Q1) != 0

outputArg1 = inputArg1;
outputArg2 = inputArg2;
end