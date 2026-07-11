function D_KL = kullback_leibler_div(q1,Q1,q2,Q2)
%kullback_leibler_div computes KullBack-Leibler Divergence from
%distribution E1 to E2. May be thought as a general score of "distance"
%and is also applicable to ellipsoids.
%   INPUT:
% q1,2 : ellipsoid center
% Q1,2 : ellipsoid shape matrix
%   OUTPUT:
% qI : resulting ellipsoid center
% QI : resulting ellipsoid shape matrix

% E(q,Q) = {x in R^n : (x-q)'*Q-1*(x-q)<=1}
% D_KL(E1,E2) = 0 if E1 and E2 are identical

n = length(q1);
D_KL = 1/2*(trace(Q2\Q1) + (q2-q1)'*inv(Q2)*(q2-q1) - n + log(det(Q2)/det(Q1)));
end