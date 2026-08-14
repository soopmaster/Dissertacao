function [w,W] = ellip2R3(q,Q)
%ellip2R3 return ellipsoid resulting from the intersect of
%ellipsoid with R3 plane
%   INPUT:
% q : ellipsoid center
% Q : ellipsoid shape matrix

% E(q,Q) = {x in R^n : (x-q)'*Q-1*(x-q)<=1}
% P(C,B) = {x in R^n : C'*x=B}
% E_H(w,W) = E(q,Q) intersect P(C,B) != 0

nx = size(Q,1);

C = eye(nx);
C = C(:,4:end);
B = zeros(size(C,2),1);

[~,W] = ellipIntersectPoly(zeros(nx,1),Q,C,B);
w = q(1:3);

end