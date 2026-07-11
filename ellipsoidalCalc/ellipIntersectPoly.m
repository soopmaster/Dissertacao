function [w,W] = ellipIntersectPoly(q,Q,C,B)
%ellipIntersectHyp return ellipsoid resulting from the intersect of
%ellipsoid with polytope
%   INPUT:
% q : ellipsoid center
% Q : ellipsoid shape matrix
% C : polytope normal vectors (C = [c1, c2, ..., cm])
% B : polytope offsets (C = [b1, b2, ..., bm])
%   OUTPUT:
% w : resulting ellipsoid center
% W : resulting ellipsoid shape matrix

% E(q,Q) = {x in R^n : (x-q)'*Q-1*(x-q)<=1}
% P(C,B) = {x in R^n : C'*x=B}
% E_H(w,W) = E(q,Q) intersect P(C,B) != 0

Ct = C;
Bt = B;
w = q;
W = Q;
while ~isempty(Bt)
    [w,W] = ellipIntersectHyp(w,W,Ct(:,end),Bt(end));
    N = length(w);
    Ct = Ct(1:N,1:end-1);
    Bt = Bt(1:end-1);
end

end