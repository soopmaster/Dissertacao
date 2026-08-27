function Q = ellipFromAlph(alph,y)
%ellipFromAlph return Shape matrix from alph and y
%   INPUT:
% alph : ellipsoid center
% y : ellipsoid shape matrix

% E(q,Q) = {x in R^n : (x-q)'*Q-1*(x-q)<=1}
% P(C,B) = {x in R^n : C'*x=B}
% E_H(w,W) = E(q,Q) intersect P(C,B) != 0

Q = y.*inv(alph'*alph);
end