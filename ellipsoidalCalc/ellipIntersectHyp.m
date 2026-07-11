function [w,W] = ellipIntersectHyp(q,Q,c,b)
%ellipIntersectHyp return ellipsoid resulting from the intersect of ellipsoid with hyperplane
%   INPUT:
% q : ellipsoid center
% Q : ellipsoid shape matrix
% c : hyperplane normal vector
% b : hyperplane offset
%   OUTPUT:
% w : resulting ellipsoid center
% W : resulting ellipsoid shape matrix

% E(q,Q) = {x in R^n : (x-q)'*Q-1*(x-q)<=1}
% H(c,b) = {x in R^n :  c'*x=b}
% E_H(w,W) = E(q,Q) intersect H(n,b) != 0

n = length(q);

% 1) Check if dist(E,H) < 0
% We won't be doing this... Assume there is guaranted intersection

% 2) Transform H into H([1 0 ... 0]',0)
v1 = c;
v2 = zeros(n,1);
v2(1) = 1;
T = orthoAlignmentMx(v1,v2);

% 3) Transform ellips into E(qt,Qt)
qt = q - b/((c'*c)^(1/2))*T*c;
Qt = T*Q*T';

% 4) Define M
M = inv(Qt);

m11 = M(1,1);
mbar = M(2:end,1);
Mbar = M(2:end,2:end);

% 5) Define E(wt,Wt)
qt1 = qt(1);
wt = qt + qt1*[-1; Mbar\mbar];
Maux = [zeros(1,size(Mbar,1)+1); [zeros(size(Mbar,1),1),inv(Mbar)]];
Wt = (1-qt1^2*(m11-mbar'*(Mbar\mbar)))*Maux;

% 6) Inverse Transform to E_H(w,W)
w = T'*wt + b/((c'*c)^(1/2))*c;
W = T'*Wt*T;

% If W is singular, reduce dimensions to not singular
row_all_zeros = 1;
while row_all_zeros <= size(W,1) && ~any(W(row_all_zeros,:)).is_zero()
    row_all_zeros = row_all_zeros+1;
end  
if row_all_zeros <= size(W,1)
    w = w(1:row_all_zeros-1);
    W = W(1:row_all_zeros-1,1:row_all_zeros-1);
end

end