function vec = m2v(M)
%M2V Vectorises Matrix M
%   Vectorises nxm matrix into mnx1 vector by stacking collumns
s = size(M);
vec = reshape(M,s(1)*s(2),1);
end