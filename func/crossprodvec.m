function vec = crossprodvec(M)
%CROSSPRODVEC map skew-symmetric matrix M to vector
%   [ 0  a  b ;     [-c;
%    -a  0  c ; ->    b;
%    -b  -c 0 ]      -a]

a = M(1,2);
b = M(1,3);
c = M(2,3);
vec = [-c;b;-a];
end