function M = v2m(vec,s1)
%UNTITLED2 converts [Nx1] vectorised matrix to [s1xs2] matrix
%   Detailed explanation goes here
s2 = length(vec)/s1;
M = reshape(vec,s1,s2);
end