function max_min = LSE(alpha,X)
%LSE LogSumExp computes a smooth approximation of maximum values in X,
% The greater alpha is, the tighter the approximation
% If alpha<0, LSE calculates the minimum
% alpha is coefficient (usually alpha == 4 is more than enough for good
% accuracy)
% X is vector of values

max_min = 1/alpha.*log(sum(exp(alpha.*X)));
end