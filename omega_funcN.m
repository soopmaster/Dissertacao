function [omega] = omega_funcN(x,N)
%UNTITLED5 get predictors of trajectory assuming 2N constraints
%   Detailed explanation goes here
if size(x,1)<=1
    x=x';
end
omega = zeros(2*N);
prod_vec = ones(1,2*N);
for k=1:N
    if k>1
        prod_vec = prod_vec.*pow_vec;
    end
    pow_vec = [zeros(1,k), 1:2*N-k];
    omega(2*k-1:2*k,:) = prod_vec.*(x.^(pow_vec));
end

end