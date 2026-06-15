function [omega] = omega_funcN(x,N)
%UNTITLED5 get predictors of trajectory assuming size(x,2)*N constraints
%   Detailed explanation goes here
if size(x,1)<=1
    x=x';
end
num = size(x,1);
omega = zeros(num*N);
prod_vec = ones(1,num*N);
for k=1:N
    if k>1
        prod_vec = prod_vec.*pow_vec;
    end
    pow_vec = [zeros(1,k), 1:num*N-k];
    omega(num*(k-1)+1:num*k,:) = prod_vec.*(x.^(pow_vec));
end

end