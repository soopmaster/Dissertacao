function [omega] = omega_func(x)
%UNTITLED5 get predictors of trajectory assuming 2n constraints
%   Detailed explanation goes here
N=length(x);
omega = zeros(2*N);
for i=1:(2*N)
    for j=1:(2*N)
        if i<=N
            omega(i,j)=x(i)^(j-1);
        else
            if j>1
                omega(i,j)=(j-1)*x(i-N)^(j-2);
            else
                omega(i,j)=0;
            end
        end
    end
end
end