function res = dpolyN(t,N,d) 
    % Compute d-th derivative of N-th order polynomial at point t
    res = [zeros(d,1); (factorial(d+(0:N-1-d))./factorial(0:N-1-d).*t.^(0:N-1-d))'];
end