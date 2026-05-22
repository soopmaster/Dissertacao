function dxdt = prmtpE(x,c)
    % Variables for ODE: x = y(k)
    %   E(xo,alph,y) = c*y
    % where
    %   c = max_i(max(eig((alph*(Mi-Df(xo))*alph^-1)'+ alph*(Mi-Df(xo))*alph^-1))
    % and
    %   Mi \in M (set of interval Mixed Jacobian Matrices)

    % Input variables:
    %   c : offset parameter

    dxdt = c*x;

end