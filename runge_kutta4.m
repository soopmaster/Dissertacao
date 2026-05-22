function yk = runge_kutta4(fun, h, y0)
    % Runge-Kutta 4 method

    k1 = fun(y0);
    k2 = fun(y0 + h/2*k1);
    k3 = fun(y0 + h/2*k2);
    k4 = fun(y0 + h*k3);

    yk = y0 + h/6*(k1+2*k2+2*k3+k4);
end