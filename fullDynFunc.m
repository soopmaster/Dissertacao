function dfull = fullDynFunc(sysDynFunc, prmtpDynFunc, xvec, nx)
%FULLDYNFUNC Generalized function to compute dynamics of system +
%parametope
% xvec = [x;
%         alph]

x = xvec(1:nx); % Extract system variables

alph = v2m(xvec(nx+1:end),nx);

dx = sysDynFunc(x);

dalph = prmtpDynFunc(alph,x);

dfull = [dx; m2v(dalph)];

end