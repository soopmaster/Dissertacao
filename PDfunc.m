function u = PDfunc(xk,ref,Kp,Kd)
%UNTITLED Summary of this function goes here
%   xk = [p; v]
%   ref = reference value

dnum = size(xk,1)/2; % number of dimensions

p = xk(1:dnum);
v = xk(dnum+1:end);

e = ref-p;
de = -v;

uD = Kd*de;
uP = Kp*e;

u = uD + uP;
end