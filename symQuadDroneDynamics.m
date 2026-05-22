function [dxsym] = symQuadDroneDynamics(symvars,I,M)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
L = 1; % arm length
R = 0.25;
m = 0.1; % area mass of drone
I = [2*m*R^2*((L+R/2)^3-(R/2)^3)/3; 2*m*R^2*((L+R/2)^3-(R/2)^3)/3; m*R*(2*L+R)*2*(R/2)^3/3]; % moment of inertia
M = m*R^2*(2*L+R); % total mass of drone
J = diag(I); % Inertia Matrix
g = 9.8;

niu_v = 0.01;
niu_w = 0.7;

aBsym = 1/M*Fsym.*[0;0;1]; %  Linear acceleration in relation to body

dpsym = vsym;
dvsym = Asym*aBsym+[0;0;-g]-niu_v.*vsym;
dAsym = crossprodmat(omegaBsym)*Asym;
domegaBsym = J\(tauBsym-crossprodmat(omegaBsym)*J*omegaBsym)-niu_w.*omegaBsym;

dAsym_vec = reshape(dAsym,[],1);

dxsym = [dpsym; dvsym; domegaBsym; dAsym_vec];
end