function [ddddp_des,dddp,ddp,ddF_des,domega_des,tau_des] = snap_controller(K1,K2,K3,K4,K5,K6,p_ref,dp_ref,ddp_ref,dddp_ref,yaw_ref,dyaw_ref,x,F,dF,I,m)
% SNAP CONTROLLER snap controller for a 6DoF Euler-Angle-Free quadrotor
% drone model
%   Based on "Quadcopter Tracking Using Euler-Angle-Free Flatness-Based
%   Control" by Aeris El Asslouj and Hossein Rastgoftar (arXiv:2212.01540v1 [cs.RO] 3 Dec 2022)
%
% INPUT:
%   K1..6 : Linear Gain Matrices
%   d..p_ref : Spatial Coordinates of Reference
%   d..yaw_ref : Yaw of Reference
%   x : Drone's State-Space vector
%   d..F : Previous Thrust Force
%   I : Inertial momentum matrix diagonal
%   m : mass of drone
% OUTPUT:
%   ddF_des : Desired Thrust second derivative (must be integrated outside)
%   domega_des : Desired Angular velocity first derivative (is used to calculate tau_des outside)

p = x(1:3);
v = x(4:6);
A_vec = x(7:15);
A = reshape(A_vec,3,[]);
omega = x(16:end);

J = diag(I);
g = 9.8;

% Inertial Frame cardinal vecs
iI = [1;0;0];
jI = [0;1;0];
kI = [0;0;1];

% Body Frame cardinal vecs to Inertial Frame
iB = A*[1;0;0];
jB = A*[0;1;0];
kB = A*[0;0;1];

% Heading vector
h_aux = (iB'*iI).*iI+(iB'*jI).*jI;
h = h_aux/norm(h_aux);
kId = crossprodmat(kI)*h;

% 1) Position Controller (obtains ddddp_des)
ddp = (F*kB-m*g*kI)/m;
dddp = (dF*kB-F*(crossprodmat(omega)*kB))/m;
ddddp_des = -K1*(dddp-dddp_ref)-K2*(ddp-ddp_ref)-K3*(v-dp_ref)-K4*(p-p_ref);

% 2) Yaw Controller (obtains ddyaw_des)
yaw = atan2(iB'*jI,iB'*iI);
dyaw = ((omega(3)*kId)'*jB-(omega(2)*kId)'*kB)/(h'*iB);
ddyaw_des = -K5*(dyaw-dyaw_ref)-K6*(yaw-yaw_ref);

% 3) Map to thrust and torque (obtains domega_des and ddF_des)
h_F = m*ddddp_des-F*(crossprodmat(omega)*(crossprodmat(omega)*kB))-2*dF*(crossprodmat(omega)*kB);
ddF_des = h_F'*kB;

h_alpha = (h_F-ddF_des*kB)/F;
V = (ddyaw_des*h+dyaw*kId)'*iB+2*dyaw*h'*(crossprodmat(omega)*iB);
alpha_i = -h_alpha'*jB;
alpha_j = h_alpha'*iB;
alpha_k = (V+alpha_j*kId'*kB-kId'*(crossprodmat(omega)*(crossprodmat(omega)*iB)))/(kId'*jB);

domega_des = [alpha_i; alpha_j; alpha_k];

% 4) Calculate torque from domega_des (obtains tau_des)

tau_des = J*domega_des+crossprodmat(omega)*(J*omega);
end