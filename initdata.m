% Called Back from Simulink mdl -> tracking_N
clear all
close all
clc
n13 = 5.0*10^-3;
n32 = 5.0*10^-3;
n30 = 1.0*10^-4;
n0 = 7.0*10^-3;
di = 0.4;
Atk = pi*(di/2)^2;

% I. Model :: SISO
% -------------------------------

A_mimo = [-n13 0 n13; 0 -n32-n0 n32; n13 n32 -(n32+n30+n13)]/Atk
B_mimo = 8*10^-2*[1 0; 0 1; 0 0]/Atk
C_mimo = [1 0 0; 0 1 0]
D_mimo = [0 0; 0 0]

sys_mimo = ss(A_mimo,B_mimo,C_mimo,D_mimo);
si = stepinfo(sys_mimo);

%A_int = [A_mimo zeros(4,1); -C_mimo zeros(2,1)]
%B_int = [B_mimo; -D_siso]
%C_int = [C_mimo 0]
%D_int = D_mimo

% II. State Feedback Control :: N
% -------------------------------
sperc = 3; % overshoot 2%

Ts_closed_loop = si.SettlingTime
Ts_design_set = 0.8*Ts_closed_loop % << Ts_design_max
Dwn = -log(1/100)/Ts_design_set;

D = sqrt(log(sperc/100)^2/(log(sperc/100)^2+pi^2));
theta_ang = acos(D);
Wn = Dwn/D;

%Ideal transfer function
%hsideal = ( Wn^2 /(s^2 + (2*D*Wn)*s + Wn^2));
%PolesDom = pole(hsideal)
%poles = [PolesDom' 8*real(PolesDom(1,1))];

lambda_1 = -Dwn;
lambda_2 = 5*lambda_1;
lambda_3 = 8*lambda_1;
%lambda_4 = 10*lambda_1;
lambda_spec = [lambda_1 lambda_2 lambda_3]
K = place(A_mimo,B_mimo,lambda_spec)
%K = [Kn(1) Kn(2) Kn(3) 0]
%Kn = Kn(4)

%Nxu = [A_mimo B_mimo;C_mimo D_mimo]^-1*[zeros(4,1); 1]
%N = Nxu(end)+ K*Nxu(1:end-1)

N = -(C_mimo*(A_mimo-B_mimo*K)^-1*B_mimo)^-1