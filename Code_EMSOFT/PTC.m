function out = dyn_func(t,x)

% Better equilibrium. Obtained numerically
% x_equilibrium = [ 0.988733871910758   1.000009438665233   1.220764191080577  -0.005929701497498]';

x_equilibrium = [  0.89877559179912198425437882605416;...
    1;...
 1.0779564350547793288065542452819;...
 -0.000000000000000026901249610632059296110104185358];


% x_equilibrium = zeros(4,1);

% Equilibrium we originally uesed (not accurate but converges to answer)
% x_equilibrium = [ 0.898775584353738  ; 0.999999788689746 ;  0.111111423310533 ];
x = x + x_equilibrium;

% System constants
c1=0.41328;
c2=200;         % Engine speed (rad/sec)
c3=-0.366;      % Coefficient for Pumping polynomial
c4=0.08979;     % Coefficient for Pumping polynomial
c5=-0.0337;     % Coefficient for Pumping polynomial
c6=0.0001;      % Coefficient for Pumping polynomial
c7=2.821;       % Coefficient for f(theta) polynomial
c8=-0.05231;    % Coefficient for f(theta) polynomial
c9=0.10299;     % Coefficient for f(theta) polynomial
c10=-0.00063;   % Coefficient for f(theta) polynomial
c11=1.0;        % Atmospheric pressure (bar)
c12=14.7;       % Stoichiometric ration of air-to-fuel
c13=0.9;        % Factor representing error between estimated and actual manifold pressure
c14=0.4;        % Proportional gain for PI controller
c15=0.4;        % Integral gain for PI controller
c16=1.0;        % lambda reference (desired A/F to stoichiometric ratio
c17=0.94513;    % Coefficient for $A/F$ polynomial
c18=-2.3981;    % Coefficient for $A/F$ polynomial
c19=1.4106;     % Coefficient for $A/F$ polynomial
c20=0.17882;    % Coefficient for $A/F$ polynomial
c21=-0.10835;   % Coefficient for $A/F$ polynomial
c22=-2.3421;    % Coefficient for square root polynomial
c23=2.7799;     % Coefficient for square root polynomial
c24=-0.3273;    % Coefficient for square root polynomial
u1 = 15.0;      % Throttle angle 15 or 25?

u1_hat = c7+c8*u1+c9*u1^2+c10*u1^3;

% Polynomial version
% out = [ c1*(2*u1_hat*(c22*x(1)^2+c23*x(1)+c24) - (c3+c4*c2*x(1)+c5*c2*x(1)^2+c6*x(1)*c2^2)) ;...
%     4*( c17+c18*((c13/c12)*(c3+c4*c2*x(3)+c5*c2*x(3)^2+c6*x(3)*c2^2)*(1+x(4)+c14*(x(2)-c16)))+c19*((c13/c12)*(c3+c4*c2*x(3)+c5*c2*x(3)^2+c6*x(3)*c2^2)*(1+x(4)+c14*(x(2)-c16)))^2+c20*(c3+c4*c2*x(1)+c5*c2*x(1)^2+c6*x(1)*c2^2)+c21*(c3+c4*c2*x(1)+c5*c2*x(1)^2+c6*x(1)*c2^2)*((c13/c12)*(c3+c4*c2*x(3)+c5*c2*x(3)^2+c6*x(3)*c2^2)*(1+x(4)+c14*(x(2)-c16))) - x(2) );...
%     c1*(2*u1_hat*(c22*x(1)^2+c23*x(1)+c24) - c13*(c3+c4*c2*x(3)+c5*c2*x(3)^2+c6*x(3)*c2^2));...
%     c15*(x(2)-c16) ];

out = [ c1*(2*u1_hat*sqrt((x(1)/c11)-(x(1)/c11)^2) - (c3 + c4*c2*x(1) + c5*c2*x(1)^2 + c6*x(1)*c2^2)); ...
          4*((c3 + c4*c2*x(1) + c5*c2*x(1)^2 + c6*x(1)*c2^2)/(c13*(c3 + c4*c2*x(3) + c5*c2*x(3)^2 + c6*x(3)*c2^2)*(1 + x(4) + c14*(x(2) - c16))) - x(2)); ...
          c1*(2*u1_hat*sqrt((x(1)/c11)-(x(1)/c11)^2) - c13*(c3 + c4*c2*x(3) + c5*c2*x(3)^2 + c6*c2^2*x(3)));...
          c15*(x(2) - c16)];


if any(imag(out)~=0)
    error(['Error: cost function returned an imaginary number.']);
end
end