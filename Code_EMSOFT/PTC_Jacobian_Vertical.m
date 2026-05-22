function Vertical = PTC_Jacobian_Vertical(x_max,x_min,J_max,J_min,index)

Vertical = zeros(4,4);
% Compute the Jacoabian at the two given points, can be directly used if
% the certain component is monotonic
x_equilibrium = [  0.89877559179912198425437882605416;...
    1;...
 1.0779564350547793288065542452819;...
 -0.000000000000000026901249610632059296110104185358];


% x_equilibrium = zeros(4,1);

% Equilibrium we originally uesed (not accurate but converges to answer)
% x_equilibrium = [ 0.898775584353738  ; 0.999999788689746 ;  0.111111423310533 ];
x_max = x_max' + x_equilibrium;
x_min = x_min' + x_equilibrium;
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


p_max = x_max(1);
r_max = x_max(2);
pe_max = x_max(3);
i_max = x_max(4);

p_min = x_min(1);
r_min = x_min(2);
pe_min = x_min(3);
i_min = x_min(4);

%if p \in [0,0.5], sqrt((p/c11)-(p/c11)^2) is monotonic, true;
%(1/c11-2*p/(c11^2))/(2*sqrt((p/c11)-(p/c11)^2)); is also monotonic
%decreasing
%if p \in [0,1.6], c3 + c4*c2*p + c5*c2*p^2 + c6*p*c2^2 is monotonic, true,
%and c4*c2 + 2*c5*c2*p + c6*(c2^2); is also monotonic decreasing
%if pe \in [0,1.6], c3 + c4*c2*pe + c5*c2*pe^2 + c6*pe*c2^2 is monotonic,
%true, and its derivative is mononial
%if r \in [-1.5,\inf], 1/(1+i+c14*(r-c_16)) is monomial, true


square_p_diff_p_max = (1/c11-2*p_max/(c11^2))/(2*sqrt((p_max/c11)-(p_max/c11)^2));
poly_p_diff_p_max = c4*c2 + 2*c5*c2*p_max + c6*(c2^2);
square_p_diff_p_min = (1/c11-2*p_min/(c11^2))/(2*sqrt((p_min/c11)-(p_min/c11)^2));
poly_p_diff_p_min = c4*c2 + 2*c5*c2*p_min + c6*(c2^2);

if index(1) == 1
    Vertical(1,1) = c1*(2*u1_hat*max(square_p_diff_p_max,square_p_diff_p_min)-min(poly_p_diff_p_max,poly_p_diff_p_min));
else
    Vertical(1,1) = c1*(2*u1_hat*min(square_p_diff_p_max,square_p_diff_p_min)-max(poly_p_diff_p_max,poly_p_diff_p_min));
end

if index(2) == 1
    Vertical(2,1) = max(J_max(2,1),J_min(2,1));
else
    Vertical(2,1) = min(J_max(2,1),J_min(2,1));
end

if index(3) == 1
    Vertical(2,2) = max(J_max(2,2),J_min(2,2));
else
    Vertical(2,2) = min(J_max(2,2),J_min(2,2));
end

if index(4) == 1
    Vertical(2,3) = max(J_max(2,3),J_min(2,3));
else
    Vertical(2,3) = min(J_max(2,3),J_min(2,3));
end

if index(5) == 1
    Vertical(2,4) = max(J_max(2,4),J_min(2,4));
else
    Vertical(2,4) = min(J_max(2,4),J_min(2,4));
end

if index(6) == 1
    Vertical(3,1) = max(J_max(3,1),J_min(3,1));
else
    Vertical(3,1) = min(J_max(3,1),J_min(3,1));
end

if index(7) == 1
    Vertical(3,3) = max(J_max(3,3),J_min(3,3));
else
    Vertical(3,3) = min(J_max(3,3),J_min(3,3));
end

Vertical(4,2) = c15;