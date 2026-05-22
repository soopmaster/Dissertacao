function J = PTC_Jacobian(x)

J = zeros(4,4);

x_equilibrium = [  0.89877559179912198425437882605416;...
    1;...
 1.0779564350547793288065542452819;...
 -0.000000000000000026901249610632059296110104185358];


% x_equilibrium = zeros(4,1);

% Equilibrium we originally uesed (not accurate but converges to answer)
% x_equilibrium = [ 0.898775584353738  ; 0.999999788689746 ;  0.111111423310533 ];
x = x' + x_equilibrium;
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


p = x(1);
r = x(2);
pe = x(3);
i = x(4);

square_p_diff_p = (1/c11-2*p/(c11^2))/(2*sqrt((p/c11)-(p/c11)^2));
poly_p_diff_p = c4*c2 + 2*c5*c2*p + c6*(c2^2);
poly_pe_diff_pe = c4*c2 + 2*c5*c2*pe + c6*(c2^2);
poly_r_diff_p = poly_p_diff_p/(c13*(c3 + c4*c2*pe + c5*c2*pe^2 + c6*pe*c2^2)*(1 + i + c14*(r - c16)));
poly_r_diff_r = -c14*(c3 + c4*c2*p + c5*c2*p^2 + c6*p*c2^2)/(c13*(c3 + c4*c2*pe + c5*c2*pe^2 + c6*pe*c2^2)*(1 + i + c14*(r - c16))^2)-1;
poly_r_diff_pe = -((c3 + c4*c2*p + c5*c2*p^2 + c6*p*c2^2)*(poly_pe_diff_pe))/(c13*(c3 + c4*c2*pe + c5*c2*pe^2 + c6*pe*c2^2)^2*(1 + i + c14*(r - c16)));
poly_r_diff_i = -(c3 + c4*c2*p + c5*c2*p^2 + c6*p*c2^2)/(c13*(c3 + c4*c2*pe + c5*c2*pe^2 + c6*pe*c2^2)*(1 + i + c14*(r - c16))^2);

J = [c1*(2*u1_hat*square_p_diff_p-poly_p_diff_p), 0, 0, 0;...
    4*poly_r_diff_p, 4*poly_r_diff_r, 4*poly_r_diff_pe, 4*poly_r_diff_i;...
    c1*2*u1_hat*square_p_diff_p, 0, -c1*c13*poly_pe_diff_pe, 0;...
    0 , c15, 0, 0];


