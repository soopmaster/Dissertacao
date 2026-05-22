function  max_21 = Vanderpol_linear_min_22(P,center)
miu = 1;


% if Box(1,1) < -1 && Box(1,2) > -1
%     A(1,1) = 1.5;
% else
%     A(1,1) = max(-3*Box(1,1)-3.0*(Box(1,1)^2)/2,-3*Box(1,2)-3.0*(Box(1,2)^2)/2);
% end

sdpvar x1 x2;
sdpvar lambda;
sdpvar fpl;
sdpvector = monolist([x1;x2], [3]);
% We want to leave off the constant term
sdpvector = sdpvector(2:3);
V = (sdpvector-center)'*P*(sdpvector-center);

F2 = sos(-fpl+miu*((1-x1^2))-lambda*(1-V));
F2 = [F2, lambda>=0];
[sol,v,Q] = solvesos(F2,-fpl,sdpsettings('verbose',0));
x_value = double(fpl);
max_21 = x_value;