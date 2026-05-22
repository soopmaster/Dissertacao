function  max_21 = inv_Vanderpol_linear_max_22(P,center)

% if Box(1,1) < -1 && Box(1,2) > -1
%     A(1,1) = 1.5;
% else
%     A(1,1) = max(-3*Box(1,1)-3.0*(Box(1,1)^2)/2,-3*Box(1,2)-3.0*(Box(1,2)^2)/2);
% end

sdpvar x1 x2;
sdpvar lambda;
sdpvar fp;
sdpvector = monolist([x1;x2], [3]);
% We want to leave off the constant term
sdpvector = sdpvector(2:3);
V = (sdpvector-center)'*P*(sdpvector-center);

F = sos(fp-((x1^2)-1)-lambda*(1-V));
F = [F, lambda>=0];
[sol,v,Q] = solvesos(F,fp);
x_value = double(fp);
max_21 = x_value;