function [gamma,Q] = QLF_restricted(A_c)
ops = sdpsettings('verbose',0);
stop_flag = 0;
delta_gamma = 0.1;
gamma = max(real(eig(A_c)))-0.1;
%gamma = -0.4;
[n,m] = size(A_c); 
while(~stop_flag)
    P = sdpvar(n,n);
    F = [P>=eye(n),P<=4*eye(n), A_c'*P+P*A_c-gamma*P <=0];
    diagnostics=optimize(F,[],ops);
    if diagnostics.problem == 0
        Q = value(P);
        clear P;
        stop_flag = 1;
    elseif gamma > 3;
        disp('gamma too large, check');
        clear P;
        stop_flag = 1;
    else
        disp(gamma);
        gamma = gamma + delta_gamma;
    end        
end
%         disp(eig(Q));