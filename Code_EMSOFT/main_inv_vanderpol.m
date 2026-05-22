clear;clc;
%%
%set up the initial conditions: total time; time step(optional); initial
%states; perturbation of the initial state;
% addpath(genpath('C:\Users\cfan10\Google Drive\LMI\yalmip'));
% addpath(genpath('C:\Users\cfan10\Google Drive\LMI\SeDuMi_1_3'));
% addpath(genpath('C:\Users\cfan10\Google Drive\LMI\SDPT3-3.02'))
%%
ops = sdpsettings('verbose',0);
enable_disp = 0;

dt = 0.002;
t_dur = 0:dt:10;
initial_state = [0.4,0.4];
%initial_state = [0.5,0.5];


Initial_P = [1 0; 0 1];
Initial_l = (0.1)^2;

delta_initial_up = [0,sqrt(Initial_l)];
%Lipschitz constant of te system, user input
Lip_constant = 1;
[t_precise,x_precise] = ode45(@inv_vanderpol,t_dur,initial_state);
[t_true,x_true] = ode45(@inv_vanderpol,t_dur,initial_state+delta_initial_up);
%[t_true_low,x_true_low] = ode45(@Engine,t_dur,initial_state+delta_initial_low);
if enable_disp
    figure;
    scatter(x_precise(:,1),x_precise(:,2));
    xlim = ylim;
    hold on;
    scatter(x_true(:,1),x_true(:,2),'r');
    scatter(initial_state(1),initial_state(2),180,'x','b','linewidth',2);
    scatter(initial_state(1)+delta_initial_up(1),initial_state(2)+delta_initial_up(2),180,'x','r','linewidth',2);
    xlabel('x1','FontSize', 24);
    ylabel('x2','FontSize', 24);
    plot_ellipsoid(x_precise(1,:)',Initial_P/Initial_l,'k');
end
step = 10;
index = 1:step:length(t_precise);
t = t_precise(index);
x = x_precise(index,:);


%
tic;
%%
%length of the simulation result
n = length(x);
%dimention of the system
dim = 2;
Reach_info = zeros(n,5);
Reach_info(1,:) = [Initial_P(1,:),Initial_P(2,:),Initial_l];
radius = zeros(length(t_precise),1);
radius(1) = sqrt(Initial_l);
%set
for i = 1:length(t)-1
    %
    M = [Reach_info(i,1:2);Reach_info(i,3:4)];
    level = Reach_info(i,5);
    %
    %     max_21 = inv_Vanderpol_linear_max_21(M/(level*exp(Lip_constant*(t(i+1)-t(i)))),x(i,:)');
    %     min_21 = inv_Vanderpol_linear_min_21(M/(level*exp(Lip_constant*(t(i+1)-t(i)))),x(i,:)');
    %     max_22 = inv_Vanderpol_linear_max_22(M/(level*exp(Lip_constant*(t(i+1)-t(i)))),x(i,:)');
    %     min_22 = inv_Vanderpol_linear_min_22(M/(level*exp(Lip_constant*(t(i+1)-t(i)))),x(i,:)');
    %     [V,D] = eig(inv_Vanderpol_Jac(mean(x_precise(index(i):index(i+1)-1,:))));
    for cnt1 = index(i):index(i+1)-1
        if mod(cnt1,10)==1
            [V,D] = eig(inv_Vanderpol_Jac(mean(x_precise(cnt1:cnt1+9,:))));
            clear D;
        end
        Box = zeros(2,2);
        rad = radius(cnt1)+ norm(x_precise(cnt1+1)-x_precise(cnt1));
        Box(1,1) = x_precise(cnt1,1)-rad*exp(Lip_constant*(t_precise(cnt1+1)-t_precise(cnt1)));
        Box(1,2) = x_precise(cnt1,1)+rad*exp(Lip_constant*(t_precise(cnt1+1)-t_precise(cnt1)));
        Box(2,1) = x_precise(cnt1,2)-rad*exp(Lip_constant*(t_precise(cnt1+1)-t_precise(cnt1)));
        Box(2,2) = x_precise(cnt1,2)+rad*exp(Lip_constant*(t_precise(cnt1+1)-t_precise(cnt1)));
        %draw_rectangle(Box);
        Jacobian = inv_Vanderpol_Jac(x_precise(cnt1,:));
        %         disp(max(real(eig(Jacobian))));
        D = V\Jacobian*V;
        %[V,D] = eig(Jacobian);
        lambda = max(eig((D+D')/2));
        %lambda = real(max(eig((D+D)/2)));
        delta = norm(inv_vanderpol_error(Box))/2;
        %disp(cnt1);
        %disp(delta);
        change_rate = lambda+delta;
        %disp(change_rate);
        radius(cnt1+1) = radius(cnt1)*exp(change_rate*(t_precise(cnt1+1)-t_precise(cnt1)));
        %circle(x_precise(cnt1,1),x_precise(cnt1,2),radius(cnt1),'r');
        M_values(cnt1,:) = inv_vanderpol_extreme(Box);
    end
    %radius(index(i):index(i+1)-1,:) = radius(index(i):index(i+1)-1,:)*cond(V);
    max_21 = max(M_values(index(i):index(i+1)-1,1));
    min_21 = min(M_values(index(i):index(i+1)-1,2));
    max_22 = max(M_values(index(i):index(i+1)-1,3));
    min_22 = min(M_values(index(i):index(i+1)-1,4));
    A1 = [0 -1;max_21 max_22];
    A2 = [0 -1;max_21 min_22];
    A3 = [0 -1;min_21 max_22];
    A4 = [0 -1;min_21 min_22];
    %disp(A2)
    eigens = [max(real(eig(A1))),max(real(eig(A2))),max(real(eig(A3))),max(real(eig(A4)))];
    if enable_disp
        disp('--------------');
        disp(i);
        disp(eigens);
    end
    
    shrink = max(eigens)-0.1;
    %if max(eig(A1'*M+M*A1))<0&& max(eig(A4'*M+M*A4))
    if max(eig(A1'*M+M*A1-shrink*M))<0&& max(eig(A4'*M+M*A4-shrink*M))<0&& max(eig(A2'*M+M*A2-shrink*M))<0
        Q = M;
        level_next = level*exp(shrink*(t(i+1)-t(i)));
        if enable_disp
            disp('remain P in last round')
        end
        P_level_original = level;
        gamma = shrink;
    else
        %Lyapunov-like function
        stop_flag = 0;
        delta_gamma = 0.1;
        %gamma = max(eigens);
        A_c = [0, -1; (min_21+max_21)/2, (min_22+max_22)/2];
        gamma = min(real(eig(A_c)));
        while(~stop_flag)
            P = sdpvar(2,2);
            F = [P>=eye(2), A_c'*P+P*A_c-gamma*P <=0];
            diagnostics=optimize(F,[],ops);
            if diagnostics.problem == 0
                Q = value(P);
                clear P;
                stop_flag = 1;
            elseif gamma > 3;
                if enable_disp
                    disp('gamma too large, check');
                end
                clear P;
                stop_flag = 1;
            else
                if enable_disp
                    disp(gamma);
                end
                gamma = gamma + delta_gamma;
            end
        end
        %         disp(eig(Q));
        if enable_disp
            disp(gamma);
        end
        Delta_matrix = [0, 0; (max_21-min_21)/2, (max_22-min_22)/2];
        SDP_Delta_matrix = Delta_matrix'*Q-Q*Delta_matrix;
        if enable_disp
            disp(Delta_matrix);
        end
        epsilon_matrix = sqrt(norm(SDP_Delta_matrix,1)*norm(SDP_Delta_matrix,Inf));
        %disp(gamma*max(eig(Q)) + epsilon_matrix);
        gamma = (gamma*min(eig(Q)) + epsilon_matrix)/min(eig(Q));
        if enable_disp
            disp(gamma);
        end
        %     disp('eigenvalue of P');disp(eig(Q));
        %     disp('eigenvalue of A1T*P+P*A1+m*P');disp(eig(A1'*Q+Q*A1+gamma*Q))
        %     disp('eigenvalue of A2T*P+P*A2+m*P ');disp(eig(A2'*Q+Q*A2+gamma*Q))
        P_level_original = min_level_set(Q,M,level);
        level_next = P_level_original*(exp(gamma*(t(i+1)-t(i))));
        if enable_disp
            plot_ellipsoid(x(i+1,:)',Q/level_next,'k');
        end
    end
    for cnt = index(i):index(i)+step-1
        if enable_disp
            plot_ellipsoid(x_precise(cnt,:)',Q/(P_level_original*(exp(gamma*(t_precise(cnt)-t(i))))),'k');
        end
%         fp_checker_inner = fix_point_checker(x_precise(cnt+1,:),x_precise(cnt,:),Q/(P_level_original*(exp(gamma*(t_precise(cnt+1)-t(i))))),Q/(P_level_original*(exp(gamma*(t_precise(cnt)-t(i))))));
%         if fp_checker_inner == 1
%             disp('find invariant set in between');
%             break;
%         end
    end
    if enable_disp
        plot_ellipsoid(x(i+1,:)',Q/level_next,'k');
    end
    Reach_info(i+1,:) = [Q(1,:),Q(2,:),level_next];
    radius(index(i+1)) = sqrt(max(eig(inv(Q/level_next))));
%     fp_checker = fix_point_checker(x(i+1,:),x(i,:),Q/level_next,M/level);
%     if fp_checker == 1 || fp_checker_inner == 1
%         disp('find invariant set');
%         break;
%     end
end

Reach_final = reshape(Reach_info(i,1:4),2,2)/Reach_info(i,5);
dia_final = sqrt(eig(inv(Reach_final)));
disp(dia_final);

toc;

