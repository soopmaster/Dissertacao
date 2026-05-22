clear;clc;
%%
%set up the initial conditions: total time; time step(optional); initial
%states; perturbation of the initial state;

ops = sdpsettings('verbose',0);

enable_disp = 0;

dt = 0.001;
t_dur = 0:dt:20;
% t_dur = [0,10];
initial_state = [0.2,0.2];

Initial_P = [2 0; 0 1];
Initial_l = 0.1^2;

delta_initial_up = [0,sqrt(Initial_l)];
%Lipschitz constant of te system, user input
Lip_constant = 1.22;
[t_precise,x_precise] = ode45(@Engine,t_dur,initial_state);
[t_true,x_true] = ode45(@Engine,t_dur,initial_state+delta_initial_up);
%[t_true_low,x_true_low] = ode45(@Engine,t_dur,initial_state+delta_initial_low);
if enable_disp
    figure;
    scatter(x_precise(:,1),x_precise(:,2));
    axis equal;
    hold on;
    scatter(x_true(:,1),x_true(:,2),'r');
    scatter(initial_state(1),initial_state(2),180,'x','b','linewidth',2);
    scatter(initial_state(1)+delta_initial_up(1),initial_state(2)+delta_initial_up(2),180,'x','r','linewidth',2);
    xlabel('x1','FontSize', 24);
    ylabel('x2','FontSize', 24);
end

%get some larger time interval for the ellipse.
index = 1:100:length(t_precise);
t = t_precise(index);
x = x_precise(index,:);
if enable_disp
    plot_ellipsoid(x(1,:)',Initial_P/Initial_l,'k');
end

tic;
%%
%length of the simulation result
n = length(x);   
%dimention of the system
dim = 2;  
Reach_info = zeros(n,5);
Reach_info(1,:) = [Initial_P(1,:),Initial_P(2,:),Initial_l];
%set 
for i = 1:n-1
    %Box is the hyper rectangle blowed by lipschitz constant
    Box = zeros(dim,2); 
    M = [Reach_info(i,1:2);Reach_info(i,3:4)];
    level = Reach_info(i,5);
    %square_extension = sqrt(sum(eig((M/level)^-1)));
    square_extension = sqrt(max(eig((M/level)^-1)));
    for j = 1:dim
            Box(j,1) = min(x(i,j),x(i+1,j)) - square_extension*exp(Lip_constant*(t(i+1)-t(i)));
            Box(j,2) = max(x(i,j),x(i+1,j)) + square_extension*exp(Lip_constant*(t(i+1)-t(i)));
    end
%      draw_rectangle(Box);
%      hold on;
    %A is the "maximized" linear dynamix matrix   
    A1 = Engine_linear_max(Box);
    A2 = Engine_linear_min(Box);
%     A1 = Engine_linear_max_ellipse(M/(level*exp(Lip_constant*(t(i+1)-t(i)))),x(i,:)');
%     A2 = Engine_linear_min_ellipse(M/(level*exp(Lip_constant*(t(i+1)-t(i)))),x(i,:)');
    if enable_disp
        disp('--------------');
        disp(i);
        disp(max(real(eig(A1))));
        disp(max(real(eig(A2))));
    end
    
    %Lyapunov-like function
    stop_flag = 0;
    delta_gamma = 0.1;
    gamma = max(max(real(eig(A1))),max(real(eig(A2))));
    P = M;
    if max(eig(A1'*P+P*A1-gamma*P)) <=0 && max(eig(A2'*P+P*A2-gamma*P)) <=0
        if enable_disp
            disp('remain P in last round')
        end
        Q = M;
        level_next = level*(exp(gamma*(t(i+1)-t(i))));
        P_level_original = level;
    else
        eigens = [real(eig(A1)); real(eig(A2))];
        gamma = (max(eigens) + max(eigens))/2;
        while(~stop_flag)
            P = sdpvar(2,2);
            F = [P>=eye(2), A1'*P+P*A1-gamma*P <=0, A2'*P+P*A2-gamma*P <=0];
            diagnostics=optimize(F,[],ops);
            if diagnostics.problem == 0
                Q = value(P);
                clear P;
                stop_flag = 1;
            elseif gamma < -3;
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
    %     disp(Q);
        if enable_disp
            disp(gamma);
        end
    %     disp('eigenvalue of P');disp(eig(Q));
    %     disp('eigenvalue of A1T*P+P*A1+m*P');disp(eig(A1'*Q+Q*A1+gamma*Q))
    %     disp('eigenvalue of A2T*P+P*A2+m*P ');disp(eig(A2'*Q+Q*A2+gamma*Q))
        P_level_original = min_level_set(Q,M,level);
        level_next = P_level_original*(exp(gamma*(t(i+1)-t(i))));
    end
    for cnt = index(i):index(i+1)-1
        if enable_disp
            plot_ellipsoid(x_precise(cnt,:)',Q/(P_level_original*(exp(-gamma*(t_precise(cnt)-t(i))))),'m');
        end
%         fp_checker_inner = fix_point_checker(x_precise(cnt+1,:),x_precise(cnt,:),Q/(P_level_original*(exp(gamma*(t_precise(cnt+1)-t(i))))),Q/(P_level_original*(exp(gamma*(t_precise(cnt)-t(i))))));
%         if fp_checker_inner == 1
%             disp('find invariant set in between');
%         break;
%         end
    end
    %plot_ellipsoid(x(i+1,:)',Q/level_next,'k');
    Reach_info(i+1,:) = [Q(1,:),Q(2,:),level_next];
%     fp_checker = fix_point_checker(x(i+1,:),x(i,:),Q/level_next,M/level);
%     if fp_checker == 1 || fp_checker_inner == 1
%         disp('find invariant set');
%     break;
%     end
end
toc;

Reach_final = reshape(Reach_info(i,1:4),2,2)/Reach_info(i,5);
dia_final = eig(inv(Reach_final));
disp(dia_final);

%plotting some result

%figure;
% scatter(x(:,1),x(:,2));
% hold on;
% scatter(x_true(:,1),x_true(:,2),'r');
% scatter(initial_state(1),initial_state(2),180,'x','b','linewidth',2);
% scatter(initial_state(1)+delta_initial_up(1),initial_state(2)+delta_initial_up(2),180,'x','r','linewidth',2);
% xlabel('x1','FontSize', 24);
% ylabel('x2','FontSize', 24);

% figure;
% subplot(1,2,1)
% scatter(t,x(:,1));
% hold on;
% scatter(t_true,x_true(:,1),'r');
% plot(t,Res_up(:,1),'k','linewidth',2);
% plot(t,Res_low(:,1),'k','linewidth',2);
% xlabel('t','FontSize', 24);
% ylabel('x1','FontSize', 24);
% 
% subplot(1,2,2)
% scatter(t,x(:,2));
% hold on;
% scatter(t_true,x_true(:,2),'r');
% plot(t,Res_up(:,2),'k','linewidth',2);
% plot(t,Res_low(:,2),'k','linewidth',2);
% xlabel('t','FontSize', 24);
% ylabel('x2','FontSize', 24);
% %
% figure;
% subplot(1,2,1)
% plot(t,Reach_dia_up(:,1));hold on;plot(t_true,x_true(:,1)-x(:,1),'r');
% plot(t,Reach_dia_low(:,1),'k');plot(t_true_low,x_true_low(:,1)-x(:,1),'r');
% subplot(1,2,2)
% plot(t,Reach_dia_up(:,2));hold on;plot(t_true,x_true(:,2)-x(:,2),'r');
% plot(t,Reach_dia_low(:,2),'k');plot(t_true,x_true_low(:,2)-x(:,2),'r');
% %
% figure;
% plot(Reach_dia_up(:,1),Reach_dia_up(:,2));hold on;plot(x_true(:,1)-x(:,1),x_true(:,2)-x(:,2),'r');
% plot(Reach_dia_low(:,1),Reach_dia_low(:,2),'k');plot(x_true_low(:,1)-x(:,1),x_true_low(:,2)-x(:,2),'r');