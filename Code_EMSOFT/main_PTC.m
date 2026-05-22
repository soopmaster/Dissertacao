clc;
clear;
close;
%
addpath(genpath('..\SeDuMi_1_3'));
addpath(genpath('..\yalmip'));

%%
%initials
inner_step = 20;
enable_disp = 0;
INDEX_CH = zeros(128,7);
for i = 0 : 127
    INDEX_CH(i+1,:) = de2bi(i,7);
end

tic;
%% set all the initial values
dt = 5e-4;
t_dur = 0:dt:5;
initial_value = 0.01*ones(4,1);
dim = 4; %dimension of the system
%time_span = [0,10];
options = odeset('RelTol',1e-3,'MaxStep',0.1);

[t_precise,x_precise] = ode45(@PTC,t_dur,initial_value,options);

Initial_P = eye(4);
Initial_l = (2e-4)^2;

if enable_disp
    figure;
    plot(x_precise(:,2),x_precise(:,4));
    hold on;
    plot_ellipsoid(x_precise(1,[2,4])',Initial_P([2,4],[2,4])/Initial_l,'k');
end

index = 1:inner_step:length(t_precise);
t = t_precise(index);
x = x_precise(index,:);

%%
n = length(x);
%dimention of the system

Reach_info = zeros(n,17);
Reach_info(1,:) = [Initial_P(1,:),Initial_P(2,:),Initial_P(3,:),Initial_P(4,:),Initial_l];

for i = 1:length(t)-1
    %for i = 1:20
    %
    if enable_disp
        disp('----------------');
        disp(['Computing ', num2str(t(i)), ' time interval ']);
    end
    M = [Reach_info(i,1:4);Reach_info(i,5:8);Reach_info(i,9:12);Reach_info(i,13:16)];
    level = Reach_info(i,17);
    
    Average_value = mean(x_precise(index(i):index(i+1)-1,:));
    A_c = PTC_Jacobian(Average_value);
    shrink = max(real(eig(A_c)))-0.1;
    %shrink = -0.01;
    
    %Local Lipschitz constant
    Lip = norm(A_c);
    %delta_theta = (abs(x(i+1,:)-x(i,:))/2+sqrt(max(eig(inv(M/level)))))*exp(Lip*(t(i+1)-t(i)));
    x_relative_max = (max(x(i+1,:),x(i,:))+sqrt(max(eig(inv(M/level)))))*exp(Lip*(t(i+1)-t(i)));
    x_relative_min = (min(x(i+1,:),x(i,:))-sqrt(max(eig(inv(M/level)))))*exp(Lip*(t(i+1)-t(i)));
    check_ch_flag = 1;
    J_max = PTC_Jacobian(x_relative_max);
    J_min = PTC_Jacobian(x_relative_min);
    for cnt_ch = 0:127
        %        index_ch = de2bi(cnt_ch,7);
        index_ch = INDEX_CH(cnt_ch+1,:);
        A_ch = PTC_Jacobian_Vertical(x_relative_max,x_relative_min,J_max,J_min,index_ch);
        if max(eig(A_ch'*M+M*A_ch-shrink*M)) > 0
            check_ch_flag = 0;
            if enable_disp
                disp('Convex hull check fail');
            end
            break;
        end
    end
    
    if check_ch_flag == 1
        Q = M;
        level_next = level*exp(shrink*(t(i+1)-t(i)));
        if enable_disp
            disp('remain P in last round');
            disp('expontial change rate:');
            disp(shrink);
        end
    else
        if enable_disp
            disp('Compute new P')
        end
        %Lyapunov-like function
        [gamma,Q] = QLF(A_c);
        % disp(gamma);
        check_new_p_flag = 1;
        for cnt_ch = 0:127
%             index_ch = de2bi(cnt_ch,7);
            index_ch = INDEX_CH(cnt_ch+1,:);
            A_ch = PTC_Jacobian_Vertical(x_relative_max,x_relative_min,J_max,J_min,index_ch);
            if max(eig(A_ch'*Q+Q*A_ch-gamma*Q)) > 0
                check_new_p_flag = 0;
                if enable_disp
                    disp('Convex hull check fail when computing new P');
                end
                break;
            end
        end
        if check_new_p_flag == 0
            delta_theta = (abs(x(i+1,:)-x(i,:))+sqrt(max(eig(inv(M/level)))))*exp(Lip*(t(i+1)-t(i)));
            Delta_matrix = PTC_Jacobian_Delta(Average_value+delta_theta,Average_value-delta_theta)/2;
            SDP_Delta_matrix = Delta_matrix'*Q+Q*Delta_matrix;
            epsilon_matrix = sqrt(norm(SDP_Delta_matrix,1)*norm(SDP_Delta_matrix,Inf));
            %disp(gamma*max(eig(Q)) + epsilon_matrix);
            gamma = (gamma*min(eig(Q)) + epsilon_matrix)/min(eig(Q));
            if enable_disp
                disp('using interval matrix disturbance to over-approximate \gamma')
                disp('expontial change rate:');
                disp(gamma);
            end
        else
            if enable_disp
                disp('Computed P satifies convex hull constrains')
                disp('expontial change rate:');
                disp(gamma);
            end
        end
        %     disp('eigenvalue of P');disp(eig(Q));
        P_level_original = min_level_set4(Q,M,level);
        level_next = P_level_original*(exp(gamma*(t(i+1)-t(i))));
        
    end
    for cnt = index(i):index(i)+inner_step-1
        %plot_ellipsoid(x_precise(cnt,[2,4])',Q([2,4],[2,4])/(P_level_original*(exp(-gamma*(t_precise(cnt)-t(i))))),'k');
        fp_checker_inner = fix_point_checker(x_precise(cnt+1,:),x_precise(cnt,:),Q/(P_level_original*(exp(gamma*(t_precise(cnt+1)-t(i))))),Q/(P_level_original*(exp(gamma*(t_precise(cnt)-t(i))))));
        if fp_checker_inner == 1
            if enable_disp
                disp('find invariant set in between');
            end
            break;
        end
    end
    if enable_disp
        plot_ellipsoid(x(i+1,[2,4])',Q([2,4],[2,4])/level_next,'k');
    end
    Reach_info(i+1,:) = [Q(1,:),Q(2,:),Q(3,:),Q(4,:),level_next];
    fp_checker = fix_point_checker(x(i+1,:),x(i,:),Q/level_next,M/level);
    if fp_checker == 1 || fp_checker_inner == 1
        disp('find invariant set');
        break;
    end
end
toc;

Reach_final = reshape(Reach_info(i,1:16),4,4)/Reach_info(i,17);
dia_final = sqrt(eig(inv(Reach_final)));
disp(dia_final);
%%
% A1 = reshape(Reach_info(3000,1:16),4,4);
% A2 = reshape(Reach_info(3001,1:16),4,4);
% l1 = Reach_info(3000,17);
% l2 = Reach_info(3001,17);
% plot_ellipsoid(x(3000,[2,4])',A1([2,4],[2,4]),'r');
% hold on;
% plot_ellipsoid(x(3001,[2,4])',A2([2,4],[2,4]),'b');
% toc;
%To do: fix point check








% [t,x] = ode45(@PTC,t_span,initial_states);
% for i = 1:length(t)-1
%     Jacobian = PTC_Jacobian(x(i,:));
%     disp(norm(Jacobian))
%y(i) = max(real(eig(Jacobian)));
%disp(max(real(eig(Jacobian))));
%     [gamma,Q] = QLF(Jacobian);
%     disp(gamma);
%     Delta = PTC_Jacobian_Delta(x(i+1,:),x(i,:));
%     disp(norm(Delta));
% end
%plot(x);