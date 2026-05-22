clear;
clc;
%
% addpath(genpath('..\..\SeDuMi_1_3'));
% addpath(genpath('..\..\yalmip'));

%%
%initials
inner_step = 10;
enable_disp = 0;
enable_plot = 0;

INDEX_CH = zeros(2^4,4);
for i = 0 : 2^4-1
    INDEX_CH(i+1,:) = de2bi(i,4);
end

tic;
%% set all the initial values
dt = 0.001;
t_dur = 0:dt:2;
initial_value = 0.1*ones(7,1);
dim = 7; %dimension of the system
%time_span = [0,10];
options = odeset('RelTol',1e-3,'MaxStep',0.1);

[t_precise,x_precise] = ode45(@biology,t_dur,initial_value,options);

Initial_P = eye(7);
Initial_l = (1e-4)^2;

if enable_plot
    figure;
    plot(x_precise(:,2),x_precise(:,4));
    hold on;
    plot_ellipsoid(x_precise(1,[2,4])',Initial_P([2,4],[2,4])/Initial_l,'k');
end

index = 1:inner_step:length(t_precise);
t = t_precise(index);
x = x_precise(index,:);

tic;
%%
n = length(x);   
%dimention of the system

Reach_info = zeros(n,50);
Reach_info(1,:) = [reshape(Initial_P,1,49),Initial_l];

for i = 1:length(t)-1
%for i = 1:1
    %
    if enable_disp
        disp('----------------');
        disp(['Computing ', num2str(t(i)), ' time interval ']);
    end
    M = reshape(Reach_info(i,1:49),7,7);
    level = Reach_info(i,50);

    Average_value = mean(x_precise(index(i):index(i+1)-1,:));
    A_c = biology_jac(Average_value);
    shrink = max(real(eig(A_c)))+0.2;
    %shrink = -0.01;
    
    %Local Lipschitz constant
    Lip = norm(A_c);   
   %delta_theta = (abs(x(i+1,:)-x(i,:))/2+sqrt(max(eig(inv(M/level)))))*exp(Lip*(t(i+1)-t(i)));
    x_relative_max = (max(x(i+1,:),x(i,:))+sqrt(max(eig(inv(M/level)))))*exp(Lip*(t(i+1)-t(i)));
    x_relative_min = (min(x(i+1,:),x(i,:))-sqrt(max(eig(inv(M/level)))))*exp(Lip*(t(i+1)-t(i)));
    %if max(eig(A_c'*M+M*A_c-shrink*M)) > 0 
    if false
        check_ch_flag = 0;
    else
        check_ch_flag = 1;
        check_all_vertices = 1;
        %[Current_A_min,Current_A_max] = LL_Jac_Vertice(x_relative_min,x_relative_max);
        %Current_Vertex = Current_A_min;
        for cnt_ver = 1:2^4
%             for cnt_inx = 1:8
%                 if INDEX_CH(cnt_ver,cnt_inx) == 1
%                     Current_Vertex(Jacbian_Function_index(cnt_inx,1),Jacbian_Function_index(cnt_inx,2)) = Current_A_max(Jacbian_Function_index(cnt_inx,1),Jacbian_Function_index(cnt_inx,2));
%                 end
%             end
            Current_Vertex = biology_Vertical(x_relative_min,x_relative_max,INDEX_CH(cnt_ver,:));
            if max(eig(Current_Vertex'*M+M*Current_Vertex - shrink*M)) > 0
                check_all_vertices = 0;
                if enable_disp
                    disp('fail using old Q and gamma')
                end
                break
            end
        end
    end  
    
    if check_ch_flag == 1 && check_all_vertices==1
        Q = M;
        if check_all_vertices == 1
            if enable_disp
                disp('Pass all vertex check for old shrink')    
            end
            gamma = shrink;
%         else
%             gamma = shrink;
%             Delta = abs(LL_Jac(x_relative_min)-LL_Jac(x_relative_max));
%             SDP_Delta_matrix = Delta'*Q+Q*Delta;
%             epsilon_matrix = sqrt(norm(SDP_Delta_matrix,1)*norm(SDP_Delta_matrix,Inf));
%             disp(epsilon_matrix);
%             disp('Estimation error!')
%             disp(gamma*min(eig(Q)) + epsilon_matrix);
%             gamma = (gamma*min(eig(Q)) + epsilon_matrix)/min(eig(Q))/2;
%             shrink = gamma;
        end
        level_next = level*exp(shrink*(t(i+1)-t(i)));
        if enable_disp
            disp('remain P in last round');
            disp('expontial change rate:');
            disp(shrink);
        end
        P_level_original = level;
    else
        if enable_disp
            disp('Compute new P')
        end
        %Lyapunov-like function
        [gamma,Q] = QLF(A_c);
        check_rest_vertices = 1;
        %[Current_A_min,Current_A_max] = LL_Jac_Vertice(x_relative_min,x_relative_max);
        %Current_Vertex = Current_A_min;
        for cnt_ver = 1:2^4
%             for cnt_inx = 1:8
%                 if INDEX_CH(cnt_ver,cnt_inx) == 1
%                     Current_Vertex(Jacbian_Function_index(cnt_inx,1),Jacbian_Function_index(cnt_inx,2)) = Current_A_max(Jacbian_Function_index(cnt_inx,1),Jacbian_Function_index(cnt_inx,2));
%                 end
%             end
            Current_Vertex = biology_Vertical(x_relative_min,x_relative_max,INDEX_CH(cnt_ver,:));
            while max(eig(Current_Vertex'*Q+Q*Current_Vertex - gamma*Q)) > 0
                if gamma <=0
                    gamma = gamma + 0.1;
                else
                    if enable_disp
                        disp('fail using new computed Q and gamma for rest check');
                    end
                    check_rest_vertices = 0;
                    break;
                end
            end
            if check_rest_vertices == 0
                break
            end
        end
       %disp(gamma);
        if check_rest_vertices == 1
            if enable_disp
                disp('pass the rest vertex check')
                disp(gamma);
            end
        else
            [gamma,Q] = QLF(A_c);
            delta_theta = (abs(x(i+1,:)-x(i,:))/2+sqrt(max(eig(inv(M/level)))))*exp(Lip*(t(i+1)-t(i)));
            Delta_matrix = abs(biology_jac(Average_value+delta_theta)-biology_jac(Average_value-delta_theta))/2;
            %disp((Delta_matrix));
            %SDP_Delta_matrix = Delta_matrix'*Q+Q*Delta_matrix;
            SDP_Delta_matrix = Delta_matrix;
            epsilon_matrix = sqrt(norm(SDP_Delta_matrix,1)*norm(SDP_Delta_matrix,Inf));
            disp('Esitimation Error!')
            disp(gamma*min(eig(Q)) + epsilon_matrix);
            gamma = (gamma*min(eig(Q)) + epsilon_matrix)/min(eig(Q))/2;
            if enable_disp
                disp('using interval matrix disturbance to over-approximate \gamma')
                disp('expontial change rate:');
                disp(gamma);
            end
        end
    %     disp('eigenvalue of P');disp(eig(Q));
        P_level_original = min_level_set7(Q,M,level);        
    end
    level_next = P_level_original*(exp(gamma*(t(i+1)-t(i))));
%     for cnt = index(i):index(i)+inner_step-1
%         %plot_ellipsoid(x_precise(cnt,[5,7])',Q([5,7],[5,7])/(P_level_original*(exp(-gamma*(t_precise(cnt)-t(i))))),'k');
%         fp_checker_inner = fix_point_checker(x_precise(cnt+1,:),x_precise(cnt,:),Q/(P_level_original*(exp(gamma*(t_precise(cnt+1)-t(i))))),Q/(P_level_original*(exp(gamma*(t_precise(cnt)-t(i))))));
%         if fp_checker_inner == 1
%             if enable_disp
%                 disp('find invariant set in between');
%             end
%             break;
%         end
%     end
    if enable_plot
        plot_ellipsoid(x(i+1,[2,4])',Q([2,4],[2,4])/level_next,'k');
    end
    Reach_info(i+1,:) = [reshape(Q,1,49),level_next];
%     fp_checker = fix_point_checker(x(i+1,:),x(i,:),Q/level_next,M/level);
%     if fp_checker == 1 || fp_checker_inner == 1
%         disp('find invariant set');
%         break;
%     end
end

toc;

Reach_final = reshape(Reach_info(i,1:49),7,7)/Reach_info(i,50);
dia_final = sqrt(eig(inv(Reach_final)));
disp(dia_final);

% %%
% n = length(x);   
% %dimention of the system
% 
% Reach_info = zeros(n,50);
% Reach_info(1,:) = [reshape(Initial_P,1,49),Initial_l];
% 
% for i = 1:length(t)-1
% %for i = 1:20
%     %
%     if enable_disp
%         disp('----------------');
%         disp(['Computing ', num2str(t(i)), ' time interval ']);
%     end
%     M = reshape(Reach_info(i,1:49),7,7);
%     level = Reach_info(i,50);
% 
%     Average_value = mean(x_precise(index(i):index(i+1)-1,:));
%     A_c = biology_jac(Average_value);
%     shrink = max(real(eig(A_c)))+5;
%     %shrink = -0.01;
%     
%     %Local Lipschitz constant
%     Lip = norm(A_c);   
%    %delta_theta = (abs(x(i+1,:)-x(i,:))/2+sqrt(max(eig(inv(M/level)))))*exp(Lip*(t(i+1)-t(i)));
%     x_relative_max = (max(x(i+1,:),x(i,:))+sqrt(max(eig(inv(M/level)))))*exp(Lip*(t(i+1)-t(i)));
%     x_relative_min = (min(x(i+1,:),x(i,:))-sqrt(max(eig(inv(M/level)))))*exp(Lip*(t(i+1)-t(i)));
%     if max(eig(A_c'*M+M*A_c-shrink*M)) > 0 
%         check_ch_flag = 0;
%     else
%         check_ch_flag = 1;
%     end   
%     if check_ch_flag == 1
%         Q = M;
%         gamma = shrink;
%         Delta = biology_delta(x_relative_min,x_relative_max);
%         SDP_Delta_matrix = Delta'*Q+Q*Delta;
%         epsilon_matrix = sqrt(norm(SDP_Delta_matrix,1)*norm(SDP_Delta_matrix,Inf));
%         gamma = (gamma*min(eig(Q)) + epsilon_matrix)/min(eig(Q));
%         shrink = gamma;
%         level_next = level*exp(shrink*(t(i+1)-t(i)));
%         if enable_disp
%             disp('remain P in last round');
%             disp('expontial change rate:');
%             disp(shrink);
%         end
%         P_level_original = level;
%     else
%         if enable_disp
%             disp('Compute new P')
%         end
%         %Lyapunov-like function
%         [gamma,Q] = QLF(A_c);
%        %disp(gamma);
%        check_new_p_flag = 0;
%         if check_new_p_flag == 0
%             delta_theta = (abs(x(i+1,:)-x(i,:))+sqrt(max(eig(inv(M/level)))))*exp(Lip*(t(i+1)-t(i)));
%             Delta_matrix = biology_delta(Average_value+delta_theta,Average_value-delta_theta)/2;
%             %disp(norm(Delta_matrix))
%             SDP_Delta_matrix = Delta_matrix'*Q+Q*Delta_matrix;
%             epsilon_matrix = sqrt(norm(SDP_Delta_matrix,1)*norm(SDP_Delta_matrix,Inf));
%             disp(gamma*max(eig(Q)) + epsilon_matrix);
%             gamma = (gamma*min(eig(Q)) + epsilon_matrix)/min(eig(Q));
%             if enable_disp
%                 disp('using interval matrix disturbance to over-approximate \gamma')
%                 disp('expontial change rate:');
%                 disp(gamma);
%             end
%         else
%             if enable_disp
%                 disp('Computed P satifies convex hull constrains')
%                 disp('expontial change rate:');
%                 disp(gamma);
%             end
%         end
%     %     disp('eigenvalue of P');disp(eig(Q));
%         P_level_original = min_level_set7(Q,M,level);
%         level_next = P_level_original*(exp(gamma*(t(i+1)-t(i))));
%     end
%     for cnt = index(i):index(i)+inner_step-1
%         %plot_ellipsoid(x_precise(cnt,[2,4])',Q([2,4],[2,4])/(P_level_original*(exp(-gamma*(t_precise(cnt)-t(i))))),'k');
%         fp_checker_inner = fix_point_checker(x_precise(cnt+1,:),x_precise(cnt,:),Q/(P_level_original*(exp(gamma*(t_precise(cnt+1)-t(i))))),Q/(P_level_original*(exp(gamma*(t_precise(cnt)-t(i))))));
%         if fp_checker_inner == 1
%             if enable_disp
%                 disp('find invariant set in between');
%             end
%             break;
%         end
%     end
%     plot_ellipsoid(x(i+1,[2,4])',Q([2,4],[2,4])/level_next,'k');
%     Reach_info(i+1,:) = [reshape(Q,1,49),level_next];
%     fp_checker = fix_point_checker(x(i+1,:),x(i,:),Q/level_next,M/level);
%     if fp_checker == 1 || fp_checker_inner == 1
%         disp('find invariant set');
%         break;
%     end
% end
% 
% toc;
% 
% Reach_final = reshape(Reach_info(i,1:49),7,7)/Reach_info(i,50);
% dia_final = sqrt(eig(inv(Reach_final)));
% disp(dia_final);

