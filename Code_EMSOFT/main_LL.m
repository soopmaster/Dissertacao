clear;clc;
%%
%set up the initial conditions: total time; time step(optional); initial
%states; perturbation of the initial state;
% addpath(genpath('C:\Users\cfan10\Google Drive\LMI\SeDuMi_1_3'));
% addpath(genpath('C:\Users\cfan10\Google Drive\LMI\yalmip'));
% addpath(genpath('C:\Users\cfan10\Google Drive\LMI\SDPT3'))
%%
ops = sdpsettings('verbose',0);
enable_disp = 0;
enable_plot = 0;

dt = 0.001;
t_dur = 0:dt:1;
inner_step = 50;
initial_state = [1.2, 1.105, 1.5, 2.4, 1, 0.1, 0.45]';
%initial_state = [0.8691    0.3680    0.5589    0.7537    0.2209    0.0841    0.2743]';
%initial_state = 1*ones(7,1);
dimension = 7;

load Initial_P;
Initial_P = reshape(Initial_P,7,7);
%Initial_P = eye(7);
%Initial_P(1,1) = (0.01);
Initial_l = (0.025)^2;

[U,V] = eig(inv(Initial_P/(Initial_l)));
%delta_initial_up = [sqrt(Initial_l),0,0,0,0,0,0]';
delta_initial_up = U(:,7)*sqrt(V(7,7))*5;
%Lipschitz constant of te system, user input
Lip_constant = 0;
[t_precise,x_precise] = ode45(@LL,t_dur,initial_state);
[t_true,x_true] = ode45(@LL,t_dur,initial_state+delta_initial_up);
%[t_true_low,x_true_low] = ode45(@Engine,t_dur,initial_state+delta_initial_low);

%% Get index of variables
INDEX_CH = zeros(2^5,5);
for i = 0 : 2^5-1
    INDEX_CH(i+1,:) = de2bi(i,5);
end
%Jacbian_Function_index = [3,2;3,3;4,3;4,4;5,4;5,5;7,2;7,7];

%%
if enable_plot
    figure;
    scatter(x_precise(:,1),x_precise(:,2));
    xlim = ylim;
    hold on;
    scatter(x_true(:,1),x_true(:,2),'r');
    scatter(initial_state(1),initial_state(2),180,'x','b','linewidth',2);
    scatter(initial_state(1)+delta_initial_up(1),initial_state(2)+delta_initial_up(2),180,'x','r','linewidth',2);
    xlabel('x1','FontSize', 24);
    ylabel('x2','FontSize', 24);
    plot_ellipsoid(x_precise(1,[1,2])',Initial_P([1,2],[1,2])/Initial_l,'k');
end



index = 1:inner_step:length(t_precise);
t = t_precise(index);
x = x_precise(index,:);

% figure
% for i = 1:size(t_precise,1)
%     lambda = max(real(eig(LL_Jac(x_precise(i,:)))));
%     tmp(i) = lambda;
% end
% plot(tmp)

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
    A_c = LL_Jac(Average_value);
    shrink = max(real(eig(A_c)));
    %shrink = -0.01;
    
    %Local Lipschitz constant
    Lip = norm(A_c);   
   %delta_theta = (abs(x(i+1,:)-x(i,:))/2+sqrt(max(eig(inv(M/level)))))*exp(Lip*(t(i+1)-t(i)));
    x_relative_max = (max(x(i+1,:),x(i,:))+sqrt(max(eig(inv(M/level)))))*exp(Lip*(t(i+1)-t(i)));
    x_relative_min = (min(x(i+1,:),x(i,:))-sqrt(max(eig(inv(M/level)))))*exp(Lip*(t(i+1)-t(i)));
    if max(eig(A_c'*M+M*A_c-shrink*M)) > 0 
        check_ch_flag = 0;
    else
        check_ch_flag = 1;
        check_all_vertices = 1;
        %[Current_A_min,Current_A_max] = LL_Jac_Vertice(x_relative_min,x_relative_max);
        %Current_Vertex = Current_A_min;
        for cnt_ver = 1:2^5
%             for cnt_inx = 1:8
%                 if INDEX_CH(cnt_ver,cnt_inx) == 1
%                     Current_Vertex(Jacbian_Function_index(cnt_inx,1),Jacbian_Function_index(cnt_inx,2)) = Current_A_max(Jacbian_Function_index(cnt_inx,1),Jacbian_Function_index(cnt_inx,2));
%                 end
%             end
            Current_Vertex = LL_Jac_Vertice(x_relative_min,x_relative_max,INDEX_CH(cnt_ver,:));
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
        for cnt_ver = 1:2^5
%             for cnt_inx = 1:8
%                 if INDEX_CH(cnt_ver,cnt_inx) == 1
%                     Current_Vertex(Jacbian_Function_index(cnt_inx,1),Jacbian_Function_index(cnt_inx,2)) = Current_A_max(Jacbian_Function_index(cnt_inx,1),Jacbian_Function_index(cnt_inx,2));
%                 end
%             end
            Current_Vertex = LL_Jac_Vertice(x_relative_min,x_relative_max,INDEX_CH(cnt_ver,:));
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
            %[gamma,Q] = QLF(A_c);
            delta_theta = (abs(x(i+1,:)-x(i,:))/2+sqrt(max(eig(inv(M/level)))))*exp(Lip*(t(i+1)-t(i)));
            Delta_matrix = abs(LL_Jac(Average_value+delta_theta)-LL_Jac(Average_value-delta_theta))/2;
            %disp((Delta_matrix));
            %SDP_Delta_matrix = Delta_matrix'*Q+Q*Delta_matrix;
            SDP_Delta_matrix = Delta_matrix;
            epsilon_matrix = sqrt(norm(SDP_Delta_matrix,1)*norm(SDP_Delta_matrix,Inf));
            if enable_disp
                disp('Esitimation Error!')
                disp(gamma*min(eig(Q)) + epsilon_matrix);
            end
            gamma = (gamma*min(eig(Q)) + epsilon_matrix)/min(eig(Q))/2;
            if enable_disp
                disp('using interval matrix disturbance to over-approximate \gamma')
                disp('expontial change rate:');
                disp(gamma);
            end
        end
        % Change the shape of ellipsoid
        %P_level_original is the changed ellipsoid
        P_level_original = min_level_set7(Q,M,level);
        if enable_plot
            plot_ellipsoid(x(i,[1,2])',M([1,2],[1,2])/level,'b');
            plot_ellipsoid(x(i,[1,2])',Q([1,2],[1,2])/P_level_original,'r');
        end               
    end
    disp(i)
    disp(gamma)
    level_next = P_level_original*(exp(gamma*(t(i+1)-t(i))));
    for cnt = index(i):index(i)+inner_step-1
%         disp(exp(gamma*(t_precise(cnt)-t(i))));
        if enable_plot
            plot_ellipsoid(x_precise(cnt,[1,2])',Q([1,2],[1,2])/(P_level_original*(exp(gamma*(t_precise(cnt)-t(i))))),'r');
        end
%         fp_checker_inner = fix_point_checker(x_precise(cnt+1,:),x_precise(cnt,:),Q/(P_level_original*(exp(gamma*(t_precise(cnt+1)-t(i))))),Q/(P_level_original*(exp(gamma*(t_precise(cnt)-t(i))))));
%         if fp_checker_inner == 1
%             if enable_disp
%                 disp('find invariant set in between');
%             end
%             break;
%         end
    end
    if enable_plot
        plot_ellipsoid(x(i+1,[1,2])',Q([1,2],[1,2])/level_next,'k');
    end
    Reach_info(i+1,:) = [reshape(Q,1,49),level_next];
%     fp_checker = fix_point_checker(x(i+1,:),x(i,:),Q/level_next,M/level);
%     if fp_checker == 1 || fp_checker_inner == 1
%         disp('find invariant set');
%         break;
%     end
end

toc;
Reach_initial = reshape(Reach_info(1,1:49),7,7)/Reach_info(1,50);
dia_initial = sqrt(eig(inv(Reach_initial)));
Reach_final = reshape(Reach_info(i,1:49),7,7)/Reach_info(i,50);
dia_final = sqrt(eig(inv(Reach_final)));
disp(prod(dia_final)/prod(dia_initial));

disp('Flow*');
final_delta_for_flow = [1.436292-1.353641, 1.171099-1.100125, 0.725349-0.652628, 1.865653-1.759446, 0.644066-0.590763, 0.140873-0.131672, 0.202599-0.171660];
initial_delta_for_flow = 0.05*ones(1,7);
disp(prod(final_delta_for_flow)/prod(initial_delta_for_flow));