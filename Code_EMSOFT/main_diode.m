

% a = 0.35;
% b = 0.2;
% 
% critial_of_d = roots([4*5*803.712 -3*4*1086.288 2*3*551.088 -2*124.548]);
% 
% x = -0.01:0.01:0.6;
% y = 5*803.712*x.^4-4*1086.288*x.^3+3*551.088*x.^2-2*124.548*x+10.656;
% [V1,V2] = diode_Vertical(a,b,critial_of_d);
% plot(x,y);
% hold on;
% scatter(a,-V1(1,1));
% scatter(b,-V2(1,1));

clc;
clear;
close;
%
% addpath(genpath('..\SeDuMi_1_3'));
% addpath(genpath('..\yalmip'));

%%
%initials
enable_disp = 0;

tic;
%% set all the initial values
dt = 5e-3;
t_dur = 0:dt:9;
initial_value = [0.475, 0.1];
dim = 2; %dimension of the system
%time_span = [0,10];
options = odeset('AbsTol',1e-1,'MaxStep',0.1);
ops_sdp = sdpsettings('verbose',0);

critial_of_d = roots([4*5*803.712 -3*4*1086.288 2*3*551.088 -2*124.548]);
delta_gamma = 0.1;

[t_precise,x_precise] = ode45(@diode,t_dur,initial_value,options);

Initial_P = eye(2);
Initial_l = (0.5e-2)^2;
stop_flag = 0;

if enable_disp
    figure;
    plot(x_precise(:,1),x_precise(:,2));
    hold on;
    plot_ellipsoid(x_precise(1,:)',Initial_P/Initial_l,'r');
end
step = 10;
index = 1:step:length(t_precise);
t = t_precise(index);
x = x_precise(index,:);

%%
n = length(x);   
%dimention of the system

Reach_info = zeros(n,5);
Reach_info(1,:) = [Initial_P(1,:),Initial_P(2,:),Initial_l];

for i = 1:length(t)-1
%for i = 1:20
    %
    if enable_disp
        disp('----------------');
        disp(['Computing ', num2str(i), ' th time interval ']);
    end
    M = [Reach_info(i,1:2);Reach_info(i,3:4)];
    level = Reach_info(i,5);


    
    %Local Lipschitz constant
    Lip = norm(diode_Jacobian(x(i,:)));   
   %delta_theta = (abs(x(i+1,:)-x(i,:))/2+sqrt(max(eig(inv(M/level)))))*exp(Lip*(t(i+1)-t(i)));
    x_relative_max = (max(x(i+1,:),x(i,:))+sqrt(max(eig(inv(M/level)))))*exp(Lip*(t(i+1)-t(i)));
    x_relative_min = (min(x(i+1,:),x(i,:))-sqrt(max(eig(inv(M/level)))))*exp(Lip*(t(i+1)-t(i)));
    [A1,A2] = diode_Vertical(x_relative_max,x_relative_min,critial_of_d);
    gamma = max(real([eig(A1);eig(A2)])) + 0.1;
    check_ch_flag = 1;
    if max(eig(A1'*M + M*A1 - gamma*M))>0 ||  max(eig(A2'*M + M*A2 - gamma*M))>0
        check_ch_flag = 0;
    end
    if check_ch_flag == 1
        Q = M;
        P_level_original = level;
        level_next = level*exp(gamma*(t(i+1)-t(i)));
        if enable_disp
            disp('remain P in last round');
            disp('expontial change rate:');
            disp(gamma);
        end
    else
        gamma = min(real([eig(A1);eig(A2)]));
        if enable_disp
            disp('Compute new P')
        end
        while(~stop_flag)
            P = sdpvar(2,2);
            F = [P>=eye(2), A1'*P+P*A1 - gamma*P <=0, A2'*P+P*A2 - gamma*P <=0];
            diagnostics=optimize(F,[],ops_sdp);
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
        if enable_disp
            disp('Exponential change rate:')
            disp(gamma);
        end
    %     disp('eigenvalue of P');disp(eig(Q));
        P_level_original = min_level_set(Q,M,level);
        level_next = P_level_original*(exp(gamma*(t(i+1)-t(i))));

    end
    for cnt = index(i):index(i)+step-1
        if enable_disp
            plot_ellipsoid(x_precise(cnt,:)',Q/(P_level_original*(exp(gamma*(t_precise(cnt)-t(i))))),'k');
        end
    end
    if enable_disp
        plot_ellipsoid(x(i+1,:)',Q/level_next,'k');
    end
    Reach_info(i+1,:) = [Q(1,:),Q(2,:),level_next];
end
toc;
Reach_final = reshape(Reach_info(i,1:4),2,2)/Reach_info(i,5);
dia_final = eig(inv(Reach_final));
disp(dia_final);


%To do: fix point check