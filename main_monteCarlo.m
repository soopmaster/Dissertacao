% Main_monteCarlo file of thesis project
% Runs all generalized functions for 1 vehichle
% Offline, initializes cell vector of vehiStruct
% Then initializes parametope propagation variables
% Online, runs specific vehichle and parametope dynamics 
% Every inner_step runs predicted trajectory and parametopes
% Initialize vehi_vec, containing vehiSt (model) and ctrlSt (controller)

% This version simulates just one vehichle and parametopes. 
% Does monte-carlo simulation on randomly sampled points in ellipsoid
% boundary to validate parametope propagation.

close all;
clearvars;
set(groot,'DefaultTextInterpreter','latex');

% vehiSimStruct : struct (struct with components of main sim)
%   .vehiSt : Struct with functions and variables of model
%   .ctrlSt : Struct with functions and variables of controller
%   .x
%   .u
%   .alph
%   .y
%   .traj
%   .pathPlot
%   .ops

% === Sim Settings === %
simT = 3; % Sim time seconds
Ts = 0.01;
knum = simT/Ts;
%calc_traj_dt = 0.1; % How often trajectory is recalculated [s]
calc_traj_dt = Ts;
recalc_step = calc_traj_dt/Ts;
pred_hor = 50;
discnum = 5; % It will actually be discnum+1
disc_step = floor(pred_hor/discnum);
disc_dt = disc_step*Ts; % Time step between Reachability discs
prmtp_ind = [1,disc_step:disc_step:pred_hor];
discnum = length(prmtp_ind);

enable_monteCarlo = 1;
enable_shadow = 1;
enable_disp = 0;
make_vid = 0;

sim.tI = 2;
sim.pA0 = [20;-20;6];
sim.pAf = [0;0;0];
% critv_ref = [-5;5;-5];
sim.jAf = [0;0;0];
% sim.aAf = [0;0;0];
dir = [0;0;0];
sim.aAf = dir;
sim.vAf = dir;
% sim.aAf = 10*dir/norm(dir);
% sim.vAf = sim.tI*sim.aAf;

% === Init vehiStruct === &
vehi_vec = init_all_vehi_local(sim,Ts,pred_hor); % Cell vec with vehiSimStruct
numVehi = length(vehi_vec);
% vehi_vec{1} is assumed to be the Invader 

for ind=1:numVehi
    vehi = vehi_vec{ind};
    nx = vehi.vehiSt.nx;
    nx_full = nx+nx^2;

    if ind == 1
        % === Define Critical Area === %
        critp_ref = [0;0;0];
        % critv_ref = [-5;5;-5];
        critv_ref = [0;0;0];
        pts_ref = [critp_ref; critv_ref; vehi.vehiSt.x0(3*2+1:end)];
        vehi.pts_ref{1} = [vehi.vehiSt.x0, repmat(pts_ref,1,pred_hor)];
    else
        pts_ref = repmat(vehi_vec{1}.vehiSt.x0,1,pred_hor+1);
        vehi.pts_ref{1} = pts_ref;
        %vehi.pts_ref = repmat(vehi.vehiSt.x0,1,pred_hor+1);
    end

    % === Init Simulation Vectors === %
    vehi.x = zeros(nx,knum);
    vehi.x(:,1) = vehi.vehiSt.x0;
    vehi.u = zeros(vehi.vehiSt.nu,knum);
    vehi.u(:,1) = vehi.vehiSt.u0;

    % vehi.alph = zeros(vehi.vehiSt.nx^2,knum/recalc_step);
    % vehi.alph(:,1) = reshape(vehi.vehiSt.alph0,[],1);
    % 
    vehi.xpred = cell(0,0);
    vehi.upred = cell(0,0);
    
    %vehi.error_track = [];

    % Define window bounds

    vehi_vec{ind} = vehi;
end

vehiAux = vehi_vec{1}; % Auxilliary vehi without plots for last prediction
vehiAux.pts_ref{1} = [vehi.vehiSt.x0, repmat(pts_ref,1,knum)];
vehiAux.ops.plot.edgecolor = 'r'; vehiAux.ops.plot.wirestyle = '-'; vehiAux.ops.plot.linewidth = 0.5;

% MSE vec for parametope. Columns are pred step, Rows sim step
MSE = zeros(pred_hor-1,0); 
x_branch = {};

textratio = (21.5-6)/21.5; % text to page size ratio in latex
textratio = 0.6;

% === Start Figure === %
fig=figure('WindowState','maximized');
% ax3D = subplot(1,2,1);
% axMSE = subplot(1,2,2);
% axSide = subplot(2,2,4);
% axes(ax3D)
hold on;
grid on
view(3)
[title_3D,subtitle_3D] = title(['Parametope Validation Sim (Monte-Carlo), N=',num2str(pred_hor)], ...
    'k = 0','Interpreter','latex');
xlabel('$x$','Interpreter','latex');
ylabel('$y$','Interpreter','latex');
zlabel('$z$','Interpreter','latex');
for ind=1:numVehi
    pts_ref = vehi_vec{ind}.pts_ref{1};
    vehi_vec{ind}.prmtpPatch = {};
    vehi_vec{ind}.monteCarloPlot = {};
    
    vehi_vec{ind}.wayptsPlot = plot3(pts_ref(1,:),pts_ref(2,:),pts_ref(3,:),'xr','DisplayName','target');
    %vehi_vec{ind}.trajPlot = plot3(vehi_vec{ind}.traj(1,:), vehi_vec{ind}.traj(2,:), vehi_vec{ind}.traj(3,:) ,'-.k','DisplayName','r_T');
    vehi_vec{ind}.bodyPlot = vehi_vec{ind}.vehiSt.bodyPlot();
    %vehi_vec{ind}.pathPlot = animatedline(vehi_vec{ind}.x(1,1), vehi_vec{ind}.x(2,1), vehi_vec{ind}.x(3,1),'LineWidth',2.5, 'Color','b','LineStyle','-','DisplayName','Path');
end
vehiAux.prmtpShadow = {}; % Parametope Shadow patch

axis equal
%axis manual
% ac = get(ax3D, 'Children');   %no hidden objects
% copyobj( ac, axTop );
% copyobj( ac, axSide );
%%
% === Simulation === %
telapse_total = tic;
k=1;
% norm(vehi_vec{1}.x(1:3,k)-critp_ref)>=0.1
stop = 0;
framek = 1;
while k<knum && ~stop
    % Update subtitle
    subtitle_3D.String = ['k = ',num2str(k)];
    disp(['=== Computing ', num2str(k), '/', num2str(knum), ' interval ===']);
    disp(['t = ', num2str((k-1)*Ts),'s']);
    telapse_k = tic;
    for kVehi=1:numVehi
        vehi = vehi_vec{kVehi};
        vehiSt = vehi.vehiSt;
        nx = vehiSt.nx;
        nx_full = nx+nx^2;
        ops = vehi.ops;

        % === Real-Time === %
        % Add disturbances!
        [xpred_aux, upred_aux, auxStaux] = calcNxt(vehi.x(:,k), vehi.pts_ref{end,1}, vehi, Ts, pred_hor);
        vehi.xpred{end+1,1} = xpred_aux(:,1:pred_hor);
        vehi.upred{end+1,1} = upred_aux(:,1:pred_hor);
        vehi.ctrlSt.auxSt = auxStaux;

        vehi.x(:,k+1) = vehi.xpred{end,1}(1:nx,2);
        vehi.u(:,k) = vehi.upred{end,1}(:,1);
        %vehi.error_track(:,k) = vehi.ctrlSt.auxSt.error_track;
        stop = norm(vehi_vec{1}.x(1:3,k)-critp_ref)<=0.1;

        if enable_shadow
        % Compute Final parametope
            [xpred_aux, upred_aux, auxStaux] = calcNxt(vehiAux.x(:,k), vehiAux.pts_ref{end,1}, vehiAux, Ts, knum);
            vehiAux.xpred{end+1,1} = xpred_aux(:,1:knum);
            vehiAux.upred{end+1,1} = upred_aux(:,1:knum);
            vehiAux.x(:,k+1) = vehiAux.xpred{end,1}(1:nx,2);
            vehiAux.u(:,k) = vehiAux.upred{end,1}(:,1);
            %vehiAux.error_track(:,k) = vehiAux.ctrlSt.auxSt.error_track;
            vehiAux.ctrlSt.auxSt = auxStaux;

             % Find prediction closest to ground
            col = [];
            tol = 1e-3;
            last_xpred = vehiAux.xpred{end,1};
            if ~isempty(last_xpred)
                while isempty(col)
                    [~,col] = find(sqrt(sum(last_xpred(3,:).^2,1))<= tol);
                    tol = tol*10;
                end
                grd_ind = col(1);
            else
                grd_ind = 0;
            end
        end

        %[alph_vec, vehi.ypred{end+1,1}] = calcParametopes(vehi.xpred{end,1}, vehi.upred{end,1}, vehi, Ts);

        if enable_monteCarlo
            % [xVal] = prmtpMonteCarloVal_surf(vehi.xpred{end,1}(1:nx,:),alph_vec,vehi.ypred{end,1},10,Ts,vehiSt);
            [x_branch{end+1,1},MSE(:,k)] = prmtpMonteCarloVal_tree(vehi.xpred{end,1}(1:nx,:),vehi.xpred{end,1}(nx+(1:nx^2),:),vehi.xpred{end,1}(end,:),1,Ts,vehiSt);
        end

        % === Update Plot === %
        plotBodyFunc = vehiSt.plotBodyFunc; 
        vehi.bodyPlot = plotBodyFunc(vehi.x(:,k), vehi.bodyPlot);
        %addpoints(vehi.pathPlot,vehi.x(1,k), vehi.x(2,k), vehi.x(3,k));
        if enable_shadow
                vehiAux.prmtpShadow = plotParametope(vehiAux.xpred{end,1}(end,grd_ind),vehiAux.xpred{end,1}(nx+1:end-1,grd_ind),vehiAux.xpred{end,1}(1:nx,grd_ind),nx,vehiAux.ops,vehiAux.prmtpShadow,2);
        end
        if enable_disp
            vehi.prmtpPatch = plotParametope(vehi.xpred{end,1}(end,prmtp_ind),vehi.xpred{end,1}(nx+1:end-1,prmtp_ind),vehi.xpred{end,1}(1:nx,prmtp_ind),nx,ops,vehi.prmtpPatch,3);
            if enable_monteCarlo
                %x_leaf_ind = x_leaf(prmtp_ind,1);
                %x_leaf_ind = x_leaf;
                vehi.monteCarloPlot = plotMonteCarlo(x_branch{end,1},prmtp_ind,vehi.monteCarloPlot);
            end
                %vehi.prmtpPatch = plotParametope(vehiSt.ysim0.*ones(1,length(prmtp_ind)),vehi.xpred(nx+1:end,prmtp_ind),vehi.xpred(1:nx,prmtp_ind),nx,ops,vehi.prmtpPatch);
            %vehi.prmtpPatch = plotParametope(vehi.ypred,vehi.alphpred,vehi.xpred(:,prmtp_ind),nx,ops,vehi.prmtpPatch);
            if 0&&k==0
                axis equal
                axis manual
            end
        end

        vehi_vec{kVehi} = vehi;
        
        %legend({'target', 'r_T', 'Path'});
    end
    disp(['Computing time for interval: ', num2str(toc(telapse_k)), 's']);
    disp('----------------');
    disp(' ');
    fontsize(11*sqrt(2)/0.9/textratio,"points");
    drawnow;
    pause(Ts);
    if make_vid == 1
        frame(k) = getframe(fig);
    end
    if k==1 || k==90 ||k==117
        exportgraphics(fig,['reachable_grd_sim',num2str(framek),'.pdf']);
        framek = framek+1;
    end
    k = k+1;
end
hold off;
disp(['Total sim computing time: ', num2str(toc(telapse_total)), 's']);

%%
if make_vid == 1
    vid=VideoWriter('test_monteCarlo.mp4','MPEG-4');
    vid.FrameRate=1/Ts;
    open(vid);
    writeVideo(vid,frame);
    close(vid);
end

if enable_monteCarlo
    % Plot MSE 
    fig_err = figure('WindowState','maximized');
    [title_MSE,subtitle_MSE] = title('Parametope Validation MSE', ...
        'Norm of distance between Monte-Carlo predictions and Ellipsoid surface','Interpreter','latex');
    grid on;
    hold on;
    xlabel('$k_{sim}$','Interpreter','latex');
    ylabel('$k_{pred}$','Interpreter','latex');
    zlabel('$MSE$','Interpreter','latex');
    surf(MSE);
    colorbar;
    fontsize(11*sqrt(2)/0.9/textratio,"points");
    hold off;
    exportgraphics(fig_err,['monteCarloMSE1_N',num2str(pred_hor),'.pdf']);
    view(3);
    exportgraphics(fig_err,['monteCarloMSE2_N',num2str(pred_hor),'.pdf']);

    % Plot 3 frames of Sim
    fig_sim = figure;
    hold on;
    grid on
    view(3)
    frame_ind = 1:floor((k-1)/4):(k-1);
    [title_3frame,subtitle_3frame] = title(['Parametope Validation Sim (Monte-Carlo), N=',num2str(pred_hor)], ...
        ['k = ',num2str(frame_ind)],'Interpreter','latex');
    xlabel('$x$','Interpreter','latex');
    ylabel('$y$','Interpreter','latex');
    zlabel('$z$','Interpreter','latex');
    axis equal;
    for ind=frame_ind
        vehi.vehiSt.plotBodyFunc(vehi.x(:,ind), vehi.vehiSt.bodyPlot());
        plotParametope(vehi.xpred{ind,1}(end,prmtp_ind),vehi.xpred{ind,1}(nx+1:end-1,prmtp_ind),vehi.xpred{ind,1}(1:nx,prmtp_ind),nx,ops,{},3);
        plotMonteCarlo(x_branch{ind,1},prmtp_ind,{});
    end
    fontsize(11*1/0.9,"points");
    exportgraphics(fig_sim,['monteCarloSim_N',num2str(pred_hor),'.pdf']);
    % exportgraphics(fig_sim,['reachable_grd_sim',num2str(framek),'.pdf']);
    
end

function vehi_vec = init_all_vehi_local(sim,Ts,pred_hor)
%INIT_ALL_VEHI outputs cell vector of vehiStruct
%   Detailed explanation goes here

vehi_vec = cell(1,1);

% Define Inital Positions
p10 = [5;-5;2];
x0 = {p10};
tI = 2;

vehi_vec{1}.vehiSt = init_QuadDroneNMPC(x0{1},Ts);

% vehi_vec{1}.ctrlSt = init_MellingerCtrl(pred_hor,2,vehi_vec{1}.vehiSt);
vehi_vec{1}.ctrlSt = init_MellingerCtrl(sim.vAf,sim.aAf,sim.jAf,pred_hor,tI,vehi_vec{1}.vehiSt);
% vehi_vec{1}.ctrlSt = init_MellingerCtrl_v3(pred_hor,Ts,vehi_vec{1}.vehiSt);

% === Define ops === %
ops = sdpsettings('verbose',0);
ops.plot.shade = 0.15;
ops.plot.linewidth=0.05;
ops.plot.wirestyle=':';
ops.plot.edgecolor='#EB4B2D';
%ops.plot.edgecolor=[235 75 45];
vehi_vec{1}.ops = ops;

% === Init Plot Obs === %
L = 0.5;
p_NB = [L;0;0];
p_SB = [-L;0;0];
p_EB = [0;L;0];
p_OB = [0;-L;0];
p_vec = [p_NB; p_EB; p_SB; p_OB];
vehi_vec{1}.p_vec = p_vec;

end

function [x_branch,MSE] = prmtpMonteCarloVal_tree(x_vec,alph_vec,y_vec,N,Ts,vehiSt)
% N : number of random inputs to propagate per tau cube face
% Monte-Carlo starts here

[nx,hor] = size(x_vec);
nx = vehiSt.nx;
num_in = 3; % Number of inputs combinations to validate
xVal = {};
x_leaf = {};

F_max = vehiSt.vbds_init(end-4+1,2);
F_min = vehiSt.vbds_init(end-4+1,1);
tau_min = vehiSt.vbds_init(end-3+1:end,1);
tau_max = vehiSt.vbds_init(end-3+1:end,2);

ux = linspace(tau_min(1),tau_max(1),num_in);
uy = linspace(tau_min(2),tau_max(2),num_in);
uz = linspace(tau_min(3),tau_max(3),num_in);
[uX,uY,uZ] = meshgrid(ux,uy,uz);

columize = @(A)reshape(permute(A, [3,2,1]), [],1,1);
tauU = [columize(uX), columize(uY), columize(uZ)]';
% Find the column where all values are 0 and remove it
tauUnorm = zeros(1,size(tauU,2));
for k=1:size(tauU,2)
    tauUnorm(k) = norm(tauU(:,k));
end
tauU = tauU(:,(tauUnorm>0));
U_branch_max = [F_max.*ones(1,size(tauU,2)); tauU];
U_branch_min = [F_min.*ones(1,size(tauU,2)); tauU];
U_branch = [U_branch_min, U_branch_max];
num_branch = size(U_branch,2);

% Squared error vec
SE = zeros(1,num_branch);
MSE = zeros(hor-1,1);

x_branch = cell(hor,1);
x_branch{1,1} = repmat(x_vec(1:nx,1),1,num_branch);
% x_leaf = cell(hor,1);
% x_leaf{1,1} = repmat(x_vec(1:nx,1),1,N*6*num_branch);

for k=1:hor-1
    for k_br=1:num_branch
        % 1. Compute points once with randomly selected boundary inputs (leaf)
        % For every face of tau cube, select random value between bounds
        % a = tau_min;
        % b = tau_max;
        % tau_curr = a + (b-a).*rand(3,N*6);
        % for k_cb=1:3
        %     % This only works when N=1 for now! Fix
        %     tau_curr(k_cb,k_cb*2-1) = repmat(tau_min(1),1,N);
        %     tau_curr(k_cb,k_cb*2) = repmat(tau_max(1),1,N);
        % end
        % U_leaf_max = [F_max.*ones(1,N*6); tau_curr];
        % U_leaf_min = [F_min.*ones(1,N*6); tau_curr];
        % U_leaf = [U_leaf_min, U_leaf_max];
        % for k_lf=1:size(U_leaf,2)
        % 
        %     x_leaf{k+1,1}(:,end+1) = runge_kutta4(@(y) vehiSt.partDynFunc(y,U_leaf(:,k_lf)), Ts, x_branch{k,1}(:,k_br));
        % end
        % 2. Propagate corners of input cube (branch)
        x_branch{k+1,1}(:,k_br) = runge_kutta4(@(y) vehiSt.partDynFunc(y,U_branch(:,k_br)), Ts, x_branch{k,1}(:,k_br));
    end
    % Cell index in x__ is samples at time k

    % 3. Calculate error to surface of next ellipsoid
    if k<=hor
        alph = v2m(alph_vec(:,k),nx);
        y = y_vec(k);
        q = x_vec(:,k);
        Q = (inv(alph'*alph)); % Shape Matrix
        [w,W] = ellip2R3(q,Q);
    
        for k_br=1:num_branch
            % Squared error vec
            SE(k_br) = ellipDistPoint(w,y,W,x_branch{k+1,1}(1:3,k_br))^2;
        end
        MSE(k) = (1/num_branch)*sum(SE);
    end
    

end

end

function [xVal] = prmtpMonteCarloVal_surf(x_vec,alph_vec,y_vec,N,Ts,vehiSt)
% Monte-Carlo starts here

[nx,hor] = size(x_vec);
nx = vehiSt.nx;
num_in = 3; % Number of inputs combinations to validate
xVal = {};

F_max = vehiSt.vbds_init(end-4+1,2);
tau_min = vehiSt.vbds_init(end-3+1:end,1);
tau_max = vehiSt.vbds_init(end-3+1:end,2);

ux = linspace(tau_min(1),tau_max(1),num_in);
uy = linspace(tau_min(2),tau_max(2),num_in);
uz = linspace(tau_min(3),tau_max(3),num_in);
[uX,uY,uZ] = meshgrid(ux,uy,uz);

columize = @(A)reshape(permute(A, [3,2,1]), [],1,1);
tauU = [columize(uX), columize(uY), columize(uZ)]';
%tauU = [v2m(uX'), v2m(uY'), v2m(uZ')]';
U = [F_max.*ones(1,size(tauU,2)); tauU];

for k=1:hor
    alph = v2m(alph_vec(:,k),nx);
    y = y_vec(k);
    q = x_vec(:,k);
    Q = y.*(inv(alph'*alph)); % Shape Matrix

    % 1. Sample points on surface of ellipsoid
    [~,~,x_opt,K] = ellipPoints(q,Q,N);

    % 2. For each point, propagate forward combinations of inputs at
    %   bounds.
    for k_x=1:N-1
        xPrv = x_opt(:,k_x);
        xNxt = [];
        for k_u=1:size(U,2)
            xNxt(:,end+1) = runge_kutta4(@(y) vehiSt.dynFunc(y,U(:,k_u)), Ts, xPrv);
        end
    end

    % 3. Calculate error to surface of next ellipsoid
    xVal{end+1,1} = xNxt;
end

end

function traj_new = calcTraj(xstart,xend,kcurr,Ts,knum,vehi)
vehiSt = vehi.vehiSt;
% Update trajt and pts_ref
tcurr = kcurr*Ts;
vehi.trajt = [tcurr, vehi.trajt(2:end)];
vehi.pts_ref = [xstart, xend];

trajFunc = vehiSt.trajFunc;
traj = trajFunc(vehi.trajt, vehi.pts_ref); % Compute Trajectory based on waypoints
%knum = size(vehi.traj,2);
%traj = [vehi.traj(:,1:(kcurr-1)), traj, traj(:,end).*ones(size(traj,1),knum-size(traj,2)-kcurr+1)]; % Pad Trajectory vec with last value
traj = [traj, traj(:,end).*ones(size(traj,1),knum-size(traj,2)-kcurr+1)]; % Pad Trajectory vec with last value
traj_new = [traj; 0.*ones(1,size(traj,2)); 0.*ones(1,size(traj,2))]; % Yaw and dYaw
            
end

function [alph_vec, y] = calcParametopes(x, u, vehi, Ts)
vehiSt = vehi.vehiSt;
nx = vehiSt.nx;

alph_vec = zeros(nx^2,size(x,2));
alph_vec(:,1) = reshape(vehiSt.alph0,[],1);

y = zeros(1,size(x,2));
y(1) = vehiSt.ysim0;

%parametope_c = zeros(1,size(x,2));

% u_max = [vehiSt.vbds_init(nx+1,2); zeros(3,1)];
u_max = vehiSt.vbds_init(nx+1:end,2);

for k=1:size(x,2)-1
    % === Propagate Parametope === %
    jacFunc = vehiSt.jacFunc;
    Dfxo = jacFunc(x(:,k), u(:,k));
    
    alph = reshape(alph_vec(:,k),nx,[]);
    alphk_1 = runge_kutta4(@(y) prmtpU(y,-Dfxo,nx),Ts,alph);
    alph_vec(:,k+1) = reshape(alphk_1,[],1);
    
    % Calculate maximum/minimum achievable value for x in next timestep
    xmin = runge_kutta4(@(y) vehiSt.partDynFunc(y,vehiSt.vbds_init(nx+1:end,1)),Ts,x(1:nx,k));
    %xmax = runge_kutta4(@(y) vehiSt.partDynFunc(y,vehiSt.vbds_init(nx+1:end,2)),Ts,x(1:nx,k)); % In case for near-perfect prediction accross the horizon (prediction uncertainties don't stack)
    if k==1
        xmax = x(1:nx,k);
    end
    xmax = runge_kutta4(@(y) vehiSt.partDynFunc(y,u_max),Ts,xmax); % (prediction uncertainties do stack)
    vbds_ext(:,1) = vehiSt.vbds_init(:,1) + [xmin; zeros(vehiSt.nu,1)];
    vbds_ext(:,2) = vehiSt.vbds_init(:,2) + [xmax; zeros(vehiSt.nu,1)];

    % vbds_ext(:,1) = 1e-1.*vehiSt.vbds_init(:,1)+[x(1:nx,k); u(:,k)];
    % vbds_ext(:,2) = 1e-1.*vehiSt.vbds_init(:,2)+[x(1:nx,k); u(:,k)];
    % vbds_ext(:,1) = vehiSt.vbds_init(:,1)+[x(1:nx,k); zeros(vehiSt.nu,1)];
    % vbds_ext(:,2) = vehiSt.vbds_init(:,2)+[x(1:nx,k); zeros(vehiSt.nu,1)];
    
    % Calculate c just once outside of runge-kutta
    % Propagate y with parametope method
    %[eig_prmtp_c_vec,eig_parametope_c] = prmtpEllips_c(vehiSt.rStrt,vehiSt.symvars,vbds_ext,-Dfxo,alph_vec(:,k),nx,vehiSt.rCode,vehiSt.rHasVarMx,vehiSt.rConstDfMx);
    % Propagate y with scaling method
    % ratio = xmax./x(1:nx,k);
    % parametope_c = max(ratio(1:3));
    
    if k>1
        % parametope_c = (norm(xmax(1:3)-x(1:3,k))/sqrt(y(k-1)))^2;
        parametope_c = (norm(xmax(1:3)-x(1:3,k))-sqrt(y(k-1)))/sqrt(y(k-1));
        y(k) = runge_kutta4(@(yr) prmtpE(yr,parametope_c/Ts), Ts, y(k-1));
        % y(k+1) = runge_kutta4(@(yr) 1/parametope_c, Ts, y(k));
    else
        y(k) = norm(xmax(1:3)-x(1:3,k))^2;
    end

end
end

function [xpred, u, auxSt] = calcNxt(xk, xd, vehi, Ts, pred_hor)
% Calculate next step for vehi from x0 to xd
% Also Compute prediction
vehiSt = vehi.vehiSt;
ctrlSt = vehi.ctrlSt;
auxSt = ctrlSt.auxSt;

nx = vehiSt.nx;

% Compute predictions whether controller returns horizon of u values or not

states = zeros(nx,pred_hor);
states(:,1) = xk(1:nx);

alphs = zeros(nx^2,pred_hor);
alphs(:,1) = m2v(vehiSt.alph0);

ys = zeros(1,pred_hor);
ys(1) = vehiSt.ysim0;

u_max = [vehiSt.vbds_init(nx+1,2); zeros(3,1)];

% === Calculate Control action === %
[u,auxSt] = ctrlSt.ctrlFunc(xk(:,1), xd, auxSt);

% === Calculate System === %
for ind=1:size(u,2)
    % Propagate States
    states(:,ind+1) = runge_kutta4(@(y) vehiSt.partDynFunc(y,u(:,ind)), Ts, states(:,ind));

    % Propagate alph
    Dfx0 = vehiSt.jacFunc(states(:,ind),u(:,ind));
    alphk_1 = runge_kutta4(@(y) prmtpU(y,-Dfx0),Ts,v2m(alphs(:,ind),nx));
    alphs(:,ind+1) = m2v(alphk_1);

    % Propagate y
    % Calculate maximum/minimum achievable value for x in next timestep
    xmin = runge_kutta4(@(y) vehiSt.partDynFunc(y,vehiSt.vbds_init(end-3:end,1)),Ts,states(:,ind));
    %xmax = runge_kutta4(@(y) vehiSt.partDynFunc(y,vehiSt.vbds_init(nx+1:end,2)),Ts,x(1:nx,k)); % In case for near-perfect prediction accross the horizon (prediction uncertainties don't stack)
    if ind==1
        xmax = states(1:nx,ind);
    end
    xmax = runge_kutta4(@(y) vehiSt.partDynFunc(y,u_max),Ts,xmax); % (prediction uncertainties do stack)
    vbds_ext(:,1) = vehiSt.vbds_init(:,1) + [xmin; zeros(vehiSt.nu,1)];
    vbds_ext(:,2) = vehiSt.vbds_init(:,2) + [xmax; zeros(vehiSt.nu,1)];

    % Calculate c just once outside of runge-kutta
    % Propagate y with parametope method
    %[eig_prmtp_c_vec,eig_parametope_c] = prmtpEllips_c(vehiSt.rStrt,vehiSt.symvars,vbds_ext,-Dfx0,alphs(:,ind),nx,vehiSt.rCode,vehiSt.rHasVarMx,vehiSt.rConstDfMx);
    % Propagate y with scaling method
    % ratio = xmax./x(1:nx,k);
    % parametope_c = max(ratio(1:3));
    
    if ind>1
        % parametope_c = (norm(xmax(1:3)-x(1:3,k))/sqrt(y(k-1)))^2;
        parametope_c = (norm(xmax(1:3)-states(1:3,ind))-sqrt(ys(ind-1)))/sqrt(ys(ind-1));
        ys(ind) = runge_kutta4(@(yr) prmtpE(yr,parametope_c/Ts), Ts, ys(ind-1));
        % y(k+1) = runge_kutta4(@(yr) 1/parametope_c, Ts, y(k));
    else
        ys(ind) = norm(xmax(1:3)-states(1:3,ind))^2;
    end
end
xpred = [states(:,1:end-1); alphs(:,1:end-1); ys];
auxSt.ys = ys;

end

function [x, u] = calcPath(x0, vehi, t, Ts)
% Run trajectory prediction of vehi from x0 throughout t
vehiSt = vehi.vehiSt;
nx = vehiSt.nx;

x = zeros(nx,length(t)+1);
x(:,1) = x0;
u = zeros(vehiSt.nu,length(t));

for k=t
    kin = k-t(1)+1;
    % === Calculate Control action === %
    ctrlFunc = vehiSt.contrlFunc;
    u(:,kin) = ctrlFunc(x(:,kin), vehi.traj(:,k));
    %u(:,kin) = ctrlFunc(x(:,kin), xd);

    % === Calculate System === %
    dynFunc = vehiSt.dynFunc;
    x(:,kin+1) = runge_kutta4(@(y) dynFunc(y,u(:,kin)), Ts, x(:,kin));

end

end

function traj = plotTraj(x, traj)
    traj.XData = x(1,:);
    traj.YData = x(2,:);
    traj.ZData = x(3,:);
end

function MCplot = plotMonteCarlo(x_branch,prmtp_ind,MCplot_old)
% MCplot has two components: The branch plots, which will be 26 lines, and
% the leaf plots which will be N scatter plots (N is number of plotted
% ellipsoids)

% x_branch{1} : contains all branch samples at k=1
% x_leaf{1} : contains all leaf samples at k=1

nx = size(x_branch{1},1);

knum_branch = size(x_branch,1);
%knum_leaf = size(x_leaf,1);
num_branch = size(x_branch{1},2);
%num_leaf = size(x_leaf{1},2);

% Reshape x_branch so that each cell contains all samples of a branch instead
x_branchaux = cell(num_branch,1);
for k_br=1:num_branch
    x_branchaux{k_br,1} = zeros(nx,knum_branch);
    for k_t=1:knum_branch
        x_branchaux{k_br,1}(:,k_t) = x_branch{k_t,1}(:,k_br);
    end
end
x_branch = x_branchaux;

% create a default color map ranging from red to light pink
color_start = [1, 0, 0];
color_end = [255, 192, 203]/255;

if isempty(MCplot_old)
    MCplot.monteCarloPlot = struct;
    MCplot.monteCarloPlot.branch = cell(num_branch,1);
    %MCplot.monteCarloPlot.leaf = cell(knum_leaf,1);
    % Create MonteCarlo branch plots
    colorMapLength = knum_branch;
    colors_p = [linspace(color_start(1),color_end(1),colorMapLength)', linspace(color_start(2),color_end(2),colorMapLength)', linspace(color_start(3),color_end(3),colorMapLength)'];
    for k_br=1:num_branch
        MCplot.monteCarloPlot.branch{k_br,1} = plot3(zeros(1,knum_branch), zeros(1,knum_branch), zeros(1,knum_branch),'LineWidth',0.05, 'Color',[21,155,81]/255,'LineStyle','-', ...
            'Marker','*','MarkerIndices',prmtp_ind+1,'MarkerSize',0.1,'MarkerEdgeColor','r','DisplayName','');
    end
    % % Create MonteCarlo leaf plots
    % colorMapLength = knum_leaf;
    % colors_p = [linspace(color_start(1),color_end(1),colorMapLength)', linspace(color_start(2),color_end(2),colorMapLength)', linspace(color_start(3),color_end(3),colorMapLength)'];
    % for k=1:knum_leaf
    %     MCplot.monteCarloPlot.leaf{k,1} = plot3(zeros(1,num_leaf), zeros(1,num_leaf), zeros(1,num_leaf),'LineWidth',0.05, 'Color','r','LineStyle','none','Marker','*','DisplayName','');
    % end
else
    MCplot = MCplot_old;
end

% Update branch plots
for k=1:num_branch
    MCplot.monteCarloPlot.branch{k,1}.XData = x_branch{k,1}(1,:);
    MCplot.monteCarloPlot.branch{k,1}.YData = x_branch{k,1}(2,:);
    MCplot.monteCarloPlot.branch{k,1}.ZData = x_branch{k,1}(3,:);
end
% % Update leaf plots
% for k=1:knum_leaf
%     MCplot.monteCarloPlot.leaf{k,1}.XData = x_leaf{k,1}(1,:);
%     MCplot.monteCarloPlot.leaf{k,1}.YData = x_leaf{k,1}(2,:);
%     MCplot.monteCarloPlot.leaf{k,1}.ZData = x_leaf{k,1}(3,:);
% end

end

function p_vec = plotParametope(y,alph_vec,x,nx,ops,patch_vec,dim)
if isempty(patch_vec)
    p_vec = cell(1,size(y,2));
else
    p_vec = patch_vec;
end

for k=1:size(p_vec,2)
    if k>size(y,2)
        %set(p_vec{k},'LineWidth',0);
        set(p_vec{k},'FaceAlpha',0);
        p_vec{k}.Vertices = [0;0;0]';
        p_vec{k}.Faces = 1;
    else
        alph = v2m(alph_vec(:,k),nx);
        Q = (y(k)).*(inv(alph'*alph)); % Shape Matrix
        % Q = ((1)^-1.*alph'*alph); % Shape Matrix

        [w,W] = ellip2R3(x(:,k),Q);
        if dim==2
            [w,W] = ellipIntersectHyp(w,W,[0;0;1],0);
        end
        
        [~,~,x_opt,K] = ellipPoints(w,W,50);
    
        if isempty(patch_vec)
            p = patch('Vertices', x_opt','Faces',K,'FaceColor', ops.plot.edgecolor);
            set(p,'LineStyle',ops.plot.wirestyle);   
            set(p,'LineWidth',ops.plot.linewidth);
            set(p,'EdgeColor',ops.plot.edgecolor);
            set(p,'Facealpha',ops.plot.shade);     
            lighting gouraud;
            %view(3);
            camlight('headlight','infinite');
            camlight('headlight','infinite');
            camlight('right','local');
            camlight('left','local');   
            %drawnow;
        
            p_vec{k} = p;
        else
            p_vec{k}.Faces = K;
            p_vec{k}.Vertices = x_opt';
        end
    end
end
end
