% Mainv7 file of thesis project
% Runs all generalized functions for any number of vehichles
% Offline, initializes cell vector of vehiStruct
% Then initializes parametope propagation variables
% Online, runs specific vehichle and parametope dynamics 
% Every inner_step runs predicted trajectory and parametopes
% Defender tries to intercept invader
% Using NMPC to align parametopes
% Initialize vehi_vec, containing vehiSt (model) and ctrlSt (controller)

% This version is more organised

close all;
clearvars;

% vehiSimStruct : struct (struct with components of main sim)
%   .vehiSt : Struct with functions and variables of model
%   .ctrlSt : Struct with functions and variables of controller
%   .x
%   .u
%   .alph
%   .y
%   .xpred : [states; vec(alph); ypred]
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
pred_hor = 20;
discnum = 5;
disc_step = floor(pred_hor/discnum);
disc_dt = disc_step*Ts; % Time step between Reachability discs
prmtp_ind = [1,disc_step:disc_step:pred_hor];
discnum = length(prmtp_ind);

enable_prmtp = 1;
enable_disp = 1;
make_vid = 1;
load_vehi_vec = 0;

% === Init vehiStruct === &
if load_vehi_vec && exist('prmtpV5_data.mat') 
    load prmtpV5_data.mat;
else
    vehi_vec = init_all_vehi_local(Ts,pred_hor); % Cell vec with vehiSimStruct
end
numVehi = length(vehi_vec);
% vehi_vec{1} is assumed to be the Invader 

for ind=1:numVehi
    vehi = vehi_vec{ind};
    nx = vehi.vehiSt.nx;
    nx_full = nx+nx^2;

    if ind == 1
        % === Define Critical Area === %
        critp_ref = [0;0;0.5];
        % critv_ref = [-5;5;-5];
        critv_ref = [0;0;0];
        pts_ref = [critp_ref; critv_ref; vehi.vehiSt.x0(3*2+1:end)];
        vehi.pts_ref{1} = [vehi.vehiSt.x0, repmat(pts_ref,1,pred_hor)];
    else
        pts_ref = repmat([vehi_vec{1}.vehiSt.x0;m2v(vehi_vec{1}.vehiSt.alph0);vehi_vec{1}.vehiSt.ysim0],1,pred_hor);
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
    %vehi.ypred = cell(0,0);

    vehi.error_track = [];

    % Define window bounds

    vehi_vec{ind} = vehi;
end

%vehiPred = vehi_vec{1}; % Auxilliary vehi for prediction

% === Start Figure === %
fig=figure('WindowState','maximized');
% ax3D = subplot(2,2,[1 3]);
% axTop = subplot(2,2,2);
% axSide = subplot(2,2,4);
% axes(ax3D)
hold on;
grid on
view(3)
[title_3D,subtitle_3D] = title('Adversary Two Drone Sim', ...
    't = 0 s','Interpreter','latex');
xlabel('$x$','Interpreter','latex');
ylabel('$y$','Interpreter','latex');
zlabel('$z$','Interpreter','latex');
for ind=1:numVehi
    pts_ref = vehi_vec{ind}.pts_ref{1};
    vehi_vec{ind}.prmtpPatch = {};
    vehi_vec{ind}.wayptsPlot = plot3(pts_ref(1,:),pts_ref(2,:),pts_ref(3,:),'xr','DisplayName','target');
    %vehi_vec{ind}.trajPlot = plot3(vehi_vec{ind}.traj(1,:), vehi_vec{ind}.traj(2,:), vehi_vec{ind}.traj(3,:) ,'-.k','DisplayName','r_T');
    vehi_vec{ind}.bodyPlot = vehi_vec{ind}.vehiSt.bodyPlot();
    vehi_vec{ind}.pathPlot = animatedline(vehi_vec{ind}.x(1,1), vehi_vec{ind}.x(2,1), vehi_vec{ind}.x(3,1),'LineWidth',2.5, 'Color','b','LineStyle','-','DisplayName','Path');
end
axis equal
axis manual
% ac = get(ax3D, 'Children');   %no hidden objects
% copyobj( ac, axTop );
% copyobj( ac, axSide );
%%
% === Simulation === %
telapse_total = tic;
k=1;
while k<knum && norm(vehi_vec{1}.x(1:3,k)-vehi_vec{2}.x(1:3,k))>=0.5 
    % Update subtitle
    subtitle_3D.String = ['t = ',num2str(k*Ts), 's'];
    disp(['=== Computing ', num2str(k), '/', num2str(knum), ' interval ===']);
    disp(['t = ', num2str((k-1)*Ts),'s']);
    telapse_k = tic;
    for kVehi=1:numVehi
        vehi = vehi_vec{kVehi};
        vehiSt = vehi.vehiSt;
        nx = vehiSt.nx;
        nx_full = nx+nx^2;
        ops = vehi.ops;

        if k>1
            % === Recalculate Trajectory (for Defender only) === %
            if kVehi~=1
                %vehi.pts_ref(:,end) = vehi_vec{1}.x(:,k);
                vehi.pts_ref{end+1,1} = vehi_vec{1}.xpred{end,1};
            end

        end
       
        % === Real-Time === %
        % Add disturbances!
        [vehi.xpred{end+1,1}, vehi.upred{end+1,1}, vehi.ctrlSt.auxSt] = calcNxt(vehi.x(:,k), vehi.pts_ref{end,1}, vehi, Ts, pred_hor);
        vehi.x(:,k+1) = vehi.xpred{end,1}(1:nx,2);
        vehi.u(:,k) = vehi.upred{end,1}(:,1);
        vehi.error_track(:,k) = vehi.ctrlSt.auxSt.error_track;

        %[alph_vec, vehi.ypred{end+1,1}] = calcParametopes(vehi.xpred{end,1}, vehi.upred{end,1}, vehi, Ts);

        % === Update Plot === %
        plotBodyFunc = vehiSt.plotBodyFunc; 
        vehi.bodyPlot = plotBodyFunc(vehi.x(:,k), vehi.bodyPlot);
        addpoints(vehi.pathPlot,vehi.x(1,k), vehi.x(2,k), vehi.x(3,k));
        if enable_disp
            vehi.prmtpPatch = plotParametope(vehi.xpred{end,1}(end,prmtp_ind),vehi.xpred{end,1}(nx+1:end-1,prmtp_ind),vehi.xpred{end,1}(1:nx,prmtp_ind),nx,ops,vehi.prmtpPatch);
            %vehi.prmtpPatch = plotParametope(vehiSt.ysim0.*ones(1,length(prmtp_ind)),vehi.xpred(nx+1:end,prmtp_ind),vehi.xpred(1:nx,prmtp_ind),nx,ops,vehi.prmtpPatch);
            %vehi.prmtpPatch = plotParametope(vehi.ypred,vehi.alphpred,vehi.xpred(:,prmtp_ind),nx,ops,vehi.prmtpPatch);
            if k==1
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
    drawnow;
    pause(Ts);
    if make_vid == 1
        frame(k) = getframe(fig);
    end
    k = k+1;
end
hold off;
disp(['Total sim computing time: ', num2str(toc(telapse_total)), 's']);

if make_vid == 1
    vid=VideoWriter('test_v7.mp4','MPEG-4');
    vid.FrameRate=1/Ts;
    open(vid);
    writeVideo(vid,frame);
    close(vid);
end

% Plot Error Graph
err = vehi_vec{2}.error_track;
figure;
subplot(3,1,1);
grid on;
hold on;
plot(Ts.*(0:k-2),err(1,:),'DisplayName','Scaled Dist Error');
plot(Ts.*(0:k-2),err(4,:),'DisplayName','UnScaled Dist Error');
legend;
hold off;
subplot(3,1,2);
grid on;
hold on;
plot(Ts.*(0:k-2),err(2,:),'DisplayName','Scaled Control Action');
plot(Ts.*(0:k-2),err(5,:),'DisplayName','UnScaled Control Action');
legend;
hold off;
subplot(3,1,3);
grid on;
hold on;
plot(Ts.*(0:k-2),err(3,:),'DisplayName','Scaled Min DKL');
plot(Ts.*(0:k-2),err(6,:),'DisplayName','UnScaled Min DKL');
xlabel('$t [s]$','Interpreter','latex');
legend;
hold off;

%save prmtpV5_data.mat vehi_vec % Save sim results (probably can't support casadi solver)

function vehi_vec = init_all_vehi_local(Ts,pred_hor)
%INIT_ALL_VEHI outputs cell vector of vehiStruct
%   Detailed explanation goes here

vehi_vec = cell(1,2);

% Define Inital Positions
p10 = [5;-5;6];
% p20 = [-5;5;0];
p20 = [-1;1;1];
x0 = {p10, p20};

vehi_vec{1}.vehiSt = init_QuadDroneNMPC(x0{1},Ts);
vehi_vec{2}.vehiSt = init_QuadDroneNMPC(x0{2},Ts);
%vehi_vec{2}.vehiSt = init_doubleInteg3D(x0{2},Ts);

% vehi_vec{1}.ctrlSt = init_NMPCctrl(vehi_vec{1}.vehiSt.nx,5,Ts,vehi_vec{1}.vehiSt);

vehi_vec{1}.ctrlSt = init_MellingerCtrl(pred_hor,2,vehi_vec{1}.vehiSt);
%vehi_vec{2}.ctrlSt = init_NMPCctrl_v2(vehi_vec{1}.vehiSt.nx,pred_hor,Ts,vehi_vec{2}.vehiSt);
vehi_vec{2}.ctrlSt = init_NMPCctrl_v4(vehi_vec{1}.vehiSt.nx,pred_hor,Ts,vehi_vec{2}.vehiSt);
%vehi_vec{2}.ctrlSt = init_NMPCctrl(vehi_vec{1}.vehiSt.nx,pred_hor,Ts,vehi_vec{2}.vehiSt);
%vehi_vec{2}.ctrlSt = init_MellingerCtrl_v2(pred_hor,Ts,vehi_vec{2}.vehiSt);
%vehi_vec{2}.ctrlSt = init_MellingerCtrl(4,vehi_vec{2}.vehiSt);

% === Define ops === %
ops = sdpsettings('verbose',0);
ops.plot.shade = 0.15;
ops.plot.linewidth=0.05;
ops.plot.wirestyle=':';
ops.plot.edgecolor='#EB4B2D';
%ops.plot.edgecolor=[235 75 45];
vehi_vec{1}.ops = ops;

ops = sdpsettings('verbose',0);
ops.plot.shade = 0.15;
ops.plot.linewidth=0.05;
ops.plot.wirestyle=':';
ops.plot.edgecolor='#11F4BF';
%ops.plot.edgecolor=[17 244 191];
vehi_vec{2}.ops = ops;

% === Init Plot Obs === %
L = 0.5;
p_NB = [L;0;0];
p_SB = [-L;0;0];
p_EB = [0;L;0];
p_OB = [0;-L;0];
p_vec = [p_NB; p_EB; p_SB; p_OB];
vehi_vec{1}.p_vec = p_vec;

L = 0.5;
p_NB = [L;0;0];
p_SB = [-L;0;0];
p_EB = [0;L;0];
p_OB = [0;-L;0];
p_vec = [p_NB; p_EB; p_SB; p_OB];
vehi_vec{2}.p_vec = p_vec;

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

u_max = [vehiSt.vbds_init(nx+1,2); zeros(3,1)];

for k=1:size(x,2)-1
    % === Propagate Parametope === %
    jacFunc = vehiSt.jacFunc;
    Dfxo = jacFunc(x(:,k), u(:,k));
    
    alph = reshape(alph_vec(:,k),nx,[]);
    alphk_1 = runge_kutta4(@(y) prmtpU(y,-Dfxo),Ts,alph);
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
    [eig_prmtp_c_vec,eig_parametope_c] = prmtpEllips_c(vehiSt.rStrt,vehiSt.symvars,vbds_ext,-Dfxo,alph_vec(:,k),nx,vehiSt.rCode,vehiSt.rHasVarMx,vehiSt.rConstDfMx);
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
    [eig_prmtp_c_vec,eig_parametope_c] = prmtpEllips_c(vehiSt.rStrt,vehiSt.symvars,vbds_ext,-Dfx0,alphs(:,ind),nx,vehiSt.rCode,vehiSt.rHasVarMx,vehiSt.rConstDfMx);
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

function p_vec = plotParametope(y,alph_vec,x,nx,ops,patch_vec)
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
