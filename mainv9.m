% Mainv9 file of thesis project
% Runs all generalized functions for any number of vehichles
% Offline, initializes cell vector of vehiStruct
% Then initializes parametope propagation variables
% Online, runs specific vehichle and parametope dynamics 
% Every inner_step runs predicted trajectory and parametopes
% Defender tries to intercept invader
% Using NMPC to align parametopes
% Initialize vehi_vec, containing vehiSt (model) and ctrlSt (controller)

% This version is for simulations and performance analysis

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
sim.simT = 3; % Sim time seconds
sim.Ts = 0.01;
sim.knum = sim.simT/sim.Ts;
%calc_traj_dt = 0.1; % How often trajectory is recalculated [s]
sim.calc_traj_dt = sim.Ts;
sim.recalc_step = sim.calc_traj_dt/sim.Ts;
sim.pred_hor = 50;
sim.discnum = 5;
sim.disc_step = floor(sim.pred_hor/sim.discnum);
sim.disc_dt = sim.disc_step*sim.Ts; % Time step between Reachability discs
sim.prmtp_ind = [1,sim.disc_step:sim.disc_step:sim.pred_hor];
sim.discnum = length(sim.prmtp_ind);
sim.vid_name = ['DvAvid_N',num2str(sim.pred_hor)];
sim.fig_name = ['DvAfig_N',num2str(sim.pred_hor)];
sim.err_name = ['DvAerr_N',num2str(sim.pred_hor)];
sim.wOn = 0;

% Critical Area Ellipse
sim.Ecrit.q = [0;0;0];
sim.Ecrit.r = 2;
sim.Ecrit.Q = sim.Ecrit.r^2.*eye(3);
ops.plot.shade = 0.5;
ops.plot.linewidth=0.5;
ops.plot.wirestyle='-';
ops.plot.edgecolor='b';
ops.plot.Rspacename = '$\mathcal{G}_C$';
%ops.plot.edgecolor=[235 75 45];
sim.Ecrit.ops = ops;

sim.enable_shadow = 1;
sim.enable_disp = 1;
sim.enable_render = 1;
sim.make_vid = 0;
sim.load_vehi_vec = 0;

sim.tI = 2; % Time till A hits C
% === Define pA0, pAf, vAf, aAf === %
% sim.pA0 = [5;-5;6];
sim.pA0 = [20;-20;6];
sim.pAf = [sim.Ecrit.q(1:2);0.1];
% critv_ref = [-5;5;-5];
sim.jAf = [0;0;0];
% sim.aAf = [0;0;0];
dir = [-1;1;-1];
sim.aAf = 10*dir/norm(dir);
sim.vAf = sim.tI*sim.aAf;

% Define Inital Positions
p20 = [2;2;1];
p30 = [-5;-2;1];
sim.x0 = {sim.pA0, p20};
sim.Dnum = length(sim.x0)-1;

simT = sim.simT; % Sim time seconds
Ts = sim.Ts;
knum = sim.knum;
%calc_traj_dt = 0.1; % How often trajectory is recalculated [s]
calc_traj_dt = sim.calc_traj_dt;
recalc_step = sim.recalc_step;
pred_hor  = sim.pred_hor;
discnum = sim.discnum;
disc_step = sim.disc_step;
disc_dt = sim.disc_dt; % Time step between Reachability discs
prmtp_ind = sim.prmtp_ind;
discnum = sim.discnum;

tI = sim.tI;
% === Define pA0, pAf, vAf, aAf === %
pAf = sim.pAf;
% critv_ref = [-5;5;-5];
aAf = sim.aAf;
vAf = sim.vAf;

% Disturbance Variables
wT = 5; % Disturb pulse width
wPulse = 0; % Signals if disturb pulse is underway
wZ = 0; % Checks whether disturb function has crossed 0
wOn = [];

% Critical Area Ellipse
Ecrit = sim.Ecrit;

enable_shadow = sim.enable_shadow;
enable_disp = sim.enable_disp;
enable_render = sim.enable_render;
make_vid = sim.make_vid;
load_vehi_vec = sim.load_vehi_vec;

% === Init vehiStruct === &
if load_vehi_vec && exist('prmtpV8_data.mat') 
    load prmtpV8_data.mat;
else
    vehi_vec = init_all_vehi_local(sim); % Cell vec with vehiSimStruct
    vehiAux_vec = vehi_vec;
end
numVehi = length(vehi_vec);
% vehi_vec{1} is assumed to be the Invader 
if sim.Dnum > 1
    % Multi-defender scenario
    sim.vid_name = [sim.vid_name,'_m'];
    sim.fig_name = [sim.fig_name,'_m'];
    sim.err_name = [sim.err_name,'_m'];
else
    sim.vid_name = [sim.vid_name,'_s'];
    sim.fig_name = [sim.fig_name,'_s'];
    sim.err_name = [sim.err_name,'_s'];
end
if sim.wOn
    % Attacker Disturb on
    sim.vid_name = [sim.vid_name,'_unp'];
    sim.fig_name = [sim.fig_name,'_unp'];
    sim.err_name = [sim.err_name,'_unp'];
end

for ind=1:numVehi
    vehi = vehi_vec{ind};
    nx = vehi.vehiSt.nx;
    nx_full = nx+nx^2;

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

    %vehi.error_track = [];
    vehi.safeArea = [];
    vehi.w = zeros(3,1);

    % === Auxilliary vehi without plots for last prediction
    vehiAux_vec{ind} = vehi;
    vehiAux_vec{ind}.grd_ind = [];
    %vehiAux_vec{ind}.grdEllip = {}; % Cell_vec with ground ellipse 
    if ind == 1
        pts_ref = [pAf; vAf; vehi.vehiSt.x0(3*2+1:end)];
        vehi.pts_ref{1} = [vehi.vehiSt.x0, repmat(pts_ref,1,pred_hor)];

        vehiAux_vec{ind}.pts_ref{1} = [vehi.vehiSt.x0, repmat(pts_ref,1,knum)];
        vehiAux_vec{ind}.ops.plot.edgecolor = 'r'; vehiAux_vec{ind}.ops.plot.wirestyle = '-'; vehiAux_vec{ind}.ops.plot.linewidth = 0.5;
        vehiAux_vec{ind}.ops.plot.Rspacename = '$\Psi$';
    else
        pts_ref = repmat([vehi_vec{1}.vehiSt.x0;m2v(vehi_vec{1}.vehiSt.alph0);vehi_vec{1}.vehiSt.ysim0],1,pred_hor);
        vehi.pts_ref{1} = pts_ref;

        vehiAux_vec{ind}.pts_ref{1} = repmat(pts_ref(:,1),1,pred_hor);
        vehiAux_vec{ind}.ops.plot.edgecolor = 'g'; vehiAux_vec{ind}.ops.plot.wirestyle = '-'; vehiAux_vec{ind}.ops.plot.linewidth = 0.5;
        vehiAux_vec{ind}.ops.plot.Rspacename = '$\Upsilon$';

        %vehi.pts_ref = repmat(vehi.vehiSt.x0,1,pred_hor+1);
    end

    % Define window bounds

    vehi_vec{ind} = vehi;
end

%vehiPred = vehi_vec{1}; % Auxilliary vehi for prediction

%%
% === Simulation === %
telapse_total = tic;
k=1;
stop = 0;

while k<knum && ~stop 
    % Update subtitle
    subtitle_3D.String = ['t = ',num2str(k*Ts), 's'];
    disp(['=== Computing ', num2str(k), '/', num2str(knum), ' interval ===']);
    disp(['t = ', num2str((k-1)*Ts),'s']);
    telapse_k = tic;
    for kVehi=1:numVehi
        vehi = vehi_vec{kVehi};
        vehiAux = vehiAux_vec{kVehi};
        vehiSt = vehi.vehiSt;
        nx = vehiSt.nx;
        nx_full = nx+nx^2;
        ops = vehi.ops;

        if k>1
            % === Recalculate Trajectory (for Defender only) === %
            if kVehi~=1
                %vehi.pts_ref(:,end) = vehi_vec{1}.x(:,k); 
                vehi.pts_ref{end+1,1} = vehi_vec{1}.xpred{end,1};
                stop = stop || norm(vehi_vec{1}.x(1:3,k)-vehi.x(1:3,k))<=3;
            else
                stop = stop || vehi.x(3,k)<=0;
            end

        end
       
        if sim.wOn && kVehi == 1 
            % Add disturbances!
            % Changes desired position
            % Disturb holds for wT timesteps
            % Frequency scales with distance to target
            wFreq = norm(vehi.x(1:3,k)-Ecrit.q)^2;
            if wZ*sin(2*pi*wFreq*(k-1)*Ts) <= 0 && ~wPulse
                wPulse = 1;
                % w = vehi.vehiSt.tauBmax*randn(3,1);
                waux = [randn(2,1);0];
                w = sim.Ecrit.r.*waux/norm(waux);
                wD = wT;
            end
            wZ = sin(2*pi*wFreq*(k-1)*Ts);

            if wPulse
                wD = wD-1;
                wPulse = logical(wD);
            else
                w = zeros(3,1);
            end
            wOn(:,k) = w(1:3);
            vehi.w = w;
            % w = [w;zeros(size(vehi.pts_ref{end,1},1)-3,1)];
            % vehi.pts_ref{end,1} = vehi.pts_ref{end,1} + [w;zeros(size(vehi.pts_ref{end,1},1)-3,1)];
        end
        % === Real-Time === %
        wFull = [vehi.w;zeros(size(vehi.pts_ref{end,1},1)-3,1)];
        [vehi.xpred{end+1,1}, vehi.upred{end+1,1}, vehi.ctrlSt.auxSt] = calcNxt(vehi.x(:,k), vehi.pts_ref{end,1} + wFull, vehi, Ts, pred_hor);
        vehi.x(:,k+1) = vehi.xpred{end,1}(1:nx,2);
        vehi.u(:,k) = vehi.upred{end,1}(:,1);
        %vehi.error_track(:,k) = vehi.ctrlSt.auxSt.error_track;

        vehi_vec{kVehi} = vehi;

        if enable_shadow
        % === Compute Final parametope
            if kVehi==1
                active_shadow = 1;
                [xpred_aux, upred_aux, auxStaux] = calcNxt(vehiAux.x(:,k), vehiAux.pts_ref{end,1}, vehiAux, Ts, knum);
                vehiAux.xpred{end+1,1} = xpred_aux(:,1:knum);
                vehiAux.upred{end+1,1} = upred_aux(:,1:knum);
                vehiAux.ctrlSt.auxSt = auxStaux;
        
                vehiAux.x(:,k+1) = vehiAux.xpred{end,1}(1:nx,2);
                vehiAux.u(:,k) = vehiAux.upred{end,1}(:,1);
                vehiAux.error_track(:,k) = vehiAux.ctrlSt.auxSt.error_track;
            else
                % 1. Use Attacker's xpred and upred
                % 2. alph0, ysim0 from ellipsoid with lowest D_KL
                DKL_min = Inf;
                qAs = vehi_vec{1}.xpred{end,1}(1:nx,:);
                qDs = vehi_vec{kVehi}.xpred{end,1}(1:nx,:);

                pred_dist = sqrt(sum(((qAs(1:3,:) - qDs(1:3,:)).^2),1));
                [min_pred_dist, min_dist_ind] = min(pred_dist);
                yA = vehi_vec{1}.xpred{end,1}(end,min_dist_ind);
                yD = vehi_vec{kVehi}.xpred{end,1}(end,min_dist_ind);

                ymin = sqrt(yA)+sqrt(yD);
                active_shadow = min_pred_dist <= ymin;

                if active_shadow
                    alphA = v2m(vehi_vec{1}.xpred{end,1}(nx+(1:nx^2),min_dist_ind),nx);
                    alphD = v2m(vehi_vec{kVehi}.xpred{end,1}(nx+(1:nx^2),min_dist_ind),nx);
                    
                    % Some future compatibility problems here: The dynamics used SHOULD be the Attacker's, AND alph used should be compatible.
                    % Right now theres no guarantee they are compatible (this only works because they so happen to have the same dynamics)
                    xtrajA = vehiAux_vec{1}.xpred{end,1}(:,pred_hor+1:end);
                    utrajA = vehiAux_vec{1}.upred{end,1}(:,pred_hor+1:end);
                    [alphs, ys] = calcParametopes(xtrajA, utrajA, alphD, yD, vehi_vec{1}, Ts);

                    offset = [qDs(1:3,min_dist_ind);zeros(nx-3,1)]-vehiAux_vec{1}.xpred{end,1}(1:nx,1);
                    vehiAux.xpred{end+1,1} = [vehiAux_vec{1}.xpred{end,1}(1:nx,:)+offset; [vehiAux_vec{1}.xpred{end,1}(nx+1:end,1:pred_hor), [alphs; ys]]];
                else 
                    vehiAux.xpred{end+1,1} = [];
                end
                
            end
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
            vehiAux.grd_ind(end+1) = grd_ind;

            if vehiAux.grd_ind(end)~=0
                % Get Safeguarded Ground Area
                xSafe = vehiAux.xpred{end,1}(1:nx,vehiAux.grd_ind(end));
                alphSafe = v2m(vehiAux.xpred{end,1}(nx+(1:nx^2),vehiAux.grd_ind(end)),nx);
                ySafe = vehiAux.xpred{end,1}(end,vehiAux.grd_ind(end));

                Qsafe = ellipFromAlph(alphSafe,ySafe);

                [wSafe,Wsafe] = ellip2R3(xSafe,Qsafe);
                [wSafe,Wsafe] = ellipIntersectHyp(wSafe,Wsafe,[0;0;1],0);
            else
                wSafe = [0;0];
                Wsafe = zeros(2);
            end
            vehiAux.grdEllip(k).w = wSafe;
            vehiAux.grdEllip(k).W = Wsafe;

            if kVehi~=1
                % Intersect Defender Safe Area with Attacker Ground area
                % [wSafe,Wsafe] = ellipIntersectEllip(wSafe,Wsafe,Ecrit.q(1:2),Ecrit.Q(1:2,1:2));
                [wI,WI] = ellipIntersectEllip(vehiAux.grdEllip(k).w,vehiAux.grdEllip(k).W,vehiAux_vec{1}.grdEllip(k).w,vehiAux_vec{1}.grdEllip(k).W);
                vehiAux.safeArea(end+1) = ellipVol(WI);
            else
                vehiAux.safeArea(end+1) = ellipVol(vehiAux.grdEllip(k).W);
            end
            
            % vehiAux.grdEllip{end+1,1} = E;
        end

        vehi_vec{kVehi} = vehi;
        vehiAux_vec{kVehi} = vehiAux;
        
        %legend({'target', 'r_T', 'Path'});
    end
    disp(['Computing time for interval: ', num2str(toc(telapse_k)), 's']);
    disp('----------------');
    disp(' ');
    k = k+1;
end
% hold off;
disp(['Total sim computing time: ', num2str(toc(telapse_total)), 's']);

% Plot Error Graph Considering multiple D
for ind=1:numVehi
    if ind < numVehi
        dist_err(ind,:) = sqrt(sum((vehi_vec{ind+1}.x(1:3,1:k)-vehi_vec{1}.x(1:3,1:k)).^2,1));
    end
    safeArea(ind,:) = vehiAux_vec{ind}.safeArea;
end
% dist_err = vehi_vec{2}.x(1:3,:)-vehi_vec{1}.x(1:3,:);
% safeAreaA = vehiAux_vec{1}.safeArea;
% safeAreaD = vehiAux_vec{2}.safeArea;
utrack = vehi_vec{2}.u;
t = Ts.*(0:k-1);

textratio = 0.6;
fig_err = figure('WindowState','maximized');
title('Performance in Area Defense',['$N_{H} = $',num2str(pred_hor)])
subplot(2,1,1);
grid on;
hold on;
title('$||p_D-p_A||_2$','Interpreter','latex')
dplot = plot(t,dist_err,'LineWidth',1);
for ind=1:numVehi-1
    dplot(ind).DisplayName = ['$D_',num2str(ind),'$'];
end
legend('Interpreter','latex');
hold off;
subplot(2,1,2);
grid on;
hold on;
title('$vol(\mathcal{G})$','Interpreter','latex')
aplot = plot(t,[zeros(size(safeArea,1),1),safeArea],'LineWidth',1);
aplot(1).DisplayName = 'A'; 
for ind=2:numVehi
    aplot(ind).DisplayName = ['$D_',num2str(ind-1),'$'];
end
legend('Interpreter','latex');
% fontsize(10,"points");
fontsize(11*sqrt(2)/0.9/textratio,"points");
hold off;
% savepdf(fig_err);
exportgraphics(fig_err,['DvAfigs\',sim.err_name,'.pdf']);
% saveas(fig_err,['DvAfigs\',sim.err_name,'.pdf']);

figure;
hold on;
title('Diturbance On')
hold off;
plot(sqrt(sum(wOn.^2,1)));

figure;
hold on;
title('Thrust force')
hold off;
plot(vehi_vec{2}.u(1,1:k))

% === Save Data
vehi_save = cell(1,numVehi);
for kVehi=1:numVehi
    vehi_save{kVehi}.x = vehi_vec{kVehi}.x;
    vehi_save{kVehi}.u = vehi_vec{kVehi}.u;
    vehi_save{kVehi}.xpred = vehi_vec{kVehi}.xpred;
    vehi_save{kVehi}.upred = vehi_vec{kVehi}.upred;
    vehi_save{kVehi}.safeArea = vehi_vec{kVehi}.safeArea;
    vehi_save{kVehi}.pts_ref = vehi_vec{kVehi}.pts_ref;
end
save prmtpV9_data.mat vehi_save % Save sim results (probably can't support casadi solver)

if enable_render
    render(vehi_vec,vehiAux_vec,sim);
end

function vehi_vec = init_all_vehi_local(sim)
%INIT_ALL_VEHI outputs cell vector of vehiStruct
%   Detailed explanation goes here
Ts = sim.Ts;
pred_hor = sim.pred_hor;
tI = sim.tI;
Dnum = sim.Dnum;
x0 = sim.x0;

vehi_vec = cell(1,1+Dnum);

vehi_vec{1}.vehiSt = init_QuadDroneNMPC(x0{1},Ts);
vehi_vec{1}.ctrlSt = init_MellingerCtrl(sim.vAf,sim.aAf,sim.jAf,pred_hor,tI,vehi_vec{1}.vehiSt);
for k=2:Dnum+1
    vehi_vec{k}.vehiSt = init_QuadDroneNMPC(x0{k},Ts);
    vehi_vec{k}.ctrlSt = init_NMPCctrl_v5(vehi_vec{1}.vehiSt.nx,pred_hor,Ts,vehi_vec{k}.vehiSt);
end
%vehi_vec{2}.vehiSt = init_doubleInteg3D(x0{2},Ts);

% vehi_vec{1}.ctrlSt = init_NMPCctrl(vehi_vec{1}.vehiSt.nx,5,Ts,vehi_vec{1}.vehiSt);

% vehi_vec{1}.ctrlSt = init_MellingerCtrl_v3(pred_hor,Ts,vehi_vec{1}.vehiSt);
%vehi_vec{2}.ctrlSt = init_NMPCctrl_v2(vehi_vec{1}.vehiSt.nx,pred_hor,Ts,vehi_vec{2}.vehiSt);
% vehi_vec{2}.ctrlSt = init_NMPCctrl_v4(vehi_vec{1}.vehiSt.nx,pred_hor,Ts,vehi_vec{2}.vehiSt);
%vehi_vec{2}.ctrlSt = init_NMPCctrl(vehi_vec{1}.vehiSt.nx,pred_hor,Ts,vehi_vec{2}.vehiSt);
%vehi_vec{2}.ctrlSt = init_MellingerCtrl_v2(pred_hor,Ts,vehi_vec{2}.vehiSt);
%vehi_vec{2}.ctrlSt = init_MellingerCtrl(4,vehi_vec{2}.vehiSt);

% === Define ops === %
ops = sdpsettings('verbose',0);
ops.plot.shade = 0.15;
ops.plot.linewidth=0.05;
ops.plot.wirestyle=':';
ops.plot.edgecolor='#EB4B2D';
ops.plot.Rspacename = '$\mathcal{A}$';
%ops.plot.edgecolor=[235 75 45];
vehi_vec{1}.ops = ops;

ops = sdpsettings('verbose',0);
ops.plot.shade = 0.15;
ops.plot.linewidth=0.05;
ops.plot.wirestyle=':';
ops.plot.edgecolor='#11F4BF';
ops.plot.Rspacename = '$\mathcal{D}$';
%ops.plot.edgecolor=[17 244 191];

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
for k=2:Dnum+1
    vehi_vec{k}.ops = ops;
    vehi_vec{k}.p_vec = p_vec;
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

function [alph_vec, y] = calcParametopes(x, u, alph0, y0, vehi, Ts)
vehiSt = vehi.vehiSt;
nx = vehiSt.nx;

alph_vec = zeros(nx^2,size(x,2));
alph_vec(:,1) = m2v(alph0);
%alph_vec(:,1) = reshape(vehiSt.alph0,[],1);

y = zeros(1,size(x,2));
y(1) = y0;
%y(1) = vehiSt.ysim0;

u_max = [vehiSt.vbds_init(nx+1,2); zeros(3,1)];

for k=1:size(x,2)-1
    % === Propagate Parametope === %
    jacFunc = vehiSt.jacFunc;
    Dfxo = jacFunc(x(1:nx,k), u(:,k));
    
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
% d : disturbance to input
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
%u(:,1) = u(:,1)+[0;d]; % Add distubances

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
        %Q = (y(k)).*(inv(alph'*alph)); % Shape Matrix
        % Q = ((1)^-1.*alph'*alph); % Shape Matrix
        Q = ellipFromAlph(alph,y(k));

        [w,W] = ellip2R3(x(:,k),Q);
        if dim==2
            [w,W] = ellipIntersectHyp(w,W,[0;0;1],0);
        end
        
        [~,~,x_opt,K] = ellipPoints(w,W,100);
    
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

function render(vehi_vec,vehiAux_vec,sim)
prmtp_ind = sim.prmtp_ind;
disc_step = sim.disc_step;
Ts = sim.Ts;
enable_shadow = sim.enable_shadow;
enable_disp = sim.enable_disp;
make_vid = sim.make_vid;
Ecrit = sim.Ecrit;

% === Start Figure === %
fig=figure('WindowState','normal');
% ax3D = subplot(2,2,[1 3]);
% axTop = subplot(2,2,2);
% axSide = subplot(2,2,4);
% axes(ax3D)
hold on;
grid on
view(3)
[title_3D,subtitle_3D] = title('Adversary Drone Sim', ...
    't = 0 s','Interpreter','latex');
xlabel('$x$','Interpreter','latex');
ylabel('$y$','Interpreter','latex');
zlabel('$z$','Interpreter','latex');
handle2display = [];
name2display = {};

Ecritpatch = plotParametope(1,Ecrit.Q,Ecrit.q,3,Ecrit.ops,{},2);
Ecritpatch = Ecritpatch{1};
Ecritpatch.DisplayName = Ecrit.ops.plot.Rspacename;
handle2display(end+1) = Ecritpatch;

for ind=1:size(vehi_vec,2)
    pts_ref = vehi_vec{ind}.pts_ref{1};

    vehiAux_vec{ind}.prmtpShadow = {}; % Parametope Shadow patch
    vehiAux_vec{ind}.prmtpPatch = {}; % Parametope Shadow patch
    vehi_vec{ind}.prmtpPatch = {};
    handle2display(end+1) = patch(NaN,NaN,NaN,'EdgeColor',vehi_vec{ind}.ops.plot.edgecolor, ...
        'FaceColor', vehi_vec{ind}.ops.plot.edgecolor, ...
        'LineStyle',vehi_vec{ind}.ops.plot.wirestyle, ...
        'LineWidth',vehi_vec{ind}.ops.plot.linewidth, ...
        'Facealpha',vehi_vec{ind}.ops.plot.shade,'DisplayName',vehi_vec{ind}.ops.plot.Rspacename);
    handle2display(end+1) = patch(NaN,NaN,NaN,'EdgeColor',vehiAux_vec{ind}.ops.plot.edgecolor, ...
        'FaceColor', vehiAux_vec{ind}.ops.plot.edgecolor, ...
        'LineStyle',vehiAux_vec{ind}.ops.plot.wirestyle, ...
        'LineWidth',vehiAux_vec{ind}.ops.plot.linewidth, ...
        'Facealpha',vehiAux_vec{ind}.ops.plot.shade,'DisplayName',vehiAux_vec{ind}.ops.plot.Rspacename);

    vehi_vec{ind}.wayptsPlot = plot3(pts_ref(1,:),pts_ref(2,:),pts_ref(3,:),'xr','DisplayName','target');
    vehi_vec{ind}.bodyPlot = vehi_vec{ind}.vehiSt.bodyPlot();
    vehi_vec{ind}.pathPlot = animatedline(vehi_vec{ind}.x(1,1), vehi_vec{ind}.x(2,1), vehi_vec{ind}.x(3,1),'LineWidth',2.5, 'Color','b','LineStyle','-','DisplayName','Path');
end
axis equal
axis manual
% ac = get(ax3D, 'Children');   %no hidden objects
% copyobj( ac, axTop );
% copyobj( ac, axSide );
%legend(handle2display,'Interpreter','latex');
knum = length(vehi_vec{1}.xpred);

%frame_ind = 1:floor((knum)/4):(knum);
frame_ind = floor(linspace(1,knum,4));
if ~make_vid
    % Save 4 frames of figure %
    sim_ind = frame_ind;
else
    sim_ind = 1:knum;
end
% frame = [];
framek = 1;
textratio = 0.6;
for k=sim_ind 
    % Update subtitle
    subtitle_3D.String = ['t = ',num2str(k*Ts), 's'];

    for kVehi=1:length(vehi_vec)
        vehi = vehi_vec{kVehi};
        vehiAux = vehiAux_vec{kVehi};
        vehiSt = vehi.vehiSt;
        nx = vehiSt.nx;
        nx_full = nx+nx^2;
        ops = vehi.ops;

        % === Update Plot === %
        plotBodyFunc = vehiSt.plotBodyFunc; 
        vehi.bodyPlot = plotBodyFunc(vehi.x(:,k), vehi.bodyPlot);
        if k==1
            path_ind = 1;
        else
            ind = find(sim_ind==k);
            path_ind = sim_ind(ind-1):sim_ind(ind);
            path_ind = path_ind(2:end);
        end
        addpoints(vehi.pathPlot,vehi.x(1,path_ind), vehi.x(2,path_ind), vehi.x(3,path_ind));
        active_shadow = vehiAux.grd_ind(k) ~= 0;
        if enable_shadow && active_shadow
                vehiAux.prmtpShadow = plotParametope(vehiAux.xpred{k,1}(end,vehiAux.grd_ind(k)),vehiAux.xpred{k,1}(nx+1:end-1,vehiAux.grd_ind(k)),vehiAux.xpred{k,1}(1:nx,vehiAux.grd_ind(k)),nx,vehiAux.ops,vehiAux.prmtpShadow,2);
        end
        if enable_disp
            vehi.prmtpPatch = plotParametope(vehi.xpred{k,1}(end,prmtp_ind),vehi.xpred{k,1}(nx+1:end-1,prmtp_ind),vehi.xpred{k,1}(1:nx,prmtp_ind),nx,ops,vehi.prmtpPatch,3);
            if 0&&enable_shadow && kVehi~=1 && active_shadow
                shadow_ind = [prmtp_ind(end):disc_step:vehiAux.grd_ind(k)];
                vehiAux.prmtpPatch = plotParametope(vehiAux.xpred{k,1}(end,shadow_ind),vehiAux.xpred{k,1}(nx+1:end-1,shadow_ind),vehiAux.xpred{k,1}(1:nx,shadow_ind),nx,vehiAux.ops,vehiAux.prmtpPatch,3);
            end
        end
        if 0&&k==1
            axis equal
            axis manual
        end
        vehi_vec{kVehi} = vehi;
        vehiAux_vec{kVehi} = vehiAux;
    end
    legend(handle2display,'Interpreter','latex');
    fontsize(24,"points");
    % fontsize(11*sqrt(2)/0.9/textratio,"points");
    drawnow;
    %pause(Ts);
 
    if ~all(frame_ind-k)
        % Save 4 frames as figs
        exportgraphics(fig,['DvAfigs\',sim.fig_name,num2str(framek),'.pdf']);
        %saveas(fig,['DvAfigs\',sim.fig_name,num2str(framek),'.pdf']);
        framek = framek + 1;
        % figs(k) = copyobj(gca,fig);
    end
    if make_vid == 1
        frames(k) = getframe(fig);
    end
end
hold off;

if make_vid == 1
    vid=VideoWriter(['DvAfigs\',sim.vid_name,'.mp4'],'MPEG-4');
    vid.FrameRate=1/Ts;
    open(vid);
    writeVideo(vid,frames);
    close(vid);
end
% Save 4 frames as figs
% for k=frame_ind
%     saveas(figs(k),['DvAfigs\DvAsim_k',num2str(k),'_single_N',num2str(sim.pred_hor),'.png']);
% end

end

