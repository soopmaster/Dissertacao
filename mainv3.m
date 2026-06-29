% Mainv3 file of thesis project
% Runs all generalized functions for any number of vehichles
% Offline, initializes cell vector of vehiStruct
% Then initializes parametope propagation variables
% Online, runs specific vehichle and parametope dynamics 
% Every inner_step runs predicted trajectory and parametopes
% Defender tries to intercept invader

close all;
clearvars;

% vehiSimStruct : struct (struct with components of main sim)
%   .vehiSt;
%   .x
%   .u
%   .alph
%   .y
%   .traj
%   .pathPlot
%   .ops

% === Sim Settings === %
simT = 5; % Sim time seconds
Ts = 0.01;
knum = simT/Ts;
%calc_traj_dt = 0.1; % How often trajectory is recalculated [s]
calc_traj_dt = Ts;
recalc_step = calc_traj_dt/Ts;
disc_dt = 0.1; % Time step between Reachability discs
%disc_dt = Ts;
disc_step = disc_dt/Ts;
discnum = 5;
enable_prmtp = 1;

% === Init vehiStruct === &
vehi_vec = local_init_all_vehi(Ts); % Cell vec with vehiSimStruct
numVehi = length(vehi_vec);
% vehi_vec{1} is assumed to be the Invader

vehi_vec{1}.trajt = [0 4];
p0_ref = vehi_vec{1}.vehiSt.x0(1:3);
v0_ref = vehi_vec{1}.vehiSt.x0(4:6);
p1_ref = [5.5;-5.5;.5];
v1_ref = [0;0;0];
a1_ref = [0;0;0];
j1_ref = [0;0;0];
p2_ref = [0;0;2];
v2_ref = [0;0;0];
a2_ref = [0;0;0];
j2_ref = [0;0;0];
pts_ref{1} = zeros(3*4,length(vehi_vec{1}.trajt));
pts_ref{1}(1:6,1) = [p0_ref;v0_ref];
%pts_ref{1}(:,2) = [p1_ref; v1_ref; a1_ref; j1_ref];
pts_ref{1}(:,end) = [p2_ref; v2_ref; a2_ref; j2_ref];
vehi_vec{1}.pts_ref = pts_ref{1};

vehi_vec{2}.trajt = [0 4];
p0_ref = vehi_vec{2}.vehiSt.x0(1:3);
v0_ref = vehi_vec{2}.vehiSt.x0(4:6);
p1_ref = [0;0;1];
v1_ref = [0;0;0];
a1_ref = [0;0;0];
j1_ref = [0;0;0];
p2_ref = [0;0;2];
v2_ref = [0;0;0];
a2_ref = [0;0;0];
j2_ref = [0;0;0];
pts_ref{2} = zeros(3*4,length(vehi_vec{2}.trajt));
pts_ref{2}(1:6,1) = [p0_ref;v0_ref];
%pts_ref{2}(:,2) = [p1_ref; v1_ref; a1_ref; j1_ref];
pts_ref{2}(:,end) = [p2_ref; v2_ref; a2_ref; j2_ref];
%pts_ref{2}(:,end) = [p2_ref; v2_ref];
vehi_vec{2}.pts_ref = pts_ref{2};
for ind=1:numVehi
    vehi = vehi_vec{ind};

    % === Init Trajectories === %
    trajFunc = vehi.vehiSt.trajFunc;
    traj = trajFunc(vehi.trajt, vehi.pts_ref); % Compute Trajectory based on waypoints
    traj = [traj, traj(:,end).*ones(size(traj,1),knum-size(traj,2))]; % Pad Trajectory vec with last value
    vehi.traj = [traj; 0.*ones(1,knum); 0.*ones(1,knum)]; % Yaw and dYaw
    % if ind~=1
    %     vehi.traj = zeros(size(vehi.traj,1),0);
    % end

    % === Init Simulation Vectors === %
    vehi.x = zeros(vehi.vehiSt.nx,knum);
    vehi.x(:,1) = [vehi.traj(1:6,1); vehi.vehiSt.x0(7:end)];
    vehi.u = zeros(vehi.vehiSt.nu,knum);
    vehi.u(:,1) = vehi.vehiSt.u0;

    vehi.alph = zeros(vehi.vehiSt.nx^2,knum/recalc_step);
    vehi.alph(:,1) = reshape(vehi.vehiSt.alph0,[],1);

    vehi.y = zeros(1,knum/recalc_step);
    vehi.y(1) = vehi.vehiSt.ysim0;

    vehi_vec{ind} = vehi;
end

% === Start Figure === %
enable_disp = 0;
fig=figure('WindowState','maximized');
hold on;
grid on
axis equal
view(3)
xlabel('$x$','Interpreter','latex');
ylabel('$y$','Interpreter','latex');
zlabel('$z$','Interpreter','latex');
for ind=1:numVehi
    pts_ref = vehi_vec{ind}.pts_ref;
    vehi_vec{ind}.prmtpPatch = {};
    vehi_vec{ind}.wayptsPlot = plot3(pts_ref(1,:),pts_ref(2,:),pts_ref(3,:),'xr','DisplayName','target');
    %vehi_vec{ind}.trajPlot = plot3(vehi_vec{ind}.traj(1,:), vehi_vec{ind}.traj(2,:), vehi_vec{ind}.traj(3,:) ,'-.k','DisplayName','r_T');
    vehi_vec{ind}.bodyPlot = vehi_vec{ind}.vehiSt.bodyPlot();
    vehi_vec{ind}.pathPlot = animatedline(vehi_vec{ind}.x(1,1), vehi_vec{ind}.x(2,1), vehi_vec{ind}.x(3,1),'LineWidth',2.5, 'Color','b','LineStyle','-','DisplayName','Path');
end

% === Simulation === %
telapse_total = tic;
for k=1:knum-1 
    disp(['=== Computing ', num2str(k*Ts), '/', num2str(simT), 's ===']);
    for kVehi=1:numVehi
        vehi = vehi_vec{kVehi};
        ops = vehi.ops;

        if ~mod(k+recalc_step-1,recalc_step)
            telapse_k = tic;
            disp('----------------');
            disp(['Computing ', num2str(k), ' time interval ']);

            % === Predict === %
            [vehi.xpred, vehi.upred] = calcPath(vehi.x(:,k), vehi, [k:knum-1], Ts);
            %prmtp_ind = 1:inner_step:size(vehi.xpred,2);
            if discnum*disc_step >= size(vehi.xpred,2)
                prmtp_indend = size(vehi.xpred,2)-1;
            else
                prmtp_indend = discnum*disc_step;
            end
            prmtp_ind = 1:disc_step:prmtp_indend;
            [vehi.alphpred, vehi.ypred] = calcParametopes(vehi.xpred(:,prmtp_ind), vehi.upred(:,prmtp_ind), vehi, Ts*disc_step);

            disp(['Computing time for interval: ', num2str(toc(telapse_k)), 's']);

            % === Recalculate Trajectory (for Defender only) === %
            if kVehi~=1
                % Assumes models have at least position and velocity, but no acceleration or snap
                % Maybe use vehiSt function to get pts
                % Set waypoint as a point in a line between the current
                % point and the invader's position, at twice the distance
                dist = -vehi_vec{1}.x(1:3,k) + vehi.x(1:3,k);
                d_norm = norm(dist);
                dist_n = 1/d_norm.*dist;

                xstart = [vehi.x(1:6,k); vehi.traj(7:end-2,k)]; % Ignore Yaw for now...
                %xend = [vehi_vec{1}.x(1:6,k); vehi_vec{1}.traj(7:end-2,k)];
                %xend = [vehi_vec{1}.x(1:3,k); zeros(size(xstart,1)-3,1)];
                %xend = [2*d_norm*dist_n; zeros(size(xstart,1)-3,1)];
                xend = [vehi_vec{1}.x(1:3,k); dist_n*10; dist_n*1; dist_n*.1];
                traj_aux = recalcTraj(xstart, xend, k, Ts, knum, vehi);
                vehi.traj = [vehi.traj(:,1:k-1), traj_aux];
            end

            % === Update Plot === %
            %vehi.trajPlot = plotTraj(vehi.xpred, vehi.trajPlot);
            if enable_disp
                vehi.prmtpPatch = plotParametope(vehi.ypred,vehi.alphpred,vehi.xpred(:,prmtp_ind),vehi.vehiSt.nx,ops,vehi.prmtpPatch);
            end

        end
       
        % === Real-Time === %
        % Add disturbances!
        [x_aux, u_aux] = calcPath(vehi.x(:,k), vehi, k, Ts);
        vehi.x(:,k+1) = x_aux(:,end);
        vehi.u(:,k) = u_aux(:,end);
        %[vehi.alph(:,k+1), vehi.y(:,k+1)] = calcParametopes(vehi.x(:,k), vehi.u(:,k), vehi, Ts);

        % === Update Plot === %
        plotBodyFunc = vehi.vehiSt.plotBodyFunc; 
        vehi.bodyPlot = plotBodyFunc(vehi.x(:,k), vehi.bodyPlot);
        addpoints(vehi.pathPlot,vehi.x(1,k), vehi.x(2,k), vehi.x(3,k));

        vehi_vec{kVehi} = vehi;

        %legend({'target', 'r_T', 'Path'});
    end
    drawnow;
    pause(Ts);
    frame(k) = getframe(fig);
end
hold off;
disp(['Total sim computing time: ', num2str(toc(telapse_total)), 's']);

vid=VideoWriter('test_v3.mp4','MPEG-4');
vid.FrameRate=1/Ts;
open(vid);
writeVideo(vid,frame);
close(vid);

function vehi_vec = local_init_all_vehi(Ts)
%INIT_ALL_VEHI outputs cell vector of vehiStruct
%   Detailed explanation goes here

vehi_vec = cell(1,2);

% Define Inital Positions
p10 = [5;-5;0];
% p20 = [-5;5;0];
p20 = [4;-4;0];
x0 = {p10, p20};

vehi_vec{1}.vehiSt = init_QuadDroneMellinger(x0{1},Ts);
vehi_vec{2}.vehiSt = init_QuadDroneMellinger(x0{2},Ts);

% vehi_vec{1}.ctrlSt = init_NMPCctrl(vehi_vec{1}.vehiSt.nx,5,Ts,vehi_vec{1}.vehiSt);

% vehi_vec{1}.ctrlSt = init_MellingerCtrl(4,vehi_vec{1}.vehiSt);
% vehi_vec{2}.ctrlSt = init_NMPCctrl(vehi_vec{1}.vehiSt.nx,30,Ts,vehi_vec{2}.vehiSt);
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

function traj_new = recalcTraj(xstart,xend,kcurr,Ts,knum,vehi)
% Update trajt and pts_ref
tcurr = kcurr*Ts;
vehi.trajt = [tcurr, vehi.trajt(2:end)];
vehi.pts_ref = [xstart, xend];

trajFunc = vehi.vehiSt.trajFunc;
traj = trajFunc(vehi.trajt, vehi.pts_ref); % Compute Trajectory based on waypoints
%knum = size(vehi.traj,2);
%traj = [vehi.traj(:,1:(kcurr-1)), traj, traj(:,end).*ones(size(traj,1),knum-size(traj,2)-kcurr+1)]; % Pad Trajectory vec with last value
traj = [traj, traj(:,end).*ones(size(traj,1),knum-size(traj,2)-kcurr+1)]; % Pad Trajectory vec with last value
traj_new = [traj; 0.*ones(1,size(traj,2)); 0.*ones(1,size(traj,2))]; % Yaw and dYaw
            
end

function [alph_vec, y] = calcParametopes(x, u, vehi, Ts)

alph_vec = zeros(vehi.vehiSt.nx^2,size(x,2));
alph_vec(:,1) = reshape(vehi.vehiSt.alph0,[],1);

y = zeros(1,size(x,2));
y(1) = vehi.vehiSt.ysim0;

for k=1:size(x,2)-1
    % === Propagate Parametope === %
    jacFunc = vehi.vehiSt.jacFunc;
    Dfxo = jacFunc(x(:,k), u(:,k));
    
    alph = reshape(alph_vec(:,k),vehi.vehiSt.nx,[]);
    alphk_1 = runge_kutta4(@(y) prmtpU(y,-Dfxo,vehi.vehiSt.nx),Ts,alph);
    alph_vec(:,k+1) = reshape(alphk_1,[],1);
    
    vbds_ext(:,1) = vehi.vehiSt.vbds_init(:,1)+[x(:,k); u(:,k)];
    vbds_ext(:,2) = vehi.vehiSt.vbds_init(:,2)+[x(:,k); u(:,k)];
    
    % Calculate c just once outside of runge-kutta
    [alph_err,parametope_c] = prmtpEllips_c(vehi.vehiSt.rStrt,vehi.vehiSt.symvars,vbds_ext,-Dfxo,alph_vec(:,k),vehi.vehiSt.nx,vehi.vehiSt.rCode,vehi.vehiSt.rHasVarMx,vehi.vehiSt.rConstDfMx);
    y(k+1) = runge_kutta4(@(yr) prmtpE(yr,parametope_c^-2), Ts, y(k));
    % y(k+1) = runge_kutta4(@(yr) prmtpE(yr,parametope_c), Ts, y(k));
end
end

function [x, u] = calcPath(x0, vehi, t, Ts)
% Run trajectory prediction of vehi from x0 throughout t

x = zeros(vehi.vehiSt.nx,length(t)+1);
x(:,1) = x0;
u = zeros(vehi.vehiSt.nu,length(t));

for k=t
    kin = k-t(1)+1;
    % === Calculate Control action === %
    ctrlFunc = vehi.vehiSt.contrlFunc;
    u(:,kin) = ctrlFunc(x(:,kin), vehi.traj(:,k));

    % === Calculate System === %
    dynFunc = vehi.vehiSt.dynFunc;
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
        alph = reshape(alph_vec(:,k),nx,[]);
        alph_pre = (1*y(k))^-1.*alph'*alph;
        X.G = inv(sqrtm(alph_pre));
        X.c = x(:,k);
        X.A = zeros(0,nx);
        X.b = zeros(0,1);
        X.idx = [nx];
        X.type = [2];
        [Fc,p] = compileCCG(X);
        %plot(Fc,p(1:3),ops.plot.edgecolor,300,ops);
        x_opt = plot(Fc,p(1:3),ops.plot.edgecolor,50,ops);
        x_opt = x_opt{end};
        K = convhulln(x_opt');
    
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

