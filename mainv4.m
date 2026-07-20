% Mainv3 file of thesis project
% Runs all generalized functions for any number of vehichles
% Offline, initializes cell vector of vehiStruct
% Then initializes parametope propagation variables
% Online, runs specific vehichle and parametope dynamics 
% Every inner_step runs predicted trajectory and parametopes
% Defender tries to intercept invader
% Using NMPC to align parametopes
% Initialize vehi_vec, containing vehiSt (model) and ctrlSt (controller)

close all;
clearvars;

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
simT = 5; % Sim time seconds
Ts = 0.01;
knum = simT/Ts;
%calc_traj_dt = 0.1; % How often trajectory is recalculated [s]
calc_traj_dt = Ts;
recalc_step = calc_traj_dt/Ts;
pred_hor = 60;
discnum = 5;
disc_step = floor(pred_hor/discnum);
disc_dt = disc_step*Ts; % Time step between Reachability discs

enable_prmtp = 1;
enable_disp = 1;
make_vid = 0;

% === Init vehiStruct === &
vehi_vec = init_all_vehi_local(Ts,pred_hor); % Cell vec with vehiSimStruct
numVehi = length(vehi_vec);
% vehi_vec{1} is assumed to be the Invader 

for ind=1:numVehi
    vehi = vehi_vec{ind};
    nx = vehi.vehiSt.nx;
    nx_full = nx+nx^2;

    if ind ~= 1
        %pts_ref = vehi_vec{1}.vehiSt.x0;
        %vehi.pts_ref = [vehi.vehiSt.x0, pts_ref];
        vehi.pts_ref = repmat(vehi.vehiSt.x0,1,pred_hor+1);
    else
        % === Define Critical Area === %
        p1_ref = [0;0;0];
        v1_ref = [0;0;0];
        pts_ref = [p1_ref; v1_ref; vehi.vehiSt.x0(3*2+1:end)];
        vehi.pts_ref = [vehi.vehiSt.x0, repmat(pts_ref,1,1)];
        %vehi.pts_ref = [vehi.vehiSt.x0, repmat(pts_ref,1,pred_hor)];
    end

    % === Init Simulation Vectors === %
    vehi.x = zeros(nx_full,knum);
    vehi.x(:,1) = vehi.pts_ref(:,1);
    vehi.u = zeros(vehi.vehiSt.nu,knum);
    vehi.u(:,1) = vehi.vehiSt.u0;

    vehi.alph = zeros(vehi.vehiSt.nx^2,knum/recalc_step);
    vehi.alph(:,1) = reshape(vehi.vehiSt.alph0,[],1);

    vehi.y = zeros(1,knum/recalc_step);
    vehi.y(1) = vehi.vehiSt.ysim0;

    vehi_vec{ind} = vehi;
end

%vehiPred = vehi_vec{1}; % Auxilliary vehi for prediction
vehi_vecPred = vehi_vec; % Auxiliary vehi_vec for prediction

% === Start Figure === %
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
k=1;
while k<knum && norm(vehi_vec{1}.x(1:3,k)-vehi_vec{2}.x(1:3,k))>=5e-1 
    disp(['=== Computing ', num2str(k), '/', num2str(knum), ' interval ===']);
    disp(['t = ', num2str((k-1)*Ts),'s']);
    telapse_k = tic;
    for kVehi=1:numVehi
        vehi = vehi_vec{kVehi};
        vehiPred = vehi_vecPred{kVehi};
        vehiSt = vehi.vehiSt;
        nx = vehiSt.nx;
        nx_full = nx+nx^2;
        ops = vehi.ops;

        if ~mod(k+recalc_step-1,recalc_step)

            if kVehi == 1
                % === Predict === %
                % pred_t = knum-k-1;
                % vehi.xpred = zeros(nx_full,pred_t+1);
                % vehi.xpred(:,1) = vehi.x(:,k);
                % vehi.upred = zeros(vehiSt.nu,pred_t);
                vehiPred.x(1:nx,k) = vehi.x(1:nx,k);
                vehiPred.x(nx+1:end,k) = vehiPred.x(nx+1:end,1);
                vehiPred.ctrlSt.auxSt = vehi.ctrlSt.auxSt;
                for pred_k=k:min(k+pred_hor,knum)
                    %[vehi.xpred(:,pred_k+1), vehi.upred(:,pred_k), vehi.ctrlStpred.auxSt] = calcNxt(vehi.xpred(:,pred_k), vehi.pts_ref(:,end) , vehi, Ts);
                    [vehiPred.x(:,pred_k+1), vehiPred.u(:,pred_k), vehiPred.ctrlSt.auxSt] = calcNxt(vehiPred.x(:,pred_k), vehiPred.pts_ref(:,2:end) , vehiPred, Ts);
                end
                % Test parametope plot
                %vehi.prmtpPatch = plotParametope(vehiSt.ysim0,vehiPred.x(nx+1:end,pred_k),vehiPred.x(1:nx,pred_k),nx,ops,vehi.prmtpPatch);
                %vehi.prmtpPatch = plotParametope(vehiSt.ysim0.*ones(1,length(k:pred_k)),vehiPred.x(nx+1:end,k:pred_k),vehiPred.x(1:nx,k:pred_k),nx,ops,vehi.prmtpPatch);
            
                % if discnum*disc_step >= size(vehi.xpred,2)
                %     prmtp_indend = size(vehi.xpred,2)-1;
                % else
                %     prmtp_indend = discnum*disc_step;
                % end
                % prmtp_ind = 1:disc_step:prmtp_indend;
                % [vehi.alphpred, vehi.ypred] = calcParametopes(vehi.xpred(:,prmtp_ind), vehi.upred(:,prmtp_ind), vehi, Ts*disc_step);
    
                % === Recalculate Trajectory (for Defender only) === %
            else
                %vehi.pts_ref(:,end) = vehi_vec{1}.x(:,k);
                vehi.pts_ref(:,2:(pred_k-k+1)) = vehi_vecPred{1}.x(:,k:pred_k-1);
            end

            % === Update Plot === %
            %vehi.trajPlot = plotTraj(vehi.xpred, vehi.trajPlot);
            if enable_disp && kVehi==1
                prmtp_ind = k+1:disc_step:pred_k;
                vehi.prmtpPatch = plotParametope(vehiSt.ysim0.*ones(1,length(prmtp_ind)),vehiPred.x(nx+1:end,prmtp_ind),vehiPred.x(1:nx,prmtp_ind),nx,ops,vehi.prmtpPatch);
                %vehi.prmtpPatch = plotParametope(vehi.ypred,vehi.alphpred,vehi.xpred(:,prmtp_ind),nx,ops,vehi.prmtpPatch);
            end

        end
       
        % === Real-Time === %
        % Add disturbances!
        [vehi.x(:,k+1), vehi.u(:,k), vehi.ctrlSt.auxSt] = calcNxt(vehi.x(:,k), vehi.pts_ref(:,2:end), vehi, Ts);

        % === Update Plot === %
        plotBodyFunc = vehiSt.plotBodyFunc; 
        vehi.bodyPlot = plotBodyFunc(vehi.x(:,k), vehi.bodyPlot);
        addpoints(vehi.pathPlot,vehi.x(1,k), vehi.x(2,k), vehi.x(3,k));

        vehi_vec{kVehi} = vehi;
        vehi_vecPred{kVehi} = vehiPred;

        %legend({'target', 'r_T', 'Path'});
    end
    aux = v2m(vehi_vec{1,1}.x(18+1:end,pred_k),18);
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
    vid=VideoWriter('test_v4.mp4','MPEG-4');
    vid.FrameRate=1/Ts;
    open(vid);
    writeVideo(vid,frame);
    close(vid);
end

function vehi_vec = init_all_vehi_local(Ts,pred_hor)
%INIT_ALL_VEHI outputs cell vector of vehiStruct
%   Detailed explanation goes here

vehi_vec = cell(1,2);

% Define Inital Positions
p10 = [5;-5;5];
% p20 = [-5;5;0];
p20 = [-1;1;0];
x0 = {p10, p20};

vehi_vec{1}.vehiSt = init_QuadDroneNMPC(x0{1},Ts);
vehi_vec{2}.vehiSt = init_QuadDroneNMPC(x0{2},Ts);

% vehi_vec{1}.ctrlSt = init_NMPCctrl(vehi_vec{1}.vehiSt.nx,5,Ts,vehi_vec{1}.vehiSt);

vehi_vec{1}.ctrlSt = init_MellingerCtrl(2,vehi_vec{1}.vehiSt);
%vehi_vec{2}.ctrlSt = init_NMPCctrl_v2(vehi_vec{1}.vehiSt.nx,pred_hor,Ts,vehi_vec{2}.vehiSt);
%vehi_vec{2}.ctrlSt = init_NMPCctrl(vehi_vec{1}.vehiSt.nx,pred_hor,Ts,vehi_vec{2}.vehiSt);
vehi_vec{2}.ctrlSt = init_MellingerCtrl_v2(pred_hor,Ts,vehi_vec{2}.vehiSt);
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

for k=1:size(x,2)-1
    % === Propagate Parametope === %
    jacFunc = vehiSt.jacFunc;
    Dfxo = jacFunc(x(:,k), u(:,k));
    
    alph = reshape(alph_vec(:,k),nx,[]);
    alphk_1 = runge_kutta4(@(y) prmtpU(y,-Dfxo,nx),Ts,alph);
    alph_vec(:,k+1) = reshape(alphk_1,[],1);
    
    vbds_ext(:,1) = vehiSt.vbds_init(:,1)+[x(:,k); u(:,k)];
    vbds_ext(:,2) = vehiSt.vbds_init(:,2)+[x(:,k); u(:,k)];
    
    % Calculate c just once outside of runge-kutta
    [alph_err,parametope_c] = prmtpEllips_c(vehiSt.rStrt,vehiSt.symvars,vbds_ext,-Dfxo,alph_vec(:,k),nx,vehiSt.rCode,vehiSt.rHasVarMx,vehiSt.rConstDfMx);
    y(k+1) = runge_kutta4(@(yr) prmtpE(yr,parametope_c^-2), Ts, y(k));
end
end

function [xk_1, u, auxSt] = calcNxt(xk, xd, vehi, Ts)
% Calculate next step for vehi from x0 to xd
vehiSt = vehi.vehiSt;
ctrlSt = vehi.ctrlSt;

% === Calculate Control action === %
[u,auxSt] = ctrlSt.ctrlFunc(xk, xd, ctrlSt.auxSt);

% === Calculate System === %
dynFunc = vehiSt.dynFunc;
xk_1 = runge_kutta4(@(y) dynFunc(y,u), Ts, xk);

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
