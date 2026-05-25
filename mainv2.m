% Main file of thesis project
% Runs all generalized functions for any number of vehichles
% Offline, initializes cell vector of vehiStruct
% Then initializes parametope propagation variables
% Online, runs specific vehichle and parametope dynamics 
% Every inner_step runs predicted trajectory and parametopes

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
inner_step = 20;
enable_prmtp = 1;

% === Init vehiStruct === &
vehi_vec = init_all_vehi(Ts); % Cell vec with vehiSimStruct
numVehi = length(vehi_vec);
% vehi_vec{1} is assumed to be the Invader

trajt{1} = [0 1 4];
p0_ref = [5;-5;0];
v0_ref = [0;0;0];
p1_ref = [5.5;-5.5;.5];
v1_ref = [0;0;0];
a1_ref = [0;0;0];
j1_ref = [0;0;0];
p2_ref = [0;0;2];
v2_ref = [0;0;0];
a2_ref = [0;0;0];
j2_ref = [0;0;0];
pts_ref{1} = zeros(3*4,length(trajt{1}));
pts_ref{1}(1:6,1) = [p0_ref;v0_ref];
pts_ref{1}(:,2) = [p1_ref; v1_ref; a1_ref; j1_ref];
pts_ref{1}(:,end) = [p2_ref; v2_ref; a2_ref; j2_ref];

trajt{2} = [0 4];
p0_ref = [-5;5;0];
v0_ref = [0;0;0];
p1_ref = [0;0;1];
v1_ref = [0;0;0];
a1_ref = [0;0;0];
j1_ref = [0;0;0];
p2_ref = [0;0;2];
v2_ref = [0;0;0];
a2_ref = [0;0;0];
j2_ref = [0;0;0];
pts_ref{2} = zeros(3*4,length(trajt{2}));
pts_ref{2}(1:6,1) = [p0_ref;v0_ref];
%pts_ref{2}(:,2) = [p1_ref; v1_ref; a1_ref; j1_ref];
pts_ref{2}(:,end) = [p2_ref; v2_ref; a2_ref; j2_ref];
for ind=1:numVehi
    vehi = vehi_vec{ind};

    % === Init Trajectories === %
    trajFunc = vehi.vehiSt.trajFunc;
    traj = trajFunc(trajt{ind}, pts_ref{ind}); % Compute Trajectory based on waypoints
    traj = [traj, traj(:,end).*ones(size(traj,1),knum-size(traj,2))]; % Pad Trajectory vec with last value
    vehi.traj = [traj; 0.*ones(1,knum); 0.*ones(1,knum)]; % Yaw and dYaw

    % === Init Simulation Vectors === %
    vehi.x = zeros(vehi.vehiSt.nx,knum);
    vehi.x(:,1) = [vehi.traj(1:6,1); vehi.vehiSt.x0(7:end)];
    vehi.u = zeros(vehi.vehiSt.nu,knum);
    vehi.u(:,1) = vehi.vehiSt.u0;

    vehi.alph = zeros(vehi.vehiSt.nx^2,knum/inner_step);
    vehi.alph(:,1) = reshape(vehi.vehiSt.alph0,[],1);

    vehi.y = zeros(1,knum/inner_step);
    vehi.y(1) = vehi.vehiSt.ysim0;

    vehi_vec{ind} = vehi;
end

% === Start Figure === %
enable_disp = 0;
fig=figure;
hold on;
grid on
axis equal
view(3)
xlabel('$x$','Interpreter','latex');
ylabel('$y$','Interpreter','latex');
zlabel('$z$','Interpreter','latex');
for ind=1:numVehi
    p_vec = vehi_vec{ind}.p_vec;

    Xvec = [p_vec(1:3:12)];
    Yvec = [p_vec(2:3:12)];
    Zvec = [p_vec(3:3:12)];
    vehi_vec{ind}.wayptsPlot = plot3(pts_ref{ind}(1,:),pts_ref{ind}(2,:),pts_ref{ind}(3,:),'xr','DisplayName','target');
    %vehi_vec{ind}.trajPlot = plot3(vehi_vec{ind}.traj(1,:), vehi_vec{ind}.traj(2,:), vehi_vec{ind}.traj(3,:) ,'r','DisplayName','r_T');
    vehi_vec{ind}.bodyPlot = patch(Xvec, Yvec, Zvec, [0 0 1],'DisplayName','Drone Body');
    vehi_vec{ind}.pathPlot = animatedline(vehi_vec{ind}.x(1,1), vehi_vec{ind}.x(2,1), vehi_vec{ind}.x(3,1),'LineWidth',2.5, 'Color','b','LineStyle','-','DisplayName','Path');
end

% === Simulation === %
telapse_total = tic;
for k=1:knum-1 
    disp(['=== Computing ', num2str(k*Ts), '/', num2str(simT), 's ===']);
    for kVehi=1:numVehi
        vehi = vehi_vec{kVehi};
        ops = vehi.ops;

        % === Predict === %
        if ~mod(k+inner_step-1,inner_step)
            telapse_k = tic;
            disp('----------------');
            disp(['Computing ', num2str(k), ' time interval ']);

            [vehi.xpred, vehi.upred] = calcPath(vehi.x(:,k), vehi, [k:knum-1], Ts, inner_step);
            prmtp_ind = 1:inner_step:size(vehi.xpred,2);
            [vehi.alphpred, vehi.ypred] = calcParametopes(vehi.xpred(:,prmtp_ind), vehi.upred(:,prmtp_ind), vehi, Ts*inner_step);

            disp(['Computing time for interval: ', num2str(toc(telapse_k)), 's']);

            if enable_disp
                plot3(vehi.xpred(1,:), vehi.xpred(2,:), vehi.xpred(3,:), '-.k');
                plotParametope(vehi.ypred,vehi.alphpred,vehi.xpred(:,prmtp_ind),vehi.vehiSt.nx,ops,fig);
            end

        end
       
        % === Real-Time === %
        [x_aux, u_aux] = calcPath(vehi.x(:,k), vehi, k, Ts, inner_step);
        vehi.x(:,k+1) = x_aux(:,end);
        vehi.u(:,k) = u_aux(:,end);

        % === Update Plot === %
        bodyFunc = vehi.vehiSt.bodyFunc; 
        p_vec = bodyFunc(vehi.x(:,k), vehi.p_vec);
        Xvec = [p_vec(1:3:12)];
        Yvec = [p_vec(2:3:12)];
        Zvec = [p_vec(3:3:12)];
        vehi.bodyPlot.XData = Xvec;
        vehi.bodyPlot.YData = Yvec;
        vehi.bodyPlot.ZData = Zvec;
        
        addpoints(vehi.pathPlot,vehi.x(1,k), vehi.x(2,k), vehi.x(3,k));

        vehi_vec{kVehi} = vehi;

        %legend({'target', 'r_T', 'Path'});
        drawnow;
        pause(Ts);
    end
end
hold off;
disp(['Total sim computing time: ', num2str(toc(telapse_total)), 's']);

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
end
end

function [x, u] = calcPath(x0, vehi, t, Ts, inner_step)
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

function plotParametope(y,alph_vec,x,nx,ops,fig)
for k=1:length(y)
    alph = reshape(alph_vec(:,k),nx,[]);
    alph_pre = (1/10*y(k))^-1.*alph;
    X.G = inv(sqrtm(alph_pre));
    X.c = x(:,k);
    X.A = zeros(0,nx);
    X.b = zeros(0,1);
    X.idx = [nx];
    X.type = [2];
    [Fc,p] = compileCCG(X);
    plot(Fc,p(1:3),ops.plot.edgecolor,300,ops);
    %drawnow;
end
end
