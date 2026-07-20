function [ctrlSt] = init_MellingerCtrl(pred_hor,tf,vehi)
%init_mellingerCtrl creates mellinger ctrlSt ( for parametope)
%
% ctrlSt : struct (struct with components of prmtp Mellinger)
% MAIN (interfaces with main.m)
%   .auxSt : contains auxiliary variables unique to this controller that
%       may need to be modified online
%   .contrlFunc
% AUXILIARY (does not interface with main.m)

% === Initialize Auxiliary variables for mellinger controller === %
K_val = 1;
KK.Kp = 1*K_val.*eye(3); % dddp
KK.Kv = 1*K_val.*eye(3); % ddp
KK.KR = 1*K_val.*eye(3);   % dp
KK.Kom = 1*K_val.*eye(3);        % p
auxSt.KK = KK;

auxSt.dr = zeros(3,1);
auxSt.ddr = zeros(3,1);
auxSt.dddr = zeros(3,1);
auxSt.t0 = 0;

auxSt.error_track = zeros(1,1);

ctrlSt.auxSt = auxSt;

% === Specific Function Handles === %
%ctrlSt.trajFunc = @(t,pts) quadDrone_traj(t,pts,ctrlSt);
ctrlSt.ctrlFunc = @(x0,xs,auxSt) quadController(x0, xs, auxSt, pred_hor, tf, vehi);
end

function [u,rNxt] = trajNxtK(x0,W,tnxt,vehi,auxSt)
nx = vehi.nx;

for k=1:4
    rNxt{k} = W*[dpolyN(tnxt,8,k-1)];
end

yaw_ref= 0;
dyaw_ref= 0;
u = mellinger_controller(auxSt.KK,rNxt{1},rNxt{2},rNxt{3},rNxt{4},yaw_ref,dyaw_ref,x0(1:nx),vehi.I,vehi.M);

% Clamp u to bounds
u(1) = min(max(u(1),vehi.Fmin),vehi.Fmax);
tauBnorm = norm(u(2:3));
if tauBnorm == 0
    u(2:4) = [0;0;0];
elseif tauBnorm > vehi.tauBmax
    scale = vehi.tauBmax/tauBnorm;
    u(2:4) = scale.*u(2:4);
end
end

function [u, auxSt] = quadController(x0, xs, auxSt, pred_hor, tf, vehi)
    % x0 : current state [x0; alph0]
    % xs : desired state [ps; vs; ~]
    % auxSt : contains acceleration and jerk
    xs = xs(:,end);

    Ts = vehi.Ts;
    nx = vehi.nx;
    nx_full = nx+nx^2;

    rd = xs(1:3);
    drd = xs(4:6);
    % In the future, let these be unknowns
    ddrd = zeros(3,1);
    dddrd = zeros(3,1);
 
    r0 = x0(1:3); dr0 = x0(4:6); ddr0 = auxSt.ddr; dddr0 = auxSt.dddr;
    t0 = auxSt.t0;
    t1 = t0+Ts;
    tf = max(tf,t1);
    W = trajectory3Dv2([t0;tf], Ts, [r0,rd], [dr0,drd], [ddr0,ddrd], [dddr0,dddrd]);
    W = W{1};

    xpred = zeros(length(x0),pred_hor);
    xpred(:,1) = x0;
    tnxt = t0;
    for ind=1:pred_hor
        tnxt = tnxt+Ts;
        [u(:,ind),rNxt] = trajNxtK(xpred(:,ind),W,tnxt,vehi,auxSt);
        xpred(:,ind+1) = runge_kutta4(@(y) vehi.dynFunc(y,u(:,ind)), Ts, xpred(:,ind));
        if ind==1
            auxSt.dr = rNxt{2}; auxSt.ddr = rNxt{3}; auxSt.dddr = rNxt{4};
            auxSt.t0 = tnxt;
        end
    end

end