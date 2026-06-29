function [ctrlSt] = init_MellingerCtrl(tf,vehi)
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

auxSt.ddr = zeros(3,1);
auxSt.dddr = zeros(3,1);
auxSt.t0 = 0;

ctrlSt.auxSt = auxSt;

% === Specific Function Handles === %
%ctrlSt.trajFunc = @(t,pts) quadDrone_traj(t,pts,ctrlSt);
ctrlSt.ctrlFunc = @(x0,xs,auxSt) quadController(x0, xs, auxSt, tf, vehi);
end

function varargout = trajNxtK(x0_cell, xs_cell, t0, tf, Ts)
    % calculate trajectory and return next point

    r0 = x0_cell{1}; rd = xs_cell{1};
    dr0 = x0_cell{2}; drd = xs_cell{2};
    ddr0 = x0_cell{3}; ddrd = xs_cell{3};
    dddr0 = x0_cell{4}; dddrd = xs_cell{4};

    % tf is some function of the distance between points
    %tf = floor(norm(rd-r0)/2/Ts)*Ts; % Assume average speed
    t1 = t0+Ts;
    tf = max(tf,t1);
    w_aux = trajectory3Dv2([t0;tf], Ts, [r0,rd], [dr0,drd], [ddr0,ddrd], [dddr0,dddrd]);
    w_aux = w_aux{1};

    traj_order = 4;
    prod_vec = ones(1,2*traj_order);
    for k=1:nargout
        pow_vec = 2*traj_order-k:-1:1;
        rx = polyval(prod_vec.*w_aux(1,end:-1:k),t1);
        ry = polyval(prod_vec.*w_aux(2,end:-1:k),t1);
        rz = polyval(prod_vec.*w_aux(3,end:-1:k),t1);
        
        prod_vec = prod_vec(1:length(pow_vec)).*pow_vec;
        r = [rx; ry; rz]; % Coordinate of trajectory position
        varargout{k} = r;
    end
end

function [u, auxSt] = quadController(x0, xs, auxSt, tf, vehi)
    % x0 : current state [x0; alph0]
    % xs : desired state [ps; vs; ~]
    % auxSt : contains acceleration and jerk

    Ts = vehi.Ts;
    nx = vehi.nx;

    r = xs(1:3);
    dr = xs(4:6);
    % In the future, let these be unknowns
    ddr = zeros(3,1);
    dddr = zeros(3,1);
    yaw_ref = 0;
    dyaw_ref = 0;

    
    x0_cell = {x0(1:3), x0(4:6), auxSt.ddr, auxSt.dddr};
    xs_cell = {r,dr,ddr,dddr};
    t0 = auxSt.t0;
    [rNxt,drNxt,ddrNxt,dddrNxt] = trajNxtK(x0_cell, xs_cell, t0, tf, Ts);
    auxSt.ddr = ddrNxt; auxSt.dddr = dddrNxt; auxSt.t0 = t0+Ts;

    u = mellinger_controller(auxSt.KK,rNxt,drNxt,ddrNxt,dddrNxt,yaw_ref,dyaw_ref,x0(1:nx),vehi.I,vehi.M);

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