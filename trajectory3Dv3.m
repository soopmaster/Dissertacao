function varargout = trajectory3Dv3(varargin)
%trajectory3D Create trajectory that passes through given waypoints at times t with
%given velocity
% Given two position/velocity waypoints, automatically compute 7th order
% trajectory (has acceleration) that tries to be as fast as possible
% 1) Compute 3rd order trajectory from waypoints
% 2) Sample 4 points from this trajectory with shriking delay from each
% other
% 3) Calculate direction of velocity at these points and set the magnitude
% as growing from v0 to vf
% 4) Use these points as parameters to compute final 7th order trajectory

% waypt = [pt1x pt2x ... ptNx; pt1y pt2y ... ptNy; pt1z pt2z ... ptNz]
% v = [v1x v2x ... vNx; v1y v2y ... vNy; v1z v2z ... vNz]

% w_star{i} = [wx0, wx1, wx2, wx3;
%              wy0, wy1, wy2, wy3;
%              wz0, wz1, wz2, wz3;

traj_order = nargin-1; % Defines whether it's got velocity, acceleration, etc.

dt = varargin{1};

t0 = 0;
tf = 1;
w_aux = trajectory3Dv2([t0;tf], dt, varargin{2}, varargin{3});

t_aux = linspace(t0,tf,4);
rx = polyval(w_aux(1,end:-1:1),t_aux);
ry = polyval(w_aux(2,end:-1:1),t_aux);
rz = polyval(w_aux(3,end:-1:1),t_aux);
r_aux = [rx;ry;rz];

pow_vec = 2*traj_order-1:-1:1;
drx = polyval(pow_vec.*w_aux(1,end:-1:2),t_aux);
dry = polyval(pow_vec.*w_aux(2,end:-1:2),t_aux);
drz = polyval(pow_vec.*w_aux(3,end:-1:2),t_aux);
dr_aux = [drx;dry;drz];

waypt_num = size(t_aux,2);

% Continue from here

time_cell = cell(waypt_num-1,1);
for i=1:waypt_num-1
    time_cell{i} = t(i):dt:t(i+1);
end

% matrix to store the wstar for each segment
w_starx= zeros(waypt_num-1,2*traj_order);
w_stary= zeros(waypt_num-1,2*traj_order);
w_starz = zeros(waypt_num-1,2*traj_order);

w_star = cell(waypt_num-1,1);

% compute Omega for each segment
for i = 1:(waypt_num-1)
    data_x = [];
    data_y = [];
    data_z = [];
    for k=1:traj_order
        data = varargin{2+k};
        data_x = [data_x; data(1,i:i+1)'];
        data_y = [data_y; data(2,i:i+1)'];
        data_z = [data_z; data(3,i:i+1)'];
    end
    Omega = omega_funcN([t(i);t(i+1)],traj_order)';
    w_starx(i,:) = (Omega*Omega') \ (Omega*data_x);
    w_stary(i,:) = (Omega*Omega') \ (Omega*data_y);
    w_starz(i,:) = (Omega*Omega') \ (Omega*data_z);

    w_star{i} = [w_starx(i,:);w_stary(i,:);w_starz(i,:)];
end

prod_vec = ones(1,2*traj_order);
for k=1:nargout
    pow_vec = 2*traj_order-k:-1:1;
    rx = polyval(prod_vec.*w_starx(1,end:-1:k),t(1));
    ry = polyval(prod_vec.*w_stary(1,end:-1:k),t(1));
    rz = polyval(prod_vec.*w_starz(1,end:-1:k),t(1));
    
    for i=1:waypt_num-1
        time = time_cell{i};
        seg_size = length(time);
    
        rx = [rx, polyval(prod_vec.*w_starx(i,end:-1:k),time(2:end))];
        ry = [ry, polyval(prod_vec.*w_stary(i,end:-1:k),time(2:end))];
        rz = [rz, polyval(prod_vec.*w_starz(i,end:-1:k),time(2:end))];
    end
    prod_vec = prod_vec(1:length(pow_vec)).*pow_vec;
    r = [rx; ry; rz]; % Coordinate of trajectory position
    varargout{k} = r;
end

end