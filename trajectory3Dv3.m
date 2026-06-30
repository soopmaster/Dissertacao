function varargout = trajectory3Dv3(varargin)
%trajectory3D Create trajectory that passes through given waypoints at times t with
%given velocity, acceleration, jerk at waypoints
% Traj_order is number of points*derivatives given
% Max traj_order is 8; if achieved break up into segments
% If only one output is asked, return trajectory polynomial values

% waypt = [pt1x pt2x ... ptNx; pt1y pt2y ... ptNy; pt1z pt2z ... ptNz]
% v = [v1x v2x ... vNx; v1y v2y ... vNy; v1z v2z ... vNz]
% t = [t1 t2 ... tN]

% w_star{i} = [wx0, wx1, wx2, wx3;
%              wy0, wy1, wy2, wy3;
%              wz0, wz1, wz2, wz3;

t = varargin{1};
dt = varargin{2};

der_num = nargin-2; % Defines whether it's got velocity, acceleration, etc.
waypt_num = length(t);
traj_order = der_num*waypt_num;
if traj_order>8
    seg_num = waypt_num-1;
    traj_order = 2*der_num;
else
    seg_num = 1;
end

time_cell = cell(seg_num,1);
for i=1:seg_num
    time_cell{i} = t(i):dt:t(i+1);
end

% matrix to store the wstar for each segment
w_starx= zeros(seg_num,traj_order);
w_stary= zeros(seg_num,traj_order);
w_starz = zeros(seg_num,traj_order);

w_star = cell(seg_num,1);

% compute Omega for each segment
for i = 1:(seg_num)
    data_x = [];
    data_y = [];
    data_z = [];
    seg_ind = waypt_num/seg_num*(i-1)+1:waypt_num/seg_num*i;
    for k=1:der_num
        data = varargin{2+k};
        data_x = [data_x; data(1,seg_ind)'];
        data_y = [data_y; data(2,seg_ind)'];
        data_z = [data_z; data(3,seg_ind)'];
    end
    Omega = omega_funcN([t(seg_ind)],der_num)';
    w_starx(i,:) = (Omega*Omega') \ (Omega*data_x);
    w_stary(i,:) = (Omega*Omega') \ (Omega*data_y);
    w_starz(i,:) = (Omega*Omega') \ (Omega*data_z);

    w_star{i} = [w_starx(i,:);w_stary(i,:);w_starz(i,:)];
end

if nargout<=1
    varargout{1} = w_star;
    return
end

prod_vec = ones(1,traj_order);
for k=1:nargout
    pow_vec = traj_order-k:-1:1;
    rx = polyval(prod_vec.*w_starx(1,end:-1:k),t(1));
    ry = polyval(prod_vec.*w_stary(1,end:-1:k),t(1));
    rz = polyval(prod_vec.*w_starz(1,end:-1:k),t(1));
    
    for i=1:seg_num
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