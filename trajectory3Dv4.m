function W = trajectory3Dv4(varargin)
%trajectory3D Create trajectory that passes through given waypoints at times t with
%given velocity, acceleration, jerk at waypoints
% Traj_order is number of constraints given
% If only one output is asked, return trajectory polynomial values
% INPUT: (t, P, V,...)
% t = [t0, t1, t2, ... tN]
% P = [p0, p1, p2, ..., pN]
% V = [v0, v1, v2, ..., vN]
% pi = [pix; piy; piz] if defined as constraint
% vi = [NaN; NaN; NaN] id undefined

% w_star{i} = [wx0, wx1, wx2, wx3;
%              wy0, wy1, wy2, wy3;
%              wz0, wz1, wz2, wz3;

t = varargin{1};
%dt = varargin{2};
way = cell(1,length(varargin)-1);
traj_order = 0;
not_nan_ind_mat = logical([]);
Px = [];
Py = [];
Pz = [];
for k=2:length(varargin)
    way{k-1} = varargin{k}; % = P, V, ...
    not_nan_ind = ~isnan(way{k-1}(1,:));
    not_nan_ind_mat = [not_nan_ind_mat; not_nan_ind];
    traj_order = traj_order + sum(double(not_nan_ind));
    Px = [Px; varargin{k}(1,not_nan_ind)'];
    Py = [Py; varargin{k}(2,not_nan_ind)'];
    Pz = [Pz; varargin{k}(3,not_nan_ind)'];
end
T = [];
for k=1:length(varargin)-1
    T = [T; dpolyN(t(not_nan_ind_mat(k,:))',traj_order,k-1)'];
end

Wx = (T'*T)\(T'*Px);
Wy = (T'*T)\(T'*Py);
Wz = (T'*T)\(T'*Pz);

W = [Wx';Wy';Wz'];

end