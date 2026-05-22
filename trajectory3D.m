function [r] = trajectory3D(waypt,v,t,res)
%trajectory3D Create trajectory that passes through given waypoints at times t with
%given velocity

% waypt = [pt1x pt2x ... ptNx; pt1y pt2y ... ptNy; pt1z pt2z ... ptNz]
% v = [v1x v2x ... vNx; v1y v2y ... vNy; v1z v2z ... vNz]
% t = [t1 t2 ... tN]
% res = num of points per segment

% w_star{i} = [wx0, wx1, wx2, wx3;
%              wy0, wy1, wy2, wy3;
%              wz0, wz1, wz2, wz3;

waypt_num = size(t,2);

time = linspace(t(1),t(end),res*(waypt_num-1));

% matrix to store the wstar for each segment
w_starx= zeros(waypt_num-1,4);
w_stary= zeros(waypt_num-1,4);
w_starz = zeros(waypt_num-1,4);

w_star = cell(waypt_num-1,1);

% compute Omega for each segment
for i = 1:(waypt_num-1)
    wayptx = waypt(1,i:(i+1))';
    waypty = waypt(2,i:(i+1))';
    wayptz = waypt(3,i:(i+1))';
    vx_vec = v(1,i:(i+1))';
    vy_vec = v(2,i:(i+1))';
    vz_vec = v(3,i:(i+1))';
    Omega = omega_func([t(i);t(i+1)])';
    w_starx(i,:) = (Omega*Omega') \ (Omega*[wayptx;vx_vec]);
    w_stary(i,:) = (Omega*Omega') \ (Omega*[waypty;vy_vec]);
    w_starz(i,:) = (Omega*Omega') \ (Omega*[wayptz;vz_vec]);

    w_star{i} = [w_starx(i,:);w_stary(i,:);w_starz(i,:)];
end


r = zeros(3,res*(waypt_num-1)); % Coordinate of trajectory position

for i=1:waypt_num-1
    r(1,(i-1)*res+1:i*res) = polyval(w_starx(i,end:-1:1),time((i-1)*res+1:res*i));
    r(2,(i-1)*res+1:i*res) = polyval(w_stary(i,end:-1:1),time((i-1)*res+1:res*i));
    r(3,(i-1)*res+1:i*res) = polyval(w_starz(i,end:-1:1),time((i-1)*res+1:res*i));
end

end