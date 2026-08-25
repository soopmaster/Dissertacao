function d = ellipDistPoint(varargin)
%ellipDistPoint return distance between ellipsoid and point
%   INPUT:
% q : ellipsoid center
% r : ellipsoid radius
% Q : ellipsoid shape matrix
% p : point

% E(q,Q) = {x in R^n : (x-q)'*Q^-1*(x-q)<=r^2}

if length(varargin) == 3
    q = varargin{1};
    Q = varargin{2};
    p = varargin{3};
    r = sqrt(max(eig(Q)));
    Q = 1/r^2.*Q;
else
    q = varargin{1};
    r = varargin{2};
    Q = varargin{3};
    p = varargin{4};
end

rp = sqrt((p-q)'*inv(Q)*(p-q));

d = norm(p*(1-r/rp));
end