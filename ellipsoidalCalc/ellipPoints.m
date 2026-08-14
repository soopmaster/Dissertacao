function [Fc,p,x_opt,K] = ellipPoints(w,W,N)
%ellipPoints samples N points from surface of ellipsoid defined by w,W
%   INPUT:
% w : ellipsoid center
% W : ellipsoid shape matrix

% Such that for ellipsoid in quadric form: 
%   E(w,W) = {x in R3 : (x-w)'*inv(W)*(x-w) <= 1}
% And in sphere form:
%   E(w,W) = {s in S(1) : sqrtm(W)*s + w}

nx = length(w);

G = sqrtm(W);
if ~isreal(G)
    G = 1e-6.*eye(nx);
end
X.G = G;
X.c = w;
X.A = zeros(0,nx);
X.b = zeros(0,1);
X.idx = [nx];
X.type = [2];
[Fc,p] = compileCCG(X);
% x_opt = plot(Fc,p(1:3),ops.plot.edgecolor,50,ops);
x_opt = plot(Fc,p(1:nx),'k',N);
x_opt = x_opt{end};
K = convhulln(x_opt');
end