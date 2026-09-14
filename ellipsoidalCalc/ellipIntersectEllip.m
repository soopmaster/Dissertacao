function [qI,QI] = ellipIntersectEllip(q1,Q1,q2,Q2)
%ellipIntersectEllip find external approximation of minimum volume
%ellipsoid of intersection
%   INPUT:
% q1,2 : ellipsoid center
% Q1,2 : ellipsoid shape matrix
%   OUTPUT:
% qI : resulting ellipsoid center
% QI : resulting ellipsoid shape matrix

% E(q,Q) = {x in R^n : (x-q)'*Q-1*(x-q)<=1}
% E_I(qI,QI) = E(q1,Q1) intersect E(q1,Q1) != 0

n = length(q1);

if ~any(any(Q1)) || ~any(any(Q2)) 
    qI = zeros(n,1);
    QI = zeros(n);
    return
end

W1 = inv(Q1); W2 = inv(Q2);

X = @(r) r*W1 + (1-r)*W2;
alpha = @(r) 1 - r*(1-r)*((q2-q1)'*(W2/X(r)*W1)*(q2-q1));
Qp = @(r) alpha(r)*inv(X(r));
qp = @(r) X(r)\(r*W1*q1+(1-r)*W2*q2);

f = @(r) alpha(r)*det(X(r))^2*trace(X(r)\(W1-W2))-...
        n*det(X(r))^2*(2*qp(r)'*(W1*q1-W2*q2)+qp(r)'*(W2-W1)*qp(r)-q1'*W1*q1+q2'*W2*q2);

try 
    [root,fval,exitflag,output] = fzero(f,[0 1]);
catch ME
    % If f has no roots, then one ellipse contains the other.
    % Return smallest ellipse
    if ellipVol(Q1)>ellipVol(Q2)
        qI = q2;
        QI = Q2;
    else
        qI = q1;
        QI = Q1;
    end
    return
end

qI = qp(root);
QI = Qp(root);
%[~,intersect_true] = chol(QI);
intersect_true = all(eig(QI)>0);

if ~intersect_true
    % If QI is impossible, there is no intersection
    qI = zeros(n,1);
    QI = zeros(n);
end

end