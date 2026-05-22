function [qdot] = qtrndot(q,angv)
%   Return derivative of quaternion q with angular velocity angv

qdot = 1/2*[eye(3)*q(4)+crossprodmat(q(1:3)); -q(1:3)']*angv;

end