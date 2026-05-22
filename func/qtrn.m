function [q] = qtrn(e,Phi)
%   Return quaternion associated to euler axis (e) and euler angle (Phi)

q = [e.*sin(Phi/2); cos(Phi/2)];

end