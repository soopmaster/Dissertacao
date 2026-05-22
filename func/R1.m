function Rx = R1(x)
%   Return DCM of rotation of x degrees around axis 1 
    Rx = [1 0 0;
          0 cos(x) -sin(x);
          0 sin(x) cos(x)];
end