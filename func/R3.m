function Rz = R3(x)
%   Return DCM of rotation of x degrees around axis 3 
    Rz = [cos(x) -sin(x) 0;
          sin(x) cos(x) 0;
          0 0 1];
end