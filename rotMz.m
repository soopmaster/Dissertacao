function rotMz = rotMz(x)
% == Euler Angles == %
% Rotation Matrix around z axis (yaw)
rotMz = [cos(x) -sin(x) 0;
          sin(x) cos(x) 0;
          0 0 1]; 
end