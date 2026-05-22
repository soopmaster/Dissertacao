function rotMy = rotMy(x)
% == Euler Angles == %
% Rotation Matrix around y axis (pitch)
rotMy = [cos(x) 0 sin(x);
          0 1 0;
          -sin(x) 0 cos(x)]; 
end