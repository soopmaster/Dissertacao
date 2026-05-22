function rotMx = rotMx(x)
% == Euler Angles == %
% Rotation Matrix around x axis (roll)
rotMx = [1 0 0;
          0 cos(x) -sin(x);
          0 sin(x) cos(x)]; 
end