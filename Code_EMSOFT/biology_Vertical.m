function Vertical = biology_Vertical(x_max,x_min,index)


if index(1) == 1
   x3 = min(x_min(3),x_max(3));
else
   x3 = max(x_min(3),x_max(3));
end

if index(2) == 1
    x4 = min(x_min(4),x_max(4));
else
    x4 = max(x_min(4),x_max(4));
end

if index(3) == 1
    x5 = min(x_min(5),x_max(5));
else
    x5 = max(x_min(5),x_max(5));
end

if index(4) == 1
    x6 = min(x_min(6),x_max(6));
else
    x6 = max(x_min(6),x_max(6));
end


Vertical = [-0.4 0 50*x4 50*x3 0 0 0;...
    0.4 -1 0 0 0 0 0;
    0 1 -50*x4 -50*x3 0 0 0;...
    0 0 -50*x4 -50*x3 50*x6 50*x5 0;...
    0 0 50*x4 50*x3 -50*x6 -50*x5 0;...
    0 0 0 0 -50*x6 -50*x5 0.5;...
    0 0 0 0 50*x6 50*x5 -0.5;];
