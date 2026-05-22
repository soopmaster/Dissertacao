function J = LL_Jac_Vertice(y_min,y_max,INDEX)
y_min = reshape(y_min,7,1);
y_max = reshape(y_max,7,1);

if INDEX(1) ==1
    x2 = min(y_min(2),y_max(2));
else
    x2 = max(y_min(2),y_max(2));
end

if INDEX(2) ==1
    x3 = min(y_min(3),y_max(3));
else
    x3 = max(y_min(3),y_max(3));
end

if INDEX(3) ==1
    x4 = min(y_min(4),y_max(4));
else
    x4 = max(y_min(4),y_max(4));
end

if INDEX(4) ==1
    x5 = min(y_min(5),y_max(5));
else
    x5 = max(y_min(5),y_max(5));
end

if INDEX(5) ==1
    x7 = min(y_min(7),y_max(7));
else
    x7 = max(y_min(7),y_max(7));
end

J = [-0.9, 0, 1.4, 0, 0, 0, 0;...
    0, -1.5, 0, 0, 2.5, 0, 0;...
    0, -0.8*x3, -0.8*x2, 0, 0, 0, 0.60;...
    0, 0, -1.3*x4, -1.3*x3, 0, 0, 0;...
    0.70, 0, 0, -1.0*x5, -1.0*x4, 0, 0;...
    0.30, 0, 0, 0, 0, -3.1, 0;...
    0, -1.5*x7, 0, 0, 0, 1.8, -1.5*x2;];