function dy = bruss(t,y)
dy = zeros(2,1);
dy(1) = 1+(y(1))^2*y(2)-2.5*y(1);
dy(2) = 1.5*y(1)-(y(1))^2*y(2);