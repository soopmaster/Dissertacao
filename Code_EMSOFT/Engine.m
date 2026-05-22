function dy = Engine(t,y)
dy = zeros(2,1);
dy(1) = -y(2)-3.0*(y(1))^2/2-(y(1))^3/2;
dy(2) = 3*y(1)-y(2);