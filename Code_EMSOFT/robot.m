function dy = robot(t,y)
dy = zeros(4,1);
m = 1;
l = 3;
kp1 = 2;
kp2 = 1;
kd1 = 2;
kd2 = 1;

u1 = kp1;
u2 = kp2;

dy(1) = y(3);
dy(2) = y(4);
dy(3) = (-2*m*y(2)*y(3)*y(4)-kp1*y(1)-kd1*y(3))/(m*y(2)*y(2)+l/3)+(kp1*u1)/(m*y(2)*y(2)+l/3);
dy(4) = y(2)*y(3)*y(3)-kp2*y(2)/m-kd2*y(4)/m+kp2*u2/m;