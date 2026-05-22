function e = inv_vanderpol_error(Box)

D = zeros(2,2);
D(1,:) = 0;

mul_max = max([Box(1,2)*Box(2,2),Box(1,2)*Box(2,1),Box(1,1)*Box(2,2),Box(1,1)*Box(2,1)]);
mul_min = min([Box(1,2)*Box(2,2),Box(1,2)*Box(2,1),Box(1,1)*Box(2,2),Box(1,1)*Box(2,1)]);

D(2,1) = 2*(mul_max-mul_min);

x1_max = max(abs(Box(1,1)),abs(Box(1,2)));
x1_min = min(abs(Box(1,1)),abs(Box(1,2)));
D(2,2) = (x1_max^2-x1_min^2);

e = norm(D+D');