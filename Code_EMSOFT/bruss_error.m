function D = bruss_error(Box)

mul_max = max([Box(1,2)*Box(2,2),Box(1,2)*Box(2,1),Box(1,1)*Box(2,2),Box(1,1)*Box(2,1)]);
mul_min = min([Box(1,2)*Box(2,2),Box(1,2)*Box(2,1),Box(1,1)*Box(2,2),Box(1,1)*Box(2,1)]);

D(1,1) = 2* (mul_max-mul_min);
D(1,2) = max(abs(Box(1,1)),abs(Box(1,2)))^2 - min(abs(Box(1,1)),abs(Box(1,2)))^2;
D(2,1) = 2* (mul_max-mul_min);
D(2,2) = max(abs(Box(1,1)),abs(Box(1,2)))^2 - min(abs(Box(1,1)),abs(Box(1,2)))^2;
