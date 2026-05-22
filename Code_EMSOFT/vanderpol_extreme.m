function M = vanderpol_extreme(Box)

miu = 1;
%return the Maximum value of J(2,1), the Minimum value of J(2,1), the
%Maximum value of J(2,2), the Minimum value of J(2,2)
M = zeros(4,1);


mul_max = max(-2*[Box(1,2)*Box(2,2),Box(1,2)*Box(2,1),Box(1,1)*Box(2,2),Box(1,1)*Box(2,1)]);
mul_min = min(-2*[Box(1,2)*Box(2,2),Box(1,2)*Box(2,1),Box(1,1)*Box(2,2),Box(1,1)*Box(2,1)]);

M(1) = -1+miu*mul_max;
M(2) = -1+miu*mul_min;

x1_max = max(abs(Box(1,1)),abs(Box(1,2)));
x1_min = min(abs(Box(1,1)),abs(Box(1,2)));

if Box(1,1)<0&&Box(1,2)>0
    M(3)= 1-x1_max^2;
    M(4)= 1;
else
    M(3) = 1-x1_max^2;
    M(4) = 1-x1_min^2;
end