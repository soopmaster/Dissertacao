function A = Vanderpol_Linear_min_min(Box)
miu = 1;
A = zeros(size(Box));
A(1,1) = 0;
A(1,2) = 1;

mul_max = max(-2*[Box(1,2)*Box(2,2),Box(1,2)*Box(2,1),Box(1,1)*Box(2,2),Box(1,1)*Box(2,1)]);
mul_min = min(-2*[Box(1,2)*Box(2,2),Box(1,2)*Box(2,1),Box(1,1)*Box(2,2),Box(1,1)*Box(2,1)]);

A(2,1) = -1+miu*mul_min;

if Box(1,1)<0 && Box(1,2)>0
    A(2,2) = min(miu,miu*(1-(max(abs(Box(1,1)),abs(Box(1,2))))^2));
else
    A(2,2) = miu*(1-(max(abs(Box(1,1)),abs(Box(1,2))))^2);
end

