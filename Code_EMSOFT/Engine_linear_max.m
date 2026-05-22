function A = Engine_linear_max(Box)

A = [0, -1;3,-1];

if Box(1,1) < -1 && Box(1,2) > -1
    A(1,1) = 1.5;
else
    A(1,1) = max(-3*Box(1,1)-3.0*(Box(1,1)^2)/2,-3*Box(1,2)-3.0*(Box(1,2)^2)/2);
end