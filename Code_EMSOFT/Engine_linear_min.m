function A = Engine_linear_min(Box)

A = [0, -1;3,-1];
A(1,1) = min(-3*Box(1,1)-3.0*(Box(1,1)^2)/2,-3*Box(1,2)-3.0*(Box(1,2)^2)/2);