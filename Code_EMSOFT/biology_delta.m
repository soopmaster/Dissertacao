function J = biology_delta(x_min,x_max)

J = abs(biology_jac(x_min) - biology_jac(x_max));