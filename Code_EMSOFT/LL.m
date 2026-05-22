function dy = LL(t,y)
    y = reshape(y,7,1);
    x1 = y(1);
    x2 = y(2);
    x3 = y(3);
    x4 = y(4);
    x5 = y(5);
    x6 = y(6);
    x7 = y(7);

    
    x1_dot = 1.4 * x3 - 0.9 * x1;
    x2_dot = 2.5 * x5 - 1.5 * x2;
    x3_dot = 0.6 * x7 - 0.8 * x3 * x2;
    x4_dot = 2.0 - 1.3 * x3 * x4;
    x5_dot = 0.7 * x1 - x4 * x5;
    x6_dot = 0.3 * x1 - 3.1 * x6;
    x7_dot = 1.8 * x6 - 1.5 * x7 * x2;
    
    dy = [x1_dot,x2_dot,x3_dot,x4_dot,x5_dot,x6_dot,x7_dot]';