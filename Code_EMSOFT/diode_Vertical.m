function [V1,V2] = diode_Vertical(x_max, x_min,critial)

    C = 1;
    L = 1;
    R = 0.2;

    V_d_max = max(x_max(1),x_min(1));
    V_d_min = min(x_max(1),x_min(1));

    for i = 1:length(critial)
        if V_d_min<=critial(i) && V_d_max >= critial(i)        
            dI_dV_d_max = max([poly11(V_d_min) poly11(V_d_max) poly11(critial(i))]);
            dI_dV_d_min = min([poly11(V_d_min) poly11(V_d_max) poly11(critial(i))]);
        else
            dI_dV_d_max = max([poly11(V_d_min) poly11(V_d_max)]);
            dI_dV_d_min = min([poly11(V_d_min) poly11(V_d_max)]);
        end
    end

    % dx(1) = (-I_d + x(2))/C;
    % dx(2) = (-x(1) - R*x(2) + V_in)/L;

    V1 = [-dI_dV_d_max/C, 1/C; -1/L, -R/L];
    V2 = [-dI_dV_d_min/C, 1/C; -1/L, -R/L];
end

function y = poly11(x)
    y = 5*803.712*x^4-4*1086.288*x^3+3*551.088*x^2-2*124.548*x+10.656;
end
