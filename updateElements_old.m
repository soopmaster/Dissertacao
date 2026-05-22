function [range] = updateElements_old_old(range)
%Update .r field of existing range struct based on bounds of vars
%   Detailed explanation goes here

% range : struct (struct with components of expr)
%   .f : sym                        (expr)
%   .op : "plus"/"times"/"power"/"cos"...  (operation between RHS and LHS)
%   .a : struct range               (LHS of expr)
%   .b : struct range               (RHS of expr)
%   .r : [a1,a2]                    (range of expr)

    oprt = range.op;
    
    switch oprt
        case {"integer","constant","rational","vpareal","complex"}
            return

        case "variable"
            range.r = findRange(range.f);
            return
            
        case "plus"
            range.a = updateElements_old(range.a);
            range.b = updateElements_old(range.b);
            xy = [range.a.r(1)+range.b.r(1), range.a.r(2)+range.b.r(2)];
            range.r = [min(xy), max(xy)];
            return

        case "times"
            range.a = updateElements_old(range.a);
            range.b = updateElements_old(range.b);
            xy = [range.a.r(1)*range.b.r(1), range.a.r(1)*range.b.r(2), range.a.r(2)*range.b.r(1), range.a.r(2)*range.b.r(2)];
            range.r = [min(xy), max(xy)];
            return

        case "power"
            % If power is uneven, it's monotonic
            % If power is even, check if critical point is between bounds
            % If power is negative,...
            if range.b<0
                range.a = updateElements_old(range.a);
                if not((range.a.r(1)<0) && (0<range.a.r(2)))
                    xy = [1/range.a.r(2), 1/range.a.r(1)];
                elseif range.a.r(2) ~= 0 && range.a.r(2) == 0
                    xy = [-inf,1/range.a.r(1)];
                elseif range.a.r(2) == 0 && range.a.r(2) ~= 0
                    xy = [1/range.a.r(2),inf];
                else
                    xy = [-inf,inf];
                end
                range.r = [min(xy),max(xy)];
            elseif is_even(range.b)
                % Critical point is at 0
                range.a = updateElements_old(range.a);
                if (range.a.r(1)<0) && (0<range.a.r(2))
                    xy = [range.a.r(1)^range.b, range.a.r(2)^range.b, 0];
                else
                    xy = [range.a.r(1)^range.b, range.a.r(2)^range.b];
                end
                range.r = [min(xy),max(xy)];
            else
                range.a = updateElements_old(range.a);
                xy = [range.a.r(1)^range.b, range.a.r(2)^range.b];
                range.r = [min(xy),max(xy)];
            end
            return

        case "cos"
            range.a = updateElements_old(range.a);
            a1 = range.a.r(1);
            a2 = range.a.r(2);
            % If range.a.r is greater than 2*pi range is between [-1,1]
            num_pi = floor((a2-a1)/pi);
            xy = [cos(a1),cos(a2)];
            if num_pi>2
                xy = [-1,1];
            else
                if a1<=0
                    n = ceil(a1/pi)+1*~mod(a1,pi);
                else
                    n = floor(a1/pi)+1*~mod(a1,pi); 
                end
                while n*pi < a2
                    xy = [xy,cos(n*pi)];
                    n = n+1;
                end

            end
            range.r = [min(xy),max(xy)];
            return

        case "sin"
            range.a = updateElements_old(range.a);
            a1 = range.a.r(1);
            a2 = range.a.r(2);
            % If range.a.r is greater than 2*pi range is between [-1,1]
            num_pi = floor((a2-a1)/pi);
            xy = [sin(a1),sin(a2)];
            if num_pi>2
                xy = [-1,1];
            else
                if a1<=0
                    n = ceil(2*a1/pi)+1*~mod(a1,pi/2);
                else
                    n = floor(2*a1/pi)+1*~mod(a1,pi/2); 
                end
                while (n+1/2)*pi < a2
                    xy = [xy,cos((n+1/2)*pi)];
                    n = n+1;
                end

            end
            range.r = [min(xy),max(xy)];
            return

        case "tan"
            % Do Something...
        
        case "exp"
            % Do Something...
    end


end

function bool = is_even(num)
bool = ~mod(num,2);
end

function elem_out = addIfNew(elem_aux,elem_out)
ind = 0;
while ind < length(elem_out)
    if isequal(elem_aux,elem_out(ind+1))
        break
    end
    ind = ind+1;
end
if ind >= length(elem_out)
    elem_out(end+1) = elem_aux;
end
end