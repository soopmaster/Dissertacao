function [range] = getElements_old_old(expr)
%Segregate multi-variable expr into a composite of single-variable or monotonic sub-expressions
%   Detailed explanation goes here

% range : struct (struct with components of expr)
%   .f : sym                        (expr)
%   .op : "plus"/"times"/"power"/"cos"...  (operation between RHS and LHS)
%   .a : struct range               (LHS of expr)
%   .b : struct range               (RHS of expr)
%   .r : [a1,a2]                    (range of expr)

    oprt = symFunType(expr);
    
    switch oprt
        case {"integer","constant","rational","vpareal","complex"}
            range.f = expr;
            range.r = [double(expr), double(expr)];
            range.a = struct([]);
            range.b = struct([]);
            range.op = "constant";
            return

        case "variable"
            range.f = expr;
            range.r = findRange(expr);
            range.a = struct([]);
            range.b = struct([]);
            range.op = "variable";
            return
            
        case "plus"
            subexpr_vec = children(expr);
            subexpr_vec = [subexpr_vec{:}].';
            a = subexpr_vec(1);
            b = fold(@plus,subexpr_vec(2:end));
            range.f = expr;
            range.op = "plus";
            range.a = getElements_old(a);
            range.b = getElements_old(b);
            xy = [range.a.r(1)+range.b.r(1), range.a.r(2)+range.b.r(2)];
            range.r = [min(xy), max(xy)];
            return

        case "times"
            subexpr_vec = children(expr);
            subexpr_vec = [subexpr_vec{:}].';
            a = subexpr_vec(1);
            b = fold(@times,subexpr_vec(2:end));
            range.f = expr;
            range.op = "times";
            range.a = getElements_old(a);
            range.b = getElements_old(b);
            xy = [range.a.r(1)*range.b.r(1), range.a.r(1)*range.b.r(2), range.a.r(2)*range.b.r(1), range.a.r(2)*range.b.r(2)];
            range.r = [min(xy), max(xy)];
            return

        case "power"
            % If power is uneven, it's monotonic
            % If power is even, check if critical point is between bounds
            % If power is negative,...
            subexpr_vec = children(expr);
            subexpr_vec = [subexpr_vec{:}].';
            a = subexpr_vec(1);
            pwr_num = double(subexpr_vec(2));
            if pwr_num<0
                range.f = expr;
                range.op = "power";
                range.a = getElements_old(a^abs(pwr_num));
                range.b = -1;
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
            elseif is_even(pwr_num)
                % Critical point is at 0
                range.f = expr;
                range.op = "power";
                range.a = getElements_old(a);
                range.b = pwr_num;
                if (range.a.r(1)<0) && (0<range.a.r(2))
                    xy = [range.a.r(1)^pwr_num, range.a.r(2)^pwr_num, 0];
                else
                    xy = [range.a.r(1)^pwr_num, range.a.r(2)^pwr_num];
                end
                range.r = [min(xy),max(xy)];
            else
                range.f = expr;
                range.op = "power";
                range.a = getElements_old(a);
                range.b = pwr_num;
                xy = [range.a.r(1)^pwr_num, range.a.r(2)^pwr_num];
                range.r = [min(xy),max(xy)];
            end
            return

        case "cos"
            subexpr_vec = children(expr);
            subexpr_vec = [subexpr_vec{:}].';
            a = subexpr_vec(1);
            range.f = expr;
            range.a = getElements_old(a);
            range.b = [];
            range.op = "cos";
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
            subexpr_vec = children(expr);
            subexpr_vec = [subexpr_vec{:}].';
            a = subexpr_vec(1);
            range.f = expr;
            range.a = getElements_old(a);
            range.b = [];
            range.op = "sin";
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