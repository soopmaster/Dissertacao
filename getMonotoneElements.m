function [elem_out] = getMonotoneElements(expr,elem_in)
%Get the list of monotone sub-expressions in expression
%   Detailed explanation goes here
    expr = expand(expr);
    oprt = symFunType(expr);
    elem_out = elem_in;
    if any(strcmp(oprt,["integer","constant","rational","vpareal","complex"]))
        return
    elseif isExprMonotone(expr,symvar(expr))
        if isequal(oprt, "times")
            expr_aux = children(expr);
            expr_aux = [expr_aux{:}].';
            num_ind = isSymType(expr_aux,"constant");
            expr_aux = expr_aux(not(num_ind)); % Strip constants
            expr = expand(fold(@times, expr_aux));
        end
        elem_out = addIfNew(expr,elem_out);
        return
    end
    % If plus, apply getMonotoneElements on first parcel then other parcel
    % If times, apply getMonotoneElements in parcels that are powers, then other
    %   parcels
    % If power, apply getMonotoneElements

    switch oprt
        case {"integer","constant","rational","vpareal","complex"}
            return
            
        case "plus"
            subexpr_vec = children(expr);
            subexpr_vec = [subexpr_vec{:}].';
            for ind=1:length(subexpr_vec)
                elem_out = getMonotoneElements(subexpr_vec(ind),elem_out);
            end

        case "times"
            %expr = expand(expr);
            subexpr = children(expr);
            subexpr = [subexpr{:}].';
            power_ind = isSymType(subexpr,"power");
            power_vec = subexpr(power_ind); % contains subexpressions of input that are powers
            rest_vec = subexpr(not(power_ind)); % contains subexpressions that are not powers
            a = expand(fold(@times, rest_vec));
            subexpr_vec = [a; power_vec];
            for ind=1:length(subexpr_vec)
                elem_out = getMonotoneElements(subexpr_vec(ind),elem_out);
            end

        case "power"
            subexpr = children(expr,1);
            elem_out = getMonotoneElements(subexpr,elem_out);
    end


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