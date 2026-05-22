function [elem_out] = getNonConstElements(expr,elem)
%Get the list of non-constant elements in expression
%   Detailed explanation goes here

    in_vec = children(expr);
    in_vec = [in_vec{:}].';
    operation_type = symFunType(expr);

    const_ind = isSymType(in_vec,"constant");
    in_vec = in_vec(not(const_ind));
    
    power_ind = isSymType(in_vec,"power");
    power_vec = in_vec(power_ind); % contains subexpressions of input that are powers
    rest_vec = in_vec(not(power_ind)); % contains subexpressions that are not powers

    elem_out = elem;

    if not(isempty(power_vec))
        power_cell=children(power_vec);
        power_cell = [power_cell{:}].';
    
        power_cell_ind = isSymType(power_cell,"constant");
    
        power = power_cell(power_cell_ind);
    
        if power < 0
            elem_aux = power_cell(not(power_cell_ind))^power;
            elem_out = addIfNew(elem,elem_aux,elem_out);
        end
        %elem_out = elem;
    end
    if operation_type == "power"
        elem_aux = expr;
        elem_out = addIfNew(elem,elem_aux,elem_out);
    end
    if operation_type == "times" && not(isempty(rest_vec))
        elem_aux = fold(@times, rest_vec); % add product of non-powers to list
        elem_out = addIfNew(elem,elem_aux,elem_out);

        % check if expression in power_vec exists in list
        % if not, add to list
        
    elseif operation_type == "plus"
        %elem_out = elem;
        for ind=1:length(in_vec)
            elem_out = getNonConstElements(in_vec(ind),elem_out);
        end
    end
end

function elem_out = addIfNew(elem,elem_aux,elem_out)
aux_is_pwr = isSymType(elem_aux,"power");
if aux_is_pwr
    pwr_expr = children(elem_aux,1);
    pwr_num = children(elem_aux,2);
end
ind = 1;
while ind < length(elem_out)+1
    if isequal(elem_aux,elem_out(ind))
        break
    end
    if aux_is_pwr
        if isSymType(elem_out(ind),"power") && has(pwr_expr,children(elem_out(ind),1))
            if abs(pwr_num) < abs(children(elem_out(ind),2))
                elem_out(ind) = elem_aux;
                return
            else
                return
            end
        end
    end
    ind = ind+1;
end
if ind >= length(elem_out)
    elem_out(end+1) = elem_aux;
end
end