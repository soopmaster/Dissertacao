function bool = isExprMonotone(expr,vars)
% Check if expr is monotonic in all axis (vars) along the entire domain
% Prove for all vars that the derivative of the function is always the same sign
% Assume expr is continuous

jac = jacobian(expr,vars);
for ind=1:length(jac)
    d_aux = jac(ind);
    % % If d_aux has vars that aren't vars(ind), consider them as constants
    % vars_daux = symvar(d_aux);
    % if not(isempty(vars_daux)) && (length(vars_daux)>1 || not(isequal(vars(ind),vars_daux)))
    %     bool = false;
    %     return
    % end

    zeroStruct = solve(d_aux,vars(ind),"ReturnConditions",true); % Zeros of Jacobian according to var(ind)
    zNames = fieldnames(zeroStruct);
    zero = getfield(zeroStruct,zNames{1});
    if hasSymType(zero,"variable")
        bool = false;
        return
    else
        zero = double(zero);
    end

    if not(isempty(zero))
        dd_aux = jacobian(d_aux,vars(ind));
        critStruct = solve(dd_aux,vars(ind),"ReturnConditions",true); % Critical points of jacobian
        cNames = fieldnames(critStruct);
        crit = sort(double(getfield(critStruct,cNames{1})));
    
        for jnd=1:length(zero)
            if zero(jnd)>0
                gt_crit_ind = zero(jnd)>crit; % Index of values in crit where zero is greater
            else
                gt_crit_ind = zero(jnd)<crit; % Index of values in crit where zero is greater
            end
            if all(gt_crit_ind)
                % If all values all values are one, zero is not between critical points and function is non-monotonic 
                bool = false;
                return
            elseif length(gt_crit_ind)~=1 && any(gt_crit_ind)
                % Else, check if critpoints on either side of zero have
                %   different signs
                k = find(gt_crit_ind==0);
                if zero(jnd)>0
                    c1 = crit(k-1);
                    c2 = crit(k);
                else
                    c1 = crit(k);
                    c2 = crit(k+1);
                end
                % Substitute other vars with constants
                other_vars = symvar(d_aux);
                fc1 = subs(d_aux,vars(ind),c1);     
                fc2 = subs(d_aux,vars(ind),c2);                
                fprod = fc1*fc2;
                fprod = subs(fprod, other_vars,ones(1,length(other_vars)));
                if  double(fprod) < 0
                    bool = false;
                    return
                end
            end
        end
    end
end

bool = true;

end