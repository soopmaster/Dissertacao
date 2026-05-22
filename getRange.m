function range_vec = getRange(expr,vars)
% Let expr be f:R^n->R continuous and diferentiable and let vars be x \in [a,b]
% Computes range (minimum and maximum) of f(x)
% If f is globally monotonic or just in [a,b], then f \in [f(a),f(b)]

nx = length(vars);
% range_vec = [f(a), a1, a2, ...;
%             [f(b), b1, b2, ...]
range_vec = zeros(2,nx+1);

x_bounds = findRange(expr);
a = x_bounds(:,1)';
b = x_bounds(:,2)';
range_vec(1,2:end) = a;
range_vec(2,2:end) = b;

f_a = subs(expr,vars,a);
f_b = subs(expr,vars,b);
range_vec(1,1) = f_a;
range_vec(2,1) = f_b;

jac = jacobian(expr,vars);
critStruct = solve(jac,vars,"ReturnConditions",true); % Critical points of function
crit_x = zeros(1,nx);
for ind=1:nx
    crit_x(ind) = getfield(critStruct,string(vars(ind)));
end
f_crit = subs(expr,vars,crit_x);
range_vec(size(range_vec,1)+1,1) = f_crit;
range_vec(size(range_vec,1),2:end) = crit_x;

range_vec = sort(range_vec)

% for ind=1:length(jac)
%     d_aux = jac(ind);
%     % % If d_aux has vars that aren't vars(ind), consider them as constants
%     % vars_daux = symvar(d_aux);
%     % if not(isempty(vars_daux)) && (length(vars_daux)>1 || not(isequal(vars(ind),vars_daux)))
%     %     bool = false;
%     %     return
%     % end
% 
%     zeroStruct = solve(d_aux,vars(ind),"ReturnConditions",true); % Zeros of Jacobian according to var(ind)
%     zNames = fieldnames(zeroStruct);
%     zero = getfield(zeroStruct,zNames{1});
%     if hasSymType(zero,"variable")
%         % Not monotonic, do further analysis
%         bool = false;
%         return
%     else
%         zero = double(zero);
%     end
%     if not(isempty(zero))
%         dd_aux = jacobian(d_aux,vars(ind));
%         critStruct = solve(dd_aux,vars(ind),"ReturnConditions",true); % Critical points of jacobian
%         cNames = fieldnames(critStruct);
%         crit = sort(double(getfield(critStruct,cNames{1})));
%         for jnd=1:length(zero)
%             if zero(jnd)>0
%                 gt_crit_ind = zero(jnd)>crit; % Index of values in crit where zero is greater
%             else
%                 gt_crit_ind = zero(jnd)<crit; % Index of values in crit where zero is greater
%             end
%             if all(gt_crit_ind)
%                 % If values are one, zero is not between critical points and function is non-monotonic 
%                 % Not monotonic, do further analysis
%                 bool = false;
%                 return
%             elseif length(gt_crit_ind)~=1 && any(gt_crit_ind)
%                 % Else, check if critpoints on either side of zero have
%                 %   different signs
%                 k = find(gt_crit_ind==0);
%                 if zero(jnd)>0
%                     c1 = crit(k-1);
%                     c2 = crit(k);
%                 else
%                     c1 = crit(k);
%                     c2 = crit(k+1);
%                 end
%                 % Substitute other vars with constants
%                 other_vars = symvar(d_aux);
%                 fc1 = subs(d_aux,vars(ind),c1);     
%                 fc2 = subs(d_aux,vars(ind),c2);                
%                 fprod = fc1*fc2;
%                 fprod = subs(fprod, other_vars,ones(1,length(other_vars)));
%                 if  double(fprod) < 0
%                     % Not monotonic, do further analysis
%                     bool = false;
%                     return
%                 end
%             end
%         end
%     end
% end
% % Monotonic, range is given by [f(a),f(b)]
% 
% bool = true;

end