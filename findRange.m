function [expr_bds_vec] = findRange(expr)
%   Find range of symbolic expression, assuming its bounded

%   First extract upper and lower bounds of vars in expr in vector form
%   Then get all unique combinations of vars assuming they can assume their
%   bounds
%   Then evaluate expr to all combinations and pick the minimum and maximum
%   Then check if local extrema are within range; if so ignore else define
%   as range

% expr_bds_vec = [x1min, x1max;
%                [x2min, x2max;
%                [...  , ...  ]

    expr_vars = symvar(expr);
    expr_bounds = assumptions(expr);
    expr_bounds = reshape(expr_bounds.',[],2).';
    
    expr_bds_vec = zeros(length(expr_vars),2);
    
    for i=1:length(expr_vars)
        for j=1:2
            expr_bds_vec(i,j) = double(findSymType(expr_bounds(j,i),"constant"));
        end
    end

    %comb = combinations(expr_bds_vec)
end
