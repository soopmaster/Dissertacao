function [range_struct,constDfMx,code,hasVarMx] = initRangeStruct(Df,dim,symvars,varbounds)
% Inital Cell matrix with range structs of Df
range_struct = cell(dim);
num_nconst=0;
constDfMx = zeros(dim); % Build matrix with constant elements of Df
hasVarMx = zeros(dim);
for j=1:size(Df,2)
    for i=1:size(Df,1)
        if hasSymType(Df(i,j),"variable")
            range_struct{i,j} = getElements(Df(i,j),symvars,varbounds);
            if range_struct{i,j}.r(2)-range_struct{i,j}.r(1) <= 10^-6
                %range_struct{i,j} = range_struct{i,j}.r(2);
                %constDfMx(i,j) = range_struct{i,j}.r(2);
                hasVarMx(i,j) = -1;
            else
                num_nconst=num_nconst+1;
                hasVarMx(i,j) = 1;
            end
        else
            range_struct{i,j} = double(Df(i,j));
            constDfMx(i,j) = double(Df(i,j));
        end
    end
end

% Build "code" variable to help with IMJ
code = zeros(num_nconst,2^num_nconst);
for k=1:2^num_nconst
    code(:,k) = str2num(dec2bin(k-1,num_nconst)')+1;
end

end