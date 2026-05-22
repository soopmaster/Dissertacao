function [coord_circ] = nsphere(n,div_num)
%   NSPHERE returns the Cartesian Coordinates of the unitary n-sphere with
% INPUT VARS:
%   n : dimension
%   div_num : number of phi subdivisions
% OUTPUT VARS:
%   coord_circ = [p1, p2, ..., pn]

phi = zeros(n-1,div_num);
for ind=1:n-2
    phi(ind,:) = linspace(0,pi,div_num);
end
phi(end,:) = linspace(0,2*pi,div_num);

coord_circ = zeros(n,div_num^(n-1)-1);
for pt=1:size(coord_circ,2)
    ang_ind = flip(str2num(dec2base(pt-1,div_num,n-1)'))+1;
    for j=1:n
        k_start = j-1;
        if j < n-1
            coord_circ(j,pt) = cos(phi(j,ang_ind(j)));
        elseif j == n-1
            coord_circ(j,pt) = sin(phi(j,ang_ind(j)));
            %coord_circ(j,pt) = 1;
        elseif j == n
            coord_circ(j,pt) = cos(phi(j-1,ang_ind(j-1)));
            k_start = n-2;
        end
        for k=k_start:-1:1
            coord_circ(j,pt) = coord_circ(j,pt)*sin(phi(k,ang_ind(k)));
        end
    end
    scatter3(coord_circ(1,pt),coord_circ(2,pt),coord_circ(3,pt));
    hold on;
end
end