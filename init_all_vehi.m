function vehi_vec = init_all_vehi(Ts)
%INIT_ALL_VEHI outputs cell vector of vehiStruct
%   Detailed explanation goes here

vehi_vec = cell(1,2);

vehi_vec{1}.vehiSt = init_QuadDroneMellinger(Ts);
vehi_vec{2}.vehiSt = init_QuadDroneMellinger(Ts);

end