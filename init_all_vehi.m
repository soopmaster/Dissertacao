function vehi_vec = init_all_vehi(Ts,pred_hor)
%INIT_ALL_VEHI outputs cell vector of vehiStruct
%   Detailed explanation goes here

vehi_vec = cell(1,2);

% Define Inital Positions
p10 = [5;-5;5];
% p20 = [-5;5;0];
p20 = [-1;1;0];
x0 = {p10, p20};

vehi_vec{1}.vehiSt = init_QuadDroneNMPC(x0{1},Ts);
vehi_vec{2}.vehiSt = init_QuadDroneNMPC(x0{2},Ts);

% vehi_vec{1}.ctrlSt = init_NMPCctrl(vehi_vec{1}.vehiSt.nx,5,Ts,vehi_vec{1}.vehiSt);

vehi_vec{1}.ctrlSt = init_MellingerCtrl(2,vehi_vec{1}.vehiSt);
vehi_vec{2}.ctrlSt = init_NMPCctrl_v2(vehi_vec{1}.vehiSt.nx,pred_hor,Ts,vehi_vec{2}.vehiSt);
%vehi_vec{2}.ctrlSt = init_NMPCctrl(vehi_vec{1}.vehiSt.nx,pred_hor,Ts,vehi_vec{2}.vehiSt);
%vehi_vec{2}.ctrlSt = init_MellingerCtrl_v2(pred_hor,Ts,vehi_vec{2}.vehiSt);
%vehi_vec{2}.ctrlSt = init_MellingerCtrl(4,vehi_vec{2}.vehiSt);

% === Define ops === %
ops = sdpsettings('verbose',0);
ops.plot.shade = 0.15;
ops.plot.linewidth=0.05;
ops.plot.wirestyle=':';
ops.plot.edgecolor='#EB4B2D';
%ops.plot.edgecolor=[235 75 45];
vehi_vec{1}.ops = ops;

ops = sdpsettings('verbose',0);
ops.plot.shade = 0.15;
ops.plot.linewidth=0.05;
ops.plot.wirestyle=':';
ops.plot.edgecolor='#11F4BF';
%ops.plot.edgecolor=[17 244 191];
vehi_vec{2}.ops = ops;

% === Init Plot Obs === %
L = 0.5;
p_NB = [L;0;0];
p_SB = [-L;0;0];
p_EB = [0;L;0];
p_OB = [0;-L;0];
p_vec = [p_NB; p_EB; p_SB; p_OB];
vehi_vec{1}.p_vec = p_vec;

L = 0.5;
p_NB = [L;0;0];
p_SB = [-L;0;0];
p_EB = [0;L;0];
p_OB = [0;-L;0];
p_vec = [p_NB; p_EB; p_SB; p_OB];
vehi_vec{2}.p_vec = p_vec;

end