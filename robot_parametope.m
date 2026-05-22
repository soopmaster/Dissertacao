% Definition of Parametope:
%   ((xo,alph,y))g := {x \in R^n : g(alph,x-xo)<=y}
% Controlled Embeddings:
%   dxo = f(xo,wo),   dalph = U(xo,alph,y), dy = E(xo,alph,y)
% Elipsoidal Reachable Sets:
%   g(alph,x-xo) = norm(alph*(x-xo))
%   U(xo,alph,y) = -alph*Df(xo)
%   E(xo,alph,y) = c*y
% where
%   c = max_i(max(eig((alph*(Mi-Df(xo))*alph^-1)'+ alph*(Mi-Df(xo))*alph^-1))
% and
%   Mi \in M (set of interval Mixed Jacobian Matrices)

close all;
clearvars;

% System is robot arm
syms q1 q2 z1 z2 
symvars=[q1, q2, z1, z2];
%symvars=[z1,z2,q1,q2];

varmin = [-0.1;-0.1;-0.01;-0.01];
varmax = [0.1;0.1;0.01;0.01];
varbounds = [varmin, varmax];
varbounds_init = varbounds;

u1 = 2;
u2 = 1;
m = 1;
M = 1;
L = sqrt(3);
kp1 = 2;
kp2 = 1;
kd1 = 2;
kd2 = 1;

f = [z1;
     z2; 
     1/(m*q2^2+M*L^2/3)*(-2*m*q2*z1*z2-kd1*z1+kp1*(u1-q1));
     q2*z1^2+1/m*(-kd2*z2+kp2*(u2-q2))];

Df = jacobian(f, [z1, z2, q1, q2]);
Df_vec = reshape(Df,1,[]);

ops = sdpsettings('verbose',0);
enable_disp = 1;
dt = 0.01;
t_dur = 0:dt:10;
inner_step = 20;
initial_state = [1.505,1.505,0.005,0.005]';
dimension = 4;

Initial_alph = [2.5063   -0.0276    0.9436    0.0280 ; -0.0276    2.4722    0.0013    0.9537; 0.9436    0.0013    2.2441    0.0151; 0.0280    0.9537 0.0151    2.2845]^(1/2);
initial_offset=0.02;

[t_precise,x_precise] = ode45(@robotArmODE,t_dur,initial_state);
x_precise = x_precise';

%[t,y] = ode45(@(t,y) robot_prmtpU(y,Dfxo_vec(x:,k)),[0 dt],alph(:,k));

%[t,y] = ode45(@(t,y) robot_prmtpE(y,alph,Df_vec,symvars,varmin,varmax,x(:,k),alph_vec(:,k),dim),[0 dt],y(:,k));

assume(varmin'+initial_state'<symvars & symvars<varmax'+initial_state');

if enable_disp
    figure;
    plot(x_precise(1,:),x_precise(2,:));
    xlim = ylim;
    hold on;
    scatter(initial_state(1),initial_state(2),180,'x','b','linewidth',2);
    xlabel('x1','FontSize', 12);
    ylabel('x2','FontSize', 12);
    drawnow;
    if enable_disp
        plot_ellipsoid_offset(x_precise([1,2],1),Initial_alph([1,2],[1,2]),initial_offset,'r');
    end
end

index = 1:inner_step:length(t_precise);
tsim = t_precise(index);
x = x_precise(:,index);
knum = length(tsim);

alph_vec = zeros(dimension^2,knum);
alph_vec(:,1) = reshape(Initial_alph,[],1);

ysim = zeros(1,knum);
ysim(1) = initial_offset;

% Inital Cell matrix with range structs of Df
range_struct = cell(size(Df));
num_nconst=0;
constDfMx = zeros(dimension); % Build matrix with constant elements of Df
for j=1:size(Df,2)
    for i=1:size(Df,1)
        if hasSymType(Df(i,j),"variable")
            range_struct{i,j} = getElements(Df(i,j),symvars,varbounds);
            num_nconst=num_nconst+1;
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

% Auxiliary matrix
hasVarMx = hasSymType(Df,"variable");

telapse_total = tic;
for k=1:knum-1
    telapse_k =tic;
    if enable_disp
        disp('----------------');
        disp(['Computing ', num2str(tsim(k)), ' time interval ']);
    end

    %Dfxo = Dfxo_vec(Df,symvars,x(:,k)');
    Dfxo = reshape(robot_jac(x(:,k)),1,[]);
    %reshape(Dfxo,dimension,[])
    % [t,y] = ode45(@(t,y) robot_prmtpU(y,Dfxo,dimension),[0 dt],alph_vec(:,k));
    % alph_vec(:,k+1) = y(end,:)';
    alph_vec(:,k+1) = runge_kutta4(@(y) prmtpU(y,Dfxo,dimension),dt*inner_step,alph_vec(:,k));
    alph = reshape(alph_vec(:,k+1),dimension,[]);

    %assume(varmin'+x(:,k)'<symvars & symvars<varmax'+x(:,k)');
    varbounds(:,1) = varbounds_init(:,1)+x(:,k);
    varbounds(:,2) = varbounds_init(:,2)+x(:,k);
    % [t,y] = ode45(@(t,y) robot_prmtpE(y,Df_vec,symvars,x(:,k)',alph_vec(:,k),dimension),[0 dt],ysim(:,k));
    % ysim(:,k+1) = y(end,:)';

    % Calculate c just once outside of runge-kutta
    parametope_c = prmtpEllips_c(range_struct,symvars,varbounds,Dfxo,alph_vec(:,k),dimension,code,hasVarMx,constDfMx);
    ysim(k+1) = runge_kutta4(@(y) prmtpE(y,parametope_c), dt*inner_step, ysim(k));
   
    disp(['Computing time for interval: ', num2str(toc(telapse_k)), 's']);
    plot_ellipsoid_offset(x([1,2],k),alph(1:2,1:2),ysim(k),'r');
    drawnow;

end
disp(['Total sim computing time: ', num2str(toc(telapse_total)), 's']);
