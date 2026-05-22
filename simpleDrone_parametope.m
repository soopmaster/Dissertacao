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

% System is simple drone
syms psym vsym ssym wsym [3 1]
syms Fsym
usym = [Fsym; wsym];
symvars=[psym; vsym; ssym];
nx = length(symvars);
symvars_ext = [symvars;usym];

varmin = [-0.*ones(3,1); -0.*ones(3,1); -0.*ones(3,1)];
varmax = -varmin;
varbounds = [varmin, varmax];
varbounds_init = varbounds;
ubounds = [2; pi/4.*ones(3,1)];
varbounds_ext = [varbounds; [-ubounds, ubounds]];

m = 0.5; % mass of drone

g = [0;0;-9.8];
niu_v = 0.01;

fx = [vsym;
      g-niu_v*vsym;
      zeros(3,1)];

gx = [zeros(3,4);
      1/m.*rotM(ssym)*[0;0;1], zeros(3,3);
      zeros(3,1), eye(3)];

dx = fx + gx*usym;

Df = jacobian(dx, symvars);
Df_vec = reshape(Df,1,[]);

ops = sdpsettings('verbose',0);
enable_disp = 1;
dt = 0.01;
t_dur = 0:dt:10;
inner_step = 20;
initial_state = [zeros(nx,1)];
u = [1;0.1;0;0];
dim = nx;

while true
  A = rand(dim, dim);
  if rank(A) == dim; break; end    %will be true nearly all the time
end
alph = A' * A;
alph=eye(dim);
Initial_alph=0.1.*alph^(1/2);
initial_offset=0.01;

[t_precise,x_precise] = ode45(@(t,x) simpleDroneODE(x,u),t_dur,initial_state);
x_precise = x_precise';

%[t,y] = ode45(@(t,y) robot_prmtpU(y,Dfxo_vec(x:,k)),[0 dt],alph(:,k));

%[t,y] = ode45(@(t,y) robot_prmtpE(y,alph,Df_vec,symvars,varmin,varmax,x(:,k),alph_vec(:,k),dim),[0 dt],y(:,k));

assume(varmin'+initial_state'<=symvars & symvars<=varmax'+initial_state');

if enable_disp
    figure;
    plot3(x_precise(1,:),x_precise(2,:),x_precise(3,:));
    grid on;
    hold on;
    scatter3(initial_state(1),initial_state(2),initial_state(3),180,'x','b','linewidth',2);
    xlabel('x','FontSize', 12);
    ylabel('y','FontSize', 12);
    zlabel('z','FontSize', 12);
    drawnow;
    if enable_disp
        plot_ellipsoid_offset(x_precise(1:3,1),Initial_alph(1:3,1:3),initial_offset,'r');
    end
end

index = 1:inner_step:length(t_precise);
tsim = t_precise(index);
x = x_precise(:,index);
knum = length(tsim);

alph_vec = zeros(dim^2,knum);
alph_vec(:,1) = reshape(Initial_alph,[],1);

ysim = zeros(1,knum);
ysim(1) = initial_offset;

[range_struct,constDfMx,code,hasVarMx] = initRangeStruct(Df,dim,symvars_ext,varbounds_ext);

telapse_total = tic;
for k=1:knum-1
    telapse_k =tic;
    if enable_disp
        disp('----------------');
        disp(['Computing ', num2str(tsim(k)), ' time interval ']);
    end

    %Dfxo = Dfxo_vec(Df,symvars,x(:,k)');
    Dfxo = reshape(simpleDrone_jac([x(:,k);u]),1,[]);
    %reshape(Dfxo,dimension,[])
    % [t,y] = ode45(@(t,y) robot_prmtpU(y,Dfxo,dimension),[0 dt],alph_vec(:,k));
    % alph_vec(:,k+1) = y(end,:)';
    alph_vec(:,k+1) = runge_kutta4(@(y) prmtpU(y,Dfxo,dim),dt*inner_step,alph_vec(:,k));
    alph = reshape(alph_vec(:,k+1),dim,[]);

    %assume(varmin'+x(:,k)'<symvars & symvars<varmax'+x(:,k)');
    varbounds(:,1) = varbounds_init(:,1)+x(:,k);
    varbounds(:,2) = varbounds_init(:,2)+x(:,k);
    % [t,y] = ode45(@(t,y) robot_prmtpE(y,Df_vec,symvars,x(:,k)',alph_vec(:,k),dimension),[0 dt],ysim(:,k));
    % ysim(:,k+1) = y(end,:)';

    % Calculate c just once outside of runge-kutta
    parametope_c = prmtpEllips_c(range_struct,symvars_ext,varbounds_ext,-Dfxo,alph_vec(:,k),dim,code,hasVarMx,constDfMx);
    ysim(k+1) = runge_kutta4(@(y) prmtpE(y,parametope_c), dt*inner_step, ysim(k));
   
    disp(['Computing time for interval: ', num2str(toc(telapse_k)), 's']);
    plot_ellipsoid_offset(x(1:3,k),alph(1:3,1:3),ysim(k),'r');
    drawnow;

end
disp(['Total sim computing time: ', num2str(toc(telapse_total)), 's']);
