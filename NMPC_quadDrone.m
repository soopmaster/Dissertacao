function X = NMPC_quadDrone(XU0)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

% The objective is to align D's parametope with A's
% This is, to maximize the inner product between the principal axis
% of A's and D's ellipsoid parametope
% While also minimizing distance to 0
% So we want both of these in the cost function
% Problem is, the parametope has its own dynamics, which depend on drone's
% dynamics

function [c,ceq] = stateconstraints(XU,Ceq,X0)
    c = []; % Compute nonlinear inequalities
    ceq = Ceq(XU,X0); % Compute nonlinear equalities
end

function J = cost(X,U,Y,Qbar,Qbarm,Rbar)
    xind_step = nx^2+nx;
    yind_step = ny^2+ny;
    Cx = [eye(3), zeros(3,nx-3)];
    Cy = [eye(3), zeros(3,ny-3)];

    Dbar = [];
    xbar = [];
    Abar = [];
    ybar = [];
    for k=0:N
        Dnxk = reshape(X(k*xind_step+(nx+1:xind_step)),nx,[]);
        Dk = proj3(Dnxk);
        Dbar = [Dbar; Dk];

        Anyk = reshape(Y(k*yind_step+(ny+1:yind_step)),ny,[]);
        Ak = proj3(Anyk);
        Abar = [Abar, Ak];
    
        xnxk = X(k*xind_step+(1:nx));
        xk = Cx*xnxk;
        xbar = [xbar; xk];

        ynyk = Y(k*yind_step+(1:ny));
        yk = Cy*ynyk;
        ybar = [ybar; yk];
    end


    J = -m2v(Abar)'*m2v((Qbarm*Dbar)') + (ybar-xbar)'*Qbar*(ybar-xbar) + U'*Rbar*U;
end


nonlcon = @(XU) stateconstraints(XU,Ceq,X0);
X = fmincon(@cost,XU0,[],[],[],[],[],[],nonlcon);
end

fk = @(x,u) runge_kutta4(@(xfull) fullDynFunc(@(x) vehiSt.dynFunc(x,u), @(alph,x) prmtpU(alph,-vehiSt.jacFunc(x,u),nx),xfull,nx),Ts,x);

F = @(X,U,X0) [X(1:nx)-X0;               % x0 = init
            fk(X(1:nx),U(1:nu))-X(1*nx+(1:nx));   % f(x0,u0) = x1
            fk(X(1*nx+(1:nx)),U(1*nu+(1:nu)))-X(2*nx+(1:nx));   % f(x1,u1) = x2
            fk(X(2*nx+(1:nx)),U(2*nu+(1:nu)))-X(3*nx+(1:nx));   % f(x2,u2) = x3
            fk(X(3*nx+(1:nx)),U(3*nu+(1:nu)))-X(4*nx+(1:nx));   % f(x3,u3) = x4
            fk(X(4*nx+(1:nx)),U(4*nu+(1:nu)))-X(5*nx+(1:nx));];  % f(x4,u4) = x5

Ceq = @(XU,X0) F(XU(1:(N+1)*nx,1),XU((N+1)*nx+(1:(N*nu)),1),X0);