% Study of Kullback-Leibler Divergence as a metric to compare ellipsoids

close all;
clearvars;

knum = 300;
Ts = 0.01;

DKL12 = zeros(1,knum);
DKL21 = zeros(1,knum);

q1 = zeros(3,knum);
Q1 = eye(3);

q2 = zeros(3,knum);
q2(:,1) = [-3;0;0];
Q2 = 0.1.*eye(3);

v2 = 2;
f = @(x) [v2;0;0];
% dq1 = [v;0;0]

[F1,p1,x_opt1,K1] = ellipPoints(q1(:,1),Q1,50);
[F2,p2,x_opt2,K2] = ellipPoints(q2(:,1),Q2,50);

% figure;
% [title_3D,subtitle_3D] = title('$D_{KL}$ Study', ...
%     'k = 0','Interpreter','latex');
% hold on
% grid on
% axis equal
% view(2)
% patch1 = patch('Vertices', x_opt1','Faces',K1,'FaceColor', 'r','FaceAlpha',0.5);
% patch2 = patch('Vertices', x_opt2','Faces',K2,'FaceColor', 'b','FaceAlpha',0.5);
% camlight('headlight','infinite');
for k=1:knum
    %subtitle_3D.String = ['k = ',num2str(k)];
    %q2 = q2 + v2*Ts.*[1;0;0];
    q2(:,k+1) = runge_kutta4(@(x)f(x), Ts, q2(:,k));
    % [F2,p2,x_opt2,K2] = ellipPoints(q2(:,k),Q2,50);
    % patch2.Faces = K2;
    % patch2.Vertices = x_opt2';

    DKL12(k) = kullback_leibler_div(q1(:,k),Q1,q2(:,k),Q2);
    DKL21(k) = kullback_leibler_div(q2(:,k),Q2,q1(:,k),Q1);
    % pause(Ts);
    % drawnow
end

for k=[1,knum/2,knum]
    % subtitle_3D.String = ['k = ',num2str(k)];
    [F2,p2,x_opt2,K2] = ellipPoints(q2(:,k),Q2,50);

    figure;
    [title_3D,subtitle_3D] = title('$D_{KL}$ Study', ...
        ['k = ',num2str(k)],'Interpreter','latex');
    hold on
    grid on
    axis equal
    view(2)
    patch1 = patch('Vertices', x_opt1','Faces',K1,'FaceColor', 'r','FaceAlpha',0.5);
    patch2 = patch('Vertices', x_opt2','Faces',K2,'FaceColor', 'b','FaceAlpha',0.5);
    legend('$\varepsilon_1$', '$\varepsilon_2$','Interpreter','latex')
    camlight('headlight','infinite');
    drawnow
end


figure
grid on;
hold on;
title('DKL Asymmetry', 'Interpreter','latex')
ylim([0 Inf])
xlabel('$k$','Interpreter','latex')
plot(DKL12,'LineWidth',1)
plot(DKL21,'LineWidth',1)
legend('$D_{KL}(\varepsilon_1 || \varepsilon_2)$', '$D_{KL}(\varepsilon_2 || \varepsilon_1)$','Interpreter','latex')