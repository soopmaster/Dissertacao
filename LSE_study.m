% Study of LSE performance for different coefs, compare to minimum

close all;
clearvars;

load prmtpV7_data.mat;

vehi = vehi_save{1,end};
DKL = vehi.error_track(7+1:end,:);

maxDKL = max(DKL,[],'all');
minDKL = min(DKL,[],'all');

DKL_min = min(DKL,[],1);
c = -0.01;
for k=1:size(DKL,2)
    DKL_LSE(k) = LSE(c,DKL(:,k)');
end
MSE = 1/size(DKL,2).*sqrt(DKL_min-DKL_LSE).^2;

fig1=figure('WindowState','maximized')
hold on;
title('min vs LSE (MSE)','Interpreter','latex')
subtitle(['$\sigma$ =',num2str(c), ', min val = ',num2str(minDKL), ', max val = ',num2str(maxDKL)],'Interpreter','latex')
xlabel('k','Interpreter','latex')
grid on;
ylim([-Inf Inf])
plot(MSE,'LineWidth',1.5)
inf_ind = find(MSE==-Inf);
plot(inf_ind,max(MSE)*0.5.*ones(1,length(inf_ind)),'rx')
legend('MSE','Inf in LSE')
fontsize(20,"points");
exportgraphics(fig1,'LSE_study1.pdf')

fig2=figure('WindowState','maximized')
hold on;
title('min vs. LSE','Interpreter','latex')
subtitle(['$\sigma$ =',num2str(c), ', min val = ',num2str(minDKL), ', max val = ',num2str(maxDKL)],'Interpreter','latex')
xlabel('k','Interpreter','latex')
grid on;
ylim([-Inf 49000])
plot(DKL_LSE,'LineWidth',1.5)
plot(DKL_min,'LineWidth',1.5)
inf_ind = find(DKL_LSE==Inf);
plot(inf_ind,49000*0.5.*ones(1,length(inf_ind)),'bx')
legend('LSE','min','Inf in LSE')
fontsize(20,"points");
exportgraphics(fig2,'LSE_study2.pdf')