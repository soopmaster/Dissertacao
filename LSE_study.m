% Study of LSE performance for different coefs, compare to minimum

close all;
clearvars;

load prmtpV7_data.mat;

vehi = vehi_save{1,end};
DKL = vehi.error_track(7+1:end,:);

maxDKL = max(DKL,[],'all');
minDKL = min(DKL,[],'all');

DKL_min = min(DKL,[],1);
c = -0.1;
for k=1:size(DKL,2)
    DKL_LSE(k) = LSE(c,DKL(:,k)');
end
MSE = 1/size(DKL,2).*sqrt(DKL_min-DKL_LSE).^2;

figure
hold on;
title('MSE(minDKL,lseDKL)','Interpreter','latex')
subtitle(['smoothing=',num2str(c), ', min=',num2str(minDKL), ', max=',num2str(maxDKL)],'Interpreter','latex')
xlabel('k','Interpreter','latex')
grid on;
ylim([-Inf Inf])
plot(MSE)
inf_ind = find(MSE==-Inf);
plot(inf_ind,max(MSE)*0.5.*ones(1,length(inf_ind)),'rx')
legend('MSE','Inf')

figure
hold on;
title('minDKL vs. lseDKL)','Interpreter','latex')
subtitle(['smoothing=',num2str(c), ', min=',num2str(minDKL), ', max=',num2str(maxDKL)],'Interpreter','latex')
xlabel('k','Interpreter','latex')
grid on;
ylim([-Inf 49000])
plot(DKL_LSE)
plot(DKL_min)
inf_ind = find(DKL_LSE==Inf);
plot(inf_ind,49000*0.5.*ones(1,length(inf_ind)),'bx')
legend('lseDKL','minDKL','Inf in lseDKL')