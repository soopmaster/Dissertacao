clc; clear;
load A;

dim = 28;
intervals = 30/0.1;

TELASPE = tic;
Reach_info = zeros(intervals,dim^2+1);

[gamma, Q] = QLF(A);
Reach_info(1,1:dim^2) = reshape(A,1,dim^2);
Reach_info(1,dim^2+1) = 0.1;
for i = 2:intervals
    Reach_info(i,1:dim^2) = reshape(A,1,dim^2);
    Reach_info(i,dim^2+1) = 0.1*exp(gamma * i * 0.1);
end
toc(TELASPE)

Reach_initial = reshape(Reach_info(1,1:dim^2),dim,dim)/Reach_info(1,dim^2+1);
dia_initial = sqrt(eig(inv(Reach_initial)));
Reach_final = reshape(Reach_info(i,1:dim^2),dim,dim)/Reach_info(i,dim^2+1);
dia_final = sqrt(eig(inv(Reach_final)));
disp('F/I:');
disp(prod(dia_final)/prod(dia_initial));

sum_volume = 0;
for i = 1:intervals
    Reach_ellip = reshape(Reach_info(i,1:dim^2),dim,dim)/Reach_info(i,dim^2+1);
    dia_now = sqrt(eig(inv(Reach_ellip)));
    current_volume = prod(dia_now);
    sum_volume = sum_volume + current_volume;
end
sum_volume = sum_volume/intervals;
disp('Average/Initial:')
disp(sum_volume/prod(dia_initial));