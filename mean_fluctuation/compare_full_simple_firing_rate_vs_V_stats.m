
clc;clear all;close all;

R=  load('~/Desktop/012-201505191310-17121_1432011996519_RYG.mat');
% %%%%%%%%%%%%% compare V_mean time serious

V_mean = R.pop_stats.V_mean{1};
V_std = R.pop_stats.V_std{1};
I_mean = R.pop_stats.I_input_mean{1} + 1.5; % 1.5 is the constant current
I_std = R.pop_stats.I_input_std{1};
num_spikes = R.num_spikes{1};

sigma = 10/0.1;
size = 5*sigma;
x = linspace(-size / 2, size / 2, size);
gaussFilter = exp(-x .^ 2 / (2 * sigma ^ 2));
gaussFilter = gaussFilter / sum (gaussFilter); %
num_spikes_sm = transpose(filter(gaussFilter,1, num_spikes'));



% load pre-calculated table for (eq 3.4)
table = load('v_CV_tables.mat');
% manually fix some (possibly) numerical errors
table.CV_mat(1:13,1) = 1;
table.CV_mat(end-2:end,1) = 0;
table.v_mat(end-2:end,1) = 0.5;
% something is still wrong here!!!
% transpose for interp2
table.miu_V_mat = table.miu_V_mat';
table.sigma_V_mat = table.sigma_V_mat';
table.CV_mat = table.CV_mat';
table.v_mat = table.v_mat';


figure(1);
subplot(1,2,1);
surf(table.miu_V_mat, table.sigma_V_mat, table.v_mat*1000)
xlim([0 20]);
ylim([0 4]);

subplot(1,2,2);
grid on;hold on;
plot3(V_mean+70, V_std, num_spikes/500/(10^-4));


