
clc;clear all;close all;

R=  load('~/Desktop/012-201505191310-17121_1432011996519_RYG.mat');
% %%%%%%%%%%%%% compare V_mean time serious

V_mean = R.pop_stats.V_mean{1}  + 70;
V_std = R.pop_stats.V_std{1};
I_mean = R.pop_stats.I_input_mean{1} + 1.5; % 1.5 is the constant current
I_std = R.pop_stats.I_input_std{1};


% Nlags = 1*10^3;
% [c, lags] = autocorr(I_mean(1:10:end), Nlags);
% figure(4);
% plot(lags/1000, c);

num_spikes = R.num_spikes{1};
sigma = 100;
size = 5*sigma;
x = linspace(-size / 2, size / 2, size);
gaussFilter = exp(-x .^ 2 / (2 * sigma ^ 2));
gaussFilter = gaussFilter / sum (gaussFilter); %
num_spikes_sm = transpose(filter(gaussFilter,1, num_spikes'));

figure(1);
seg_1 = 1:40*10^4;
axm(1) = subplot(5,1,1);
plot( (1:40*10^3)*10^-3, num_spikes_sm(1:40*10^3) ); ylabel('spike count');
axm(2) = subplot(5,1,2);
plot( seg_1*10^-4, I_mean(seg_1) ); ylabel('I_{mean}');
axm(3) =subplot(5,1,3);
plot( seg_1*10^-4, I_std(seg_1) );ylabel('I_{std}');
axm(4) = subplot(5,1,4);
plot( seg_1*10^-4, V_mean(seg_1) );ylabel('V_{mean}');
axm(5) = subplot(5,1,5);
plot( seg_1*10^-4, V_std(seg_1) );ylabel('V_{std}');
linkaxes(axm, 'x');
%%%
xData = V_mean;
yData = V_std;
zData = num_spikes;
cData = num_spikes_sm;
clear I_mean V_meam I_std;

figure(2);
set(gcf,'position', [680   349   663   726]);
% subplot(3,1,1);
% plot( (1:40*10^3)*10^-3, cData(1:40*10^3) );
% lh = line([1 1]*10^-4, [0  1], 'color', 'g');

clear R;

axmat(1) =subplot(2,1,1);
grid on;
i_begin = 1*10^4;
i_step = 1;
i_end = 3*10^4;
seg = i_begin:i_step:i_end;
hold on;
plot3(xData(seg), yData(seg), zData(seg),'-b');
xlabel('V_{mean}')
ylabel('V_{std}')
zlabel('Spikes')
xlim([0 20]);
ylim([1 3]);

axmat(2) =subplot(2,1,2);
grid on;
hold on;
i_begin = 8*10^4;
i_step = 1;
i_end = 10*10^4;
seg = i_begin:i_step:i_end;
plot3(xData(seg), yData(seg), zData( seg ),'-r');
xlabel('V_{mean}')
ylabel('V_{std}')
zlabel('Spikes')
xlim([0 20]);
ylim([1 3]);

% figure(3);
% window_length =  0.5*10^4; % 0.5 sec
% dt = 0.1;
% freq_range = 1:0.2:100; % Hz
% fs = 1/dt*1000; % sampling frequency in Hz
% [~,ff,tt,pp] = spectrogram(zData(seg_1), window_length,round(window_length*0.9), freq_range, fs, 'yaxis');  %
% % Setting 'yaxis' to display frequency on the y-axis and time on the x-axis
% % pp is a matrix representing the Power Spectral Density (PSD) of each segment
% 
% tt_jump = 1;
% tt = tt(1:tt_jump:end);
% pp_db = transpose(10*log10(abs(pp(:,1:tt_jump:end))));
% 
% % % % waterfall plot too slow
% % subplot(2,1,1)
% % waterfall(ff, tt, pp_db);
% % xlabel('Hz');
% % ylabel('Time (sec)');
% % zlabel('PSD (dB)');
% % set(gca,'xscale','log');
% 
% imagesc(ff, tt, pp_db);
% xlabel('Hz');
% ylabel('Time (sec)');
