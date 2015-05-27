
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

xData = V_mean;
yData = V_std;
zData = num_spikes;
cData = num_spikes_sm;



%%%
clear I_mean V_meam I_std;


i_begin = 5*10^4;
i_step = 10;
i_end = 8*10^4;

seg = i_begin:i_step:i_end;


figure(1);
set(gcf,'position', [680   349   663   726]);

subplot(3,1,1);
%raster_plot_all(R,[500, zeros(1,8)])
plot( (1:40*10^4)*10^-4, cData(1:40*10^4) );
lh = line([1 1]*10^-4, [0  1], 'color', 'g');

clear R;

subplot(3,1,2:3);
grid on;


xlim(minmax(xData(seg)));
ylim(minmax(yData(seg)));
zlim(minmax(zData(seg)));

cmin = min(cData(seg))-1;
cmax = max(cData(seg));
crange = cmax - cmin;
creso = 1000;
cmap = cool(creso);

hold on;
i_pre = i_begin;
for i = seg
    
%     if mod(i, 2000) == 0
%         cla
%     end
%     if mod(i, 1000) == 0
%     subplot(3,1,1);
%     set(lh,'XData', [i i]*10^-4);
%     end
    
    subplot(3,1,2:3);
    plot3([xData(i_pre) xData(i)], [yData(i_pre), yData(i)], [zData(i_pre) zData(i)],...
        '-','color',  cmap( ceil((cData(i)-cmin)/crange*creso), :) );
    title(sprintf('%.1f ms', i*0.1))
    % pause(0.001);
    
    i_pre = i;
end


