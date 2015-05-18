
clc;clear all;close all;
R = load('001-201505181500-15181_1431927741273_RYG.mat');


% %%%%%%%%%%%%% compare V_mean time serious

V_mean = R.popV.V_mean{1};
V_std = R.popV.V_std{1};

I_mean = ones(size(V_mean))*1.5; % constant current
I_std2 = zeros(size(V_mean));
for i = 1:9
    I_mean = I_mean + R.synI{i}.I_mean;
    I_std2 = I_std2 + (R.synI{i}.I_std).^2; % this part is probably wrong!!!s
end


% % subsampling
% DT = 2;
% I_mean = I_mean(1:DT:end);
% V_mean = V_mean(1:DT:end);
% V_std = V_std(1:DT:end);
% I_std2 = I_std2(1:DT:end);


xData = I_mean;
yData = V_mean;
zData = I_std2.^0.5;

i_begin = 4.8*10^4;
i_step = 2;
i_end = 5.5*10^4;


seg = i_begin:i_step:i_end;

figure(1);
set(gcf,'position', [680   349   663   726]);
subplot(3,1,1);
raster_plot_all(R,[500, zeros(1,8)])
lh = line([1 1]*10^-4, [1 500], 'color', 'g');

subplot(3,1,2:3);
xlim(minmax(xData(seg)));
ylim(minmax(yData(seg)));
zlim(minmax(zData(seg)));

hold on;
;
i_pre = i_begin;
for i = seg
    
    if mod(i, 100) == 0
    subplot(3,1,1);
    set(lh,'XData', [i i]*10^-4);
    end
    
    subplot(3,1,2:3);
    plot3([xData(i_pre) xData(i)], [yData(i_pre), yData(i)], [zData(i_pre) zData(i)], '-');
    title(sprintf('%.1f ms', i*0.1))
    pause(0.001);
    
    i_pre = i;
end



