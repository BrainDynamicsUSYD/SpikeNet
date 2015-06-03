clc;clear;close all;



R=  load('~/Desktop/012-201505191310-17121_1432011996519_RYG.mat');
V_mean = R.pop_stats.V_mean{1};
V_std = R.pop_stats.V_std{1};
I_mean = R.pop_stats.I_input_mean{1} + 1.5; % 1.5 is the constant current
I_std = R.pop_stats.I_input_std{1};
num_spikes = R.num_spikes{1};



V_simu = zeros(size(V_mean));

V_simu(1) = V_mean(1);

dt = 0.1;
t_m = 15;
scale = 0.8;
for i = 2:length(V_mean)
    dV = -V_simu(i-1)/t_m + I_mean(i-1)*scale;
    
    V_simu(i) = V_simu(i-1) + dV*dt;
end

plot(V_simu)
%plot(V_mean, V_simu)