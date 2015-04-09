clc;clear;close all;

load('example_constant_current_400sec_metaresults_improved_detection_08_200.mat');
% % deterministic system

% load('example_400sec_metaresults.mat'); % noisy system

theta_i = 2;

switch_seq = switch_seq_cell{theta_i};
switch_seq_loop = switch_seq_loop_num_cell{theta_i};

% min_length = 30;
% loop_8 = [];
% for loop  = unique(switch_seq_loop)
%     seq = switch_seq(switch_seq_loop == loop);
%     if length(seq) >= min_length && length( unique(seq) ) >= 7
%         loop_8 = [loop_8 loop];
%     end
% end

low_du = low_du_cell{theta_i};
low_du_loop = switch_level_loop_num_cell{theta_i};

high_du = high_du_cell{theta_i};
high_du_loop = high_du_loop_num_cell{theta_i};

switch_seq = switch_seq_cell{theta_i};
switch_seq_loop = switch_seq_loop_num_cell{theta_i};

cluster_low_du = [];

for loop = unique(high_du_loop)
    low_du_tmp = low_du(low_du_loop == loop);
    high_du_tmp = high_du(high_du_loop == loop);
    seq_tmp = switch_seq(switch_seq_loop == loop);
    
    % recover high_start
    low_du_tmp = low_du_tmp(:)';
    high_du_tmp = high_du_tmp(:)';
    
    high_start = cumsum([1 (low_du_tmp + high_du_tmp(1:end-1))]);
    
    for c = unique(seq_tmp)
        high_du_c = high_du_tmp(seq_tmp == c);
        high_start_c = high_start(seq_tmp == c);
        if length(high_du_c) >= 2
            low_du_c = high_start_c(2:end) - high_start_c(1:end-1) - high_du_c(1:end-1);
            cluster_low_du = [cluster_low_du low_du_c];
        end
    end
    
end

bin_num = 200;

figure(1);set(gcf, 'color', 'w');hold on;

[y,x] = hist(cluster_low_du, bin_num);
y = y/sum(y); % normalize it for comparison
plot(x,y,'r');

% fit data
f = fit(x(x>10^3.5)', y(x>10^3.5)', 'power1') % non-linear regression, use GUI toolbox for more details!
alpha_low = f.b; % y=a*x^b
m_low = f.a;
plot( x, m_low*x.^alpha_low, 'k--' )


% %%  rigorous power law test!!
% x = x( randperm(length(x), 10000) ); % sub-sampling
% power_data = x;
% [alpha, xmin, L] = plfit(cluster_low_du);
% [alpha, xmin, ntail] = plvar(power_data); % this would take a minute
% % what's the relationship between results given by plfit and plvar??
% plplot(power_data, xmin, alpha)
% %%% check p value: this would take a minute
% plpva(power_data,xmin) 
% disp('p value should be larger than 0.1 for powerlaw to be a good fit');
% alpha, xmin


% cd /import/yossarian1/yifan/Project1/double_well_data/
% low_du_model = CollectVectorYG('low_du','low_du','*mat');
% [y,x] = hist(low_du_model, bin_num);
% y = y/sum(y); % normalize it for comparison
% plot(x,y,'b:');

set(gca, 'xscale', 'log')
set(gca, 'yscale', 'log')



