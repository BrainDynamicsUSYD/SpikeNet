clc;clear;close all;

R=  load('~/Desktop/012-201505191310-17121_1432011996519_RYG.mat');
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

% seed the matlab rand function! The seed is global.
% Be very careful about that!!!!!!!!!!!
loop_num = 1;
date_now = datestr(now,'yyyymmddHHMM-SSFFF');
scan_temp = textscan(date_now,'%s','Delimiter','-');
rand_seed = loop_num*10^5+eval(scan_temp{1}{2}); % scan_temp{1}{2} is the SSFFF part
rng(rand_seed,'twister'); %  Note that this effect is global!
fprintf('Random number generator seed is: %f\n', rand_seed );

% time
dt = 0.1; % ms
step_tot = 1*10^4; % 100 bio-second took 70mins to simulate
% T = (1:step_tot)*dt/1000;% sec

% %% define the parameters
t_m = 10; %ms    <--- scale it carefully to 15ms !?

% simulation
N_pop = 1;

%
mu_ext_range = minmax(I_mean)
sigma_ext_range = minmax(I_std)

% external inputs
%%%%%%%% try two correlated noises for mu_ext and sigma_ext

load('tot_mean_std_distribution.mat')
% tau_1 = 2e+04; % what's it in the full system? 3.8971e+04, 4 sec
figure(1);
subplot(1,2,1);
hist(exp_mean/10^4, 50);
%set(gca,'xscale','log');
xlabel('Current correlation time (sec)')
subplot(1,2,2);
hist(exp_std/10^4, 50);
xlabel('Current fluctuation correlation time (sec)')
%set(gca,'xscale','log');


tau_1 = mean(exp_mean)
r1 = exp_corr_gaussian_noise(step_tot,tau_1);

% tau_2 = 5e+03; % 4.8828e+03, 500 ms
tau_2 = mean(exp_std)
r2 = exp_corr_gaussian_noise(step_tot,tau_2);


% compensate for 15ms memberane time constant instead of 10ms
% and for many other shits
I_mean = I_mean * 8; % 8
I_std = I_std * 4; % 4

mu_ext = r1*std(I_mean) + mean(I_mean); %
sigma_ext = r2*std(I_std) + mean(I_std); %
        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%mu_ext = 5*rand(1, step_tot);

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

% memberane time constant
t_m = 10; %ms

% bookkeeping
miu_V = zeros(N_pop, step_tot);
sigma2_V = zeros(N_pop, step_tot);
firing_rate = zeros(N_pop, step_tot);
CV_hist =  zeros(N_pop, step_tot);

% initial condition
miu_V(:,1) = rand(N_pop,1)*0;
sigma2_V(:,1) = rand(N_pop,1)*1;

% external inputs (identical to all populations)
% mu_ext = 0.5*ones(1, step_tot);
if nargin == 6
    sigma_ext = zeros(1, step_tot);
end

% 4e4 steps take about 1 min to finish
for t = 2:step_tot
    miu_I = mu_ext(t)/t_m;
    sigma2_I = sigma_ext(t)^2/(t_m/2);
    % integrate the ODEs
    d_miu_V = - miu_V(:,t-1)/t_m + miu_I;
    d_sigma2_V = - sigma2_V(:,t-1)/(t_m/2) + sigma2_I;
    
    miu_V(:,t) = miu_V(:,t-1) + d_miu_V*dt;
    sigma2_V(:,t) = sigma2_V(:,t-1) + d_sigma2_V*dt;
end

for t = 2:step_tot
    miu_V_tmp = miu_V(t);
    sigma_V_tmp =  sigma2_V(t)^0.5;
    v_pre = get_eq3_4_table_lookup(miu_V_tmp, sigma_V_tmp, table,  'v');
    if isnan(v_pre)
        t
        %v_pre = 0; % this is not right
        %[out] = get_eq3_4(miu_V_tmp, sigma_V_tmp, 'v');
        t_m = 10; %ms
        t_ref = 2; %ms
        V_th = 20; %mV
        V_res = 10; %mV is V_res the same as V_r ?? Yes.
        upper = (V_th-miu_V_tmp)/(sqrt(2)*sigma_V_tmp);
        lower = (V_res-miu_V_tmp)/(sqrt(2)*sigma_V_tmp);
        % % symbolic functions are not very reliable !!!
        % syms x y;
        % f1 = exp(x^2)*(1+erf(x));
        % v = 1/(t_ref + t_m*sqrt(pi)*double(int(f1, x, lower, upper)));  % double() does numerical integration
        %%%%%%%%%%%%%%%%%%%%
        % another approach: 2D integral over some non-rectangle area, also faster
        fx = @(x) exp(x.^2).*(1+erf(x));
        v_pre = 1/(t_ref + t_m*sqrt(pi)*integral(fx, lower, upper));
    end
    firing_rate(t) = v_pre;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%R.firing_rate becomes NaN because the value is too close to zero (about < 1e-57)

figure(2);
for i = 1:N_pop
    subplot(N_pop, 1, i);
    plot(dt*(1:step_tot), firing_rate(i,:)*1000);
end
ylabel('Hz');
xlabel('ms');




