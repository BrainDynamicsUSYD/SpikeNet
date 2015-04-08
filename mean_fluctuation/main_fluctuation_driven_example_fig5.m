

% this script reproduces the results shown in figure 5 of the paper
% Mean-driven and fluctuation-driven persistent activity in recurrent networks
% That is, an example of fluctuation driven bistability.

clc;close;clear all;

%% load pre-calculated table for (eq 3.4)
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

% show plots of table for (eq 3.4)
show_figure = 0;
if show_figure == 1
figure(1);
set(gcf,'color','w');
subplot(1,2,2);
surf(table.miu_V_mat,table.sigma_V_mat,table.CV_mat );
xlabel('\mu_V');
ylabel('\sigma_V');
zlabel('CV');

subplot(1,2,1);
surf(table.miu_V_mat,table.sigma_V_mat,table.v_mat*1000 );
xlabel('\mu_V');
ylabel('\sigma_V');
zlabel('Firing Rate (Hz)');
end


%% define the parameters
t_m = 10; %ms
N_pop = 2; 
% c is number of connection
c = [1 1;    % [EE  IE
    1 1];   %  EI  II]
% j is the connection strength
j = [16.56 16.56;    % [EE  IE
    -11.56 -11.56];   %  EI  II]
% the above c and j values give c_miu = 5 and c_sigma = 20.2
% time
dt = 0.1; % ms
step_tot = 35000;
T = (1:step_tot)*dt/1000;% sec
% bookkeeping
miu_V = zeros(N_pop, step_tot);
sigma2_V = zeros(N_pop, step_tot);
firing_rate = zeros(N_pop, step_tot);
CV_hist =  zeros(N_pop, step_tot);
% initial condition
miu_V(:,1) = ones(N_pop,1)*5.124;
sigma2_V(:,1) = ones(N_pop,1)*5.496^2;
% external inputs
mu_ext  = ones(1,step_tot)*5; %mV
sigma_ext = ones(1,step_tot)*5;  %mV
sigma_ext(5000:10000) = 7; %mV, elevated from 5 to 7


%% simulation
tic;  % about 1 min to finish
for t = 2:step_tot
    
    miu_I = ones(N_pop,1)*mu_ext(t)/t_m;
    sigma2_I = ones(N_pop,1)*sigma_ext(t)^2/(t_m/2);
    % calculate the mean and std of currents
    for i_pre = 1:N_pop
        % solve the mean and CV of first passage time of the O-U process, i.e, (eq 3.4)
        miu_V_pre = miu_V(i_pre, t-1);
        sigma_V_pre =  sqrt(sigma2_V(i_pre, t-1));
        [ out ] = get_eq3_4_table_lookup(miu_V_pre, sigma_V_pre, table,  'both');
        v_pre = out(1);
        CV_pre = out(2);
        firing_rate(i_pre, t) = v_pre;
        CV_hist(i_pre, t) = CV_pre;
        for j_post = 1:N_pop
            c_tmp = c(i_pre, j_post);
            j_tmp = j(i_pre, j_post);
            miu_I(j_post) = miu_I(j_post) + c_tmp*j_tmp*v_pre;
            sigma2_I(j_post) = sigma2_I(j_post) + c_tmp*j_tmp^2*v_pre*CV_pre^2;
        end
    end
    % integrate the ODEs
    d_miu_V = - miu_V(:,t-1)/t_m + miu_I;
    d_sigma2_V = - sigma2_V(:,t-1)/(t_m/2) + sigma2_I;
    
    miu_V(:,t) = miu_V(:,t-1) + d_miu_V*dt;
    sigma2_V(:,t) = sigma2_V(:,t-1) + d_sigma2_V*dt;
    
end
toc;


%% reproduce the figure 5
figure(1);
set(gcf,'color','w');
%
subplot(5,1,1);hold on;box on;
for i = 1:N_pop
    plot(T, CV_hist(i,:),'k');
end
plot([0.5 0.5],[0 1.5],'--');
plot([1 1],[0 1.5],'--');
ylabel('CV');
%
subplot(5,1,2:3);hold on;box on;
for i = 1:N_pop
    plot(T,firing_rate(i,:)*1000,'k');
end
plot([0.5 0.5],[0 90],'--');
plot([1 1],[0 90],'--');
ylabel('Firing Rate (Hz)');
xlabel('Time (s)');
%
subplot(5,1,4:5);hold on;box on;
mu_ext = 5; %mV
sigma_ext = 5; %mV
c_mu = 5; %mV
c_sigma = 20.2; %mV
fh1 = @(mu_V, sigma_V) get_eq3_4_table_lookup(mu_V, sigma_V, table,  'v') - (mu_V-mu_ext)./(t_m*c_mu);
fh2 = @(mu_V, sigma_V) (sigma_V.^2-sigma_ext^2)/(t_m/2*c_sigma^2) - get_eq3_4_table_lookup(mu_V, sigma_V, table,  'v*CV^2');
% ezplot may take several minutes
h1 = ezplot(fh1,[2 , 25, 2, 20]);
h2 = ezplot(fh2,[2 , 25, 2, 20]);
set(h1, 'color', 'r','linestyle','-');
set(h2, 'color', 'k','linestyle','--');
title([]);


