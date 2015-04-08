% nullclines
clc;clear;close all;




%%%%%%%%%%% table should be loaded as this
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
%%%%%%%%%%% table should be loaded as above



% %%%%%%%%%%% table should be loaded as this
% table_aug = load('v_CV_tables_aug.mat');
% % manually fix some (possibly) numerical errors 
% x_up = table_aug.CV_mat;
% x_up(isnan(x_up)) = 1;
% x_up(isinf(x_up)) = 1;
% table_aug.CV_mat(1:30,:) = x_up(1:30,:);
% 
% x_down = table_aug.CV_mat;
% x_down(isnan(x_down)) = 0;
% table_aug.CV_mat(31:60,:) = x_down(31:60,:);
% table_aug.v_mat(isnan(table_aug.v_mat)) = 0.5;
% % something is still wrong here!!!
% 
% % transpose for interp2
% table_aug.miu_V_mat = table_aug.miu_V_mat'; 
% table_aug.sigma_V_mat = table_aug.sigma_V_mat'; 
% table_aug.CV_mat = table_aug.CV_mat'; 
% table_aug.v_mat = table_aug.v_mat';
% %%%%%%%%%%% table should be loaded as above



%%%%%%%%%%% table should be loaded as this
table_fig4 = load('v_CV_tables_fig4.mat');
% transpose for interp2
table_fig4.miu_V_mat = table_fig4.miu_V_mat'; 
table_fig4.sigma_V_mat = table_fig4.sigma_V_mat'; 
table_fig4.CV_mat = table_fig4.CV_mat'; 
table_fig4.v_mat = table_fig4.v_mat';
%%%%%%%%%%% table should be loaded as above


case_num = 2 ;

t_m = 10;

if case_num == 1
    mu_ext = 18; %mV
    sigma_ext = 0.65; %mV
    c_mu = 7.2; %7.2mV
    c_sigma = 1; %mV
    figure(1);
    set(gcf,'color','w');
    hold on;
    fh1 = @(mu_V, sigma_V) get_eq3_4_table_lookup(mu_V, sigma_V, table_fig4,  'v') - (mu_V-mu_ext)./(t_m*c_mu);
    fh2 = @(mu_V, sigma_V) (sigma_V.^2-sigma_ext^2)/(t_m/2*c_sigma^2) - get_eq3_4_table_lookup(mu_V, sigma_V, table_fig4,  'v*CV^2');
    % ezplot may take several minutes 
    h1 = ezplot(fh1,[17 ,22, 0.5, 0.8]);
    h2 =  ezplot(fh2,[17 ,22, 0.5, 0.8]);

    
elseif  case_num == 2
    mu_ext = 5; %mV
    sigma_ext = 5; %mV
    c_mu = 5; %mV
    c_sigma = 20.2; %mV
    
    figure(1);
    set(gcf,'color','w');
    hold on;
    fh1 = @(mu_V, sigma_V) get_eq3_4_table_lookup(mu_V, sigma_V, table,  'v') - (mu_V-mu_ext)./(t_m*c_mu);
    fh2 = @(mu_V, sigma_V) (sigma_V.^2-sigma_ext^2)/(t_m/2*c_sigma^2) - get_eq3_4_table_lookup(mu_V, sigma_V, table,  'v*CV^2');
    % ezplot may take several minutes 
    h1 = ezplot(fh1,[2 , 25, 2, 20]);
    h2 = ezplot(fh2,[2 , 25, 2, 20]);
    
    
end

set(h1, 'color', 'r','linestyle','-');
set(h2, 'color', 'k','linestyle','--');