%

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


% %% define the parameters
t_m = 10; %ms



figure(1);
hold on;box on;

c_mu = 5; %mV
c_sigma = 20.2; %mV

for mu_ext = 3:7 %5; %mV

fh1 = @(mu_V, sigma_V) get_eq3_4_table_lookup(mu_V, sigma_V, table,  'v') - (mu_V-mu_ext)./(t_m*c_mu);
% ezplot may take several minutes
h1 = ezplot(fh1,[2 , 25, 2, 20]);
set(h1, 'color', 'r','linestyle','-');
title([]);
end


for sigma_ext = 3:7 %5; %mV
 
fh2 = @(mu_V, sigma_V) (sigma_V.^2-sigma_ext^2)/(t_m/2*c_sigma^2) - get_eq3_4_table_lookup(mu_V, sigma_V, table,  'v*CV^2');
h2 = ezplot(fh2,[2 , 25, 2, 20]);
set(h2, 'color', 'k','linestyle','--');

title([]);
end

set(gcf,'color','w')


