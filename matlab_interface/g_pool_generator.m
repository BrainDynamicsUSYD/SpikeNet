function [ g, g_info ] = g_pool_generator( N, mu_p, s_p )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
% % pulished pdf curve
% mu_p = -0.702; % published data
% s_p = 0.9355; % published data
% subplot(2,1,1); hold on;
% w = 10.^(-2:0.01:1);
% f_d = 1/(s_p *sqrt(2*pi));
% f = f_d*exp( -(log(w) - mu_p).^2/(2*s_p ^2))./w;
% plot( w, f, 'b')
% set(gca,'xscale','log');
% % % reverse-engineer the  distribution parameters to match the published pdf
% % % curve
% EPSP_mu_norm = mu_p-s_p^2; % an analytical relationship!!!!! How to prove it???
% EPSP_sigma_norm = s_p;
% EPSP_pool_generator = @(N)lognrnd(EPSP_mu_norm,EPSP_sigma_norm,[1,N]); 
% EPSP_pool = EPSP_pool_generator(1000000);
% counts = histc(EPSP_pool, logspace(-2,1,100)); % log-scale is very different from linear-scale!!!
% counts = counts/max(counts)*max(f);
% plot(logspace(-2,1,100), counts,'rx');
% set(gca,'xscale','log')
% xlim([10^-2,10^1]);
% subplot(2,1,2);
% hist(log10(EPSP_pool),100);
% xlim([-2,1]);
% 
% [ ~, fit_EPSP_2_g ] = g_EPSP_conversion( 0.01:0.01:1 );
% g_pool_generator = @(N)fit_EPSP_2_g( EPSP_pool_generator(N) );



EPSP_max = 20; % mV

% pulished pdf curve for log-normally distributed EPSP
% mu_p = -0.702; % published data
% s_p = 0.9355; % published data

% mu_p = -0.2; % 
% s_p = 0.5; %

disp('Note that log-bin correction is being used!')
EPSP_mu_norm = mu_p-s_p^2; % log-bin correction to parameter mu!!
EPSP_sigma_norm = s_p;

[EPSP_mu, EPSP_sigma] = lognstat(EPSP_mu_norm, EPSP_sigma_norm) %#ok<NOPRT>

EPSP_pool_generator = @(N)lognrnd(EPSP_mu_norm,EPSP_sigma_norm,[1,N]); 


[ ~, fit_EPSP_2_g ] = g_EPSP_conversion( );

EPSP = EPSP_pool_generator(N);
while max(EPSP) > EPSP_max
    EPSP(EPSP > EPSP_max) = EPSP_pool_generator( sum(EPSP > EPSP_max) );
end

g = fit_EPSP_2_g( EPSP );
g = g(:)';

g_info.EPSP_max = EPSP_max;
g_info.mu_p = mu_p;
g_info.s_p = s_p;
g_info.fit_EPSP_2_g = fit_EPSP_2_g;

end

