% try to generate two cross-correlated noises, which themselves are
% auto-correlated
close all;
N = 10^6;

tau_1 = 1000;
r1 = exp_corr_gaussian_noise(N,tau_1);

tau_2 = 500;
r2 = exp_corr_gaussian_noise(N,tau_2);

alpha = 0.4; % zero time-lag correlation of the two signals
r3 = alpha*r1 + sqrt(1-alpha^2)*r2; 
% [ref] Generating partially correlated noise—A comparison of methods
% by Hartmann and Cho

figure(1);

[acf, lags] = autocorr(r1, 3000);
plot(lags, acf,'b');
hold on;

[acf, lags] = autocorr(r3, 3000);
plot(lags, acf,'k');
hold on;

[acf, lags] = crosscorr(r1, r3, 3000);
plot(lags, acf,'r')

corrcoef(r1, r3)