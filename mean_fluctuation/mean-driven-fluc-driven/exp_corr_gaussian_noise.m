function r = exp_corr_gaussian_noise(N,tau)
% % Gaussian Markov process?
% % Reference:
% % How to generate exponentially correlated Gaussian random numbers
% % by Markus Deserno
% % URL: http://www.cmu.edu/biolphys/deserno/pdf/corr_gaussian_random.pdf
% e.g.
% N = 10000;
% tau = 100; % autocorrelation time constant in time steps
% 
% another useful reference is 
% Generating partially correlated noise—A comparison of methods
% by Hartmann and Cho

f = exp(-1/tau);
g = randn(1,N);

r = zeros(1,N);
r(1) = g(1);
for i = 2:N
    r(i) = f*r(i-1) + sqrt(1-f^2)*g(i);
end

r = (r - mean(r))/std(r); % is this necessary?
end

