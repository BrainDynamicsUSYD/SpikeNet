function [X, pmf] = pmf_yg(pmf_name, para, x_range, N)
% [X] = pmf_yg(pmf_name, para, x_range, N)
% generate random numbers given the discretized x_range
% subject to the specified discretized probabilty mass function 
% defined within that range
%
% pmf_name='logn'
% para = [mu, sigma]
% mu and sigma: mean and std of the lognormal distribution
%
% pmf_name='power'
% para = [alpha]
% alpha: exponent (absolute value)
%
% pmf_name='exp'
% para = [ tau ]
% tau: "time"-constant
%
% x_range: the discretized values that the random number can take
% N: number of samples in the output X
%
% Note that the actual mean and std of the generated sample may be
% different from the give values because of the finite range of x_range
%
% performance:
% with N ~ 10^4, it takes 0.005 sec to generate X, which is OK.

x_range = sort(x_range); % must be ascending

switch lower(pmf_name)
    case 'logn'
        mu = para(1);
        sigma = para(2);
        mu_n = 0.5*log(mu^4/(mu^2+sigma^2));
        sigma_n = sqrt(log((mu^2+sigma^2)/mu^2));
        pmf_unnorm = 1./x_range.*exp( -(log(x_range)-mu_n).^2/(2*sigma_n^2)); % un-normalized
        pmf = pmf_unnorm./(sum(pmf_unnorm));
        cmf = cumsum(pmf);
    case 'power'
        alpha = para(1);
        pmf_unnorm = x_range.^(-alpha);
        pmf = pmf_unnorm./(sum(pmf_unnorm));
        cmf = cumsum(pmf);
    case 'exp'
        tau = para(1);
        pmf_unnorm = exp(-x_range./tau);
        pmf = pmf_unnorm./(sum(pmf_unnorm));
        cmf = cumsum(pmf);
end


u = sort(rand(1, N)); % must be ascending
X = zeros(1,N);
% Inverse transform sampling (discretized version)
i_cmf = 1;
for i_u = 1:N
    while u(i_u) > cmf(i_cmf)
        i_cmf = i_cmf + 1;
    end
    X(i_u) = x_range(i_cmf);
end
X = X(randperm(N)); % shuffle the order

end
