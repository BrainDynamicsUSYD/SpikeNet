function R = get_lagged_cov(R)
% this functions calculates the squared Frobenius norm of M, the matrix of
% lagged covariances among neurons.
% see Fast sampling-based Inference in Balanced Neuronal Networks

sample_size = 400;
sample_neurons = randperm(R.N(1), sample_size);
lagged_cov_win_ms = 50; % ms
len = round(lagged_cov_win_ms/R.reduced.dt);
% r = movmean(transpose(R.reduced.spike_hist{1}), len, 'Endpoints','discard');
r = movingmean_glen(transpose(R.reduced.spike_hist{1}(sample_neurons,:)), len, 1, 1);
r = r(ceil(len/2):end-ceil(len/2),:);

lags_ms = linspace(0, 200, 10);
lagged_cov = [];

for lag_tmp = lags_ms
    lag_step = round(lag_tmp/R.reduced.dt);
    lagged_cov_tmp = cov([r(lag_step+1:end,:) r(1:end-lag_step,:)]);
    lagged_cov = [lagged_cov sum(sum(lagged_cov_tmp(1:end/2,end/2:end).^2))]; %#ok<AGROW>
    lag_tmp/max(lags_ms) %#ok<NOPRT>
end
lagged_cov = lagged_cov/lagged_cov(1);

R.Analysis.lagged_cov = lagged_cov;
R.Analysis.lagged_cov_win_ms = lagged_cov_win_ms;
R.Analysis.lagged_cov_lags_ms = lags_ms;

end