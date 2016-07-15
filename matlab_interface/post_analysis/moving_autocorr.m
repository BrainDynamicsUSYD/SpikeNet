function [mid_points, lags, acc_mat] = moving_autocorr(x, window_size, lagNum, window_number)
% does this make sense? Yes!

% lagNum = round(window_size/4);

start_points = round(linspace(1, length(x)-window_size-1, window_number));
start_points = unique(start_points);
end_points = start_points + window_size;
mid_points = round((start_points + end_points)/2);

acc_mat = [];
for i = 1:length(end_points)
    [acc_tmp, lags] = autocorr( x(start_points(i):end_points(i)), lagNum );
    acc_tmp = acc_tmp(:);
    acc_mat = [acc_mat acc_tmp];
end

end