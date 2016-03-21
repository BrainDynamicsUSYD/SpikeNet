function [ degree_in, degree_out ] = logn_in_out_degree( N, deg_mean, deg_std, in_out_corrcoef  )
%[ degree_in, degree_out ] = logn_in_out_degree( N, deg_mean, deg_std, in_out_corrcoef  )
%   deg_mean = [deg_mean_in, deg_mean_out]
%   deg_std = [deg_std_in, deg_std_out]

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
mismatch_tol = min(deg_mean)/2;
mismatch = mismatch_tol + 1;
while mismatch > mismatch_tol
    [ degree_in_out ] = ceil(my_logn_rand( deg_mean, deg_std, ...
        [1 in_out_corrcoef; in_out_corrcoef, 1], N, 'mu_sigma_log' ));
    mismatch = abs(sum(degree_in_out(1:end-1,1)) - sum(degree_in_out(1:end-1,2)));
end
% Regenerate the last pair to best match sum(in) == sum(out)
% If there is still any mismatch, manually correct it.
% A lousy technique, cannot believe it works fine
[ degree_in_out_re ] = ceil(my_logn_rand( deg_mean, deg_std, ...
    [1 in_out_corrcoef; in_out_corrcoef, 1], N, 'mu_sigma_log' ));
mismatch = sum(degree_in_out(1:end-1,1)) - sum(degree_in_out(1:end-1,2));
corr = degree_in_out_re(:,2) - degree_in_out_re(:,1);
[resi, corr_best_ind] = min( abs(mismatch - corr) );
last_pair = degree_in_out_re(corr_best_ind,:);
corr_best = last_pair(2) - last_pair(1);
if resi > 0 && mismatch > corr_best
    last_pair(2) =  last_pair(2) + (mismatch - corr_best);
    degree_in_out(end,:) =  last_pair;
elseif resi > 0 && mismatch < corr_best
    last_pair(1) =  last_pair(1) + (corr_best - mismatch);
    degree_in_out(end,:) =  last_pair;
else
    degree_in_out(end,:) =  last_pair;
    if sum(degree_in_out(:,1) - degree_in_out(:,2)) ~= 0
        degree_in_out(end,:) =  flipud(last_pair);
    end
end


degree_in = ceil(degree_in_out(:,1));
degree_out = ceil(degree_in_out(:,2));

end

