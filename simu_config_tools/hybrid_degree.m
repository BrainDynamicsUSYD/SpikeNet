 function [ degree_in, degree_out ] = hybrid_degree( N, deg_mean, deg_std, r, hybrid )
%[ degree_in, degree_out ] = logn_in_out_degree( N, deg_mean, deg_std, in_out_corrcoef  )
%   deg_mean = [deg_mean_in, deg_mean_out]
%   deg_std = [deg_std_in, deg_std_out]

mismatch_tol = deg_mean/2;
mismatch = mismatch_tol + 1;

while abs(mismatch) > mismatch_tol
    deg = hybrid_pois_logn(deg_mean, deg_std, r, N, hybrid);
    mismatch = sum(deg(:,1)) - sum(deg(:,2));
end

% randomly add 1's to compensate the mismatch 
% (this should have minimal effect on the correlation)
adjust_ind = randperm(N, abs(mismatch));
if mismatch > 0
    deg(adjust_ind, 2) = deg(adjust_ind ,2) + 1;
elseif mismatch < 0
    deg(adjust_ind, 1) = deg(adjust_ind, 1) + 1;
end


% % Regenerate the last pair to best match sum(in) == sum(out)
% % If there is still any mismatch, manually correct it.
% % A lousy technique, cannot believe it works fine
% [ deg_re ] = hybrid_pois_logn(deg_mean, deg_std, r, N, hybrid);
% mismatch = sum(deg(1:end-1,1)) - sum(deg(1:end-1,2));
% corr = deg_re(:,2) - deg_re(:,1);
% [resi, corr_best_ind] = min( abs(mismatch - corr) );
% last_pair = deg_re(corr_best_ind,:);
% corr_best = last_pair(2) - last_pair(1);
% if resi > 0 && mismatch > corr_best
%     last_pair(2) =  last_pair(2) + (mismatch - corr_best);
%     deg(end,:) =  last_pair;
% elseif resi > 0 && mismatch < corr_best
%     last_pair(1) =  last_pair(1) + (corr_best - mismatch);
%     deg(end,:) =  last_pair;
% else
%     deg(end,:) =  last_pair;
%     if sum(deg(:,1) - deg(:,2)) ~= 0
%         deg(end,:) =  flipud(last_pair);
%     end
% end

degree_in = ceil(deg(:,1));
degree_out = ceil(deg(:,2));

 end

 
 function deg = hybrid_pois_logn(deg_mean, deg_std, r, N, hybrid)
     
 deg_logn = ceil(my_logn_rand( [deg_mean deg_mean], deg_std, [1 r; r, 1], N, 'mu_sigma_log' ));
 deg_pois = poissrnd_2D( [deg_mean deg_mean], N, r);
 
 logn_ind = randperm(N, round(N*hybrid));
 
 deg = deg_pois;
 if hybrid > 0
    deg(logn_ind, :) = deg_logn(logn_ind, :);
 end
 
 
 end

