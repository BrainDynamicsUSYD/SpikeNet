function [ sorted_rate, cluster_sequence, cluster_rank ] = cluster_rate_sort( rate )
    % sort the rate and each time point
    [sorted_rate, cluster_sequence] = sort(rate, 'descend');
    
    % cluster rank
    if nargout == 3
        [a_ind, b_ind] = size(cluster_sequence);
        cluster_rank = zeros(a_ind,b_ind);
        for i = 1:b_ind
            cluster_rank(cluster_sequence(:,i),i) = (1:a_ind)';
        end
    end
    
end

