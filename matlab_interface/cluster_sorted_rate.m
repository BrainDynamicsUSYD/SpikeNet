function R = cluster_sorted_rate( R )
    % get cluster rate
    R = get_cluster_rate( R );
    dt = R.reduced.dt;
    C_rate = R.cluster.rate;
        
    % sort the rate and each time point
    [sorted_rate, cluster_sequence] = sort(C_rate, 'descend');
    
    % cluster rank
    [a_ind, b_ind] = size(cluster_sequence);
    cluster_rank = zeros(a_ind,b_ind);
    for i = 1:b_ind
        cluster_rank(cluster_sequence(:,i),i) = (1:a_ind)';
    end

    % threshold the 1st symbolic sequence based on 1st sorted rate
    theta = 5:1:15; % Hz, threshold
    lt = length(theta);
    seq_1st = repmat(cluster_sequence(1,:), lt, 1); 
    for i = 1:lt
        seq_1st(i, sorted_rate(1,:) < theta(i)  ) = 0; % or NaN?
    end
    
    
    % statistical analysis on the thresholded 1st symbolic sequence
    switch_seq = cell(1,lt);
    switch_freq = zeros(1,lt);
    high_du = cell(1,lt);
    low_du = cell(1,lt);
    switch_level = cell(1,lt);
    for i = 1:lt
        [seq_tmp, high_du_tmp, low_du_tmp] = seq_postprocess( seq_1st(i,:), dt );
        
        switch_level{i} = switch_level_detect(seq_tmp);
        switch_seq{i} = seq_tmp;
        high_du{i} = high_du_tmp;
        low_du{i} = low_du_tmp;
        
        % switch frequency (order parameter?)
        switch_freq(i) = length(seq_tmp)/(sum([high_du_tmp low_du_tmp]));
    end
    
    % output results
    R.cluster.threshold = theta;
    R.cluster.switch_freq = switch_freq;
    R.cluster.switch_seq = switch_seq;
    R.cluster.switch_level = switch_level;
    R.cluster.high_du = high_du;
    R.cluster.low_du = low_du;
    R.cluster.sorted_rate = sorted_rate;
    R.cluster.sym_seq = cluster_sequence;
    R.cluster.sym_seq_1st_theta = seq_1st;
    R.cluster.rate_rank = cluster_rank;
    
end



function [switch_level] = switch_level_detect(seq)

level_num = ceil(log2(max(seq))) + 1; % total number of hierarchy levels
seq = seq - 1; % zero based index

% convert to binary strings
seq_bin = cell(size(seq));
for i = 1:length(seq)
    seq_bin{i} = dec2bin(seq(i), level_num);
end

% find level based on binary strings
switch_level = zeros(1, length(seq)-1);
for i = 1:(length(seq)-1)
    cmp = (seq_bin{i} == seq_bin{i+1});
    for j = 2:level_num
        if cmp(j) == 1 && cmp(j-1) == 0
            cmp(j) = 0;
        end
    end
    switch_level(i) = (level_num+1) - find(cmp,1,'last'); 
end

end
