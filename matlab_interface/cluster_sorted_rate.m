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
    R.cluster.rate_rank = cluster_rank;
    
end



function [switch_seq, high_du, low_du] = seq_postprocess(seq, dt)
% symbolic sequence postprocesing
% start from simple solutions!

% cut head and tail
if ~isempty(seq)
    head = seq(1);
    h = 0;
    for i = 1:length(seq)
        if seq(i) == head
            h = i;
        else
            break;
        end
    end
    seq(1:h) = [];
end

if ~isempty(seq)
    tail = seq(end);
    t = 0;
    for i = length(seq):-1:1
        if seq(i) == tail
            t = i;
        else
            break;
        end
    end
    seq(t:end) = [];
end

% cut low state from head and tail
if ~isempty(seq)
    head = seq(1);
    if head == 0 % low state
        h = 0;
        for i = 1:length(seq)
            if seq(i) == head
                h = i;
            else
                break;
            end
        end
        seq(1:h) = [];
    end
end
if ~isempty(seq)
    tail = seq(end);
    if tail == 0 % low state
        t = 0;
        for i = length(seq):-1:1
            if seq(i) == tail
                t = i;
            else
                break;
            end
        end
        seq(t:end) = [];
    end
end

% contract sequence and find durations
if isempty(seq)
    switch_seq = [];
    du = [];
else
    switch_seq = seq(1);
    du = 1;
    for i = 2:length(seq)
        if seq(i) == seq(i-1)
            du(end) = du(end)+1;
        elseif seq(i)*seq(i-1) == 0 % a normal transition (high-low or low-high)
            du = [du 1]; % new duration counter
            switch_seq = [switch_seq seq(i)]; % new entry in switch sequence
        else % a direction transition (take-over: high-high)
            du = [du 0 1]; % insert a low state but with zero length
            switch_seq = [switch_seq 0 seq(i)]; % insert a low state but with zero length
        end
    end
end

% differentiate high and low states
high_du = du(switch_seq > 0)*dt;
low_du = du(switch_seq == 0)*dt;
switch_seq(switch_seq == 0) = []; % only high states

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
