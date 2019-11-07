function R = SpikeNet2SpikeGrid_spikehist(R)
% transfer the spike history's data format from SpikeNet to SpikeGrid
disp('SpikeNet2SpikeGrid, spike history...');
tic;
Num_neu = R.N;
gsize = sqrt(sum(Num_neu));%sqrt(Num_neu);
if floor(gsize)==gsize
    spike_hist = R.spike_hist;
    MemP = false(gsize);
    MemP(2:2:gsize,2:2:gsize) = 1;
    ind_E = MemP == 0;
    ind_I = MemP == 1;
    [~,end_step] = size(spike_hist{1});
    S = false(gsize^2,end_step);
    for i=1:1:end_step
        MemP(ind_E) = spike_hist{1}(:,i);
        MemP(ind_I) = spike_hist{2}(:,i);
        S(:,i) = MemP(:);
    end
    % if one want save memory, S can be convert to sparse matrix but slower
    R.Analysis.spike_hist_combined = sparse(S);
    toc;
else % no-square, should be specified manually
    spike_hist = R.spike_hist;
    gsizex=30;
    gsizey=180;
    MemP = false(gsizex,gsizey);
    MemP(2:2:gsizex,2:2:gsizey) = 1;
    ind_E = MemP == 0;
    ind_I = MemP == 1;
    [~,end_step] = size(spike_hist{1});
    S = false(gsizex*gsizey,end_step);
    for i=1:1:end_step
        MemP(ind_E) = spike_hist{1}(:,i);
        MemP(ind_I) = spike_hist{2}(:,i);
        S(:,i) = MemP(:);
    end
    % if one want save memory, S can be convert to sparse matrix but slower
    R.Analysis.spike_hist_combined = sparse(S);
    toc;
end
end

