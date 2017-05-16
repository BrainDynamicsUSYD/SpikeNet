function  spike_hist_shuffle = raster_marginals_shuffling(spike_hist)
%%%%%%%%%%%%%% shuffling
% Implementation of the Raster Marginals Model, as introduced in
% "Population rate dynamics and multineuron firing patterns in
%  sensory cortex", Journal of Neuroscience.
%
% Basically, the shuffling preserves the average firing rate of each neuron
% the temporal population firing rate.
% In other words, the so-called raster marginals model is a 0/1 matrices
% with given marginals.
%
sample_num = length(spike_hist(:,1));
if sample_num > 100
    error('Sample is too large for shuffling!')
end
spike_hist_shuffle = spike_hist';
shuffle_num = 100*nchoosek(sample_num,2);
c = ceil(rand(shuffle_num,2)*sample_num); % two randomly selected columns
for i = 1:shuffle_num
    I = spike_hist_shuffle(:,c(i,1)) + spike_hist_shuffle(:,c(i,2)) == 1; % where the 2 columns don't coincide
    cA = spike_hist_shuffle(I,[c(i,1) c(i,2)]); % a copy of the part that matters, to make it run faster
    i01 = find(cA(:,1) == 0);
    i10 = find(cA(:,1) == 1);
    toFlip = ceil(min(length(i01), length(i10))/2); % how many 01s & 10s to flip
    i01 = i01(randperm(length(i01)));
    i01 = i01(1:toFlip);
    i10 = i10(randperm(length(i10)));
    i10 = i10(1:toFlip);
    % the flip itself:
    cA(i01,1) = true; cA(i01,2) = false;
    cA(i10,1) = false; cA(i10,2) = true;
    spike_hist_shuffle(I,[c(i,1) c(i,2)]) = cA;
end;
spike_hist_shuffle = spike_hist_shuffle'; % the shuffled spike history
end
