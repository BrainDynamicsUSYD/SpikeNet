function [R] = get_SWR_spike_phase_lock(R)


% loop through electrodes
hil_phase = [];
spike_count = [];
spike_sort_range = R.LFP.ripple_event.spike_sort_range;
% spike_sort_range = 12;


for i = 1:length(R.LFP.LFP_ripple(:,1))
    no_std = R.LFP.wavelet.peak.rp_amp_no_std{i};
    s_tmp = R.ExplVar.LFP_range_sigma;
    spike_sort_neurons = R.LFP.LFP_neurons{1}(i,:) >= 1/(s_tmp*sqrt(2*pi))*exp(-0.5*(spike_sort_range/s_tmp)^2);
    spike_count_sort = full(sum(R.spike_hist{1}(spike_sort_neurons,:)));
    
    hil = hilbert(R.LFP.LFP_ripple(i,:));
    hil_phase0 = atan2(imag(hil), real(hil));
    
    a = R.LFP.ripple_event.ripple_start_steps{i};
    b = R.LFP.ripple_event.ripple_start_steps{i} + R.LFP.ripple_event.ripple_du_steps{i} - 1;

    
    for j = 1:length(a)
        if no_std(j) > 7
            hil_phase = [hil_phase hil_phase0(a(j):b(j))]; %#ok<AGROW>
            spike_count = [spike_count spike_count_sort(a(j):b(j))]; %#ok<AGROW>
        end
    end

end

R.SWR_spike_phase_lock.hil_phase = hil_phase;
R.SWR_spike_phase_lock.spike_count  = spike_count;

end