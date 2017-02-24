function [R] = get_SWR_spike_phase_lock(R)

lag_ms = 10; %ms
lag_step = round(lag_ms/R.reduced.dt);
% loop through electrodes
hil_phase = [];
spike_count = [];
rp_amp_no_std  = [];
spike_sort_range = 6;
try
spike_sort_range = R.LFP.ripple_event.spike_sort_range;
catch
end


dt_conv = R.dt/R.reduced.dt;
for i = 1:length(R.LFP.LFP_ripple(:,1))
    no_std = R.LFP.wavelet.peak.rp_amp_no_std{i};
    s_tmp = R.ExplVar.LFP_range_sigma;
    spike_sort_neurons = R.LFP.LFP_neurons{1}(i,:) >= 1/(s_tmp*sqrt(2*pi))*exp(-0.5*(spike_sort_range/s_tmp)^2);
    spike_count_sort = full(sum(R.reduced.spike_hist{1}(spike_sort_neurons,:)));
    
    hil = hilbert(R.LFP.LFP_ripple(i,:));
    hil_phase0 = atan2(imag(hil),real(hil));
    
    a = R.LFP.ripple_event.ripple_start_steps{i};
    b = R.LFP.ripple_event.ripple_start_steps{i} + R.LFP.ripple_event.ripple_du_steps{i} - 1;
    
    
    for j = 1:length(a)
        s_ab = round(a(j)*dt_conv):round(b(j)*dt_conv);
        xx = hil_phase0( round(s_ab/dt_conv));
        yy = spike_count_sort(s_ab-lag_step );
        zz = R.LFP.LFP_ripple(i, round(s_ab/dt_conv));
        %             figure(1);
        %             plot(zz,'b')
        %             hold on;
        %             plot(yy,'r')
        %             plot(xx,'g')
        %             pause;
        %             close all
        hil_phase = [hil_phase xx]; %#ok<AGROW>
        spike_count = [spike_count yy]; %#ok<AGROW>
        rp_amp_no_std = [rp_amp_no_std no_std(j)*ones(size(yy))]; %#ok<AGROW>
    end
    
    
end

R.SWR_spike_phase_lock.hil_phase = hil_phase;
R.SWR_spike_phase_lock.spike_count  = spike_count;
R.SWR_spike_phase_lock.rp_amp_no_std  = rp_amp_no_std;
R.SWR_spike_phase_lock.lag_ms  = lag_ms;

end