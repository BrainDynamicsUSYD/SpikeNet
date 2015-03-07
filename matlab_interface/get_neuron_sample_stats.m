function R = get_neuron_sample_stats( R )
% autocorrelation is calculated instead of power spectral density
% note that Fourier transform of the autocorrelation function of a signal is the power spectrum of the signal

dt = R.dt;

disp('Getting neuron sample data statistics');

% XCC_EI = [];
for pop = 1:2
    n = length(R.neuron_sample.neuron_ind{pop});
    
    for i = 1:n
       
        E = R.neuron_sample.I_AMPA{pop}(i,:) + R.neuron_sample.I_ext{pop}(i,:); % I_ext is always excitatory
        I = -(R.neuron_sample.I_GABA{pop}(i,:));
        V = R.neuron_sample.V{pop}(i,:);

        E = E(1:10:end); % down-sampling
        I = I(1:10:end);
        V = V(1:10:end);
        
        [xcc_EI,lags] = crosscorr(E,I, round(100/dt)); % do NOT use xcorr!!
        [acc_E,~] = autocorr(E, round(100/dt) );
        [acc_I,~] = autocorr(I, round(100/dt) );
        [acc_V,~] = autocorr(V, round(100/dt) );
        
        if i == 1 && pop == 1
            XCC_EI = xcc_EI;
            ACC_E = acc_E;
            ACC_I = acc_I;
            ACC_V = acc_V;
            EI_ratio = mean(E)/mean(I);
        else
            XCC_EI = [XCC_EI; xcc_EI];
            ACC_E = [ACC_E; acc_E];
            ACC_I = [ACC_I; acc_I];
            ACC_V = [ACC_V; acc_V];
            EI_ratio = [EI_ratio mean(E)/mean(I)];
        end
        
    end
end

R.neuron_sample_stats.lag = lags;
R.neuron_sample_stats.EI_crosscorr = XCC_EI;
R.neuron_sample_stats.E_autocorr = ACC_E;
R.neuron_sample_stats.I_autocorr = ACC_I;
R.neuron_sample_stats.V_autocorr = ACC_V;
R.neuron_sample_stats.EI_ratio = EI_ratio;



end
