function [lags, xcc] = plot_EI_current_corr(R)
% plot E-I current correlation

text_fontsize = 12;


% Dump fields
dt = R.dt;
neuron_sample = R.neuron_sample;


XCC = [];
for pop = 1:length(neuron_sample.neuron_ind)
    n = length(neuron_sample.neuron_ind{pop});
    
    for i = 1:n
        x = neuron_sample.I_AMPA{pop}(i,:);
        y = -(neuron_sample.I_GABA{pop}(i,:));
        x = x(1:10:end); % down-sampling
        y = y(1:10:end);
        [xcc,lags] = crosscorr(x,y,round(100/dt)); % do NOT use xcorr!!
        %plot(lags*dt,xcc);
        
        if i == 1 && pop == 1
            XCC = xcc;
        else
            XCC = [XCC; xcc];
        end
        
    end
end

shadedErrorBar(lags*dt, mean(XCC), std(XCC) );
xlabel('lag (ms)','fontsize',text_fontsize);
ylabel('cross-correlation','fontsize',text_fontsize);

end


