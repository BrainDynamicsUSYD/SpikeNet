function R = get_EI_current_crosscorr( R )
disp('Getting EI current crosscorrelation');

dt = R.dt;
XCC = [];

for pop = 1:2
    n = length(R.neuron_sample.neuron_ind{pop});
    
    for i = 1:n
        x = R.neuron_sample.I_AMPA{pop}(i,:);
        y = -(R.neuron_sample.I_GABA{pop}(i,:));
        x = x(1:10:end); % down-sampling
        y = y(1:10:end);
        [xcc,lags] = crosscorr(x,y,round(100/dt)); % do NOT use xcorr!!
        
        if i == 1 && pop == 1
            XCC = xcc;
        else
            XCC = [XCC; xcc];
        end
        
    end
 
end

R.Balance.EI_crosscorr = XCC;
R.Balance.lag = lags;

end
