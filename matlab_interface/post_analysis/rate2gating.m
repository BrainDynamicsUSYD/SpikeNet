function func = rate2gating(tau_d, tau_r)
% covert firing rate to gating variable


% figure(1); hold on;
n_trial = 40;
x = [];
y = [];
for rate_Hz = 1:40
    rate = rate_Hz*10^-3; % KHz
    s_avg = [];
    for n = 1:n_trial
        dt = 0.1;
        step_tot = 10^5;
        s1 = zeros(1,step_tot);

        spike_counter = 0;
        for i=1:step_tot-1
            if rand < rate
                spike_counter = (tau_r/dt);
            end
            % Yifan's model
            if spike_counter > 0
                s1(i+1) = s1(i) + (1-s1(i))/(tau_r/dt);
                s1(i+1) = s1(i+1)*exp(-dt/tau_d);
            else
                s1(i+1) = s1(i)*exp(-dt/tau_d);
            end
            spike_counter = max(0,spike_counter-1);           
        end
        s_avg = [s_avg mean(s1)];
    end
    % plot(rate_Hz, mean(s_avg), 'o');
    x = [x rate_Hz]; %#ok<*AGROW>
    y = [y mean(s_avg)];
end

func = @(b,x)  (b*x./(1+b*x));
nlm = fitnlm(x',y',func,1);
b = nlm.Coefficients.Estimate;
% plot(sort(x), func(b,sort(x)),'b-');

func = @(x)  (b*x./(1+b*x));

end

