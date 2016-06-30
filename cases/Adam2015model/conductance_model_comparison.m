close all;
tau_decay = 2.0;
tau_rise = 0.5;

dt = 0.05;
step_tot = 1000;
s1 = zeros(1,step_tot);
s2 = zeros(1,step_tot);
s3 = zeros(1,step_tot);

% Calculate the peak time in Adam's model
% Note that in Yifan's model the peak time is simply tau_rise
tau_peak = log(tau_decay/tau_rise)*(tau_decay*tau_rise)/(tau_decay-tau_rise) %#ok<*NOPTS>

for i=1:step_tot-1
    t = i*dt;
    % Yifan's model
    if t <= tau_rise
        s1(i+1) = s1(i) + (1-s1(i))/(tau_rise/dt);
        s1(i+1) = s1(i+1)*exp(-dt/tau_decay);
    else
        s1(i+1) = s1(i)*exp(-dt/tau_decay);
    end
    % Adam's model
    s2(i) = (exp(-t/tau_decay)-exp(-t/tau_rise))/(tau_decay-tau_rise);
        
    % Yifan's model after conversion
    if t <= tau_peak
        s3(i+1) = s3(i) + (1-s3(i))/(tau_peak/dt); % use tau_peak instead!!!
        s3(i+1) = s3(i+1)*exp(-dt/tau_decay);
    else
        s3(i+1) = s3(i)*exp(-dt/tau_decay);
    end
end
T = (1:step_tot)*dt;
weight_factor = trapz(T, s2)/trapz(T, s3) % keep the integral the same 
figure(1);hold on;
set(gcf,'color','w');
plot(T, s1, 'r', T, s2, 'b', T, s3*weight_factor, 'r:'); % scale the coupling weights with this factor!!!
xlim([0 30])
legend('Yifan''s','Adam''s','Yifan''s after conversion');


