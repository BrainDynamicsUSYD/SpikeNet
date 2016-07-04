close all;
tau_decay = 7.0;
tau_rise = 0.5;
K = 1.5;

% Calculate the peak time in Adam's model
% Note that in Yifan's model the peak time is simply tau_rise
tau_peak = log(tau_decay/tau_rise)*(tau_decay*tau_rise)/(tau_decay-tau_rise) %#ok<*NOPTS>
% f2 = 1/(-exp(-tau_peak/tau_rise) + exp(-tau_peak/tau_decay)) % this makes
% the peak 1
f1 = 1/(tau_decay - tau_rise) % this makes the integral 1

dt = 0.1;
step_tot = 10000;
s1 = zeros(1,step_tot);
s2 = zeros(1,step_tot);
s3 = zeros(1,step_tot);
s4 = zeros(1,step_tot);
s4_trans =  zeros(1,step_tot);
s5 = zeros(1,step_tot);
s4_trans(1) = 1; % scales the peak magnitude

s6 = zeros(1,step_tot);
s7 = zeros(1,step_tot);
s8 = zeros(1,step_tot);
s6(1) = 1*K;
s7(1) = 1*K;

I1 = zeros(1,step_tot);
I3 = zeros(1,step_tot);
I4 = zeros(1,step_tot);


for i=1:step_tot-1
    t = i*dt;
    % Yifan's model
    if t <= tau_rise
        s1(i+1) = s1(i) + (1-s1(i))/(tau_rise/dt);
        I1(i) = s1(i+1);
        s1(i+1) = s1(i+1)*exp(-dt/tau_decay);
    else
        s1(i+1) = s1(i)*exp(-dt/tau_decay);
        I1(i) = s1(i);
    end
    % Adam's model (alpha function)
    s2(i) = (exp(-t/tau_decay)-exp(-t/tau_rise))*f1;
        
    % Yifan's model after conversion
    if t <= tau_peak
        s3(i+1) = s3(i) + (1-s3(i))/(tau_peak/dt); % use tau_peak instead!!!
        I3(i) = s3(i+1);
        s3(i+1) = s3(i+1)*exp(-dt/tau_decay);
    else
        s3(i+1) = s3(i)*exp(-dt/tau_decay);
        I3(i) = s3(i);
    end
    
    % Adam's model with kinetic equations (check numerical error!!!)
    s4(i+1) = s4(i) + s4_trans(i)*dt;
    I4(i) = s4(i+1);
    s4(i+1) = s4(i+1)*exp(-dt/tau_decay);
    s4_trans(i+1) = s4_trans(i)*exp(-dt/tau_rise);
    
    dsdt = -s5(i)/tau_decay +  s4_trans(i);
    s5(i+1) =  s5(i) + dsdt*dt;
    s4_trans(i+1) = s4_trans(i)*exp(-dt/tau_rise);
    
    % James' idea
    % s8(i) = s6(i)-s7(i);
    s6(i+1) = s6(i)*exp(-dt/tau_decay);
    s7(i+1) = s7(i)*exp(-dt/tau_rise);
    s8(i+1) = s6(i+1)-s7(i+1);
end
T = (1:step_tot)*dt;
weight_factor = trapz(T, s2)/trapz(T, s3) % keep the integral the same 
figure(1);hold on;
set(gcf,'color','w');
plot(T, s1, 'r', T, s2, 'b', T, s3*weight_factor, 'r:', ...  % scale the coupling weights with this factor!!!
    T, s4,'g', ...
    T, s5,'c', ...
    T, I3*weight_factor, 'k:', ...
    T, s8/K*f1,'y'); 
xlim([0 30])
legend('Yifan''s','Adam''s','Yifan''s after conversion', 'Adam''s with kinetics');

numerical_error = trapz(T, I3)/trapz(T, s3)-1


% max_s4 = max(s4)

area_I4 = trapz(T, I4)

