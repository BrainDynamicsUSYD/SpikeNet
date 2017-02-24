function Hz = get_Hz_from_Eq_34(sigma_I, tau_m, gL, Vss, tau_syn, tau_ref, Vr, Vth)
% warning::
% this function has not been thoroughly testes against known theoretical
% results!
% 
% this function solves the firing rate given by the system described by
% Eq.34 in [1] Gu & Gong 2016, The dynamics of memory retrieval in
% hierarchical networks: a modeling study 
%
% Also see [2] Harmonic oscillator in heat bath: exact simulation of
% time-lapse-recorded data and exact analytical benchmark statistics 

nspikes = 5000;
dt = 1; % 1ms instead of 0.1!!
fpt = zeros(1, nspikes); % first passage time

m = tau_m; % mass
gamma = 1+tau_m/tau_syn; % friction coefficient
k = 1/tau_syn; % Hooke's constant
w0 = sqrt(k/m);
tau = m/gamma;

% damping_ratio = gamma/(2*sqrt(m*k))
% damping_ratio = 0.5*(sqrt(tau_syn/tau_m) + sqrt(tau_m/tau_syn)); % always >= 1
% note that according to the above expression
% the harmonic oscillator is always overdamped, 
% except when tau_m = tau_syn, which gives a critically damped case.

% some constants defined by [2]
w =  sqrt( w0^2 - 1/(4*tau^2) ); % cyclic freq of the damped oscillator
D = sigma_I^2/(2*gL^2*(tau_syn+tau_m)^2); % Einstein's constant

% % note that when tau_syn and tau_m becomes similar, w blows up everything
% % in the following expressions
% J = [1/(2*w*tau)  1/w;
%     -w0^2/w  -1/(2*w*tau)]; % eq (10)
% e_Mdt = exp(-dt/(2*tau))*( cos(w*dt)*eye(2) + sin(w*dt)*J ); % eq (9)
% sigma_xx_2 = D/(4*w^2*w0^2*tau^3)*...
%     (4*w^2*tau^2 + exp(-dt/tau) * (cos(2*w*dt) - 2*w*tau*sin(2*w*dt) - 4*w0^2*tau^2)); % eq (15)
% sigma_xx = sqrt(sigma_xx_2);
% sigma_vv_2 = D/(4*w^2*tau^3)*...
%     (4*w^2*tau^2 + exp(-dt/tau) * (cos(2*w*dt) + 2*w*tau*sin(2*w*dt) - 4*w0^2*tau^2)); % eq (16)
% sigma_vv = sqrt(sigma_vv_2);
% sigma_xv = sqrt(  D/(w^2*tau^2)*exp(-dt/tau)*sin(w*dt)^2 ); % eq (17)


% To avoid blowing up near tau_syn = tau_m (damping ratio = 1), use the 
% following expressions which make use of sinc() function in matlab 
% for more info, type help sinc
J = [dt/(2*tau)  dt;
    -w0^2*dt  -dt/(2*tau)];
e_Mdt = exp(-dt/(2*tau))*( cos(w*dt)*eye(2) + sinc(w*dt/pi)*J );
sigma_xx_2 = D/(w0^2*tau) - D*exp(-dt/tau)/(w0^2*tau^3)*...
    (0.5*dt^2*(sinc(w*dt/pi))^2 + tau*dt*sinc(2*w*dt/pi) + tau^2);
sigma_xx = sqrt(sigma_xx_2);
sigma_vv_2 = D/tau - D*exp(-dt/tau)/(tau^3)*...
    (0.5*dt^2*(sinc(w*dt/pi))^2 - tau*dt*sinc(2*w*dt/pi) + tau^2);
sigma_vv = sqrt(sigma_vv_2);
sigma_xv = sqrt(  dt^2*D*exp(-dt/tau)/(tau^2)*(sinc(w*dt/pi))^2 ); 

sigma_Y_inf = sqrt(  sigma_I^2/(2*gL^2*(tau_syn+tau_m)*tau_syn*tau_m)   );

for i = 1:nspikes
    X = [Vr - Vss;  % X = V - Vss, simple variable substitution by me
         randn()*sigma_Y_inf];
    T = tau_ref;
    while X(1) + Vss < Vth
        xi = randn();
        zeta = randn();
        dX = [ sigma_xx*xi; % eq (13)
            sigma_xv^2/sigma_xx*xi + sqrt(sigma_vv^2 - sigma_xv^4/sigma_xx^2)*zeta]; % eq (14)
        X = mtimes(e_Mdt, X) + dX; % eq (7)
        T = T + dt; % record time 
        if 1000/T < 0.01 % less than 0.01Hz
            T = nan;
            break
        end
    end
    if isnan(T)
        fpt = nan;
        break;
    end
    fpt(i) = T/1000;  % sec
end

if isnan(fpt)
    Hz = 0;
else
    Hz = length(fpt)/sum(fpt); % Hz: number of spikes over total time
end

end