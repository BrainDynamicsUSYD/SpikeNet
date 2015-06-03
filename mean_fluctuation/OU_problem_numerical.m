%clc;clear;

% constants
Vr = 10;
Vth = 20;
tau_ref = 2;
tau_syn = 4; % this could be problematic 
% because there are actually two synaptic time constants, 
% that is, tau_I = 3 and tau_E = 5


% parameters (4 is too many?)
tau_m = 10;
Vss = 12;
gL = 0.1;
sigma_c = 1.4;


%%%%%%%%%%%%%% Calculate the firing rate predicted by 
% reference : Firing Frequency of Leaky Integrate-and-fire Neurons with Synaptic Current Dynamics

% get sigma_V_eff
sigma_V = sigma_c./((tau_m.^0.5).*gL); % I added this equation because I believe there is a mistake in (15.62)

% get V_th_eff & V_r_eff
% This first order correction is in good agreement with the results from
% numerical simulations for tau_syn < 0.1*tau_m
alpha = 1.03;
k = (tau_syn ./ tau_m).^0.5; % this is my own approximation, which can be problematic; 
Vth_eff_1st = (Vth - Vss)./sigma_V + alpha .* k;
Vr_eff_1st = (Vr - Vss)./sigma_V + alpha .* k;
% Second order correction (should be good for 0<k<1)
beta0 = -0.35;
beta1 = 0.27;
beta = beta0 + beta1*(Vth - Vss)./sigma_V;
Vth_eff_2nd = (Vth - Vss)./sigma_V + (alpha.*k  + beta.*k.^2);
Vr_eff_2nd = (Vr - Vss)./sigma_V + (alpha.*k  + beta.*k.^2);
% get upper & lower limits of the integral
upper = Vth_eff_2nd;
lower = Vr_eff_2nd;
% calculate firing rate 
fx = @(x) exp(x.^2).*(1+erf(x));
v_new = 1/(tau_ref + tau_m*sqrt(pi)*integral(fx, lower, upper)); % can matlab integral() be counted on?!
rate_theory = v_new * 1000 % Hz





%%%%%%%%%%%%%%%% numerical results
% % the problem can be perfectly cast into a stochastic damped harmonic oscillator
% % reference:
% % Harmonic oscillator in heat bath: Exact simulation of time-lapse-recorded data
% % and exact analytical benchmark statistics

ntrails = 2000;
dt = 0.1;
fpt = zeros(1, ntrails); % first passage time

m = tau_m; % mass
gamma = 1+tau_m/tau_syn; % friction coefficient
k = 1/tau_syn; % Hooke's constant
w0 = sqrt(k/m);
tau = m/gamma;

% damping_ratio = gamma/(2*sqrt(m*k))
damping_ratio = 0.5*(sqrt(tau_syn/tau_m) + sqrt(tau_m/tau_syn)); % always >= 1
% note that according to the above expression
% the harmonic oscillator is always overdamped, 
% except when tau_m = tau_syn, which gives a critically damped case.



w =  sqrt( w0^2 - 1/(4*tau^2) ); % cyclic freq of the damped oscillator
D = sigma_c^2/(2*gL^2*(tau_syn+tau_m)^2); % simple algebra by me

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










h = waitbar(0,'Please wait...');
for i = 1:ntrails
    X = [Vr - Vss;  % X = V - Vss, simple variable substitution by me
         randn()*sigma_c];
    T = tau_ref;
    while X(1) + Vss < Vth
        xi = randn();
        zeta = randn();
        dX = [ sigma_xx*xi; % eq (13)
            sigma_xv^2/sigma_xx*xi + sqrt(sigma_vv^2 - sigma_xv^4/sigma_xx^2)*zeta]; % eq (14)
        X = mtimes(e_Mdt, X) + dX; % eq (7)
        T = T + dt; % record time 
    end
    fpt(i) = T/1000;  % sec
    if mod(i,10) == 0
        waitbar(i/1000,h)
    end
end
% close(h);
delete(h);

rate_simu = length(fpt)/sum(fpt) % Hz: number of spikes over total time


