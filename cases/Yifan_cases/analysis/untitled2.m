
% input: 
mu_input % see Eq 20 gL_eff
sigma_input


% define synapse
td = 3; % ms; synaptic decay time constants
tr = 1; % ms; synaptic rise time constants
Vrev = 0;% ms 






% define parameters
Vth = -50; % mV; threshold
Vr = -60; % mV; reset 
V_L = -70;
tau_ref = 2; % ms;
gL = 0.0167; % miuS
Cm = 0.25; %nF

% conversion factor  from rectangular pulse to 1st order dynamics for the
% neurotransmitter
A = tr/(td+tr) + td^2/(td+tr)^2.* (1-exp(-1-tr/td));


% get sbar (below Eq.29)
sbar = (v > 0) .* tau_decay .* v .* A;  % this linear approximation could be problematic; consider use emperical results
if sum(sbar > 0.5) > 0
    warning('Linear sbar approximation could be inaccurate due to large firing rate')
end




% get Vbar
Vbar = Vss - (Vth-Vr).*tau_m_eff.*v - (Vss-Vr).*tau_ref.*v; % (15.52) ?

% get Vss
Vss = (gL*V_L + mtimes(C.*g, sbar.*Vrev))./gL_eff; 

% get gL_eff
gL_eff = gL + mu_input;

% get tau_m_eff
tau_m_eff = Cm./gL_eff;

Hz = get_Hz_from_Eq_34(sigma_I, tau_m_eff, gL_eff, Vss, tau_syn_eff, tau_ref, Vr, Vth);













