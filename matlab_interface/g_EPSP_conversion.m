function [ fit_g_2_EPSP, fit_EPSP_2_g ] = g_EPSP_conversion( varargin )
%[ fit_g_2_EPSP, fit_EPSP_2_g ] = g_EPSP_conversion( g ,varargin)
%   This function finds the EPSP amplitude given the coupling strength
%   (conductance)

% visualisze_trace = 1;

% default parameters
dt = 0.1; %ms

Cm = 0.25; % nF
V_lk = -70.0;   % Leak reversal, -70.0
g_lk = 0.0167; % muS

V_rev = 0; % for AMPA
tau_r = 1.0; %ms
tau_d = 5; %ms

V_th = -50; % threshold;

g_min = 0.001;
g_max = 0.010;

% I_K_on = 1;
% V_K = -85; % mV
% dg_K = 0.010; % muS
% g_K = dg_K;
% tau_K = 80; %ms

% read parameters
for i = 1:(length(varargin)/2)
    eval([varargin{i*2-1}, '=', num2str(varargin{i*2}), ';' ]);
end


g = linspace(g_min, g_max, 10);
g = g(:);

EPSP = zeros(size(g));
for i = 1:length(g)
    
    V_old = V_lk; % resting potential
    peaked = false;
    tau_r_step = round(tau_r/dt);
    tau_r_step_left = tau_r_step; % spike at time zero
    s = 0;
    while ~peaked && V_old < V_th
        if tau_r_step_left > 0
            s = s + 1/tau_r_step*(1-s);
            tau_r_step_left = tau_r_step_left - 1;
        end
        I = -g(i)*s*(V_old-V_rev);
        I_lk = -g_lk*(V_old-V_lk);
        
%         if I_K_on == 1 % This is not right because I_K is based on
%         % post-synaptic firing!!!s
%             I_K = -g_K*(V_old - V_K);
%             g_K = g_K * exp(-dt/tau_K);
%         else
%             I_K = 0;
%         end
        
        V_dot = (I + I_lk )/Cm;
        V_new = V_old+V_dot*dt;
        
        peaked = V_new < V_old;
        V_old = V_new;
        s = s*exp(-dt/tau_d);
        if V_old >= V_th
            V_old = NaN;
            break;
        end
    end
    
    EPSP(i) = V_old - V_lk;
    
end

if sum(isnan(EPSP)) > 0
    warning('Some given g''s are unrealistically large!')
end

fit_g_2_EPSP = fit(g,EPSP,'poly1');
fit_EPSP_2_g = fit(EPSP,g,'poly1');

end

