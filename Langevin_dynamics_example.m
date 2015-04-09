% Langevin dynamics
% d^2(x)/dt^2 + beta*dx/dt + dV/dt = noise(t)
% where V is the potential well
%
% For a harmonic well, dV/dt = k*x, that is, V = k*x^2/2, a quadractic
% potential well
%
% Note that the velocity autocorrelation function (ACF) can be solved
% analytically for the above Langevin equation. The volcity ACF is
% oscillatory with an exponentially dacaying envelope.
% The larger the noise, the quick the decay, a phenomena called "phase
% diffusion".
% That is, noise causes the lose of phase memory of the oscillation.
%
% Reference: The Langevin Equation: with Application to Stochastic
% Problems in Physics, Chemisty and Electrical Engineering (2nd), 2004,
% Chp. 3, p241-243
%
% Thanks for Yang Qi.
%
% I should do more stochastic simulations like this!
% The exercises help to built intuitions and insights!

clc;clear; close all;


period = 20/3; % 2*pi/period = sqrt(k/m)
noise_strength = 0.1;

beta = 0.1; % damping is essential: recall fluctuation-dissipation theorem.
% without the damping term, the system will diverge (numerical or dynamical divergence??)

v0 = 0;
x0 = 0;
% initial condition does not matter since the system is purely
% noise-driven.

dt = 0.01;
Ttot = 1000*period;
T = 0:dt:Ttot;

x = zeros(size(T));
v = zeros(size(T));

noise = noise_strength * randn(size(T));

x(1) = x0;
v(1) = v0;
for i = 2:length(T);
    a = -(2*pi/period)^2*x(i-1) - beta*v(i-1) + noise(i)*sqrt(dt); % be careful about the sqrt(dt)!!!
    v(i) = v(i-1) + a*dt  ;
    x(i) = x(i-1) + v(i-1)*dt ;
end


[cc, lags] = autocorr(v, ceil(100/dt));
lags = lags*dt;


figure(1);
set(gcf,'color', 'w');
t_sample = round(400*period/dt):round(420*period/dt);

subplot(3,1,1);
plot(T(t_sample), v(t_sample));
xlabel('time');
ylabel('displacement');

subplot(3,1,2);
plot(x(t_sample), v(t_sample));

subplot(3,1,3);
plot(lags,cc);
xlabel('time');
ylabel('velocity ACF');

