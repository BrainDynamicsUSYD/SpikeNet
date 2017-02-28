function plot_neuron_stats(R)

figure(1)

%%%%%%%%%%%%
subplot(1,2,1);
mu = R.neuron_stats.I_tot_time_mean{1};
sig = sqrt(R.neuron_stats.I_tot_time_var{1});

V_th = -50;
V_lk = -70;
g_lk = 0.0167;
thre = (V_th - V_lk)*g_lk;

box on
hold on;
m = linspace(0, thre, 100);
s = linspace(0, 1, 100);
plot(mu(:), sig(:),'.');

plot(thre*ones(size(m)), s,'r--');
plot(m, thre-m,'r--');
xlim([0 0.4]);
ylim([0 1]);
xlabel('\mu_I (nA)')
ylabel('\sigma_I (nA)')


%%%%%%%%%%%%
subplot(1,2,2);hold on;box on;
r = R.neuron_stats.IE_ratio{1};
e = linspace(min(r), max(r), 30);
h = histc(r,e);
bh = bar(e,h);
set(bh,'FaceColor','None')
plot(mean(r),1.2*max(h),'*');
xlabel('IE-ratio')

end