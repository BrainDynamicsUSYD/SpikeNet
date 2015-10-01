% visual the spiking pattern of the spatially embedded network

N = R.N(1);
hw = (sqrt(N)-1)/2;
[Lattice, ~] = lattice_nD(2, hw);

spike_hist = R.reduced.spike_hist{1};
dt = R.reduced.dt;
step_tot = R.reduced.step_tot;

figure(1);
set(gcf, 'color','w');
xlim([-hw hw])
ylim([-hw hw])
for t = 1:step_tot
    spikes_t = find(spike_hist(:,t));
    
    plot(Lattice(spikes_t,1), Lattice(spikes_t,2), 'rx')
    pause(0.02);
end