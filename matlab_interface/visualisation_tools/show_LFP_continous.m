function show_LFP_continous(R, choice)
dt = R.dt;
stamp = R.stamp;
samp_file = [stamp(1:end-3) '0_neurosamp'];


down_sample = 20;

if nargin == 1
    choice = 1;
end

if choice == 3
    load(samp_file, 'V');
    [N, steps] = size(V);
    X = reshape(V, [sqrt(N) sqrt(N) steps]);
    clear V;
elseif choice == 2
    load(samp_file, 'ripple_power_grid');
    X = ripple_power_grid;
    clear ripple_power_grid;
elseif choice == 1
    load(samp_file, 'LFP_power_grid');
    X = ripple_power_grid;
    clear LFP_power_grid;
end


[N_s, ~, steps] = size(X);
N = N_s^2;
if N ~= R.N(1)
    error('Not all of the excitatory neurons are sampled!')
end
w = sqrt(N);
hw = (w-1)/2;
[Lattice, ~] = lattice_nD(2, hw);
load(samp_file, 'peak');


figure('NumberTitle','off','Name','Ripple Hilbert Power','color', 'w');
mm = minmax(reshape(X(:),1,[]));
xlim([-hw hw]);
ylim([-hw hw]);
hold on;
h1 = imagesc(-hw:hw, -hw:hw,X(:,:,1), mm);
% h1 = imagesc(-hw:hw, -hw:hw,X(:,:,1));
s = find(R.spike_hist{1}(:,1));
h2 = plot(Lattice(s,2),Lattice(s,1),'r.','MarkerSize',15);

t = find(R.neuron_sample.t_ind{1});

for i = 1:down_sample:steps
    s = find(R.spike_hist{1}(:,t(i)));
    pause(0.01)
    set(h1,'CData',X(:,:,i));
    set(h2,'XData',Lattice(s,2),'YData',Lattice(s,1));
    xlabel([num2str(t(i)*dt ),' ms']);
end

end