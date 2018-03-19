% visualize the spiking pattern of the spatially embedded network
function embed_network_spike_movie(R, make_video, manual)

if nargin == 1
    make_video = 0;
end

if make_video == 1
    writerObj = VideoWriter('default_system.avi');
    open(writerObj)
end

figure('NumberTitle','off','Name','Spike Movie','color','w')
   N = R.N(1);
hw = (sqrt(N)-1)/2;
[Lattice, ~] = lattice_nD(2, hw);

spike_hist = R.reduced.spike_hist{1};
dt = R.reduced.dt;

step_tot = min(R.reduced.step_tot, 10^3);

if nargin < 3
    manual = 0;
end

if isempty(findall(0,'Type','Figure'))
    figure('NumberTitle','Off','Name','Movie');
    set(gcf, 'color','w');
end

xlim([-hw hw])
ylim([-hw hw])
box on;
hold on;
h1 = plot(0, 0, 'rx');

t_bin = 1;

for t = 1:t_bin:(step_tot-t_bin)
    [spikes_t, ~, ~] = find(spike_hist(:,t:(t+t_bin)));
    
    delete(h1);
    h1 = plot(Lattice(spikes_t,1), Lattice(spikes_t,2), 'rx');
    if manual == 0
        pause(0.01);
    else
        pause
    end
    
    xlabel([num2str(t*dt/1000),' sec']);
    
    % M(t) = getframe(gcf);
    if make_video == 1
        writeVideo(writerObj, getframe(gcf));
    end
    
end
if make_video == 1
    close(writerObj);
end

end
