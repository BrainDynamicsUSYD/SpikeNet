
function visualize_grid_and_SWR(R, LFP_centre_x, LFP_centre_y)
close all;
clc;

dt = R.dt;
ind_a_vec = R.grid.raw.ind_ab(1,:);
ind_b_vec = R.grid.raw.ind_ab(2,:);

hw = 31;
fw = 2*hw+1;
x_mean_chosen = R.grid.centre(1,:);
y_mean_chosen = R.grid.centre(2,:);
dist_std_chosen = R.grid.radius;
t_mid = R.grid.raw.t_mid;
t_mid_chosen = R.grid.t_mid;
jump_dist = R.grid.jump_dist;

[Lattice, ~] = lattice_nD(2, hw);
x_pos_o = Lattice(R.spike_hist_compressed{1}, 1);
y_pos_o = Lattice(R.spike_hist_compressed{1}, 2);


figure('Name','Vis','color','w','NumberTitle','off');
axis equal;
box on;
set(gca,'xtick',[],'ytick',[]);

xlim([-hw hw]);
ylim([-hw hw]);
hold on;


%
cm = colormap;
ms_transient = 200;
data_tmp = R.LFP.LFP_ripple_hilbert(:, round(ms_transient/dt)+1:end);
cm_range = minmax(data_tmp(:)');

for ii = 1:length(LFP_centre_x)
    elec(ii) = plot(LFP_centre_x(ii), LFP_centre_y(ii), 'b.', 'MarkerSize',69);
end


j = 1;

ang=0:0.01:2*pi;
x_shift_vs = [0 fw fw -fw -fw fw -fw 0 0 ];
y_shift_vs = [0 fw -fw fw -fw 0  0   fw -fw];

for tt = unique(sort([1:10:t_mid(end) t_mid(:)']))
    
    ripple_power = R.LFP.LFP_ripple_hilbert(:,tt);
    for ii = 1:length(LFP_centre_x)
        color_tmp = cm(min(ceil((ripple_power(ii)-cm_range(1))/diff(cm_range)*64),64),:);
        set(elec(ii), 'Color',color_tmp);
    end
    
    % ripple dectection
    detceded = sum(R.LFP.ripple_event.is_SWR(:, tt), 2) > 0;
    h4 = plot(LFP_centre_x(detceded), LFP_centre_y(detceded), '.k', 'MarkerSize',15);
    
    
    
    t = find(tt == t_mid);
    
    if ~isempty(t)
        if exist('h1','var')
            delete(h1);clear h1;
        end
        if exist('h2','var')
            delete(h2);clear h2;
        end
        if exist('h3','var')
            delete(h3);clear h3;
        end
        %t_range_tmp = (t_mid(t)-win_len/2):(t_mid(t) + win_len/2 -1);
        
        ind_range_tmp = ind_a_vec(t):ind_b_vec(t);
        h1 = plot(x_pos_o(ind_range_tmp), y_pos_o(ind_range_tmp), 'bo');
        
        if sum(t_mid_chosen == t_mid(t)) == 1
            
            x_tmp = x_mean_chosen(j);
            y_tmp = y_mean_chosen(j);
            r_cos = x_tmp+dist_std_chosen(j)*cos(ang);
            r_sin = y_tmp+dist_std_chosen(j)*sin(ang);
            
            h3 = plot(r_cos - x_shift_vs(1),r_sin - y_shift_vs(1),'b', ...
                r_cos - x_shift_vs(2),r_sin - y_shift_vs(2),'b',...
                r_cos - x_shift_vs(3),r_sin - y_shift_vs(3),'b',...
                r_cos - x_shift_vs(4),r_sin - y_shift_vs(4),'b',...
                r_cos - x_shift_vs(5),r_sin - y_shift_vs(5),'b', ...
                r_cos - x_shift_vs(6),r_sin - y_shift_vs(6),'b', ...
                r_cos - x_shift_vs(7),r_sin - y_shift_vs(7),'b', ...
                r_cos - x_shift_vs(8),r_sin - y_shift_vs(8),'b', ...
                r_cos - x_shift_vs(9),r_sin - y_shift_vs(9),'b');
            
            if j == 1
                h2 = plot( x_tmp, y_tmp, 'r>', 'MarkerSize', 8);
            else
                h2 = plot( x_tmp, y_tmp, 'rx', 'MarkerSize', 16);
                xs = sprintf('jump_size = %4.2f, radius = %4.2f', jump_dist(j-1), dist_std_chosen(j));
                xlabel(xs);
            end
            
            
            j = j + 1;
        end
        
    end
    drawnow;
    pause(0.02);
    
    delete(h4);
    ts = sprintf('time = %8.1f ms', tt*dt);
    title(ts);
    
end


end
