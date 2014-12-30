function plot_cluster_detection(R, seg, varargin)


%
shadingMode = 2; % 1 for only shading the corresponding cluster, 2 for the the entire population
i_theta = 2;
text_fontsize = 12;

% Input check and default values
if nargin < 2
    seg = 1;
end

seg_size = 4*10^4; % 2*10^4 for 2-pop, segmentation size for each plot


for i = 1:(length(varargin)/2)
    eval([varargin{i*2-1}, '=', num2str(varargin{i*2}), ';' ]);
end


% Dump fields
dt = R.reduced.dt;
step_tot = R.reduced.step_tot;
Mnum = R.ExplVar.Mnum;
% 
dt = dt/1000;

% Segmetation
seg_ind = get_seg(step_tot, seg_size, seg);


% R = cluster_sorted_rate(R);


pop = 1;


set(gca, 'xlim', [min(seg_ind), max(seg_ind)]*dt, 'XAxisLocation', 'bottom', ...
    'ylim', [1 500]);


seq_1st = R.cluster.sym_seq_1st_theta(i_theta,seg_ind);
cutoff = 0;
dt_tmp = 1; % since we are only interested in the indices
[switch_seq, high_du, ~, high_start, ~] = seq_postprocess(seq_1st, dt_tmp, cutoff);
my_ylim = get(gca,'ylim');
ymin = my_ylim(1);
ymax = my_ylim(2);

% color shading for detected high states



xc = seg_ind*dt;
yc = linspace(ymin+1,ymax, 500);
yc_dM = round((ymax-ymin)/Mnum);


[xc, yc] = meshgrid(xc, yc);
zc = zeros(size(xc));
for i = 1:length(high_du)
    t1 = seg_ind(high_start(i));
    t2 = t1 + high_du(i)-1;
    c_tmp = switch_seq(i);
    if shadingMode == 1
        zc(xc >= t1*dt & xc <= t2*dt  & yc > ymin+(c_tmp-1)*yc_dM & yc <= ymin+c_tmp*yc_dM ) = 1;
    elseif shadingMode == 2
        zc(xc >= t1*dt & xc <= t2*dt ) = 1;
    end
end
contourf(xc,yc,zc,[0 1],'LineStyle','none')
colormap([1 1 1;0 1 0]);

hold on;


sample_color = []; %ceil((1:500)/round(500/8))/8;
raster_plot(R, pop, seg, sample_color,'text_fontsize',text_fontsize);
set(gca, 'xlim', [min(seg_ind), max(seg_ind)]*dt, 'XAxisLocation', 'bottom');
xlabel('Time (sec)','fontsize',text_fontsize)%,'FontWeight','Bold')
ylabel('Neurons','fontsize',text_fontsize)%,'FontWeight','Bold');
box on;
set(gca, 'ytick', [])
hold on;



end
