function C_rate_plot( R, seg, seg_size )
%UNTITLED6 Summary of this function goes here
%   Detailed explanation goes here

% Input check and default values
if nargin < 3
    seg_size = 4*10^4; % 2*10^4 for 2-pop, segmentation size for each plot
end
if nargin < 2
    seg = 1;
end

% Dump fields
dt = R.dt;
step_tot = R.step_tot;
C_rate = R.C_rate;

% Segmetation
seg_ind = get_seg(step_tot, seg_size, seg);

T = seg_ind*dt;


% find right shift
shift = 0.1;
max_rate = max(max(C_rate));
rate_shift = max_rate*(1+shift);

% plot
hold on;
for i = 1:length(C_rate(:,1));
    plot(T, C_rate(i,seg_ind)+(i-1)*rate_shift);
end

axis off;

end
