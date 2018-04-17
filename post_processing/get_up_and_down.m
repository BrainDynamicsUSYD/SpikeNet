function [a, b] = get_up_and_down(sh, dt, varargin )
% sh: spike history with a dt = 1
% a: up-state onset time step
% b: down-state onset time step
%
% Ref: Sequential structure of neocortical spontaneous activity in vivo
% (Supplementary material)

% default parameters
if dt ~= 1
    error('dt should be 1!')
end
sh_tot = full(sum(sh));
bin_down = 30;
bin_up = 60;
up_du = 100;
up_sp_tot = 50;
down_no = 1;
up_no = 15;

% read non-default parameters
for i = 1:(length(varargin)/2)
    eval([varargin{i*2-1}, '=', num2str(varargin{i*2}),';' ]);
end


% UP state onset was defined as the time of the first spike marking a
% transition from a period of global silence (30 ms with at most one spike
% of any cell) to a period of activity (60 ms with at least 15 spikes from
% any cells). 
sh_downsum = movsum(sh_tot, bin_down,'Endpoints','discard');
sh_downsum = [NaN(1, bin_down) sh_downsum]; % sh_downsum(N) is the sum from N-bin to N-1
sh_upsum = movsum(sh_tot, bin_up,'Endpoints','discard');
sh_upsum = sh_upsum(1:end-up_du);
sh_downsum = sh_downsum(1:length(sh_upsum));
% % debug
% figure(1);hold on;
% t_sample = 100:600;
% plot(sh_downsum,'r')
% plot(sh_upsum,'b');
% [Y,X,~] = find(sh(:,t_sample));
% xdata = ( [X(:)'; X(:)']-1+t_sample(1))*dt; % sec
% ydata = [Y(:)'-1;Y(:)']-1;
% line(xdata, ydata,'color','k','linewidth',0.1);
% ylim([0 50])
% xlim(minmax(t_sample))

up_onset = find(sh_downsum <= down_no & sh_upsum >= up_no);
up_onset_diff = [diff(up_onset) 1];
up_onset = up_onset(up_onset_diff ~= 1);

% UP state offset was defined as the time of the last spike
% marking a transition to at least 30 ms of silence. 
sh_downsum = movsum(sh_tot, bin_down,'Endpoints','discard');
sh_downsum = [sh_downsum(2:end) NaN(1, length(sh_tot)-length(sh_downsum)+1)];
up_offset = find(sh_downsum <= down_no & sh_tot >= 1) ;

% UP states shorter than 200 ms, or containing <50 spikes from all neurons,
% were excluded from further analysis      
a = [];
b = [];
for at = up_onset
    bt = up_offset( find( (up_offset - at) > 0, 1) ); % find the first offset after the onset
    if bt - at - 1 >= up_du
        if sum(sh_tot(at:bt)) >= up_sp_tot
            a = [a at]; %#ok<*AGROW>
            b = [b bt];
            
%             % check at and bt
%             if sum(sh_tot(at-bin_down:at-1)) > down_no || sum(sh_tot(at:at+bin_up-1)) < up_no || ...
%                 sum(sh_tot(bt-bin_up+1:bt)) < up_no || sum(sh_tot(bt+1:bt+bin_down)) > down_no
%                 error('e')
%             end
        end
    end
end

% trim overlapping on state
i_trim = [];
for i = 2:length(a)
    if a(i) - a(i-1) < up_du
        i_trim = [i_trim i];
    end
end
a(i_trim) = [];
b(i_trim) = [];



end