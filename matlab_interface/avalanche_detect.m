function avalanche = avalanche_detect(R, plot_figure)
% Let's do it!!!
% Reference:
% [1] Benayoun et al, 2010, Avalanches in a Stochastic Model of Spiking
% Neurons


% wframe = 4; % (ms) frame length

if nargin == 1
    plot_figure = 0;
end

% dump fields
num_spikes = R.num_spikes{1};
dt = R.dt;
step_tot = R.step_tot;



% average spike intervals in the network (use population maybe?)
p_rate = sum(num_spikes)/step_tot;

step_sep = ceil(p_rate);

% segementation based on step_sep

if step_sep == 1
    warning('step_sep is only one!')
end

[A1, B1, D1, D0] = cut_binary( (num_spikes >= 1), step_sep );
ava_duration = D1(:)'; % avalanche durations
IAI = D0(:)'; % inter-avalance intervals
ava_num = length(A1); % number of avalanches
ava_tt = cell(1,ava_num); % avalanche time-course
ava_size = zeros(1,ava_num);
for i = 1:ava_num
    ava_tt{i} = num_spikes( A1(i):B1(i) );
    ava_size(i) = sum(ava_tt{i});
end



avalanche.size = ava_size;
avalanche.dt_sep = step_sep*dt;
avalanche.duration = ava_duration;
avalanche.interval = IAI;
avalanche.time_course = ava_tt;


if plot_figure == 1
%     figure('NumberTitle','off','Name','Avalanche Detection','color','w',...
%         'Position', [680   824   376   271], 'units','pixels');
    
    bin_size = 1;
    bin_edge = 1:bin_size:max(ava_size);
    bin_count = histc(ava_size,bin_edge);
    bin_count = bin_count/sum(bin_count);
    
    
    poisson_size_dist = (1-exp(p_rate)).^(1:1:max(ava_size))*exp(p_rate);
    poisson_size_dist = poisson_size_dist(poisson_size_dist >= min(bin_count(bin_count>0)) );

    plot(bin_edge,bin_count,'.');
    set(gca,'xscale','log','yscale','log');
    ylim([10^floor(log10(min(bin_count))), 1]);
    hold on;
    plot(1:length(poisson_size_dist),poisson_size_dist,'r--');

    avalanche.poisson_size_dist = poisson_size_dist;
    
end





end







function [A1, B1, D1, D0, A0, B0] = cut_binary(bin, D0_min)
% cut 1001110100000000000000000111011 using consecutive zeros with minimum
% length D0_min
% The heading and tailing parts of data are discarded
% A: start, B: end, D: duration

padding = false(1,1);
% begining-ending points (A&B points) detection
A1 = find([padding, bin(1:end-1) == 0 & bin(2:end) == 1]);  % if ...01111..., detect first 1
B1 = find([bin(1:end-1) == 1 & bin(2:end) == 0, padding]); % if ...11110...., detect last 1
[A1, B1] = pruneAB(A1,B1);
D1 = B1-A1;

% apply D0_min
for i = 1:length(A1)-1
    D0_tmp = A1(i+1) - B1(i) - 1;
    if D0_tmp < D0_min
        A1(i+1) = [];
        B1(i) = [];
        D1(i) = D1(i) + D1(i+1) + D0_tmp;
        D1(i+1) = [];
    end
end

% find A0,B0,D0
A0 = B1(1:end-1)+1;
B0 = A1(2:end)-1;
D0 = B0-A0;

end



function[A_pair, B_pair] = pruneAB(A,B)
% find A-B pairs
% cut off dangling A or B

if ~isempty(A) && ~isempty(B)
    % pairing of point A's and point B's
    A_pair = A;
    B_pair = B;
    if length(B_pair) > length(A_pair)
        % B A-B A-B A-B
        B_pair(1) = [];
    elseif length(B_pair) < length(A_pair)
        % A-B A-B A-B A
        A_pair(end) = [];
        
    elseif B_pair(1) < A_pair(1) && length(B_pair) == length(A_pair) % the latter is always true but kept for readability
        % B A-B A-B A-B A
        B_pair(1) = [];
        A_pair(end) = [];
    end
else
    disp('A and B points are empty!');
    A_pair = [];
    B_pair = [];
end

end


