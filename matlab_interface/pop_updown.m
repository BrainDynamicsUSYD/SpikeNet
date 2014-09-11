

function result = pop_updown(R, pop_ind, theta)
% A: state beginning time vector index
% B: state endding time vector index

    % down parameter
    dt = R.reduced.dt;
    spike_hist = R.reduced.spike_hist{pop_ind};
    N = R.N;
    
    % data smoothing and thresholding
    kernel_width = 50; % ms, kernel length % kernel for rate estimation
    kernel = spike_train_kernel_YG(kernel_width, dt, 'gaussian');
    pop_rate = SpikeTrainConvolve(sum(spike_hist, 1)/N(pop_ind), kernel);
    result.theta = theta;
    c_bin = pop_rate > theta; % bin: binary
    
    % begining-ending points (A&B points) detection
    padding = false(1,1);
    beginning = [padding, c_bin(1:end-1) == 0 & c_bin(2:end) == 1];  % if ...01111..., detect first 1
    ending = [c_bin(1:end-1) == 1 & c_bin(2:end) == 0, padding]; % if ...11110...., detect last 1
    
    % up and down detection
    A_up = find(beginning); % beginning point of up state
    A_up = A_up(:)'; % row vector
    B_up = find(ending); % ending point of up state
    B_up = B_up(:)'; % row vector
    A_down = B_up;
    B_down = A_up;

    [result.up_duration, result.upA, result.upB] = duration_from_AB(A_up,B_up,dt);
    [result.down_duration, result.downA, result.downB] = duration_from_AB(A_down,B_down,dt);

end




function [duration, A_pair, B_pair] = duration_from_AB(A,B,dt)
duration = [];
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
    % calculate durations
    if ~isempty(A_pair) && ~isempty(B_pair)
        duration = (B_pair - A_pair)*dt;
    end
else
    disp('A and B points for up state are empty!');
    A_pair = [];
    B_pair = [];
end
end
