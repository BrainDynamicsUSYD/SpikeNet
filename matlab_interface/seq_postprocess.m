
function [switch_seq, high_du, low_du, high_start, low_start] = seq_postprocess(seq, dt, cutoff)
% symbolic sequence postprocesing
% start from simple solutions!

if nargin == 2
    cutoff = 1;
end
h_fin = 0;

if cutoff == 1
    % cut head and tail
    if ~isempty(seq)
        head = seq(1);
        h = 0;
        for i = 1:length(seq)
            if seq(i) == head
                h = i;
            else
                break;
            end
        end
        seq(1:h) = [];
        h_fin = h_fin + h;
    end
    
    if ~isempty(seq)
        tail = seq(end);
        t = 0;
        for i = length(seq):-1:1
            if seq(i) == tail
                t = i;
            else
                break;
            end
        end
        seq(t:end) = [];
    end
    
    % cut low state from head and tail
    % why did I do this????
    % oh, because I'm only interested in the transitions between high states
    if ~isempty(seq)
        head = seq(1);
        if head == 0 % low state
            h = 0;
            for i = 1:length(seq)
                if seq(i) == head
                    h = i;
                else
                    break;
                end
            end
            seq(1:h) = [];
            h_fin = h_fin + h;
        end

    end
    if ~isempty(seq)
        tail = seq(end);
        if tail == 0 % low state
            t = 0;
            for i = length(seq):-1:1
                if seq(i) == tail
                    t = i;
                else
                    break;
                end
            end
            seq(t:end) = [];
        end
    end
end


% contract sequence and find durations
if isempty(seq)
    switch_seq = [];
    du = [];
    start = [];
else
    switch_seq = seq(1);
    du = 1;
    start = 1;
    for i = 2:length(seq)
        if seq(i) == seq(i-1)
            du(end) = du(end)+1;
        elseif seq(i)*seq(i-1) == 0 % a normal transition (high-low or low-high)
            start = [start i];
            du = [du 1]; % new duration counter
            switch_seq = [switch_seq seq(i)]; % new entry in switch sequence
        else % a direction transition (take-over: high-high)
            start = [start 0 i];
            du = [du 0 1]; % insert a low state but with zero length
            switch_seq = [switch_seq 0 seq(i)]; % insert a low state but with zero length
        end
    end
end

start = start + h_fin; % account for the cut-off in the head

% differentiate high and low states
high_du = du(switch_seq > 0)*dt;
high_start = start(switch_seq > 0)*dt;
low_du = du(switch_seq == 0)*dt;
low_start = start(switch_seq == 0)*dt;
switch_seq(switch_seq == 0) = []; % only high states

end
