function [ seq ] = persistence_requirement( seq, min_steps )
% consecutive zeros are low states
% consecutive other numbers are high states
% the persistence requirement is that the high states cannot be shorter
% than some min steps, otherwise the high states should be regarded as low
% states.
% to apply this condition, simply make those high states not long enough to
% be all zeros

n = length(seq);

i_first = find(seq(1:end-1), 1,'first'); % first non-zero element

if ~isempty(i_first) && i_first < n
    %  initialize length counter and starting index for high state
    high_du_tmp = 1;
    i_begin = i_first;
    
    for i = (i_first+1):n
        % in the middle of a high state
       if seq(i) == seq(i-1) && seq(i) ~= 0 
           high_du_tmp = high_du_tmp + 1; % keep counting the length of it
           
       % at the end of a high state
       elseif seq(i) == 0 && seq(i-1) ~= 0 
           if high_du_tmp < min_steps % if not meet the persistence requirement
               seq( i_begin:(i-1) ) = 0; % delete this high state
           end
           
       % at the beginning of a high state
       elseif seq(i) ~=0 && seq(i-1) == 0 
           high_du_tmp = 1; % reset length counter 
           i_begin = i; % reset the starting index
           
       % direct switch between two high states
       elseif seq(i) ~= seq(i-1) && seq(i) ~= 0  && seq(i-1) ~= 0  
            if high_du_tmp < min_steps % if not meet the persistence requirement
               seq( i_begin:(i-1) ) = 0; % delete the previous high state
            end
            high_du_tmp = 1; % reset length counter 
            i_begin = i; % reset the starting index
           
       % in the middle of a low state
       elseif seq(i) == 0 && seq(i-1) == 0 
           continue; % do nothing
       end
           
    end
    
    
end

end