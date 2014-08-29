function visualise_state_on_layout(p, state_now, state_old)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

% reset previous state
for i = 1:length(state_old)
    set(p{round(state_old(i))}, 'FaceAlpha', 0.03);
end
% show current state
for i = 1:length(state_now)
    set(p{round(state_now(i))}, 'FaceAlpha', 0.5);
end

end

