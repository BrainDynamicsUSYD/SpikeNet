function [new_mean, new_var] = Welford_online(new_data, old_mean, old_var, step_now, is_end)

% online mean and var calculation: Welford's method (1962, Technometrixcs)
new_mean = old_mean + (new_data - old_mean) ./ step_now;
new_var = old_var + (new_data - old_mean) .* (new_data - old_mean);

if is_end == true
    new_var = new_var / (step_now - 1);
end

end