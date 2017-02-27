function [new_mean, new_var, step_new] = nan_Welford_online(new_data, old_mean, old_var, step_old, is_end, op)

switch lower(op)
    case 'nan'
        g = ~isnan(new_data);
    case 'zero'
        g = new_data > 0;
end

new_mean = old_mean;
new_var = old_var;
% online mean and var calculation: Welford's method (1962, Technometrixcs)
new_mean(g) = old_mean(g) + (new_data(g) - old_mean(g)) ./ step_old(g);
new_var(g) = old_var(g) + (new_data(g) - old_mean(g)) .* (new_data(g) - old_mean(g));

step_new = step_old;
step_new(g) = step_new(g) + 1;

if is_end == true
    new_var = new_var ./ (step_new - 1);
end

end