function [expl, obsv_mean, obsv_std] = collect_expl_obsv_vector(expl_in, obsv)

expl = unique(expl_in);
obsv_mean = [];
obsv_std = [];
for e = expl
    obsv_mean = [obsv_mean nanmean(obsv(e == expl_in))]; %#ok<AGROW>
    obsv_std = [obsv_std nanstd(obsv(e == expl_in))]; %#ok<AGROW>
end


end



