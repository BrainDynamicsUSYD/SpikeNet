function [expl1, expl2, obsv_hist] = collect_expl_obsv_hist_array(expl_in1,expl_in2, obsv)

expl1 = unique(expl_in1);
expl2 = unique(expl_in2);
obsv_hist = cell(length(expl1), length(expl2));
for i = 1:length(expl1)
    for j = 1:length(expl2)
        L = find(expl1(i) == expl_in1 & expl2(j) == expl_in2);
        for l = L
            data_tmp = obsv{l};
            obsv_hist{i,j} = [obsv_hist{i,j} data_tmp(:)'];
        end
    end
end



