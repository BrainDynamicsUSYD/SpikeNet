function [expl1, expl2, obsv_mean, obsv_std] = collect_expl_obsv_mat(expl_in1,expl_in2, obsv)
disp('Use imagesc( expl1, expl2, transpose(obsv_mean))!')
expl1 = unique(expl_in1);
expl2 = unique(expl_in2);
obsv_mean = zeros(length(expl1), length(expl2));
obsv_std = zeros(length(expl1), length(expl2));
for i = 1:length(expl1)
    for j = 1:length(expl2)
        obsv_mean(i,j) =  nanmean(obsv(expl1(i) == expl_in1 & expl2(j) == expl_in2 ));
        obsv_std(i,j) =  nanstd(obsv(expl1(i) == expl_in1 & expl2(j) == expl_in2 )) ;
    end
end



