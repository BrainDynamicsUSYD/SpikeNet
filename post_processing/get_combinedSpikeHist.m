function get_combinedSpikeHist(Result_cell)
% save combined spike history for Adam's model
disp('Saving combined spike history for Adam ELIF model...');
Result_num = length(Result_cell);
for r_num = 1:Result_num
    R = Result_cell{r_num};
    spikehist = sparse(R.Analysis.spike_hist_combined);
    save([R.stamp(1:end-3),'ComSpikeHis.mat'],'spikehist','-v7.3')
end