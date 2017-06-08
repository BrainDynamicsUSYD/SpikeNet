function R = get_weight_vs_CC( R ) 

cc = R.Analysis.CC_pop{1};
pairs = R.Analysis.CC_pop_pairs{1};

% get connection weight matrix
s = strsplit(R.stamp,'_');
FID = [s{1},'_in.h5'];
[syn_cell] = read_syn_h5(FID, 1, 1, 1);

I = syn_cell{1}.I;
J = syn_cell{1}.J;
K = syn_cell{1}.K;
clear syn_cell;
W = sparse(I,J, K, 63^2, 63^2); % binary
clear I J K;

% collect spike-train correlation coefficients versus connection weights
cc_comp = cell(1,20);
w_comp = cell(1,20);
bin_a = -1.0:0.1:0.9;
bin_b = bin_a + 0.1;
for i = 1:length(cc)
    ind = find(cc(i) >= bin_a & cc(i) < bin_b);
    if ~isempty(ind)
        cc_comp{ind} = [cc_comp{ind} cc(i) cc(i)];
        w_comp{ind} = [w_comp{ind} W(pairs(1,i),pairs(2,i)) W(pairs(2,i),pairs(1,i))];
    end
end

% output resultss
R.w_CC.cc = cc_comp;
R.w_CC.w = w_comp;
R.w_CC.bin_ab = {bin_a, bin_b};

end


