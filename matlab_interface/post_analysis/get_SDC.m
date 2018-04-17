function R = get_SDC(R)

disp('Getting SDC...')


n_sample = 2000;

L = lattice_nD(2, 31);
I_sample = find(  L(:,1).^2+L(:,2).^2 <= n_sample/pi );
N_sample = length(I_sample);



s = strsplit(R.stamp,'_');
FID = [s{1},'_in.h5'];
[syn_cell] = read_syn_h5(FID, 1, 1, 1);

I = syn_cell{1}.I;
J = syn_cell{1}.J;
K = syn_cell{1}.K;
clear syn_cell;
W = full(sparse(I,J,ones(size(K)), 63^2, 63^2)); % binary
W = W(I_sample, I_sample); % downsample
clear I J ;

SDC_m = 1000; % sample size

SDC_mu = [];
SDC_sig = [];
for SDC_n = 3:12
    s = [];
    for j = 1:SDC_m
        ind = randperm(N_sample, SDC_n);
        s_tmp = get_sdc(W(ind,ind));
        s = [s s_tmp]; %#ok<*AGROW>
    end
    SDC_mu = [SDC_mu nanmean(s)];
    SDC_sig = [SDC_sig nanstd(s)];
end

R.SDC.n = 3:12;
R.SDC.mu = SDC_mu;
R.SDC.sig = SDC_sig;

end



function s = get_sdc(W)
k_in = sum(W,1);
k_out = sum(W,2);

s = corrcoef(k_in(:),k_out(:));
s = s(1,2);

end


