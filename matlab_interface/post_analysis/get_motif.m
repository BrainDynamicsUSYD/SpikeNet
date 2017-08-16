function R = get_motif(R)

disp('Getting motif...')
warning('Lots of hard coded numbers here.')

n_sub = 100; % number of substitute network
n_sample = 1000;

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


[p, bi_p, uni_p, z_bi, z_uni] = uni_bi_connect_prob(W);

f = get_motif_quadruple(W);

f_sub = zeros(13,1);
for i = 1:n_sub
    
    W_sub = MyRandomGraphGenerator('E_R_uni_bi','N', N_sample, 'uni_p', uni_p, 'bi_p', bi_p);
    f_sub_tmp = get_motif_quadruple(W_sub);
    f_sub = f_sub + f_sub_tmp;
end
f_sub = f_sub / n_sub;

R.motif.bi_p = bi_p;
R.motif.uni_p = uni_p;
R.motif.z_bi = z_bi;
R.motif.z_uni = z_uni;
R.motif.f = f;
R.motif.f_sub = f_sub;



end

function f = get_motif_quadruple(W)
n_sample = 1000;
[~, n] = size(W);
f = zeros(13,1);
for i = 1:n_sample
    s = randperm(n, 4);
    f_tmp = motif3struct_bin(W(s,s));
    f = f + f_tmp;
end
f = f / n_sample;
end


