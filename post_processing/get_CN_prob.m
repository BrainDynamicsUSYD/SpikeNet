function [R] = get_CN_prob(R,varargin)

n_sample = 2000;
% 12-neuron
n_samp = 2000;
s_samp = 12;


for i = 1:length(varargin)/2
    var_name = varargin{2*i-1};
    var_value = varargin{2*i};
     if isnumeric(var_value)
        eval([var_name, '=', num2str(var_value), ';']);
     else
         eval([var_name, '=''', var_value, ''';']);
     end
end

s = strsplit(R.stamp,'_');
FID = [s{1},'_in.h5'];
[syn_cell] = read_syn_h5(FID, 1, 1, 1);

L = lattice_nD(2, 31);
I_sample = find(  L(:,1).^2+L(:,2).^2 <= n_sample/pi );

I = syn_cell{1}.I;
J = syn_cell{1}.J;
K = syn_cell{1}.K;
clear syn_cell;
W = full(sparse(I,J,K, 63^2, 63^2));
A = W(I_sample, I_sample);
A = double(A>0);
clear I J K W;


cn_tot = [];
A_tot = [];
cn_p = zeros(1,5);
for n = 1:n_samp
    i_samp = randperm(length(I_sample), s_samp);
    A_samp = A(i_samp, i_samp);
    cn = A_samp'* A_samp; % A * A' gives common post-synaptic neighbours, note that cn is symmetric
    
    A_samp(logical(eye(size(A_samp)))) = NaN;
    A_samp = A_samp(:);
    A_samp(isnan(A_samp)) = [];
    A_tot = [A_tot; A_samp(:)];
    
    cn(logical(eye(size(cn)))) = NaN;
    cn = cn(:);
    cn(isnan(cn)) = [];
    
    cn_tot = [cn_tot; cn(:)]; %#ok<*AGROW>
end
[c,bin] = histc(cn_tot,-0.5:1:4.5);

for n = 1:5
    cn_p(n) = cn_p(n) + sum(A_tot(bin == n));
end
cn_p = cn_p./c(1:5)';

R.CN_prob.cn_p = cn_p;
R.CN_prob.I_sample = I_sample;
R.CN_prob.n_samp = n_samp;
R.CN_prob.s_samp = s_samp;

end