function R = get_stru_var_decomp(R)
%
warning('Lots of hard coded numbers here.')

V_rev = 0;
fw = 63;
N = fw^2;

load('stru_var_metadata.mat','stru_var')

stru_var_tmp = stru_var;
a_mean = stru_var_tmp.a_mean;
a_var = stru_var_tmp.a_var;

J_mean = mean(stru_var_tmp.J_mean(:));
J_var = mean(stru_var_tmp.J_var(:));

a_mean_mat = zeros(N, N);
a_var_mat = zeros(N, N);
for ind = 1:N
    [i,j] = ind2sub([fw fw], ind);
    a_tmp = circshift(reshape(a_mean,fw,fw), [i-round(fw/2), j-round(fw/2)]);
    %     imagesc(reshape(a_tmp,fw,fw));
    %     pause(0.1)
    a_mean_mat(ind,:) = a_tmp(:)';
    a_tmp = circshift(reshape(a_var,fw,fw), [i-round(fw/2), j-round(fw/2)]);
    a_var_mat(ind,:) = a_tmp(:)';
end


s_mean = R.syn_stats{3}.s_mean;
s_cov = R.syn_stats{3}.s_cov;
s_var = s_cov(logical(eye(length(s_cov(1,:)))));
V_mean = R.neuron_stats.V_time_mean{1} - V_rev;
V_cov = R.neuron_stats.V_time_cov{1};
V_var = s_cov(logical(eye(length(V_cov(1,:)))));


a_mean_s_cov = zeros(N, 1);
for jj = 1:N
    %%%%
    c_tmp = a_mean_mat(jj,:)';
    c_tmp = c_tmp*c_tmp';
    a_mean_s_cov(jj,1) = nansum(nansum(c_tmp.*s_cov));
    %%%%
    disp(jj/N)
end

sig2_comp_i = zeros(N,0);

pie_label = cell(1,15);
for i = 1:15 % skip 0
    pie_label{i} = '';
    vb = de2bi(i,4);
    if vb(1) == 1 % 1 => var, 0 => mean
        a = sqrt(a_var_mat);
        pie_label{i} = [pie_label{i}, 'a'];
    else
        a = a_mean_mat;
    end
    if vb(2) == 1
        J = sqrt(J_var);
        pie_label{i} = [pie_label{i}, 'J'];
    else
        J = J_mean;
    end
    if vb(3) == 1
        s = sqrt(s_var(:));
        pie_label{i} = [pie_label{i}, 's'];
    else
        s = s_mean(:);
    end
    if vb(4) == 1
        V = sqrt(V_var(:));
        pie_label{i} = [pie_label{i}, 'V'];
    else
        V = V_mean(:);
    end
    
    sig2_comp_i(:,end+1)  =  V.^2*J^2 .* ((a.^2)*(s.^2));  %#ok<AGROW,*SAGROW>
    
    
    if vb(1) == 0 && vb(3) == 1 && vb(2) == 0
        sig2_comp_i(:,end) =  V.^2*J^2 .* a_mean_s_cov;
    end
    
end

I_var_comp = mean(sig2_comp_i);
I_std_p = sqrt(sum(sig2_comp_i,2));
I_mean_p = - V_mean(:).*(a_mean_mat*J_mean*s_mean(:));

I_mean = R.syn_stats{3}.I_time_mean;
I_std = sqrt(R.syn_stats{3}.I_time_var);

% record
R.stru_var_decomp.I_var_comp_p = I_var_comp;
R.stru_var_decomp.I_std_p = I_std_p';
R.stru_var_decomp.I_mean_p = I_mean_p';
R.stru_var_decomp.I_std = I_std;
R.stru_var_decomp.I_mean = I_mean;
R.stru_var_decomp.pie_label = pie_label;

end




