function stru_var = get_structural_variance( target_in_deg, deg_tol, in_file, config_file)

N = 63^2;
fw = sqrt(N);
hw = (fw-1)/2;

a_mean = zeros(fw,fw);
a_var = zeros(fw, fw);
k_mean = zeros(fw,fw);
k_var = zeros(fw, fw);
k_steps = ones(fw, fw);
a_steps = ones(fw, fw);

i_pre = 1;
j_post = 1;

if nargin == 2
    in_file = '*in.h5';
    config_file = '*config_data.mat';
end

dir_strut = dir(in_file);
num_files = length(dir_strut);
in_files = cell(1,num_files);
for id_out = 1:num_files
    in_files{id_out} = dir_strut(id_out).name;
end

dir_strut = dir(config_file);
num_files = length(dir_strut);
conf_files = cell(1,num_files);
for id_out = 1:num_files
    conf_files{id_out} = dir_strut(id_out).name;
end

for f = 1:num_files
    
    fprintf('Loading %s and ', in_files{f});
    fprintf('%s...\n', conf_files{f});
    
    
    [syn_cell] = read_syn_h5(in_files{f}, i_pre, j_post, 1);
    load( conf_files{f},'in_degree');
    
    for i = find(abs(in_degree - target_in_deg) < deg_tol)
        J_tmp = syn_cell{1}.J( syn_cell{1}.I == i);
        K_tmp = syn_cell{1}.K( syn_cell{1}.I == i);
        [x_tmp, y_tmp] = ind2sub([fw,fw], J_tmp);
        Wi = sparse(x_tmp, y_tmp, K_tmp, fw, fw);
        [xc_tmp, yc_tmp] = ind2sub([fw,fw], i);
        Wi = circshift(Wi, [hw+1-xc_tmp, hw+1-yc_tmp]);
        
        Ai = double(Wi > 0);
        [a_mean, a_var, a_steps] = nan_Welford_online(Ai, a_mean, a_var, a_steps, false,'nan');
        
        [k_mean, k_var, k_steps] = nan_Welford_online(Wi, k_mean, k_var, k_steps, false, 'zero');
    end
    
end

a_var = a_var./(a_steps - 1);
k_var = k_var./(k_steps - 1);

stru_var.a_mean = a_mean;
stru_var.a_var = a_var;
stru_var.J_mean = k_mean;
stru_var.J_var = k_var;
stru_var.i_pre = i_pre;
stru_var.j_post = j_post;
stru_var.degree = target_in_deg;
stru_var.sample_size = a_steps;
end