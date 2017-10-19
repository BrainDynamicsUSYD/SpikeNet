function main_test_seed_HDF5(varargin)

% clc;clear;close all;
% addpath(genpath(cd));
% cd ~/tmp_data
    
loop_num = 2;

% seed the matlab rand function! The seed is global.
[FID ] = new_ygin_files_and_randseedHDF5(loop_num);


k = 1*2.4e-3; % miuSiemens


% Basic parameters
dt = 0.1;
step_tot = 1000;
N = [10; 12; ];
writeBasicParaHDF5(FID, dt, step_tot, N)
Num_pop = length(N);
discard_transient = 0; % ms

% write pop para
writePopParaHDF5(FID, 1,  'tau_ref', 3.1);
writePopParaHDF5(FID, 2,  'tau_ref', 3.2);
% write synapse para
writeSynParaHDF5(FID, 'tau_decay_AMPA', 2.0, 'Dt_trans_AMPA', 0.5);


g_ext_strength = 0.5;
writeExtConductanceSettingsHDF5(FID, 2, g_ext_strength*ones(1,N(2)), ones(1,N(2)));



%%%%%%%%%%%%%%%%%%% Chemical Connections %%%%%%%%%%%%%%%%%%%%%%%
% type(1:AMAP, 2:GABAa, 3:NMDA)

pop_type = [1 2];
for i_pre = 1:Num_pop
    for j_post = 1:Num_pop
        [I, J, ~] = find(MyRandomGraphGenerator('E_R_pre_post','N_pre',N(i_pre),'N_post',N(j_post),'p',rand));
        K = ones(size(I))*k;
        D = ones(size(I))*1;
        writeChemicalConnectionHDF5(FID,  pop_type(i_pre), i_pre,  j_post,  I,J,K,D); % (FID, type, i_pre, j_post, I, J, K, D)
    end
end


end


