function writeInhSTDPHDF5(FID, pop_ind_pre, pop_ind_post, inh_STDP_step)
% add inhibitory STDP to connection
%     FID: file id for writing data
% pop_ind: 

pop_ind_pre = pop_ind_pre - 1; % from matlab to c++ index
pop_ind_post = pop_ind_post - 1;

syn_type = 1; % GABA

n_syns = h5read(FID,'/config/syns/n_syns');
n_match = NaN; 
for n = (1:n_syns)-1
    try
        type = hdf5read(FID,['/config/syns/syn',num2str(n),'/INIT006/type']);
        i_pre = hdf5read(FID,['/config/syns/syn',num2str(n),'/INIT006/i_pre']);
        j_post = hdf5read(FID,['/config/syns/syn',num2str(n),'/INIT006/j_post']);
        
        if pop_ind_pre == i_pre && pop_ind_post == j_post &&  syn_type == type
            n_match = n;
        end
    catch ME
    end
end


if isnan(n_match)
    error('Cannot find syn with identical pop_ind_pre, pop_ind_post and syn_type = GABA!')
else
    hdf5write(FID,['/config/syns/syn',num2str(n_match),'/INIT009/inh_STDP_on_step'], inh_STDP_step,'WriteMode','append');
end


end

