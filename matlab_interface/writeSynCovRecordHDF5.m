function writeSynCovRecordHDF5(FID, pop_ind_pre, pop_ind_post, syn_type,  time_start, time_end)

% for C/C++ index convetion
pop_ind_pre = pop_ind_pre-1;
pop_ind_post = pop_ind_post-1;
syn_type = syn_type-1;

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
    error('Cannot find syn with identical pop_ind_pre, pop_ind_post and syn_type!')
else
    time_start = time_start - 1;
    time_end = time_end - 1;
    % write
    hdf5write(FID,['/config/syns/syn',num2str(n_match),'/SAMP104/time_start'],time_start,'WriteMode','append');
    hdf5write(FID,['/config/syns/syn',num2str(n_match),'/SAMP104/time_end'],time_end,'WriteMode','append');
    
end

end