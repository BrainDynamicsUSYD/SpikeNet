function writeJHLearnDirectionHDF5(FID, pop_ind_pre, pop_ind_post, direction)
% add inhibitory STDP to connection
%     FID: file id for writing data
% pop_ind: 

pop_ind_pre = pop_ind_pre - 1; % from matlab to c++ index
pop_ind_post = pop_ind_post - 1;

n_syns = h5read(FID,'/config/syns/n_syns');
n_match = NaN; 
for n = (1:n_syns)-1
    try
        type = hdf5read(FID,['/config/syns/syn',num2str(n),'/INIT006/type']);
        i_pre = hdf5read(FID,['/config/syns/syn',num2str(n),'/INIT006/i_pre']);
        j_post = hdf5read(FID,['/config/syns/syn',num2str(n),'/INIT006/j_post']);
        
        if pop_ind_pre == i_pre && pop_ind_post == j_post 
            n_match = n;
        end
    catch ME
    end
end


if isnan(n_match)
    error('Cannot find syn with identical pop_ind_pre, pop_ind_post!')
else
    hdf5write(FID,['/config/syns/syn',num2str(n_match),'/INIT016/direction'],int32(direction),'WriteMode','append');

end



end

