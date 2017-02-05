function [syn_cell] = read_syn_h5(FID, pop_ind_pre, pop_ind_post, syn_type)
% refer to the function writeChemicalConnectionHDF5


% find match
pop_ind_pre = pop_ind_pre - 1; % from matlab to c++ index
pop_ind_post = pop_ind_post - 1;
syn_type = syn_type - 1;
 
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

% read data
if isnan(n_match)
    error('Cannot find syn with identical pop_ind_pre, pop_ind_post and syn_type!')
else
    syn.I = try_h5read(FID,['/config/syns/syn',num2str(n_match),'/INIT006/I']) + 1;
    syn.J = try_h5read(FID,['/config/syns/syn',num2str(n_match),'/INIT006/J']) + 1;
    syn.K = try_h5read(FID,['/config/syns/syn',num2str(n_match),'/INIT006/K']);
    syn.D = try_h5read(FID,['/config/syns/syn',num2str(n_match),'/INIT006/D']);
end

syn_cell = {syn};

end