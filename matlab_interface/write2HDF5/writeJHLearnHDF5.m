function writeJHLearnHDF5(FID, pop_ind_pre, pop_ind_post, infsteps,learn_rate,learn_rate_all,tau, inf_scale,noise, ntype_pre,ntype_post, dropout)
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
        
        if pop_ind_pre == i_pre && pop_ind_post == j_post 
            n_match = n;
        end
    catch ME
    end
end


if isnan(n_match)
    error('Cannot find syn with identical pop_ind_pre, pop_ind_post!')
else
    hdf5write(FID,['/config/syns/syn',num2str(n_match),'/INIT016/infsteps'],int32(infsteps),'WriteMode','append');
    hdf5write(FID,['/config/syns/syn',num2str(n_match),'/INIT016/learnrate'],learn_rate,'WriteMode','append');
    hdf5write(FID,['/config/syns/syn',num2str(n_match),'/INIT016/learnrateall'],learn_rate_all,'WriteMode','append');
    hdf5write(FID,['/config/syns/syn',num2str(n_match),'/INIT016/tau'],int32(tau),'WriteMode','append');
    hdf5write(FID,['/config/syns/syn',num2str(n_match),'/INIT016/infscale'],inf_scale,'WriteMode','append');
    hdf5write(FID,['/config/syns/syn',num2str(n_match),'/INIT016/noise'],noise,'WriteMode','append');
    hdf5write(FID,['/config/syns/syn',num2str(n_match),'/INIT016/ntype_pre'],ntype_pre,'WriteMode','append');
    hdf5write(FID,['/config/syns/syn',num2str(n_match),'/INIT016/ntype_post'],ntype_post,'WriteMode','append');
    hdf5write(FID,['/config/syns/syn',num2str(n_match),'/INIT016/dropout'],dropout,'WriteMode','append');

    
%     hdf5write(FID,['/config/syns/syn',num2str(n_match),'/INIT016/infsteps'],infsteps,'WriteMode','append');
%     hdf5write(FID,['/config/syns/syn',num2str(n_match),'/INIT016/learnrate_E'],learn_rate_E,'WriteMode','append');
%     hdf5write(FID,['/config/syns/syn',num2str(n_match),'/INIT016/learnrateall_E'],learn_rate_all_E,'WriteMode','append');
%     hdf5write(FID,['/config/syns/syn',num2str(n_match),'/INIT016/learnrate_I'],learn_rate_I,'WriteMode','append');
%     hdf5write(FID,['/config/syns/syn',num2str(n_match),'/INIT016/learnrateall_I'],learn_rate_all_I,'WriteMode','append');
%     hdf5write(FID,['/config/syns/syn',num2str(n_match),'/INIT016/tau'],int32(tau),'WriteMode','append');
%     hdf5write(FID,['/config/syns/syn',num2str(n_match),'/INIT016/infscaleI'],inf_scaleI,'WriteMode','append');
%     hdf5write(FID,['/config/syns/syn',num2str(n_match),'/INIT016/infscaleE'],inf_scaleE,'WriteMode','append');
%     hdf5write(FID,['/config/syns/syn',num2str(n_match),'/INIT016/noise'],noise,'WriteMode','append');
%     hdf5write(FID,['/config/syns/syn',num2str(n_match),'/INIT016/ntype_pre'],ntype_pre,'WriteMode','append');
%     hdf5write(FID,['/config/syns/syn',num2str(n_match),'/INIT016/ntype_post'],ntype_post,'WriteMode','append');

end



end

