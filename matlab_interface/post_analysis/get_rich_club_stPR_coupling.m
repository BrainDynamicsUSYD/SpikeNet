cd ../SWR_reference_100sec_10runs/

rc_rate_z = cell(1,0);
rc_stPR_z = cell(1,0);
k = 300:25:600;
for kk = k
    rc_rate_z{end+1} = []; %#ok<*SAGROW>
    rc_stPR_z{end+1} = [];
    
    for i = 1:10
        loop = sprintf('%04g',i);
        
        n = dir([loop '*RYG.mat']);
        R = load(n.name);
        
        s = strsplit(R.stamp,'_');
        FID = [s{1},'_in.h5'];
        [syn_cell] = read_syn_h5(FID, 1, 1, 1);
        
        n_sample = 2000;
        
        L = lattice_nD(2, 31);
        I_sample = find(  L(:,1).^2+L(:,2).^2 <= n_sample/pi );
        
        I = syn_cell{1}.I;
        J = syn_cell{1}.J;
        K = syn_cell{1}.K;
        clear syn_cell;
        W = full(sparse(I,J,K, 63^2, 63^2));
        CIJ = W(I_sample, I_sample);
        clear I J K W;
        
        % definition of "degree" as used for RC coefficients
        % degree is taken to be the sum of incoming and outgoing connectons
        [~,~,degree] = degrees_dir(CIJ);
        
        SmallNodes=find(degree<=kk);       %get 'small nodes' with degree <=k
        %
        I_rc = I_sample;
        I_rc(SmallNodes) = [];
        rc_rate_z{end} = [rc_rate_z{end}(:); (R.Analysis.rate{1}(I_rc) - mean(R.Analysis.rate{1})) /  std(R.Analysis.rate{1})]; %#ok<*NOPTS>
        
        for ii = 1:R.stPR.n_trials
            ind_trial = (ii-1)*66+(1:66);
            [I_check_stPR, ~, I_stPR] = intersect(I_rc, R.stPR.sample_ind(ind_trial));
            c_tmp = R.stPR.c_norm(ind_trial);
            rc_stPR_z{end} =[ rc_stPR_z{end}(:); (c_tmp(I_stPR) - mean(c_tmp)) /  std(c_tmp)];
        end
        
    end
end

save('rich_club_stPR_analysis.mat','rc_rate_z','rc_stPR_z','k')