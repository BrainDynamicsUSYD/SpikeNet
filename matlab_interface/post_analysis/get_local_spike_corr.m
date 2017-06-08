
function R = get_local_spike_corr(R)

sh0 = R.reduced.spike_hist{1};
dt = R.reduced.dt;

hw = 31;
L = lattice_nD(2, hw);
if (hw*2+1)^2 ~=R.N(1)
    error('Not a hw=31 grid!')
end
corrcoef_sample_num = 10^3; % number of sampling pairs
CC_kernel_width = 40; % ms, kernel length

% Define kernel
kernel_type = 'square_Hz';
CC_kernel = spike_train_kernel_YG(CC_kernel_width, dt, kernel_type);


CC_local_c = cell(1,0);
for r = 10:5:30

    I_sample = find(  L(:,1).^2 + L(:,2).^2 <= r^2 ); % rectangular method
    sh = sh0(I_sample,:);

    CC_local = zeros(1, corrcoef_sample_num);
    pairs = rand_unique_pairs( length(I_sample), corrcoef_sample_num );
    
    % Calculate pair-wise corrcoef
    for i = 1:corrcoef_sample_num
        CC_local(1,i) = CorrCoefYG(sh(pairs(1,i),:),sh(pairs(2,i),:), CC_kernel);
        % display progress
    end
    
    CC_local_c{end+1} = CC_local; %#ok<AGROW>
    
end

R.local_cc.cc = CC_local_c;
R.local_cc.r = 10:5:30;

end
