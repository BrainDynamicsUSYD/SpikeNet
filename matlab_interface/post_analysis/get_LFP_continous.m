function get_LFP_continous(R)
stamp = R.stamp;
samp_file = [stamp(1:end-3) '0_neurosamp'];
fw = sqrt(R.N(1));
hw = (fw-1)/2;
try
    load(samp_file, 'ripple_power_grid');
    [N_s, ~, steps] = size(ripple_power_grid);
    N = N_s^2;
    if N ~= R.N(1)
        error('Not all of the excitatory neurons are sampled!')
    end
    fw = sqrt(N);
   
catch
end

load(samp_file, 'LFP_grid');
if ~exist('LFP_grid','var')
    disp('Creating LFP_grid...')
    load(samp_file,'I_AMPA','I_GABA');
    [N, steps] = size(I_AMPA);
    if N ~= R.N(1)
        error('Not all of the excitatory neurons are sampled!')
    end
    [Lattice, ~] = lattice_nD(2, hw);
    
    dist = lattice_nD_find_dist(Lattice, hw,  0, 0);
    LFP_range_sigma = R.ExplVar.LFP_range_sigma;
    gaus_tmp = 1/(LFP_range_sigma*sqrt(2*pi))*exp(-0.5*(dist/LFP_range_sigma).^2) .* double(dist <= LFP_range_sigma*2.5);
    m = reshape(gaus_tmp, fw, fw);
    LFP_grid = zeros(fw,fw,steps);
    for i = 1:steps
        x = reshape(abs(I_AMPA(:,i)) + abs(I_GABA(:,i)), fw, fw);
        y = convolve2(x, m, 'wrap');
        
        LFP_grid(:,:,i) = y;
        %i
    end
    clear I_AMPA I_GABA;
    save(samp_file, 'LFP_grid', '-append')
    clear LFP_grid;
end


load(samp_file, 'ripple_grid');
if ~exist('ripple_grid','var')
    disp('Creating ripple_grid...')
    load(samp_file,'LFP_grid');
    
    %%%% Ripple power
    fs = 1/(R.dt*1e-3); % sampling frequency (Hz)
    % Butterworth filter
    order = 4; % 4th order
    lowFreq = 100; % ripple band (default values for this function are 150-250 Hz)
    hiFreq = 250;
    Wn = [ lowFreq  hiFreq]/(fs/2);
    [b,a] = butter(order/2,Wn,'bandpass'); %The resulting bandpass and bandstop designs are of order 2n.
    ripple_grid = zeros(size(LFP_grid));
    for i = 1:fw
        for j= 1:fw
            ripple_tmp = filter(b,a,LFP_grid(i,j,:));
            ripple_tmp = ripple_tmp(:);
            ripple_grid(i,j,:) = ripple_tmp;
        end
        %w*(i-1)+j
    end
    save(samp_file, 'ripple_grid', '-append')
end


load(samp_file, 'ripple_power_grid');
if ~exist('ripple_power_grid','var')
    disp('Creating ripple_power_grid...')
    load(samp_file,'ripple_grid');
    
    gaus_width = 5; %ms
    [ Kernel ] = spike_train_kernel_YG( gaus_width, R.dt, 'gaussian_unit' );
    ripple_power_grid = zeros(size(ripple_grid));
    for i = 1:fw
        for j= 1:fw
            ripple_tmp = ripple_grid(i,j,:);
            hil_tmp = abs(hilbert( ripple_tmp));
            ripple_power_grid(i,j,:) = conv(hil_tmp(:)', Kernel,'same');
        end
    end
    save(samp_file, 'ripple_power_grid', '-append')
end

load(samp_file, 'ripple_fit');
if ~exist('ripple_fit','var')
    
    % get peak and ripple fit
    disp('Fitting ripple power...')
    [A_row, A_col, steps] = size(ripple_power_grid);
    [x_grid, y_grid] = meshgrid(1:A_row, 1:A_col);
    circ_gauss = fittype( @(h, sigma, x_c, y_c,  x, y) h*exp(-((x-x_c).^2 + (y-y_c).^2)/(2*sigma^2) ),...
        'independent', {'x', 'y'},...
        'dependent', 'z');
    peak = [];
    ripple_fit = cell(1,steps);
    ripple_fit_goodness = cell(1,steps);
    for i = 1:steps
        A =  ripple_power_grid(:,:,i);
        [peak_mag,I] = max(A(:));
        
        [I_row, I_col] = ind2sub([A_row, A_col],I);
        
        A_c = circshift(A, [round(A_row/2)-I_row  round(A_col/2)-I_col]);
        [A_fit, G_fit] = fit([x_grid(:), y_grid(:)], A_c(:), circ_gauss, ...
            'StartPoint', [peak_mag 5 A_row/2 A_col/2 ], ...
            'Lower', [0 0 0 0 ], ...
            'Upper', [Inf 10*max(I_row, I_col) A_row A_col]);
        peak = [peak; I_row, I_col, peak_mag]; %#ok<AGROW>
        ripple_fit{i} = A_fit;
        ripple_fit_goodness{i} = G_fit;
        % plot(A_fit, [x_grid(:), y_grid(:)], A_c(:))
        %             i/steps
        %             A_fit.sigma
    end
    
    save(samp_file,  'peak', '-append')
    save(samp_file,  'ripple_fit', 'ripple_fit_goodness', '-append')
end


load(samp_file, 'LFP_power_grid');
if ~exist('LFP_power_grid','var')
    disp('Creating LFP_power_grid...')
    
    
    load(samp_file,'LFP_grid');
    fs = 1/(R.dt*1e-3); % sampling frequency (Hz)
    % Butterworth filter
    order = 4; % 4th order
    lowFreq = 1; % broad-band
    hiFreq = 1000;
    Wn = [ lowFreq  hiFreq]/(fs/2);
    [b,a] = butter(order/2,Wn,'bandpass'); %The resulting bandpass and bandstop designs are of order 2n.
    gaus_width = 5; %ms
    [ Kernel ] = spike_train_kernel_YG( gaus_width, R.dt, 'gaussian_unit' );
    LFP_power_grid = zeros(size(LFP_grid));
    for i = 1:fw
        for j= 1:fw
            LFP_tmp = filter(b,a,LFP_grid(i,j,:));
            LFP_tmp = LFP_tmp(:);
            hil_tmp = abs(hilbert( LFP_tmp));
            LFP_power_grid(i,j,:) = conv(hil_tmp, Kernel,'same');
        end
        %w*(i-1)+j
    end
    save(samp_file, 'LFP_power_grid', '-append')
    clear LFP_power_grid LFP_grid;
end

load(samp_file, 'peak_LFP');
if ~exist('peak_LFP','var')
    disp('Creating RPI and peak_LFP...')
    % get ripper-power-index
    load(samp_file, 'ripple_power_grid','LFP_power_grid');
    peak_LFP = [];
    for i = 1:length(LFP_power_grid(1,1,:))
        RPI(1,i) = sum(sum(ripple_power_grid(:,:,i)))/sum(sum(LFP_power_grid(:,:,i))); %#ok<NASGU,AGROW>
        % LFP peak
        A = LFP_power_grid(:,:,i);
        [peak_mag,I] = max(A(:));
        [I_row, I_col] = ind2sub([fw, fw],I);
        peak_LFP = [peak_LFP; I_row, I_col, peak_mag]; %#ok<AGROW>
    end
    save(samp_file, 'RPI', 'peak_LFP', '-append')
    clear 'ripple_power_grid','LFP_power_grid';
end

load(samp_file, 'ripple_power_tot');
if ~exist('ripple_power_tot','var')
    disp('Creating ripple_power_tot and LFP_power_tot...')
    load(samp_file, 'ripple_power_grid','LFP_power_grid');
    for i = 1:length(LFP_power_grid(1,1,:))
        ripple_power_tot(1,i) = mean(mean(ripple_power_grid(:,:,i))); %#ok<NASGU,AGROW>
        LFP_power_tot(1,i) = mean(mean(LFP_power_grid(:,:,i))); %#ok<NASGU,AGROW>
    end
    save(samp_file, 'ripple_power_tot', 'LFP_power_tot', '-append')
    clear 'ripple_power_grid','LFP_power_grid';
end

end

