function [ results ] = Kreuz_SPIKE_dissimilarity( spikes, t, plot )
%UNTITLED2 Summary of this function goes here

% Backgroud knowledge:
%
% Time-scale dependent measures
% The most widely used are (cf A comparison of binless spike train measures, 2010)
%        a) Victor–Purpura distance (Victor and Purpura, 1996, 1997), cost of spike insertion, spike deletion and shifting
%        b) van Rossum distance (van Rossum, 2001), convolution with exponential kernel
%        c) Schreiber et al. similarity measure (Schreiber et al., 2003), convolution with Gaussian filter
%
% Here SPIKE-distance is time-scale independent measures, see <T. Kreuz, et al., 2013>
%
%   spikes: (number of neurons x spike timing), fill with zeros if not the same length, units in (ms)
%        t: time vector, units in (ms)
%
%   Tested with N=50 neurons, M=300 spikes each, T=40,000 time steps
%   Finishd in 10 seconds.
%   No more than 50 neurons should be used!
%   
%   Ref: http://wwwold.fi.isc.cnr.it/users/thomas.kreuz/Source-Code/Monitoring.html

if nargin == 2
    plot = 0;
end

max_sample_size = 50;
if (length(spikes(:,1))) > max_sample_size
    disp('Too many neurons!! No more than 50!');
    results = [];
else
    
    % Set up d_para and read inputs
    d_para_default = struct('tmin',[],'tmax',[],'dts',0,...
        'all_train_group_names',[],'all_train_group_sizes',[],'select_train_mode',1,'select_train_groups',[],'select_trains',[],'separators',[],'separators2',[],...
        'select_averages',[],'trigger_averages',[],'markers',[],'markers2',[],'interval_separators',[],'interval_strings',[],'comment_string','');
    d_para = d_para_default;
    d_para.tmin = t(1);
    d_para.tmax = t(end);
    d_para.dts = t(2)-t(1);
    
    % Set up f_para
    % ------------------------------------------------------------------------------------------------------------------------------
    % f_para.subplot_posi:
    %
    % 1           2         3     4        5     6       7     8         9   10       % Stimulus, Spike trains, Classic, Rate, ISI
    % Stimulus    Spikes    PSTH  GPSTH    ISI  <ISI>    Rate  <Rate>    Ia  I^A
    %
    % 11  12     13    14       15    16          17    18      19     20        21     22
    % Sa  S^a    Sa_r  S_r^a    Sa_f  S_f^a       Sta   St^a    Sta_r  St_r^a    Sta_f  St_f^a
    f_para_default = struct('imagespath',['.',filesep],'moviespath',['.',filesep],...    % Default values
        'matpath',['.',filesep],'filename','SPIKE-distance','title_string','',...
        'saving',0,'print_mode',0,'publication',0,'comment_string','','num_fig',1,'pos_fig',[130 120 1320 880],'font_size',14,'multi_figure',1,...
        'timeunit_string','[ms]','xfact',1,'ma_mode',0,'spike_mao',20,'time_mao',10,'dtm',1,...
        'mov_step',0,'mov_frames_per_second',1,'mov_num_average_frames',1,'mov_frames',[],...
        'plot_mode',1,'norm_mode',1,'color_norm_mode',1,'block_matrices',0,'dendrograms',0,'dendro_color_mode',0,'subplot_size',[],...
        'subplot_posi',[0 1  0 0  0 0 0 0  0 0    0 0 0 0 0 0     0 2 0 0 0 0]); % show spike raster plot and dissimilarity measures
    f_para_default.subplot_size = ones(1,length(f_para_default.subplot_posi(f_para_default.subplot_posi>0)));
    f_para = f_para_default;
    f_para.dtm = d_para.dts;
    if plot == 0
        f_para.subplot_posi = [0 0  0 0  0 0 0 0  0 0    0 0 0 0 0 0     0 1 0 0 0 0]; % to save the plotting time
        f_para.subplot_size = ones(1,length(f_para_default.subplot_posi(f_para_default.subplot_posi>0)));
    end
    
    % Start calculation
    results = f_distances_MEX(spikes,d_para,f_para);

    % Close figure
    if plot == 0
        close(gcf());
    end
    
end


end

