% Distances_Main_Demo © Thomas Kreuz, Nebojsa Bozanic;  Version 1.3, 29th of April 2013
%
% Calculates ISI- and SPIKE-Distance from multiple spike trains "spikes"    
% (matrix with different spike trains as rows)
%
% Information can be found under "http://www.fi.isc.cnr.it/users/thomas.kreuz/Source-Code/Monitoring.html" and/or in
%
% Kreuz T, Chicharro D, Houghton C, Andrzejak RG, Mormann F: Monitoring spike train synchrony. Submitted (2012).
% Kreuz T: Measures of spike train synchrony. Scholarpedia, 6(10):11934 (2011).
% Kreuz T, Chicharro D, Greschner M, Andrzejak RG: Time-resolved and time-scale adaptive measures of spike train synchrony. J Neurosci Methods 195, 92 (2011).
% Kreuz T, Chicharro D, Andrzejak RG, Haas JS, Abarbanel HDI: Measuring multiple spike train synchrony. J Neurosci Methods 183, 287 (2009).
% Kreuz T, Haas JS, Morelli A, Abarbanel HDI, Politi A: Measuring spike train synchrony. J Neurosci Methods 165, 151 (2007)
%
% For questions and comments please contact me at "thomaskreuz (at) cnr.it".

clear all; close all
disp(' '); disp(' ')

fixed=1;            % use always the same random numbers (0-no,1-yes)
if fixed==0 rand('twister',sum(100*clock)); randn('state',sum(100*clock)); else rand('twister',5489); randn('state',5489); end %#ok<SEPEX>

% structure 'd_para': parameters that describe the data
%
% tmin: Beginning of recording
% tmax: End of recording
% dts: Sampling interval, precision of spike times
% all_train_group_names: Names of spike train groups
% all_train_group_sizes: Sizes of respective spike train groups
% select_train_mode: Selection of spike trains (1-all,2-selected groups,3-selected trains)
% select_train_groups: Selected spike train groups (if 'select_train_mode==2')
% select_trains: Selected spike trains (if 'select_train_mode==3')
% select_averages: One or more continuous intervals for selective temporal averaging
% trigger_averages: Non-continuous time-instants for triggered temporal averaging, external (e.g. certain stimulus properties) but also internal (e.g. certain event times)
% markers: Relevant time instants
% markers2: Even more relevant time instants
% separators: Relevant seperations between groups of spike trains
% separators2: Even more relevant seperations between groups of spike trains
% interval_separators: Edges of subsections
% interval_strings: Captions for subsections
% comment_string: Additional comment on the example, will be used in file and figure names

d_para_default=struct('tmin',[],'tmax',[],'dts',1,...
    'all_train_group_names',[],'all_train_group_sizes',[],'select_train_mode',1,'select_train_groups',[],'select_trains',[],'separators',[],'separators2',[],...
    'select_averages',[],'trigger_averages',[],'markers',[],'markers2',[],'interval_separators',[],'interval_strings',[],'comment_string','');
d_para=d_para_default;


% structure 'f_para': parameters that determine the appearance of the figures (and the movie)
%
% imagespath: Path where images (postscript) will be stored
% moviespath: Path where movies (avi) will be stored
% matpath: Path where Matlab files (mat) will be stored
% filename: Name under which images, movies and Matlab files will be stored
% title_string: Appears in the figure titles
% saving: Saving to Mat-file? (0-no,1-yes)
% print_mode: Saving to postscript file? (0-no,1-yes)
% publication: Omits otherwise helpful information, such as additional comments (0-no,1-yes)
% comment_string: Additional comment on the example, will be used in file and figure names
% num_fig: Number of figure
% pos_fig: Position of figure
% font_size: Font size of labels (and headlines)
% multi_figure: Open many figures (0-no,1-yes)
% timeunit_string: Time unit, used in labels
% xfact: Conversion of time unit
% ma_mode: Moving average mode: (0-no,1-only,2-both)
% spike_mao: Order of moving average (pooled ISI)
% time_mao: Order of moving average (time)
% dtm: Sampling of measure profile, downsampling possible to facilitate memory management
% mov_step: Step size for movie frames
% mov_frames_per_second: Well, frames per second
% mov_num_average_frames: Number of frames the averages are shown at the end of the movie (if this is too small they are hardly visible)
% mov_frames: Selection of individual frames (e.g. in the movie, replaces use of mov_step as step size if non-empty)
% plot_mode: +1:vs time,+2:different measures and cuts,+4:different measures,+8:different cuts,+16:different cuts-Movie (binary addition allows all combinations)
% norm_mode: normalization of averaged bivariate measure profiles (1-Absolute maximum value 'one',2-Overall maximum value,3-Individual maximum value)
% color_norm_mode: normalization of pairwise color matrices (1-Absolute maximum value,2-Overall maximum value,3-Each one individually)
% block_matrices: Allows tracing the overall synchronization among groups of spike trains (0-no,1-yes)
% dendrograms: Cluster trees from distance matrices (0-no,1-yes)
% dendro_color_mode: coloring of dendrograms (0-no,1-groups,2-trains)
% subplot_posi: Vector with order of subplots
% subplot_size: Vector with relative size of subplots

f_para_default=struct('imagespath',['.',filesep],'moviespath',['.',filesep],...    % Default values
    'matpath',['.',filesep],'filename','Demo_','title_string','',...
    'saving',0,'print_mode',1,'publication',1,'comment_string','','num_fig',1,'pos_fig',[130 120 1320 880],'font_size',14,'multi_figure',1,...
    'timeunit_string','[ms]','xfact',1,'ma_mode',0,'spike_mao',20,'time_mao',10,'dtm',1,...
    'mov_step',0,'mov_frames_per_second',1,'mov_num_average_frames',1,'mov_frames',[],...
    'plot_mode',1,'norm_mode',1,'color_norm_mode',1,'block_matrices',0,'dendrograms',0,'dendro_color_mode',0,'subplot_size',[],...
    'subplot_posi',[0 1  0 0  0 0 0 0  0 0    0 0 0 0 0 0     0 2 0 0 0 0]);
f_para_default.subplot_size=ones(1,length(f_para_default.subplot_posi(f_para_default.subplot_posi>0)));



% Choose order and relative size of the subplots
% Use "subplot_posi" to select the order of the subplots. Use 0 if no subplot is needed.
% Make sure that the range from 1 to the desired number of subplots is covered
% Use "subplot_size" to select the size of the selected subplots.
% Its order refers to the subplots selected in subplot_posi, not to the reference-order
% (i.e., its length should be equal to the number of subplot selected in 'subplot_posi')
%
% 1           2         3     4        5     6       7     8         9   10       % Stimulus, Spike trains, Classic, Rate, ISI
% Stimulus    Spikes    PSTH  GPSTH    ISI  <ISI>    Rate  <Rate>    Ia  I^A
%
% 11  12     13    14       15    16          17    18      19     20        21     22
% Sa  S^a    Sa_r  S_r^a    Sa_f  S_f^a       Sta   St^a    Sta_r  St_r^a    Sta_f  St_f^a
%


% #######################################  Loop over examples (from simple to complicated) #######################################

% Choose the examples that you would like to see

% 21:Fig2a,22:Fig2b,31:Fig3a,32:Fig3b,41:Fig4,51:Fig5,61:Fig6,71:Fig7,81:Fig8,91:Fig9,101:Movie,121:Poisson-Scalefree,131:Poisson-Expectation

examples=[21 22 31 32 41 51 61 71 81 91 101];   
%examples=[101];   

num_examples=length(examples);
for esc=1:num_examples
    example=examples(esc)
    dataset=example;
    f_para=f_para_default;                              % Reset
    
    switch example
        case 21
            [spikes,d_para]=f_get_data_demo(dataset,d_para_default);
            f_para.subplot_posi=[0 1  0 0  0 0 0 0  0 0    0 2 0 0 0 0     0 0 0 0 0 0]; % Vector with order of subplots
            f_para.plot_mode=0;                 % +1:vs time,+2-different measures and cuts,+4-different measures,+8-different cuts,+16:different cuts-Movie
            f_para.comment_string='Fig2a';    % Additional comment on the example, will be used in file and figure names
        case 22
            [spikes,d_para]=f_get_data_demo(dataset,d_para_default);
            f_para.subplot_posi=[0 1  0 0  0 0 0 0  0 0    0 2 0 0 0 0     0 3 0 0 0 0]; % Vector with order of subplots
            f_para.plot_mode=1;                 % +1:vs time,+2-different measures and cuts,+4-different measures,+8-different cuts,+16:different cuts-Movie
            f_para.comment_string='Fig2b';    % Additional comment on the example, will be used in file and figure names
        case 31
            dataset=21;
            [spikes,d_para]=f_get_data_demo(dataset,d_para_default);
            f_para.subplot_posi=[0 1  0 0  0 0 0 0  0 0    0 0 0 0 0 0     0 0 0 2 0 0]; % Vector with order of subplots
            f_para.plot_mode=1;                      % +1:vs time,+2-different measures and cuts,+4-different measures,+8-different cuts,+16:different cuts-Movie
            f_para.ma_mode=2;                        % Moving average mode: (0-no,1-only,2-both)
            f_para.time_mao=120;                     % Order of moving average (time)
            f_para.comment_string='Fig3a_Realtime';% Additional comment on the example, will be used in file and figure names
        case 32
            dataset=22;
            [spikes,d_para]=f_get_data_demo(dataset,d_para_default);
            f_para.subplot_posi=[0 1  0 0  0 0 0 0  0 0    0 0 0 0 0 0     0 0 0 2 0 0]; % Vector with order of subplots
            f_para.plot_mode=1;                      % +1:vs time,+2-different measures and cuts,+4-different measures,+8-different cuts,+16:different cuts-Movie
            f_para.ma_mode=2;                        % Moving average mode: (0-no,1-only,2-both)
            f_para.time_mao=190;                     % Order of moving average (time)
            f_para.comment_string='Fig3b_Realtime';% Additional comment on the example, will be used in file and figure names
        case 41
            [spikes,d_para]=f_get_data_demo(dataset,d_para_default);
            f_para.subplot_posi=[0 1  0 0  0 0 0 0  0 0    0 0 0 0 0 0     0 0 0 2 0 0]; % Vector with order of subplots
            f_para.plot_mode=1;                      % +1:vs time,+2-different measures and cuts,+4-different measures,+8-different cuts,+16:different cuts-Movie
            f_para.comment_string='Fig4';          % Additional comment on the example, will be used in file and figure names
        case 51
            [spikes,d_para]=f_get_data_demo(dataset,d_para_default);
            d_para.interval_separators=500:500:d_para.tmax-500; % Edges of subsections
            f_para.subplot_posi=[0 1  0 0  0 0 0 0  0 0    0 0 0 0 0 0     0 2 0 3 0 0]; % Vector with order of subplots
            f_para.plot_mode=2;                                 % +1:vs time,+2-different measures and cuts,+4-different measures,+8-different cuts,+16:different cuts-Movie
            f_para.mov_frames=[250 750 1190 1510]; % Selection of individual frames
            f_para.xfact=1000;                  % Conversion of time unit
            f_para.timeunit_string='[s]';       % Time unit, used in labels
            f_para.comment_string='Fig5';     % Additional comment on the example, will be used in file and figure names
        case 61
            [spikes,d_para]=f_get_data_demo(dataset,d_para_default);
            d_para.interval_separators=500:500:d_para.tmax-500; % Edges of subsections
            f_para.subplot_posi=[0 1  0 0  0 0 0 0  0 0    0 0 0 0 0 0     0 2 0 3 0 0]; % Vector with order of subplots
            f_para.plot_mode=2;                    % +1:vs time,+2-different measures and cuts,+4-different measures,+8-different cuts,+16:different cuts-Movie
            f_para.mov_frames=[250 750 1250 1750]; % Selection of individual frames (e.g. in the movie, replaces use of mov_step as step size if non-empty)
            f_para.xfact=1000;                  % Conversion of time unit
            f_para.timeunit_string='[s]';       % Time unit, used in labels
            f_para.comment_string='Fig6';     % Additional comment on the example, will be used in file and figure names
        case 71
            [spikes,d_para]=f_get_data_demo(dataset,d_para_default);
            d_para.interval_separators=500:500:d_para.tmax-500; % Edges of subsections
            d_para.select_averages={[0 500];[1000 2000];[2500 3000 3500 4000]-500;[0 4000]}; % One or more continuous intervals for selective temporal averaging
            f_para.subplot_posi=[0 1  0 0  0 0 0 0  0 0    0 0 0 0 0 0     0 2 0 3 0 0]; % Vector with order of subplots
            f_para.plot_mode=2;                 % +1:vs time,+2-different measures and cuts,+4-different measures,+8-different cuts,+16:different cuts-Movie
            f_para.xfact=1000;                  % Conversion of time unit
            f_para.timeunit_string='[s]';       % Time unit, used in labels
            f_para.comment_string='Fig7';     % Additional comment on the example, will be used in file and figure names
        case 81
            [spikes,d_para]=f_get_data_demo(dataset,d_para_default);
            trig_trac1=1;
            num_spikes=find(spikes(trig_trac1,:),1,'last');
            d_para.select_averages={[d_para.tmin d_para.tmax]}; % One or more continuous intervals for selective temporal averaging
            d_para.trigger_averages{1}=round(sort(spikes(trig_trac1,1:num_spikes))/d_para.dts)*d_para.dts;       % Triggered averaging over all time instants when a certain neuron fires
            f_para.subplot_posi=[0 1  0 0  0 0 0 0  0 0    0 0 0 0 0 0     0 2 0 0 0 0]; % Vector with order of subplots
            f_para.plot_mode=8;                % +1:vs time,+2-different measures and cuts,+4-different measures,+8-different cuts,+16:different cuts-Movie
            f_para.dendrograms=1;              % Cluster trees from distance matrices (0-no,1-yes)
            f_para.comment_string='Fig8';    % Additional comment on the example, will be used in file and figure names
        case 91
            dataset=71;
            [spikes,d_para]=f_get_data_demo(dataset,d_para_default);
            d_para.markers=[500:500:d_para.tmax-500];           % Relevant time instants
            d_para.all_train_group_names={'G1';'G2';'G3';'G4'}; % Names of spike train groups
            d_para.all_train_group_sizes=[10 10 10 10];         % Sizes of respective spike train groups
            d_para.select_averages={500+[0 500 1000 1500]};     % Selected average over different intervals
            f_para.subplot_posi=[0 1  0 0  0 0 0 0  0 0    0 0 0 0 0 0     0 2 0 3 0 0]; % Vector with order of subplots
            f_para.plot_mode=8;                                 % +1:vs time,+2-different measures and cuts,+4-different measures,+8-different cuts,+16:different cuts-Movie
            f_para.block_matrices=1;                            % Allows tracing the overall synchronization among groups of spike trains (0-no,1-yes)
            f_para.dendrograms=1;                               % Cluster trees from distance matrices (0-no,1-yes)
            f_para.dendro_color_mode=1;                         % Coloring of dendrograms (0-no,1-groups,2-trains)
            f_para.comment_string='Fig9';                     % Additional comment on the example, will be used in file and figure names
        case 101
            dataset=71;
            [spikes,d_para]=f_get_data_demo(dataset,d_para_default);
            num_all_trains=size(spikes,1);
            d_para.markers=[500:500:d_para.tmax-500];           % Relevant time instants
            d_para.all_train_group_names={'G1';'G2';'G3';'G4'}; % Names of spike train groups
            d_para.all_train_group_sizes=num_all_trains/4*ones(1,4);         % Sizes of respective spike train groups
            d_para.interval_strings={'2 Cluster - AABB';'2 Cluster - ABBA';'2 Cluster - ABAB';'2 Cluster - Random association';...
                '3 Cluster - ABBC';'4 Cluster - ABCD';'8 Cluster - ABCDEFGH';'Random Spiking'}; % Captions for subsections
            f_para.subplot_posi=[0 1  0 0  0 0 0 0  0 0    0 2 0 0 0 0     0 3 0 0 0 0]; % Vector with order of subplots
            f_para.plot_mode=16;                                % +1:vs time,+2-different measures and cuts,+4-different measures,+8-different cuts,+16:different cuts-Movie
            f_para.publication=0;                               % Omits otherwise helpful information, such as additional comments (0-no,1-yes)
            f_para.block_matrices=1;                            % Allows tracing the overall synchronization among groups of spike trains (0-no,1-yes)
            f_para.dendrograms=1;                               % Cluster trees from distance matrices (0-no,1-yes)
            f_para.dendro_color_mode=1;                         % Coloring of dendrograms (0-no,1-groups,2-trains)
            full_short=2;  % 1-full,2-short
            if full_short==1
                d_para.select_averages={[d_para.tmin d_para.tmax];[0 500];[500 1000];[1000 1500];[1500 2000];[2000 2500];[2500 3000];[3000 3500];[3500 d_para.tmax];...
                    [0 1000];[500 1500];[1000 2000];[1500 2500];[2000 3000];[2500 3500];[3000 4000];[0 500 1000 1500];500+[0 500 1000 1500];1000+[0 500 1000 1500];...
                    1500+[0 500 1000 1500];2000+[0 500 1000 1500];2500+[0 500 1000 1500];[0 500 1000 1500 2000 2500 3000 3500];...
                    [500 1000 1500 2000 2500 3000 3500 4000];[d_para.tmin d_para.tmax]};     % Selected average over different intervals
                for trac=1:size(spikes,1)
                    d_para.trigger_averages{trac}=round(spikes(trac,1:size(spikes,2))/d_para.dts)*d_para.dts;       % Triggered averaging over all time instants when a certain neuron fires
                end
                f_para.mov_step=1;                       % Step size for movie frames
                f_para.mov_frames_per_second=30;         % Well, frames per second
                f_para.mov_num_average_frames=30;        % Number of frames the averages are shown at the end of the movie (if this is too small they are hardly visible)
                f_para.comment_string='Movie_Full';    % Additional comment on the example, will be used in file and figure names
            else
                d_para.select_averages={[d_para.tmin d_para.tmax];[0 500];[0 1000];[0 500 1000 1500];[0 500 1000 1500 2000 2500 3000 3500]};     % Selected average over different intervals
                tracs=1+[0 1 2 3]*num_all_trains/4;
                for trac=1:length(tracs)
                    d_para.trigger_averages{trac}=round(spikes(tracs(trac),1:size(spikes,2))/d_para.dts)*d_para.dts;       % Triggered averaging over all time instants when a certain neuron fires
                end
                f_para.mov_step=100;                   % Step size for movie frames
                f_para.mov_frames_per_second=1;         % Well, frames per second
                f_para.comment_string='Movie-Short';  % Additional comment on the example, will be used in file and figure names
            end
        case 121
            [spikes,d_para]=f_get_data_demo(dataset,d_para_default);
            d_para.select_averages={[d_para.tmin d_para.tmax]}; % Continuous interval for selective temporal averaging
            f_para.subplot_posi=[0 1  0 0  0 0 0 0  0 0    0 0 0 0 0 0     0 2 0 3 0 0]; % Vector with order of subplots
            f_para.plot_mode=8;                % +1:vs time,+2-different measures and cuts,+4-different measures,+8-different cuts,+16:different cuts-Movie
            f_para.color_norm_mode=3;          % normalization of pairwise color matrices (1-Absolute maximum value,2-Overall maximum value,3-Each one individually)
        case 131  % (roughly 1.000.000 spikes, takes quite long)
            [spikes,d_para]=f_get_data_demo(dataset,d_para_default);
            f_para.subplot_posi=[0 1  0 0  0 0 0 0  0 0    0 2 0 3 0 0     0 0 0 0 0 0]; % Vector with order of subplots
            f_para.plot_mode=1;                 % +1:vs time,+2-different measures and cuts,+4-different measures,+8-different cuts,+16:different cuts-Movie
            f_para.publication=0;               % Omits otherwise helpful information, such as additional comments (0-no,1-yes)
            %f_para.ma_mode=2;                   % Moving average mode: (0-no,1-only,2-both)
            %f_para.spike_mao=100;               % Order of moving average (isi)
    end
    f_para.dtm=d_para.dts;
    f_para.num_fig=esc*10;

    disp(' '); disp(' ')
    %f_distances_MEX
    results=f_distances_MEX(spikes,d_para,f_para)
    %results=f_distances(spikes,d_para,f_para)
    disp(' '); disp(' ')
end



