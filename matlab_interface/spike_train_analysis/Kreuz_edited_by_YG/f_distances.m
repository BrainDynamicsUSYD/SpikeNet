% f_distances © Thomas Kreuz, October 2012
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

function mean_values=f_distances(spikes,d_para,f_para)

max_memo_init=10000000;      % memory management, should be big enough to hold the basic matrices but small enough to not run out of memory


% ################################# structure 'd_para': parameters that describe the data ###################################

if isfield(d_para,'all_train_group_names') && isfield(d_para,'all_train_group_sizes') && ~isempty(d_para.all_train_group_sizes) && length(d_para.all_train_group_names)==length(d_para.all_train_group_sizes)

    num_all_groups=length(d_para.all_train_group_names);
    num_all_trains=size(spikes,1);

    cum_group=[0 cumsum(d_para.all_train_group_sizes)];
    group_vect=zeros(1,num_all_trains);
    for gc=1:num_all_groups
        group_vect(cum_group(gc)+(1:d_para.all_train_group_sizes(gc)))=gc;
    end

    if d_para.select_train_mode==1                                       % All
        select_trains=1:num_all_trains;
    elseif d_para.select_train_mode==2                                   % Selected groups
        select_trains = [];
        for gc=d_para.select_train_groups
            select_trains=[select_trains find(group_vect==gc)];
        end
    elseif d_para.select_train_mode==3                                   % Selected trains
        select_trains=d_para.select_trains;
    end
    %num_select_trains=length(select_trains);

    select_group_vect=group_vect(select_trains);
    select_groups=unique(select_group_vect);
    num_select_train_groups=length(select_groups);
    num_select_group_trains=zeros(1,num_select_train_groups);
    for scc=1:num_select_train_groups
        num_select_group_trains(scc)=length(find(select_group_vect==select_groups(scc)));
    end
    cum_num_select_group_trains=cumsum(num_select_group_trains);
    select_group_center=[0 cum_num_select_group_trains(1:num_select_train_groups-1)]+diff([0 cum_num_select_group_trains(1:num_select_train_groups)])/2;
    select_group_names=d_para.all_train_group_names(select_groups);
    spikes=spikes(select_trains,:);
else
    num_select_train_groups=1;
end

if isfield(d_para,'select_averages') && ~isempty(d_para.select_averages)
    num_select_averages=size(d_para.select_averages,1);
else
    num_select_averages=0;
end

if isfield(d_para,'trigger_averages') && ~isempty(d_para.trigger_averages)
    num_trigger_averages=size(d_para.trigger_averages,2);
else
    num_trigger_averages=0;
end

if ~isfield(d_para,'dts') || isempty(d_para.dts)
    d_para.dts=f_get_dt(spikes);  % or get it automatically from the precision of the spikes
end

if (~isfield(d_para,'tmax')  || isempty(d_para.tmax)) || (~isfield(d_para,'tmin')  || isempty(d_para.tmin))
    dmin=min(min(spikes));
    dmax=max(max(spikes));
    drange=dmax-dmin;
    d_para.tmin=dmin-0.02*drange;
    d_para.tmax=dmax+0.02*drange;
end

d_para.tmin=floor(d_para.tmin/d_para.dts)*d_para.dts;
d_para.tmax=floor(d_para.tmax/d_para.dts)*d_para.dts;

num_trains=size(spikes,1);
if num_trains<2
    error('It does not make any sense to quantify variability for just one spike train!')
end
num_pairs=(num_trains*num_trains-num_trains)/2;
num_pspikes=zeros(1,num_trains);
for trac=1:num_trains
    if any(spikes(trac,:))
        num_pspikes(trac)=find(spikes(trac,:),1,'last');
    end
end
max_num_pspikes=max(num_pspikes);

pspikes=zeros(num_trains,max_num_pspikes);          % original spikes used for plotting (in units of d_para.dts)
for trac=1:num_trains
    pspikes(trac,1:num_pspikes(trac))=round(sort(spikes(trac,1:num_pspikes(trac)))/d_para.dts)*d_para.dts;
end

% structure 's_para': parameters that describe the appearance of the individual subplots (measure time profiles)
%
% window_mode: time interval selection (1-all (recording limits),2:outer spikes,3-inner spikes,4-smaller analysis window)
% colors: Colors for 1:Stim/STs/ISI,2:Stim/STs/ISI-outside,3:Mean values,4:Mean values-outside,5:Moving averages,6:Moving averages-outside
% font_size: Font size of labels (and headlines), equals f_para.font_size
% line_width: Line width
% line_style: Line style (1-,2:,3-.,4--,5:none)
% nma: Selected moving averages (depends on f_para.ma_mode)
% spike_mao: Order of moving average (pooled ISI)
% time_mao: Order of moving average (time)
% causal: determines kind of moving average, set automatically for each measure (0-no,1-yes)
% itmin: Beginning of recording
% itmax: End of recording
% wmin: Beginning of selected window (if window_mode==4)
% wmax: End of selected window (if window_mode==4)
% num_subplots: depends on f_para.subplot_posi, set automatically
% xl: x-limits for plotting, set automatically
% yl: y-limits for plotting, set automatically
% dtm: equals f_para.dtm
% plot_mode: equals f_para.plot_mode

s_para_default=struct('window_mode',1,'colors','krkcgm','font_size',14,'line_width',1,'line_style','-',...
    'nma',1,'spike_mao',1,'time_mao',1,'causal',1,'itmin',[],'itmax',[],'wmin',[],'wmax',[],...
    'num_subplots',[],'xl',[],'yl',[],'dtm',[],'plot_mode',1);
s_para=s_para_default;
s_para.font_size=f_para.font_size;

if isfield(f_para,'ma_mode') && ~isempty(f_para.ma_mode)
    if f_para.ma_mode==0 s_para.nma=1; elseif f_para.ma_mode==1 s_para.nma=2; else s_para.nma=[1 2]; end %#ok<SEPEX>
end
if isfield(f_para,'spike_mao') && ~isempty(f_para.spike_mao)
    s_para.spike_mao=f_para.spike_mao;    % order of the moving average (pooled ISI)
end
if isfield(f_para,'time_mao') && ~isempty(f_para.time_mao)
    s_para.time_mao=f_para.time_mao;    % order of the moving average (time)
end

num_bins=100;    % number of bins for the PSTH
gs_window=5;     % window for Gaussian smoothing of PSTH


% 10:ISI,12:SPIKE-pico,14:rSPIKE-pico,16:SPIKE-time,18:rSPIKE-time,20:crSPIKE-time
mat_select2=intersect([10 12 14 16 18 20],find(f_para.subplot_posi));
[dummy,ms_indy]=sort(f_para.subplot_posi(mat_select2));
mat_select=mat_select2(ms_indy);

num_sel_measures=length(mat_select);

num_all_measures=length(f_para.subplot_posi);
psth_measures=[3 4];   int_measures=[5 6];   rate_measures=[7 8];   isi_measures=[9 10];
spike_pico_measures=[11 12];   realtime_spike_pico_measures=[13 14];
bi_pico_measures=[spike_pico_measures realtime_spike_pico_measures];
spike_time_measures=[15 16];   realtime_spike_time_measures=[17 18];  crealtime_spike_time_measures=[19 20];
realtime_measures=[realtime_spike_pico_measures realtime_spike_time_measures];
bi_measures=[isi_measures(2) spike_pico_measures(2) realtime_spike_pico_measures(2) spike_time_measures(2) ...
    realtime_spike_time_measures(2) crealtime_spike_time_measures(2)];
time_measures=[spike_time_measures realtime_spike_time_measures crealtime_spike_time_measures];

mat_pico_select=intersect([isi_measures bi_pico_measures],mat_select);
num_pico_mats=length(mat_pico_select);
mat_time_select=intersect(time_measures,mat_select);
num_time_mats=length(mat_time_select);


s_para.itmin=d_para.tmin;
s_para.itmax=d_para.tmax;
if s_para.window_mode==1                                                               % all (recording limits)
    s_para.wmin=s_para.itmin; s_para.wmax=s_para.itmax;
else
    if s_para.window_mode==2                                                           % outer spikes (overall)
        s_para.wmin=min(min(min(pspikes(pspikes~=0))));
        s_para.wmax=max(max(max(pspikes(pspikes~=0))));
    elseif s_para.window_mode==3                                                       % inner spikes (overall)
        s_para.wmin=max(pspikes(1:num_trains,1));
        s_para.wmax=inf;
        for trac=1:num_trains
            if num_pspikes(trac)>0
                if pspikes(trac,num_pspikes(trac))<s_para.itmax
                    s_para.wmax=pspikes(trac,num_pspikes(trac));
                end
            end
        end
    elseif s_para.window_mode==4 && isempty(s_para.wmin) && isempty(s_para.wmax)       % select
        s_para.wmin=d_para.tmin+0.25*(d_para.tmax-d_para.tmin);
        s_para.wmax=d_para.tmax-0.25*(d_para.tmax-d_para.tmin);
    end

    s_para.wmin=floor(s_para.wmin/d_para.dts)*d_para.dts;
    s_para.wmax=floor(s_para.wmax/d_para.dts)*d_para.dts;
    if s_para.wmin>=s_para.wmax
        disp(' '); disp(' '); error ('********** Error: W_max must be larger than W_min !!!!!')
    end
    if s_para.wmin<d_para.tmin
        disp(' '); disp(' '); error ('********** Error: W_min must not be smaller than T_min !!!!!')
    end
    if s_para.wmax>d_para.tmax
        disp(' '); disp(' '); error ('********** Error: W_max must not be larger than T_max !!!!!')
    end
end
s_para.itmin=floor(s_para.itmin/d_para.dts)*d_para.dts;
s_para.itmax=floor(s_para.itmax/d_para.dts)*d_para.dts;
itrange=s_para.itmax-s_para.itmin;
if s_para.window_mode==2
    s_para.wmin=max([s_para.itmin s_para.wmin]);
    s_para.wmax=min([s_para.itmax s_para.wmax]);
end


num_ispikes=zeros(1,num_trains);
ispikes=zeros(num_trains,max_num_pspikes+2);
for trac=1:num_trains
    dspikes=pspikes(trac,1:num_pspikes(trac));
    dspikes=dspikes(dspikes>=s_para.itmin & dspikes<=s_para.itmax);
    if any(dspikes)
        num_ispikes(trac)=find(dspikes,1,'last')+2;
        if isempty(find(dspikes(1:num_ispikes(trac)-2)==s_para.itmin,1))
            ispikes(trac,1:num_ispikes(trac))=[s_para.itmin dspikes(1:num_ispikes(trac)-2) s_para.itmax];
        else
            dummy=dspikes(1:num_ispikes(trac)-2);
            dummy2=[s_para.itmin dummy(dummy~=s_para.itmin) s_para.itmax];
            num_ispikes(trac)=length(dummy2);
            ispikes(trac,1:num_ispikes(trac))=dummy2;
        end
    else
        num_ispikes(trac)=2;
        ispikes(trac,1:num_ispikes(trac))=[s_para.itmin s_para.itmax];
    end
end
max_num_ispikes=max(num_ispikes);
ispikes=ispikes(1:num_trains,1:max_num_ispikes);

num_coins=zeros(num_trains,max_num_ispikes);
num_uspikes=zeros(1,num_trains);
uspikes=zeros(num_trains,max_num_ispikes);
for trac=1:num_trains
    dummy=unique(ispikes(trac,1:num_ispikes(trac)));
    for uic=1:length(dummy)-1
        num_coins(trac,uic)=sum(ispikes(trac,1:num_ispikes(trac))==dummy(uic));
    end
    num_uspikes(trac)=length(dummy);
    uspikes(trac,1:num_uspikes(trac))=dummy;
end
max_num_uspikes=max(num_uspikes);
uspikes=uspikes(1:num_trains,1:max_num_uspikes);
for trac=1:num_trains
    uspikes(trac,num_uspikes(trac)+1:max_num_uspikes)=987654321.123456789;
end
clear ispikes

if any(f_para.subplot_posi(psth_measures))                                              % PSTH
    bin_width=(d_para.tmax-d_para.tmin)/num_bins;
    bins=d_para.tmin:bin_width:d_para.tmax;
    bin_centers=bins(1:num_bins)+(bins(2)-bins(1))/2; %#ok<NASGU>
    dpsth=zeros(1,num_bins+1);
    for stc=1:num_trains
        dpsth=dpsth+histc(pspikes(stc,1:num_pspikes(stc)),bins);
    end
    psth=dpsth(1:num_bins)/bin_width;
    maxpsth=max(psth);
    if f_para.subplot_posi(psth_measures(2))                                                    % Gaussian smoothed PSTH
        gs_psth=f_compute_gauss_smooth(psth,gs_window)';
        maxgspsth=max(gs_psth);
    else
        maxgspsth=0;
    end
else
    maxpsth=0;
    maxgspsth=0;
end


test_diffs_pc1=0; test_diffs_pc2=0; test_diffs_pc3=0; test_diffs_pc4=0; test_diffs_t1=0; test_diffs_t2=0; test_diffs_t3=0; test_diffs_t4=0; % #########
mean_values=[];
result_str='';
if any(f_para.subplot_posi(5:num_all_measures))                                                        % ISI-SPIKE
    [all_spikes,all_indy]=sort(reshape(uspikes',1,numel(uspikes)));
    all_trains=fix((all_indy-1)/max_num_uspikes)+1;
    all_trains=all_trains(all_spikes~=987654321.123456789);
    all_spikes=all_spikes(all_spikes~=987654321.123456789);
    all_trains(1:num_trains)=0;
    all_trains(end-num_trains+1:end)=0;
    uspikes(uspikes==987654321.123456789)=0;
    clear all_indy

    all_spikes=all_spikes(num_trains:end-num_trains+1);
    all_trains=all_trains(num_trains:end-num_trains+1);
    all_isi=diff(all_spikes);
    num_all_isi=length(all_isi);
    isi=all_isi(all_isi>0);
    num_isi=length(isi);

    isi_pos=cumsum([s_para.itmin isi(1:end-1)])+isi/2;
    cum_isi=all_spikes(1)+[0 cumsum(isi)];
    clear all_spikes

    isis=zeros(num_trains,max_num_uspikes-1);
    for trac=1:num_trains
        isis(trac,1:num_uspikes(trac)-1)=diff(uspikes(trac,1:num_uspikes(trac)));
        isis(trac,isis(trac,:)~=0)=isis(trac,isis(trac,:)~=0)./num_coins(trac,isis(trac,:)~=0);
    end
    clear num_coins

    ints=zeros(num_trains,num_all_isi);
    for trac=1:num_trains
        ivs=[1 find(all_trains==trac)];
        ive=[ivs(2:length(ivs))-1 num_all_isi];
        for ic=1:num_uspikes(trac)-1
            ints(trac,ivs(ic):ive(ic))=isis(trac,ic);
        end
    end
    ints=ints(:,all_isi>0);
    mean_ints=mean(ints,1);
    %clear isis

    if any(f_para.subplot_posi(rate_measures))
        rates=1./ints;
        if f_para.subplot_posi(rate_measures(2))
            mean_rates=mean(rates,1);
        end
    end

    extra_mem=any(f_para.subplot_posi([spike_pico_measures spike_time_measures realtime_spike_pico_measures realtime_spike_time_measures crealtime_spike_time_measures]));
    if extra_mem                             % Pre
        udists=cell(num_trains);
        for trac1=1:num_trains
            if num_trains>=100 && max_num_uspikes>1000
                disp(['udistc = ',num2str([1 trac1])])
            end
            for trac2=setdiff(1:num_trains,trac1)
                udists{trac1,trac2}=zeros(1,num_uspikes(trac1));
            end
        end
        for trac1=1:num_trains
            if num_trains>=100 && max_num_uspikes>1000
                disp(['udistc = ',num2str([2 trac1])])
            end
            for trac2=setdiff(1:num_trains,trac1)
                for spc=1:num_uspikes(trac1)
                    udists{trac1,trac2}(spc)=min(abs(uspikes(trac1,spc)-uspikes(trac2,1:num_uspikes(trac2))));
                end
            end
        end
    end

    time=round((s_para.itmin+f_para.dtm)/f_para.dtm)*f_para.dtm:f_para.dtm:round(s_para.itmax/f_para.dtm)*f_para.dtm;
    len=length(time);
    if any(f_para.subplot_posi(time_measures)) || mod(f_para.plot_mode,32)>7
        start=zeros(num_trains,max(num_uspikes)-1);  % integers realtive to time vector
        ende=zeros(num_trains,max(num_uspikes)-1);
        for trac=1:num_trains
            for sc=1:num_uspikes(trac)-1
                start(trac,sc)=fix((uspikes(trac,sc)+f_para.dtm/10-s_para.itmin)/f_para.dtm)+1;
                ende(trac,sc)=fix((uspikes(trac,sc+1)+f_para.dtm/10-s_para.itmin)/f_para.dtm);
            end
        end
    end

    % ###########################################################################################################################################
    % ##################################################################### Pico-Quantities #####################################################
    % ###########################################################################################################################################

    if any(f_para.subplot_posi([spike_pico_measures realtime_spike_pico_measures]))  % SPIKE-Pre-Pico
        previs=zeros(num_trains,max_num_uspikes-1);
        for trac=1:num_trains
            previs(trac,1:num_uspikes(trac)-1)=uspikes(trac,1:num_uspikes(trac)-1);
        end
        prev_spikes=zeros(num_trains,num_all_isi);
        for trac=1:num_trains
            ivs=[1 find(all_trains==trac)];
            ive=[ivs(2:length(ivs))-1 num_all_isi];
            for ic=1:num_uspikes(trac)-1
                prev_spikes(trac,ivs(ic):ive(ic))=previs(trac,ic);
            end
        end
        prev_spikes=prev_spikes(:,all_isi>0);
        clear previs
    end

    if any(f_para.subplot_posi([spike_pico_measures]))                                                 % SPIKE-Pico
        follis=zeros(num_trains,max_num_uspikes-1);
        for trac=1:num_trains
            follis(trac,1:num_uspikes(trac)-1)=uspikes(trac,2:num_uspikes(trac));
        end
        foll_spikes=zeros(num_trains,num_all_isi);
        for trac=1:num_trains
            ivs=[1 find(all_trains==trac)];
            ive=[ivs(2:length(ivs))-1 num_all_isi];
            for ic=1:num_uspikes(trac)-1
                foll_spikes(trac,ivs(ic):ive(ic))=follis(trac,ic);
            end
        end
        foll_spikes=foll_spikes(:,all_isi>0);
        clear follis
    end
    %clear uspikes

    if any(f_para.subplot_posi([spike_pico_measures realtime_spike_pico_measures]))                                         % SPIKE-Pre-Pico
        previs_indy=zeros(num_trains,max_num_uspikes-1);
        for trac=1:num_trains
            previs_indy(trac,1:num_uspikes(trac)-1)=1:num_uspikes(trac)-1;
        end
        prev_spikes_indy=zeros(num_trains,num_all_isi);
        for trac=1:num_trains
            ivs=[1 find(all_trains==trac)];
            ive=[ivs(2:length(ivs))-1 num_all_isi];
            for ic=1:num_uspikes(trac)-1
                prev_spikes_indy(trac,ivs(ic):ive(ic))=previs_indy(trac,ic);
            end
        end
        prev_spikes_indy=prev_spikes_indy(:,all_isi>0);
        clear previs_indy
    end

    if any(f_para.subplot_posi(spike_pico_measures))                                                                       % SPIKE-Pico
        follis_indy=zeros(num_trains,max_num_uspikes-1);
        for trac=1:num_trains
            follis_indy(trac,1:num_uspikes(trac)-1)=2:num_uspikes(trac);
        end
        foll_spikes_indy=zeros(num_trains,num_all_isi);
        for trac=1:num_trains
            ivs=[1 find(all_trains==trac)];
            ive=[ivs(2:length(ivs))-1 num_all_isi];
            for ic=1:num_uspikes(trac)-1
                foll_spikes_indy(trac,ivs(ic):ive(ic))=follis_indy(trac,ic);
            end
        end
        foll_spikes_indy=foll_spikes_indy(:,all_isi>0);
        clear follis_indy
    end
    clear all_isi; clear all_trains

    % ###########################################################################################################################################
    % ##################################################################### Time-Quantities #####################################################
    % ###########################################################################################################################################

    if any(f_para.subplot_posi([spike_time_measures realtime_spike_time_measures crealtime_spike_time_measures]))
        previ=zeros(num_trains,len);
        for trac=1:num_trains
            for sc=1:num_uspikes(trac)-1
                previ(trac,start(trac,sc):ende(trac,sc))=(1+(0:(ende(trac,sc)-start(trac,sc))))*f_para.dtm;            % distance to previous spike
            end
        end

        if any(f_para.subplot_posi([spike_time_measures realtime_spike_time_measures crealtime_spike_time_measures]))
            prev_indy=zeros(num_trains,len);
            for trac=1:num_trains
                for sc=1:num_uspikes(trac)-1
                    prev_indy(trac,start(trac,sc):ende(trac,sc))=sc;                                            % index of preceding spike
                end
            end
        end

        if any(f_para.subplot_posi(spike_time_measures))
            folli=zeros(num_trains,len);
            for trac=1:num_trains
                for sc=1:num_uspikes(trac)-1
                    folli(trac,start(trac,sc):ende(trac,sc))=(ende(trac,sc)-start(trac,sc))*f_para.dtm:-f_para.dtm:0;           % distance to following spike
                end
            end

            ustart=unique(start(start>0))';
            uende=unique(ende(ende>0))';
            isi_indy=zeros(1,len);
            for ic=1:length(ustart)
                isi_indy(ustart(ic):uende(ic))=ic;
            end
        end
    end

    % ###########################################################################################################################################
    % ############################################################## Memory management ##########################################################
    % ###########################################################################################################################################

    memo=num_pico_mats*num_pairs*num_isi+num_time_mats*num_pairs*len;
    extra_memo=extra_mem*sum(num_uspikes*num_trains);
    if max_memo_init<extra_memo
        error(['Error 1: Please increase the value of the variable ''max_memo'' !!!'])
    end
    max_memo=max_memo_init-extra_memo;
    if memo>max_memo-extra_mem*sum(num_uspikes*(num_trains-1))
        if num_pico_mats>0
            max_num_isi=fix((max_memo-num_time_mats*num_pairs*len)/(num_pico_mats*num_pairs));
        else
            max_num_isi=inf;
        end
        if num_time_mats>0
            max_len=fix((max_memo-num_pico_mats*num_pairs*num_isi)/(num_time_mats*num_pairs));
        else
            max_len=inf;
        end
        [num_runs,run_indy]=max([ceil(num_isi/max_num_isi) ceil(len/max_len)]);

        if num_runs==0
            disp(' '); disp(' ');
            error(['Error 2: Please increase the value of the variable ''max_memo'' !!!'])
        end

        if run_indy==1
            run_isi_lengths=[max_num_isi*ones(1,num_runs-1) num_isi-max_num_isi*(num_runs-1)];
            run_isi_ends=cumsum(run_isi_lengths);
            run_isi_starts=min([run_isi_ends; [1 run_isi_ends(1:end-1)+1]]);

            run_time_ends=zeros(1,num_runs);
            for ruc=1:num_runs
                run_time_ends(ruc)=find(time<=d_para.tmin+cum_isi(run_isi_ends(ruc)+1),1,'last');
            end
            run_time_starts=[1 run_time_ends(1:end-1)+1];
            run_time_lengths=run_time_ends-run_time_starts+1;
        else
            run_time_lengths_init=fix([max_len*ones(1,num_runs-1) len-max_len*(num_runs-1)]);   % check f_para.dtm
            run_time_ends_init=cumsum(run_time_lengths_init);
            run_time_starts_init=[1 run_time_ends_init(1:end-1)+1];

            run_isi_ends=zeros(1,num_runs);
            for ruc=1:num_runs
                run_isi_ends(ruc)=find(d_para.tmin+cum_isi<=run_time_ends_init(ruc)*f_para.dtm,1,'last')-1;
            end
            run_isi_starts=[1 run_isi_ends(1:end-1)+1];
            run_isi_lengths=run_isi_ends-run_isi_starts+1;

            run_time_ends=zeros(1,num_runs);
            for ruc=1:num_runs
                run_time_ends(ruc)=find(time<=d_para.tmin+cum_isi(run_isi_ends(ruc)+1),1,'last');
            end
            run_time_starts=[1 run_time_ends(1:end-1)+1];
            run_time_lengths=run_time_ends-run_time_starts+1;
        end
    else
        num_runs=1;
        run_isi_lengths=num_isi;
        run_isi_starts=1;
        run_isi_ends=num_isi;
        run_time_lengths=len;
        run_time_starts=1;
        run_time_ends=len;
    end

    if any(f_para.subplot_posi(isi_measures))                                                         % ISI
        isi_ratio=zeros(1,num_isi);
    end
    if f_para.subplot_posi(spike_pico_measures(2))                                                    % SPIKE-PICO
        spike_diffs=zeros(1,num_isi);
    end
    if f_para.subplot_posi(realtime_spike_pico_measures(2))                                           % REALTIME-SPIKE-PICO
        spike_diffs_realtime=zeros(1,num_isi);
    end
    if f_para.subplot_posi(spike_time_measures(2))                                                    % SPIKE-Time
        spike_diffs_t=zeros(1,len);
    end
    if f_para.subplot_posi(realtime_spike_time_measures(2))                                           % REALTIME-SPIKE-Time
        spike_diffs_realtime_t=zeros(1,len);
    end
    if f_para.subplot_posi(crealtime_spike_time_measures(2))                                           % cREALTIME-SPIKE-Time
        spike_diffs_crealtime_t=zeros(1,len);
    end
    for ruc=1:num_runs
        run_isi_range=run_isi_starts(ruc):run_isi_ends(ruc);
        run_time_range=run_time_starts(ruc):run_time_ends(ruc);
        if num_runs>1
            disp(['run_info = ',num2str(ruc),'  (',num2str(num_runs),')'])
        end

        if any(f_para.subplot_posi(isi_measures))                                                         % ISI
            bi_isi_ratio=zeros(num_pairs,run_isi_lengths(ruc));
            pac=0;
            for trac1=1:num_trains-1
                for trac2=trac1+1:num_trains
                    pac=pac+1;
                    dummy1=find(ints(trac1,run_isi_range)<ints(trac2,run_isi_range));
                    bi_isi_ratio(pac,dummy1)=ints(trac1,run_isi_range(dummy1))./ints(trac2,run_isi_range(dummy1))-1;
                    dummy2=find(ints(trac1,run_isi_range)>=ints(trac2,run_isi_range) & ints(trac1,run_isi_range)~=0);
                    bi_isi_ratio(pac,dummy2)=-(ints(trac2,run_isi_range(dummy2))./ints(trac1,run_isi_range(dummy2))-1);
                end
            end
            if f_para.subplot_posi(isi_measures(2))
                isi_ratio(run_isi_range)=mean(abs(bi_isi_ratio),1);
            end
            clear dummy1 dummy2
            if ~f_para.subplot_posi(isi_measures(1)) && (mod(f_para.plot_mode,32)<2 || ~any(ismember(isi_measures,mat_select)))    % ############
                clear bi_isi_ratio
            end
        end

        if any(f_para.subplot_posi([spike_pico_measures]))                                         % SPIKE-Pico
            bi_spike_diffs=zeros(num_pairs,run_isi_lengths(ruc));
            pac=0;
            for trac1=1:num_trains-1
                for trac2=trac1+1:num_trains
                    pac=pac+1;
                    bi_spike_diffs(pac,1:run_isi_lengths(ruc))=...             % weighting
                        ((udists{trac1,trac2}(prev_spikes_indy(trac1,run_isi_range)).*(foll_spikes(trac1,run_isi_range)-isi_pos(run_isi_range))+...
                        udists{trac1,trac2}(foll_spikes_indy(trac1,run_isi_range)).*(isi_pos(run_isi_range)-prev_spikes(trac1,run_isi_range)))./ints(trac1,run_isi_range).*ints(trac2,run_isi_range)+...
                        (udists{trac2,trac1}(prev_spikes_indy(trac2,run_isi_range)).*(foll_spikes(trac2,run_isi_range)-isi_pos(run_isi_range))+...
                        udists{trac2,trac1}(foll_spikes_indy(trac2,run_isi_range)).*(isi_pos(run_isi_range)-prev_spikes(trac2,run_isi_range)))./ints(trac2,run_isi_range).*ints(trac1,run_isi_range))./...
                        ((ints(trac1,run_isi_range)+ints(trac2,run_isi_range)).^2/2);
                end
            end
            if f_para.subplot_posi(spike_pico_measures(2))
                spike_diffs(run_isi_range)=mean(bi_spike_diffs,1);
            end
        end

        if any(f_para.subplot_posi(spike_time_measures))                                         % SPIKE-Time
            bi_spike_diffs_t=zeros(num_pairs,run_time_lengths(ruc));
            pac=0;
            for trac1=1:num_trains-1
                for trac2=trac1+1:num_trains
                    pac=pac+1;
                    dummy=isi_indy(run_time_range)-run_isi_starts(ruc)+1;
                    bi_spike_diffs_t(pac,1:run_time_lengths(ruc))=...             % weighting
                        ((udists{trac1,trac2}(prev_indy(trac1,run_time_range)).*folli(trac1,run_time_range)+...
                        udists{trac1,trac2}(prev_indy(trac1,run_time_range)+1).*previ(trac1,run_time_range))./ints(trac1,isi_indy(run_time_range)).*ints(trac2,isi_indy(run_time_range))+...
                        (udists{trac2,trac1}(prev_indy(trac2,run_time_range)).*folli(trac2,run_time_range)+...
                        udists{trac2,trac1}(prev_indy(trac2,run_time_range)+1).*previ(trac2,run_time_range))./ints(trac2,isi_indy(run_time_range)).*ints(trac1,isi_indy(run_time_range)))./...
                        ((ints(trac1,run_isi_range(dummy))+ints(trac2,run_isi_range(dummy))).^2/2);
                end
            end
            if f_para.subplot_posi(spike_time_measures(2))
                spike_diffs_t(run_time_range)=mean(bi_spike_diffs_t,1);
            end
            if ~f_para.subplot_posi(spike_time_measures(1)) && (mod(f_para.plot_mode,32)<2 || ~any(ismember(spike_time_measures,mat_select)))
                clear bi_spike_diffs_t
            end
        end

        % ######################################################################################

        if any(f_para.subplot_posi(realtime_spike_pico_measures))                                              % REALTIME-PICO
            bi_spike_diffs_realtime=zeros(num_pairs,run_isi_lengths(ruc));
            pac=0;
            for trac1=1:num_trains-1
                for trac2=trac1+1:num_trains
                    pac=pac+1;

                    run_pre_diffs=abs(prev_spikes(trac1,run_isi_range)-prev_spikes(trac2,run_isi_range));
                    run_max_pre_spikes=max(prev_spikes([trac1 trac2],run_isi_range));

                    dummy=(prev_spikes(trac1,run_isi_range)<prev_spikes(trac2,run_isi_range)-0.00000001);
                    dummy1=trac1*dummy+trac2*(1-dummy);   % index of spike train with earlier spike
                    dummy2=trac2*dummy+trac1*(1-dummy);   % index of spike train with later spike

                    for ic=1:run_isi_lengths(ruc)
                        dudists=udists{dummy1(ic),dummy2(ic)}(prev_spikes_indy(dummy1(ic),run_isi_range(ic)));
                        bi_spike_diffs_realtime(pac,ic)=(run_pre_diffs(ic)+dudists)/4*(log(run_pre_diffs(ic)/2+cum_isi(run_isi_range(ic)+1)-run_max_pre_spikes(ic)+...
                            (run_pre_diffs(ic)/2 + cum_isi(run_isi_range(ic)+1)==run_max_pre_spikes(ic))) - ...
                            log(run_pre_diffs(ic)/2 + cum_isi(run_isi_range(ic))-run_max_pre_spikes(ic) + ...
                            (run_pre_diffs(ic)/2 + cum_isi(run_isi_range(ic))==run_max_pre_spikes(ic))) ) / isi(run_isi_range(ic));
                    end
                end
            end
            spike_diffs_realtime(run_isi_range)=mean(bi_spike_diffs_realtime,1);
            clear run_pre_diffs run_max_prev_spikes;
            if ~f_para.subplot_posi(realtime_spike_pico_measures(1)) && (mod(f_para.plot_mode,32)<2 || ~any(ismember(realtime_spike_pico_measures,mat_select)))
                clear bi_spike_diffs_realtime
            end
        end


        if any(f_para.subplot_posi(realtime_spike_time_measures))                                              % REALTIME-TIME

            bi_spike_diffs_realtime_t=zeros(num_pairs,run_time_lengths(ruc));
            pac=0;
            for trac1=1:num_trains-1
                for trac2=trac1+1:num_trains
                    pac=pac+1;
                    normy=(previ(trac1,run_time_range)+previ(trac2,run_time_range));
                    dummy=(previ(trac1,run_time_range)<previ(trac2,run_time_range)-0.00000001);
                    dummy1=trac1*(1-dummy)+trac2*dummy;   % index of spike train with earlier spike
                    dummy2=trac2*(1-dummy)+trac1*dummy;   % index of spike train with later spike
                    for sc=1:run_time_lengths(ruc)
                        bi_spike_diffs_realtime_t(pac,sc)=(abs(previ(trac1,run_time_range(sc))-previ(trac2,run_time_range(sc)))+...   % later spike
                            udists{dummy1(sc),dummy2(sc)}(prev_indy(dummy1(sc),run_time_range(sc))))...      % earlier spike
                            /(2*normy(sc)+(normy(sc)==0));
                    end
                end
            end
            clear dummy; clear dummy1; clear dummy2
            if f_para.subplot_posi(realtime_spike_time_measures(2))
                spike_diffs_realtime_t(run_time_range)=mean(bi_spike_diffs_realtime_t,1);
            end
            if ~f_para.subplot_posi(realtime_spike_time_measures(1)) && (mod(f_para.plot_mode,32)<2 || ~any(ismember(realtime_spike_time_measures,mat_select)))
                clear bi_spike_diffs_realtime_t
            end
        end

        if any(f_para.subplot_posi(crealtime_spike_time_measures))                                              % cREALTIME-TIME

            xisip=zeros(num_trains,run_time_lengths(ruc));
            cxisi=zeros(num_trains,run_time_lengths(ruc));
            run_indy=run_time_range-run_time_starts(ruc)+1;
            for trac=1:num_trains
                xisip(trac,prev_indy(trac,run_time_range)==1)=0;
                xisip(trac,prev_indy(trac,run_time_range)>1)=isis(trac,prev_indy(trac,run_time_range(prev_indy(trac,run_time_range)>1))-1);
                cxisi(trac,1:run_time_lengths(ruc))=(2*previ(trac,run_time_range) + heaviside(xisip(trac,run_indy)-2*previ(trac,run_time_range)).*xisip(trac,run_indy) )./...
                    (1+heaviside(xisip(trac,run_indy)-2*previ(trac,run_time_range)));
            end

            bi_spike_diffs_crealtime_t=zeros(num_pairs,run_time_lengths(ruc));
            pac=0;
            for trac1=1:num_trains-1
                for trac2=trac1+1:num_trains
                    pac=pac+1;
                    normy=(previ(trac1,run_time_range)+previ(trac2,run_time_range));
                    normy3=(cxisi(trac1,run_indy)+cxisi(trac2,run_indy))/2;

                    dummy=(previ(trac1,run_time_range)<previ(trac2,run_time_range)-0.00000001);
                    dummy1=trac1*(1-dummy)+trac2*dummy;   % index of spike train with earlier spike
                    dummy2=trac2*(1-dummy)+trac1*dummy;   % index of spike train with later spike
                    for sc=1:run_time_lengths(ruc)
                        bi_spike_diffs_crealtime_t(pac,sc)=(abs(previ(trac1,run_time_range(sc))-previ(trac2,run_time_range(sc)))+...   % later spike
                            udists{dummy1(sc),dummy2(sc)}(prev_indy(dummy1(sc),run_time_range(sc))))...      % earlier spike
                            ./(2*normy3(sc)+(normy3(sc)==0));
                    end
                end
            end
            clear dummy; clear dummy1; clear dummy2
            if f_para.subplot_posi(crealtime_spike_time_measures(2))
                spike_diffs_crealtime_t(run_time_range)=mean(bi_spike_diffs_crealtime_t,1);
            end
            if ~f_para.subplot_posi(crealtime_spike_time_measures(1)) && (mod(f_para.plot_mode,32)<2 || ~any(ismember(crealtime_spike_time_measures,mat_select)))
                clear bi_spike_diffs_crealtime_t
            end
        end

        % ##################################################################################################################################
        % ############################################################# Movie ##############################################################
        % ##################################################################################################################################


        if mod(f_para.plot_mode,32)>1                                                                     % ############################################

            len_isi_indy=zeros(1,run_time_lengths(ruc));
            for sc=1:run_time_lengths(ruc)
                [dummyval,dummypos]=min(abs(time(run_time_range(sc))-isi_pos(run_isi_range)));
                len_isi_indy(sc)=dummypos(1);
            end

            if isfield(f_para,'mov_frames') && ~isempty(f_para.mov_frames)
                frame_select=f_para.mov_frames;
            elseif isfield(f_para,'mov_step') && ~isempty(f_para.mov_step)
                dummy=unique([f_para.mov_step:f_para.mov_step:len]);
                if isempty(intersect(mat_select,time_measures)) && num_isi<length(dummy)              % ############################################
                    frame_select=round(unique(isi_pos));
                    frame_select(frame_select==0)=1;
                else
                    frame_select=dummy;
                end
            else
                frame_select=[];
            end
            frame_select=frame_select(frame_select>=run_time_starts(ruc) & frame_select<=run_time_ends(ruc));
            run_frame_select=frame_select-run_time_starts(ruc)+1;
            num_frame_select=length(frame_select);

            mat_str=[];
            if ruc==1
                mat_indy=nchoosek(1:num_trains,2);
            end
            mov_mat=zeros(num_sel_measures,num_frame_select+num_select_averages+num_trigger_averages,num_trains,num_trains);
            for matc=1:num_sel_measures
                switch mat_select(matc)
                    case 10,
                        mat_str{matc}='ISI';
                        dist_mat=abs(bi_isi_ratio);
                    case 12,
                        mat_str{matc}='S';
                        dist_mat=bi_spike_diffs;
                    case 14,
                        mat_str{matc}='S_r';
                        dist_mat=bi_spike_diffs_realtime;
                    case 16,
                        mat_str{matc}='S';
                        dist_mat=bi_spike_diffs_t;
                    case 18,
                        mat_str{matc}='S_r';
                        dist_mat=bi_spike_diffs_realtime_t;
                    case 20,
                        mat_str{matc}='S_{cr}';
                        dist_mat=bi_spike_diffs_crealtime_t;
                end
                if num_select_averages>0
                    if matc==1 && ruc==1
                        mov_mat_sa=zeros(num_sel_measures,num_select_averages,num_trains,num_trains);
                        mov_mat_sa_weight=zeros(num_sel_measures,num_select_averages,num_runs);
                    end
                    num_sel_ave=zeros(1,num_select_averages);
                    for sac=1:num_select_averages
                        num_sel_ave(sac)=length(d_para.select_averages{sac})/2;
                    end
                end
                if num_trigger_averages>0
                    if matc==1 && ruc==1
                        mov_mat_ta=zeros(num_sel_measures,num_trigger_averages,num_trains,num_trains);
                        mov_mat_ta_weight=zeros(num_sel_measures,num_trigger_averages,num_runs);
                    end
                end
                if mat_select(matc)<time_measures(1)                                            % pico
                    if ~isempty(run_frame_select)
                        mov_mat(matc,1:num_frame_select,sub2ind([num_trains num_trains],mat_indy(:,1),mat_indy(:,2)))=dist_mat(:,len_isi_indy(run_frame_select))';
                        mov_mat(matc,1:num_frame_select,sub2ind([num_trains num_trains],mat_indy(:,2),mat_indy(:,1)))=dist_mat(:,len_isi_indy(run_frame_select))';
                    end
                    if num_select_averages>0
                        sel_ave_weight=zeros(num_select_averages,max(num_sel_ave));
                        for sac=1:num_select_averages
                            for selc=1:num_sel_ave(sac)
                                sel_ave_weight(sac,selc)=diff(d_para.select_averages{sac}(2*selc-1:2*selc));
                            end
                        end
                        for sac=1:num_select_averages
                            for selc=1:num_sel_ave(sac)
                                pfx=unique([cum_isi(run_isi_starts(ruc):run_isi_ends(ruc)+1) cum_isi(run_isi_starts(ruc)+1:run_isi_ends(ruc))-f_para.dtm/100]);
                                pfy=reshape([dist_mat; dist_mat],num_pairs,run_isi_lengths(ruc)*2);
                                first_winspike_ind=find(pfx>=d_para.select_averages{sac}(2*selc-1),1,'first');
                                last_winspike_ind=find(pfx<=d_para.select_averages{sac}(2*selc),1,'last');
                                if ~isempty(first_winspike_ind) && ~isempty(last_winspike_ind)
                                    pfx=pfx(first_winspike_ind:last_winspike_ind);
                                    if d_para.select_averages{sac}(2*selc-1)<pfx(1) && d_para.select_averages{sac}(2*selc-1)>=cum_isi(run_isi_starts(ruc))  % interval to first spike
                                        pfy_m=pfy(:,first_winspike_ind-1);
                                    end
                                    if d_para.select_averages{sac}(2*selc)>pfx(end) && d_para.select_averages{sac}(2*selc)<=cum_isi(run_isi_ends(ruc)+1)   % interval after last spike
                                        pfy_p=pfy(:,last_winspike_ind+1);
                                    end
                                    pfy=pfy(:,first_winspike_ind:last_winspike_ind);
                                    if d_para.select_averages{sac}(2*selc-1)<pfx(1) && d_para.select_averages{sac}(2*selc-1)>=cum_isi(run_isi_starts(ruc))  % interval to first spike
                                        pfx=[d_para.select_averages{sac}(2*selc-1) pfx];
                                        pfy=[pfy_m pfy];
                                    end
                                    if d_para.select_averages{sac}(2*selc)>pfx(end) && d_para.select_averages{sac}(2*selc)<=cum_isi(run_isi_ends(ruc)+1)   % interval after last spike
                                        pfx=[pfx d_para.select_averages{sac}(2*selc)];
                                        pfy=[pfy pfy_p];
                                    end
                                    if length(pfx)>1
                                        pfx=diff(pfx);
                                        pfy=pfy(:,1:end-1);
                                        ave=sum(pfy.*repmat(pfx,size(pfy,1),1),2)/sum(pfx);
                                        mov_mat_sa_weight(matc,sac,ruc)=sum(pfx);
                                        mov_mat_sa(matc,sac,sub2ind([num_trains num_trains],mat_indy(:,1),mat_indy(:,2)))=...
                                            mov_mat_sa(matc,sac,sub2ind([num_trains num_trains],mat_indy(:,1),mat_indy(:,2)))+shiftdim(ave'*mov_mat_sa_weight(matc,sac,ruc)*sel_ave_weight(sac,selc),-1);
                                        mov_mat_sa(matc,sac,sub2ind([num_trains num_trains],mat_indy(:,2),mat_indy(:,1)))=...
                                            mov_mat_sa(matc,sac,sub2ind([num_trains num_trains],mat_indy(:,2),mat_indy(:,1)))+shiftdim(ave'*mov_mat_sa_weight(matc,sac,ruc)*sel_ave_weight(sac,selc),-1);
                                    end
                                end
                            end
                        end
                        clear pfx; clear pfy;
                    end
                    if num_trigger_averages>0
                        for tac=1:num_trigger_averages
                            indy_trigger=intersect(ceil((d_para.trigger_averages{tac}-d_para.tmin)/f_para.dtm),run_time_range)-run_time_starts(ruc)+1;  % test f_para.dtm>1
                            if ~isempty(indy_trigger)
                                mov_mat_ta_weight(matc,tac,ruc)=length(indy_trigger);
                                mov_mat_ta(matc,tac,sub2ind([num_trains num_trains],mat_indy(:,1),mat_indy(:,2)))=...
                                    mov_mat_ta(matc,tac,sub2ind([num_trains num_trains],mat_indy(:,1),mat_indy(:,2)))+shiftdim(mean(dist_mat(:,len_isi_indy(indy_trigger)),2)'*mov_mat_ta_weight(matc,tac,ruc),-1);
                                mov_mat_ta(matc,tac,sub2ind([num_trains num_trains],mat_indy(:,2),mat_indy(:,1)))=...
                                    mov_mat_ta(matc,tac,sub2ind([num_trains num_trains],mat_indy(:,2),mat_indy(:,1)))+shiftdim(mean(dist_mat(:,len_isi_indy(indy_trigger)),2)'*mov_mat_ta_weight(matc,tac,ruc),-1);
                            end
                        end

                    end
                else                                                                             % time
                    if ~isempty(run_frame_select)
                        mov_mat(matc,1:num_frame_select,sub2ind([num_trains num_trains],mat_indy(:,1),mat_indy(:,2)))=dist_mat(:,run_frame_select)';
                        mov_mat(matc,1:num_frame_select,sub2ind([num_trains num_trains],mat_indy(:,2),mat_indy(:,1)))=dist_mat(:,run_frame_select)';
                    end
                    if num_select_averages>0
                        sel_ave_weight=zeros(num_select_averages,max(num_sel_ave));
                        for sac=1:num_select_averages
                            for selc=1:num_sel_ave(sac)
                                sel_ave_weight(sac,selc)=diff(d_para.select_averages{sac}(2*selc-1:2*selc));
                                indy_select=find(time(run_time_range)>=d_para.select_averages{sac}(2*selc-1) & time(run_time_range)<=d_para.select_averages{sac}(2*selc));
                                if ~isempty(indy_select)
                                    mov_mat_sa_weight(matc,sac,ruc)=length(indy_select);
                                    mov_mat_sa(matc,sac,sub2ind([num_trains num_trains],mat_indy(:,1),mat_indy(:,2)))=...
                                        mov_mat_sa(matc,sac,sub2ind([num_trains num_trains],mat_indy(:,1),mat_indy(:,2)))+shiftdim(mean(dist_mat(:,indy_select),2)'*mov_mat_sa_weight(matc,sac,ruc)*sel_ave_weight(sac,selc),-1);
                                    mov_mat_sa(matc,sac,sub2ind([num_trains num_trains],mat_indy(:,2),mat_indy(:,1)))=...
                                        mov_mat_sa(matc,sac,sub2ind([num_trains num_trains],mat_indy(:,2),mat_indy(:,1)))+shiftdim(mean(dist_mat(:,indy_select),2)'*mov_mat_sa_weight(matc,sac,ruc)*sel_ave_weight(sac,selc),-1);
                                end
                            end
                        end
                    end
                    if num_trigger_averages>0
                        for tac=1:num_trigger_averages
                            indy_trigger=intersect(ceil((d_para.trigger_averages{tac}-d_para.tmin)/f_para.dtm),run_time_range)-run_time_starts(ruc)+1;  % test f_para.dtm>1
                            if ~isempty(indy_trigger)
                                mov_mat_ta_weight(matc,tac,ruc)=length(indy_trigger);
                                mov_mat_ta(matc,tac,sub2ind([num_trains num_trains],mat_indy(:,1),mat_indy(:,2)))=...
                                    mov_mat_ta(matc,tac,sub2ind([num_trains num_trains],mat_indy(:,1),mat_indy(:,2)))+shiftdim(mean(dist_mat(:,indy_trigger),2)'*mov_mat_ta_weight(matc,tac,ruc),-1);
                                mov_mat_ta(matc,tac,sub2ind([num_trains num_trains],mat_indy(:,2),mat_indy(:,1)))=...
                                    mov_mat_ta(matc,tac,sub2ind([num_trains num_trains],mat_indy(:,2),mat_indy(:,1)))+shiftdim(mean(dist_mat(:,indy_trigger),2)'*mov_mat_ta_weight(matc,tac,ruc),-1);
                            end
                        end
                    end
                end
            end
            if ruc==num_runs
                if num_select_averages>0
                    mov_mat(1:num_sel_measures,num_frame_select+(1:num_select_averages),1:num_trains,1:num_trains)=mov_mat_sa(1:num_sel_measures,1:num_select_averages,1:num_trains,1:num_trains)./...
                        permute(repmat(shiftdim(sum(mov_mat_sa_weight(1:num_sel_measures,1:num_select_averages,1:num_runs),3),-2),num_trains,num_trains),[3 4 1 2]);
                    for sac=1:num_select_averages
                        mov_mat(1:num_sel_measures,num_frame_select+sac,1:num_trains,1:num_trains)=...
                            mov_mat(1:num_sel_measures,num_frame_select+sac,1:num_trains,1:num_trains)/sum(sel_ave_weight(sac,1:num_sel_ave(sac)));
                    end
                    clear mov_mat_sa; clear mov_mat_sa_weight;
                end
                if num_trigger_averages>0
                    mov_mat(1:num_sel_measures,num_frame_select+num_select_averages+(1:num_trigger_averages),1:num_trains,1:num_trains)=mov_mat_ta(1:num_sel_measures,1:num_trigger_averages,1:num_trains,1:num_trains)./...
                        permute(repmat(shiftdim(sum(mov_mat_ta_weight(1:num_sel_measures,1:num_trigger_averages,1:num_runs),3),-2),num_trains,num_trains),[3 4 1 2]);
                    clear mov_mat_ta; clear mov_mat_ta_weight;
                end
            end
            clear dist_mat;

            if f_para.block_matrices && num_select_train_groups>1 && num_select_train_groups<num_trains % && ruc==num_runs
                block_mov_mat=zeros(num_sel_measures,num_frame_select+num_select_averages+num_trigger_averages,num_select_train_groups,num_select_train_groups);
                for sgc1=1:num_select_train_groups
                    for sgc2=1:num_select_train_groups
                        block_mov_mat(1:num_sel_measures,1:num_frame_select+num_select_averages+num_trigger_averages,sgc1,sgc2)=...
                            mean(mean(mov_mat(1:num_sel_measures,1:num_frame_select+num_select_averages+num_trigger_averages,select_group_vect==select_groups(sgc1),...
                            select_group_vect==select_groups(sgc2)),4),3);
                    end
                end
                for matc=1:num_sel_measures
                    mat_str{num_sel_measures+matc}=['< ',mat_str{matc},' >_G'];
                end
            end
            num_mat_subplots=num_sel_measures*(1+(f_para.block_matrices && num_select_train_groups>1 && num_select_train_groups<num_trains));
            num_all_subplots=num_mat_subplots*(1+f_para.dendrograms);         % all subplots

            if mod(f_para.plot_mode,4)>1                                                       % Snapshots, different measures and cuts
                Distances_Comp
            end

            if mod(f_para.plot_mode,8)>3                                                       % Snapshots, different cuts (subplots) for each measure (figure)
                Distances_Cuts
            end

            if mod(f_para.plot_mode,32)>7                                                       % Movie
                Distances_Movie
            end
        end
    end
end

if ~any(f_para.subplot_posi==1)
    f_para.subplot_posi=(f_para.subplot_posi-min(f_para.subplot_posi(f_para.subplot_posi>0))+1).*(f_para.subplot_posi>0);
end
s_para.num_subplots=max(f_para.subplot_posi);

if mod(f_para.plot_mode,2)>0
    if s_para.window_mode==1                                                              % limits for plotting
        pmin=s_para.itmin-0.02*itrange;
        pmax=s_para.itmax+0.02*itrange;
    else
        pmax2=-inf;
        for trac=1:num_trains
            if num_pspikes(trac)>0
                pmax2=max([pmax2 pspikes(trac,num_pspikes(trac))]);
            end
        end
        pmin2=min(pspikes(1:num_trains,1));
        prange2=pmax2-pmin2;
        pmin=pmin2-0.02*prange2;
        pmax=pmax2+0.02*prange2;
    end


    if any(f_para.subplot_posi(int_measures))
        if f_para.subplot_posi(int_measures(1))>0
            maxintval=max(max(ints));
        end
        maxmeanintval=max(mean_ints);
    else
        maxintval=0;
        maxmeanintval=0;
    end
    if any(f_para.subplot_posi(rate_measures))
        if f_para.subplot_posi(rate_measures(1))>0
            maxrateval=max(max(rates));
        end
        maxmeanrateval=max(mean_rates);
    else
        maxrateval=0;
        maxmeanrateval=0;
    end

    if any(f_para.subplot_posi(bi_measures))
        if f_para.norm_mode==1
            maxbivals=ones(1,length(bi_measures));
        else
            maxbivals=zeros(1,length(bi_measures));
            if f_para.subplot_posi(isi_measures(2))
                maxbivals(1)=max(abs(isi_ratio));
            end
            if f_para.subplot_posi(spike_pico_measures(2))
                maxbivals(2)=max(spike_diffs);
            end
            if f_para.subplot_posi(realtime_spike_pico_measures(2))
                maxbivals(3)=max(spike_diffs_realtime);
            end
            if f_para.subplot_posi(spike_time_measures(2))
                maxbivals(4)=max(spike_diffs_t);
            end
            if f_para.subplot_posi(realtime_spike_time_measures(2))
                maxbivals(5)=max(spike_diffs_realtime_t);
            end
            if f_para.subplot_posi(crealtime_spike_time_measures(2))
                maxbivals(6)=max(spike_diffs_crealtime_t);
            end
            if f_para.norm_mode==2
                maxbivals=max(maxbivals)*ones(1,length(bi_measures));
            end
        end
    else
        maxbivals=zeros(1,length(bi_measures));
    end


    subplot_numbers=zeros(1,s_para.num_subplots);
    singles=zeros(1,s_para.num_subplots);
    for suc=1:s_para.num_subplots
        singles(suc)=find(f_para.subplot_posi==suc,1,'first');                               % new subplot
        subplot_numbers(suc)=length(find(f_para.subplot_posi==suc));
    end
    doubles=intersect(setxor(1:length(f_para.subplot_posi),singles),find(f_para.subplot_posi>0));     % repeated subplot
    doublesref=singles(f_para.subplot_posi(doubles));                                                % corresponding new subplot

    relsubplot_size=f_para.subplot_size(1:find(f_para.subplot_size>0,1,'last'));
    if length(relsubplot_size)~=max(f_para.subplot_posi)
        relsubplot_size=ones(1,max(f_para.subplot_posi));
    end

    regplotsize=s_para.num_subplots*1.1;
    normsubplot_size=relsubplot_size/sum(relsubplot_size);
    dumsubplot_size=normsubplot_size*regplotsize;
    subplot_start2=cumsum(dumsubplot_size);
    subplot_size2=diff([0 subplot_start2]);

    f_para.subplot_size=zeros(1,num_all_measures);               % normalized size of subplots
    f_para.subplot_size(singles)=subplot_size2;
    f_para.subplot_size(doubles)=f_para.subplot_size(doublesref);

    subplot_start=zeros(1,num_all_measures);               % normalized position of subplots
    subplot_start(singles)=subplot_start2;
    subplot_start(doubles)=subplot_start(doublesref);

    subplot_index=[1 zeros(1,length(subplot_start)-1)];
    for spc=2:num_all_measures
        if f_para.subplot_posi(spc)>0
            subplot_index(spc)=length(find(f_para.subplot_posi(1:spc-1)==f_para.subplot_posi(spc)))+1;
        end
    end

    subplot_paras=[f_para.subplot_posi',f_para.subplot_size',subplot_start',subplot_index'];

    % #####################################################################

    figure(f_para.num_fig); clf; hold on;
    set(gcf,'Position',f_para.pos_fig)
    set(gcf,'Name',[f_para.filename,'  ',d_para.comment_string,'  ',f_para.comment_string,'  ',f_para.title_string])

    xlim([pmin pmax])
    ylim([0 regplotsize])
    s_para.xl=xlim; s_para.yl=ylim;
    for spc=1:num_all_measures
        if f_para.subplot_posi(spc)>0 && abs(subplot_start(spc)-regplotsize)>0.00001 && subplot_start(spc)~=0
            line(s_para.xl,(s_para.yl(2)-subplot_start(spc))*ones(1,2),'Color','k','LineStyle','-','LineWidth',2)
        end
    end
    yt=[]; ytl=[];

    if f_para.subplot_posi(1)>0                                                                          % Stimulus
        % here you can plot the stimulus. The example below (if uncommented) shows a sine wave.
        stim=0.5+sin((0:1/((s_para.itmax-s_para.itmin)/d_para.dts):1)*2*pi)/2;
        xvals=s_para.itmin:d_para.dts:s_para.itmax;
        if s_para.window_mode==1
            plot(xvals,s_para.yl(2)-subplot_start(1)+0.05+stim/1.1*f_para.subplot_size(1),['.-',s_para.colors(1)])
        else
            lxvals=intersect(find(xvals>=s_para.wmin),find(xvals<=s_para.wmax));
            plot(xvals,s_para.yl(2)-subplot_start(1)+0.05+stim/1.1*f_para.subplot_size(1),['.-',s_para.colors(2)])
            plot(xvals(lxvals),s_para.yl(2)-subplot_start(1)+0.05+stim(lxvals)/1.1*f_para.subplot_size(1),['.-',s_para.colors(1)])
        end
        max_stim_val=1;
        stim_lab=[-1 0 1];
        stim_pos=[0 0.5 1];
        yt=[yt s_para.yl(2)-subplot_start(1)+(0.05+stim_pos/max_stim_val)/1.1*f_para.subplot_size(1)];
        ytl=[ytl stim_lab];
        text(s_para.xl(1)-0.095*(s_para.xl(2)-s_para.xl(1)),s_para.yl(2)-subplot_start(1)+0.75/1.1*f_para.subplot_size(1),'Stimulus','Color','k','FontSize',font_size+1,'FontWeight','bold')
        line(s_para.xl,s_para.yl(2)-subplot_start(1)+0.05/1.1*f_para.subplot_size(1)*ones(1,2),'Color','k','LineStyle',':')
        line(s_para.xl,s_para.yl(2)-subplot_start(1)+1.05/1.1*f_para.subplot_size(1)*ones(1,2),'Color','k','LineStyle',':')
    end

    if f_para.subplot_posi(2)>0                                                                          % Spikes
        for trac=1:num_trains
            for sc=1:num_pspikes(trac)
                if pspikes(trac,sc)>=s_para.wmin && pspikes(trac,sc)<=s_para.wmax
                    line(pspikes(trac,sc)*ones(1,2),s_para.yl(2)-subplot_start(2)+(0.05+(num_trains-1-(trac-1)+[0.05 0.95])/num_trains)/1.1*f_para.subplot_size(2),...
                        'Color',s_para.colors(1))
                else
                    line(pspikes(trac,sc)*ones(1,2),s_para.yl(2)-subplot_start(2)+(0.05+(num_trains-1-(trac-1)+[0.05 0.95])/num_trains)/1.1*f_para.subplot_size(2),...
                        'Color',s_para.colors(2))
                end
            end
        end
        if isfield(d_para,'separators')
            for sec=1:length(d_para.separators)
                line([s_para.itmin s_para.itmax],(s_para.yl(2)-subplot_start(2)+f_para.subplot_size(2)-(0.05+d_para.separators(sec)/num_trains)/1.1*f_para.subplot_size(2))*ones(1,2),'Color','c','LineStyle','--','LineWidth',1)
            end
        end
        if isfield(d_para,'separators2')
            for sec=1:length(d_para.separators2)
                line([s_para.itmin s_para.itmax],(s_para.yl(2)-subplot_start(2)+f_para.subplot_size(2)-(0.05+d_para.separators2(sec)/num_trains)/1.1*f_para.subplot_size(2))*ones(1,2),'Color','b','LineStyle','-','LineWidth',2)
            end
        end
        if num_select_train_groups==1
            text(s_para.xl(1)-0.13*(s_para.xl(2)-s_para.xl(1)),s_para.yl(2)-subplot_start(2)+0.75/1.1*f_para.subplot_size(2),'Spike','Color','k','FontSize',s_para.font_size+1)
            text(s_para.xl(1)-0.13*(s_para.xl(2)-s_para.xl(1)),s_para.yl(2)-subplot_start(2)+0.35/1.1*f_para.subplot_size(2),'trains','Color','k','FontSize',s_para.font_size+1)
            spikesvals=[0.5 1];
            if mod(num_trains,2)==0
                spikeslab=spikesvals*num_trains;
            else
                spikeslab=fix(spikesvals*num_trains);
            end
            yt=[yt s_para.yl(2)-subplot_start(2)+(1.05-[spikeslab-0.5]/num_trains)/1.1*f_para.subplot_size(2)];
            ytl=[ytl spikeslab];
        else
            for sgc=1:num_select_train_groups
                if sgc<num_select_train_groups
                    line([s_para.itmin s_para.itmax],s_para.yl(2)-subplot_start(2)+(1.05-cum_num_select_group_trains(sgc)/num_trains/1.1*f_para.subplot_size(2))*ones(1,2),'Color','k','LineStyle',':')
                end
                text(s_para.xl(1)-0.1*(s_para.xl(2)-s_para.xl(1)),s_para.yl(2)-subplot_start(2)+f_para.subplot_size(2)-(0.05+select_group_center(sgc)/num_trains)/1.1*f_para.subplot_size(2),select_group_names{sgc},'Color','k','FontSize',f_para.font_size-1)
                %text(xl(1)-0.06*(xl(2)-xl(1)),1.05-select_group_center(sgc)/num_trains,select_group_names{sgc},'Color','k','FontSize',f_para.font_size,'FontWeight','bold')
            end
            yt=[yt fliplr(s_para.yl(2)-subplot_start(2)+(1.05-cum_num_select_group_trains/num_trains)/1.1*f_para.subplot_size(2))];
            ytl=[ytl fliplr(cum_num_select_group_trains)];
        end
        line(s_para.xl,s_para.yl(2)-subplot_start(2)+0.05/1.1*f_para.subplot_size(2)*ones(1,2),'Color','k','LineStyle',':')
        line(s_para.xl,s_para.yl(2)-subplot_start(2)+1.05/1.1*f_para.subplot_size(2)*ones(1,2),'Color','k','LineStyle',':')
    end
else
    f_para.subplot_size=zeros(1,num_all_measures);
    subplot_start=zeros(1,num_all_measures);
    subplot_index=f_para.subplot_posi;
    subplot_numbers=f_para.subplot_posi;
    subplot_paras=[f_para.subplot_posi',f_para.subplot_size',subplot_start',subplot_index'];
    maxintval=0; maxmeanintval=0; maxrateval=0; maxmeanrateval=0; maxbivals=zeros(1,length(bi_measures));
    yt=[]; ytl=[]; s_para.xl=[]; s_para.yl=[];
end

% ######################################################################################################################################
% ######################################################################################################################################
% ######################################################################################################################################

if num_trains==2
    ab_str='';
else
    ab_str='^a';
end

% 1           2         3     4        5     6       7     8         9   10       % Stimulus, Spike trains, Classic, Rate, ISI
% Stimulus    Spikes    PSTH  GPSTH    ISI  <ISI>    Rate  <Rate>    Ia  I^A
%
% 11  12     13    14       15   16      17     18        19      20                % SPIKE
% Sa  S^a    Sa_r  S_r^a    Sta  St^a    Sta_r  St_r^a    Sta_cr  St_cr^a
%
%     Output                        X-input         Y-input                    Label               Max value       Datatype
measure_paras={...
    {'';                            '';             '';                         '';                 0;              0};...    %  1  % Stimulus
    {'';                            '';             '';                         '';                 0;              0};...    %  2  % Spike train
    {'mean_psth';                   'bin_centers';  'psth';                     'PSTH';             maxpsth;        3};...    %  3  % Classic
    {'mean_psth';                   'bin_centers';  'gs_psth';                  'GPSTH';            maxgspsth;      3};...    %  4
    {'single_isi';                  'isi';          'ints';                     'ISI';              maxintval;      1};...    %  5  % Ints, Rate
    {'mean_isi';                    'isi';          'mean_ints';                '<ISI>';            maxmeanintval;  1};...    %  6
    {'single_rate';                 'isi';          'rates';                    'R';                maxrateval;     1};...    %  7
    {'mean_rate';                   'isi';          'mean_rates';               '<R>';              maxmeanrateval; 1};...    %  8
    {'xxx';                         'isi';          'abs(bi_isi_ratio)';        'Ia';               1;              1};...    %  9  % piecewise constant
    {'bi_isi_dist';                 'isi';          'isi_ratio';                ['I',ab_str];       maxbivals(1);   1};...    % 10

    {'xxx';                         'isi';          'bi_spike_diffs';           'Sa';               1;              1};...    % 11  % piecewise constant
    {'bi_spike_dist';               'isi';          'spike_diffs';              ['S',ab_str];       maxbivals(2);   1};...    % 12
    {'xxx';                         'isi';          'bi_spike_diffs_realtime';  'Sa_r';             1;              1};...    % 13
    {'bi_spike_dist_realtime';      'isi';          'spike_diffs_realtime';     ['S_r',ab_str];     maxbivals(3);   1};...    % 14
    {'xxx';                         'time';         'bi_spike_diffs_t';         'Sa';               1;              2};...    % 15  % continuous time-varying
    {'bi_spike_dist_t';             'time';         'spike_diffs_t';            ['S',ab_str];       maxbivals(4);   2};...    % 16
    {'xxx';                         'time';         'bi_spike_diffs_realtime_t';'Sa_r';             1;              2};...    % 17
    {'bi_spike_dist_realtime_t';    'time';         'spike_diffs_realtime_t';   ['S_r',ab_str];     maxbivals(5);   2};...    % 18
    {'xxx';                         'time';         'bi_spike_diffs_crealtime_t';'Sa_r';            1;              2};...    % 19
    {'bi_spike_dist_crealtime_t';   'time';         'spike_diffs_crealtime_t';  ['S_r',ab_str];     maxbivals(6);   2};...    % 20
    };
%     Output                        X-input         Y-input                 Label               Max value       Datatype

s_para.dtm=f_para.dtm; s_para.plot_mode=f_para.plot_mode;
subplot_tick=zeros(1,s_para.num_subplots);

if length(f_para.subplot_posi(f_para.subplot_posi>0))~=length(unique(f_para.subplot_posi(f_para.subplot_posi>0))) % same normalization in case of double subplots
    for supc=unique(f_para.subplot_posi(f_para.subplot_posi>0))
        if length(find(f_para.subplot_posi==supc))>1
            maxval=0;
            for supc2=find(f_para.subplot_posi==supc)
                maxval=max([maxval measure_paras{supc2}{5}]);
            end
            for supc2=find(f_para.subplot_posi==supc)
                measure_paras{supc2}{5}=maxval;
            end

        end
    end
end

for supc=3:num_all_measures
    if f_para.subplot_posi(supc)>0

        if ismember(supc,realtime_measures)
            s_para.causal=1;
        else
            s_para.causal=0;
        end
        s_para.line_width=max(subplot_index)+1-subplot_index(supc);
        %dcolors='kbrgmc'; s_para.colors=repmat(dcolors(subplot_index(supc)),1,6);

        subplot_tick(f_para.subplot_posi(supc))=subplot_tick(f_para.subplot_posi(supc))+1;
        if subplot_tick(f_para.subplot_posi(supc))<subplot_numbers(f_para.subplot_posi(supc))           % last trace (maybe only one)
            dyt=yt; dytl=ytl;
            eval(     [ '[ ' measure_paras{supc}{1} ', ' measure_paras{supc}{1} '_ave, mean_values, result_str, yt, ytl ] = f_measure_profiles( mean_values, result_str, yt, ytl, '...
                measure_paras{supc}{2} ', ' measure_paras{supc}{3} ', '''  ''', subplot_paras(' num2str(supc) ',:), ' ...
                num2str(measure_paras{supc}{5}) ', ' num2str(measure_paras{supc}{6}) ', s_para );' ]     );
            yt=dyt; ytl=dytl;
        else
            eval(     [ '[ ' measure_paras{supc}{1} ', ' measure_paras{supc}{1} '_ave, mean_values, result_str, yt, ytl ] = f_measure_profiles( mean_values, result_str, yt, ytl, '...
                measure_paras{supc}{2} ', ' measure_paras{supc}{3} ', ''' measure_paras{supc}{4} ''', subplot_paras(' num2str(supc) ',:), ' ...
                num2str(measure_paras{supc}{5}) ', ' num2str(measure_paras{supc}{6}) ', s_para );' ]     );
        end
        %         if supc==int_measures(1) && f_para.subplot_posi(int_measures(2))==f_para.subplot_posi(int_measures(1))
        %             maxmeanintval=maxintval;
        %         end
        %         if supc==rate_measures(1) && f_para.subplot_posi(rate_measures(2))==f_para.subplot_posi(rate_measures(1))
        %             maxmeanrateval=maxrateval;
        %         end

    end
end

if mod(s_para.plot_mode,2)>0

    if any(f_para.subplot_size(11:num_all_measures))                           % mark perfectly synchronous events (SPIKE-Distances only)
        common=pspikes(1,1:num_pspikes(1));
        for trc=2:num_trains
            spik2=pspikes(trc,1:num_pspikes(trc));
            common=intersect(common,spik2);
        end
        num_commons=length(common);
        for spc=11:num_all_measures
            if f_para.subplot_size(spc)>0
                for cc=1:num_commons
                    if common(cc)>=s_para.wmin && common(cc)<=s_para.wmax
                        line(common(cc)*ones(1,2),s_para.yl(2)-subplot_start(spc)+[0.05 1.05]/1.1*f_para.subplot_size(spc),...
                            'Color',s_para.colors(3),'LineStyle',':')
                    else
                        line(common(cc)*ones(1,2),s_para.yl(2)-subplot_start(spc)+[0.05 1.05]/1.1*f_para.subplot_size(spc),...
                            'Color',s_para.colors(4),'LineStyle',':')
                    end
                end
            end
        end
    end

    if s_para.window_mode>1
        line(s_para.wmin*ones(1,2),s_para.yl,'Color',s_para.colors(1),'LineStyle','-.')
        line(s_para.wmax*ones(1,2),s_para.yl,'Color',s_para.colors(1),'LineStyle','--')
    end
    line(s_para.itmin*ones(1,2),s_para.yl,'Color','k','LineStyle',':')
    line(s_para.itmax*ones(1,2),s_para.yl,'Color','k','LineStyle',':')
    line(d_para.tmin*ones(1,2),s_para.yl,'Color','k','LineStyle','-.')
    line(d_para.tmax*ones(1,2),s_para.yl,'Color','k','LineStyle','-.')
    xlabel(['Time  ',f_para.timeunit_string],'Color','k','FontSize',f_para.font_size+2)


    mark_extreme_spikes=0;
    if mark_extreme_spikes
        min_last_spike=[];
        for stc=1:num_trains
            if num_pspikes(stc)>0
                min_last_spike=min([min_last_spike spikes(stc,num_pspikes(stc))]);
            end
        end
        max_first_spike=max(spikes(:,1));
        line(max_first_spike*ones(1,2),[0.05 sum(f_para.subplot_size([1 3:num_all_measures]))-0.05],'Color','k','LineStyle','-.')
        line(min_last_spike*ones(1,2),[0.05 sum(f_para.subplot_size([1 3:num_all_measures]))-0.05],'Color','k','LineStyle','-.')
    end


    if f_para.publication==0
        title([f_para.title_string,'   ',f_para.filename,'   ',d_para.comment_string,'   ',f_para.comment_string,'   ',result_str],...
            'Color','k','FontSize',f_para.font_size+1,'FontWeight','bold')
    end
    [syt,syti]=sort(yt);
    sytl=ytl(syti);

    if strcmp(f_para.timeunit_string,'[min]')
        set(gca,'XTickLabel',get(gca,'XTick')/60000)
    elseif strcmp(f_para.timeunit_string,'[s]')
        set(gca,'XTickLabel',get(gca,'XTick')/1000)
    else
        set(gca,'XTickLabel',get(gca,'XTick'))
    end

    if isfield(f_para,'xfact')
        if f_para.xfact>1
            set(gca,'XTick',[0:f_para.xfact:s_para.itmax],'XTickLabel',0:fix(s_para.itmax/f_para.xfact))
        end
    end

    set(gca,'YTick',syt,'YTickLabel',sytl,'FontSize',f_para.font_size-1);
    set(gcf,'Color','w'); set(gca,'Color','w','Box','on');

    if isfield(d_para,'markers')
        for mac=1:length(d_para.markers)
            line(d_para.markers(mac)*ones(1,2),sum(f_para.subplot_size)-[0.05 1.05],'Color','c','LineStyle',':','LineWidth',1)
        end
    end
    if isfield(d_para,'markers2')
        for mac=1:length(d_para.markers2)
            line(d_para.markers2(mac)*ones(1,2),sum(f_para.subplot_size)-[0.05 1.05],'Color','r','LineStyle','-','LineWidth',2)
        end
    end
    if f_para.print_mode                                                                          % Create postscript file
        if f_para.publication==1
            set(gcf,'PaperOrientation','Portrait');set(gcf,'PaperType', 'A4');
            set(gcf,'PaperUnits','Normalized','PaperPosition',[0 0.45 1.0 0.5]);
        else
            set(gcf,'PaperOrientation','Landscape');set(gcf,'PaperType', 'A4');
            set(gcf,'PaperUnits','Normalized','PaperPosition', [0 0 1.0 1.0]);
        end
        psname=[f_para.imagespath,f_para.filename,d_para.comment_string,f_para.comment_string,'.ps']
        print(gcf,'-dpsc',psname);
    end
end

if f_para.saving
    save ([f_para.matpath,f_para.filename,d_para.comment_string,f_para.comment_string])
end

