function [spikes,d_para]=f_get_data_demo(dataset,d_para)

% dataset: 10-Clustering,30-splay,40-Frequency mismatch,50-Multi,60-Normalization,70-Poisson,80-Reliable/Random,90-Metric


if dataset==21                                % Bi: Frequency mismatch

    num_all_trains=2; num_spikes=12;
    spikes=zeros(num_all_trains,num_spikes);
    spikes(1,1:num_spikes)=(100:100:1200); %+rand(1,num_spikes)*50;
    spikes(2,1:num_spikes-1)=100:110:1200; %+rand(1,num_spikes-1)*50;
    d_para.dts=1; d_para.tmin=0; d_para.tmax=1300;

    spikes=spikes(:,1:find(spikes(1,:),1,'last'));
    spikes(spikes>d_para.tmax)=0;

    d_para.select_averages=[d_para.tmin d_para.tmax];
    d_para.trigger_averages=[];
    %d_para.trigger_averages{1}=[250 750];
    d_para.comment_string='Bi-Frequency-Mismatch';

elseif dataset==22                                % Multi: decreasing noise + 5 events

    d_para.dts=1; d_para.tmin=0; d_para.tmax=4000;
    num_trains=50; num_spikes=40;
    noise=[1:-1/(num_spikes/4-1):0 0];
    num_noises=length(noise);

    spikes=zeros(num_trains,num_spikes);
    spikes(1:num_trains,1:num_spikes/2)=sort(rand(num_trains,num_spikes/2),2)*d_para.tmax/2;
    num_events=5;
    for nc=1:num_events
        spikes(1:num_trains,ceil(num_spikes/2*nc/(num_events+1)))=nc*d_para.tmax/2/num_events+50*noise(ceil(num_noises-(nc-1)*num_noises/num_events)).*randn(1,num_trains)';
    end

    for trc=1:num_trains
        spikes(trc,num_spikes/4*3-1+(1:num_spikes/4+1))=100*(num_spikes/2-1)+200*(1:num_spikes/4+1)+50*noise.*randn(1,num_spikes/4+1);
    end

    spikes(spikes<d_para.tmin | spikes>d_para.tmax)=0;
    for trc=1:num_trains
        dummy=spikes(trc,1:num_spikes);
        spikes(trc,1:num_spikes)=[dummy(dummy>0) dummy(dummy==0)];
    end
    spikes=spikes(:,1:find(spikes(1,:),1,'last'));
    d_para.comment_string='Multi-Events';

elseif dataset==41                                % Non-spurious events

    d_para.dts=0.1; d_para.tmin=0; d_para.tmax=400;
    num_trains=50; num_spikes=7;

    thr1=100;

    balance=0;

    while balance~=1

        randy=rand(1,num_trains);
        spikes=zeros(num_trains,num_spikes);
        spikes(1:num_trains,1)=0.001;
        spikes(1:num_trains,2)=1+randy*80;
        spikes(1:num_trains,3)=repmat(thr1,num_trains,1)+10*(rand(num_trains,1)-0.5);
        spikes(1:num_trains,4)=200;
        spikes(1:num_trains,5)=201+randy*80;
        spikes(1:num_trains,num_spikes)=d_para.tmax; %-0.001;

        spikes(:,num_spikes-1)=spikes(:,3)+200;
        spikes((spikes(:,3)>thr1),num_spikes-1)=0;

        larger=sum(spikes(:,num_spikes-1)>0);
        balance=larger*2/num_trains;

        for trc=1:num_trains
            dummy=sort(spikes(trc,1:num_spikes));
            spikes(trc,1:num_spikes)=[dummy(dummy>0) dummy(dummy==0)];
        end
        if balance==1
            break
        end
    end
    d_para.comment_string='Non-spurious';

elseif dataset==51                                      % Clustering I

    d_para.dts=1; d_para.tmin=0; d_para.tmax=2000;
    num_trains=40; num_spikes=8;
    noise=[0.1 0.15 0.2 0.25 0.2 0.15 0.1] ;

    spikes=zeros(num_trains,num_spikes);

    for nc=1:num_spikes/4
        spikes(1:num_trains/2,nc)=(nc-0.5)/num_spikes*d_para.tmax+50*noise(1).*randn(1,num_trains/2)';
        spikes(num_trains/2+(1:num_trains/2),nc)=nc/num_spikes*d_para.tmax+50*noise(1).*randn(1,num_trains/2)';
    end

    for nc=num_spikes/4+(1:num_spikes/4)
        spikes(num_trains/4+(1:num_trains/2),nc)=(nc-0.5)/num_spikes*d_para.tmax+50*noise(2).*randn(1,num_trains/2)';
        spikes([1:num_trains/4 num_trains*3/4+(1:num_trains/4)],nc)=nc/num_spikes*d_para.tmax+50*noise(2).*randn(1,num_trains/2)';
    end

    for nc=num_spikes/2+(1:num_spikes/4)
        spikes([1:num_trains/4 num_trains/2+(1:num_trains/4)],nc)=(nc-0.5)/num_spikes*d_para.tmax+50*noise(3).*randn(1,num_trains/2)';
        spikes([num_trains/4+(1:num_trains/4) num_trains*3/4+(1:num_trains/4)],nc)=nc/num_spikes*d_para.tmax+50*noise(3).*randn(1,num_trains/2)';
    end

    rand_st=randperm(num_trains);
    for nc=num_spikes*3/4+(1:num_spikes/4)
        spikes(rand_st(1:num_trains/2),nc)=(nc-0.5)/num_spikes*d_para.tmax+50*noise(4).*randn(1,num_trains/2)';
        spikes(rand_st(num_trains/2+(1:num_trains/2)),nc)=nc/num_spikes*d_para.tmax+50*noise(4).*randn(1,num_trains/2)';
    end

    spikes(spikes>0)=spikes(spikes>0)-60;

    spikes(spikes<d_para.tmin | spikes>d_para.tmax)=0;
    for trc=1:num_trains
        dummy=spikes(trc,1:num_spikes);
        spikes(trc,1:num_spikes)=[dummy(dummy>0) dummy(dummy==0)];
    end

    spikes(10,3)=spikes(10,3)-15;
    spikes(4,3)=spikes(4,3)+25;
    spikes(4,4)=spikes(4,4)+30;

    spikes(37,3)=spikes(37,3)-20;
    spikes(39,3)=spikes(39,3)+15;

    d_para.interval_separators=500:500:d_para.tmax-500; % Edges of subsections
    d_para.interval_strings={'2 Cluster - AABB';'2 Cluster - ABBA';'2 Cluster - ABAB';'2 Cluster - Random association'}; % Captions for subsections
    d_para.markers=d_para.interval_separators; % Edges of subsections
    d_para.comment_string='Clustering-1';

elseif dataset==61                                      % Clustering II

    d_para.dts=1; d_para.tmin=0; d_para.tmax=2000;
    num_trains=40; num_spikes=8;
    noise=[0.1 0.15 0.2 0.25 0.2 0.15 0.1] ;

    spikes=zeros(num_trains,num_spikes);

    for nc=1:num_spikes/4
        spikes(1:num_trains/4,nc)=(nc-0.25)/num_spikes*d_para.tmax+50*noise(5).*randn(1,num_trains/4)'-60;
        spikes(num_trains/4+(1:num_trains/2),nc)=(nc-0.5)/num_spikes*d_para.tmax+50*noise(5).*randn(1,num_trains/2)'-60;
        spikes(num_trains*3/4+(1:num_trains/4),nc)=nc/num_spikes*d_para.tmax+50*noise(5).*randn(1,num_trains/4)'-60;
    end

    for nc=num_spikes/4+(1:num_spikes/4)
        spikes(1:num_trains/4,nc)=nc/num_spikes*d_para.tmax+50*noise(6).*randn(1,num_trains/4)'-30;
        spikes(num_trains/4+(1:num_trains/4),nc)=(nc-0.25)/num_spikes*d_para.tmax+50*noise(6).*randn(1,num_trains/4)'-30;
        spikes(num_trains/2+(1:num_trains/4),nc)=(nc-0.5)/num_spikes*d_para.tmax+50*noise(6).*randn(1,num_trains/4)'-30;
        spikes(num_trains*3/4+(1:num_trains/4),nc)=(nc-0.75)/num_spikes*d_para.tmax+50*noise(6).*randn(1,num_trains/4)'-30;
    end

    for nc=num_spikes/2+(1:num_spikes/4)
        spikes(1:num_trains/8,nc)=(nc-0.11)/num_spikes*d_para.tmax+50*noise(7)/2.*randn(1,num_trains/8)';
        spikes(num_trains/8+(1:num_trains/8),nc)=(nc-0.22)/num_spikes*d_para.tmax+50*noise(7)/2.*randn(1,num_trains/8)';
        spikes(num_trains/4+(1:num_trains/8),nc)=(nc-0.33)/num_spikes*d_para.tmax+50*noise(7)/2.*randn(1,num_trains/8)';
        spikes(num_trains*3/8+(1:num_trains/8),nc)=(nc-0.44)/num_spikes*d_para.tmax+50*noise(7)/2.*randn(1,num_trains/8)';
        spikes(num_trains/2+(1:num_trains/8),nc)=(nc-0.55)/num_spikes*d_para.tmax+50*noise(7)/2.*randn(1,num_trains/8)';
        spikes(num_trains*5/8+(1:num_trains/8),nc)=(nc-0.66)/num_spikes*d_para.tmax+50*noise(7)/2.*randn(1,num_trains/8)';
        spikes(num_trains*3/4+(1:num_trains/8),nc)=(nc-0.77)/num_spikes*d_para.tmax+50*noise(7)/2.*randn(1,num_trains/8)';
        spikes(num_trains*7/8+(1:num_trains/8),nc)=(nc-0.88)/num_spikes*d_para.tmax+50*noise(7)/2.*randn(1,num_trains/8)';
    end

    for nc=num_spikes*3/4+(1:num_spikes/4)
        spikes(1:num_trains,nc)=nc/num_spikes*d_para.tmax-d_para.tmax/num_spikes.*rand(1,num_trains)';
    end

    %spikes(spikes>0)=spikes(spikes>0)-60;

    spikes(spikes<d_para.tmin | spikes>d_para.tmax)=0;
    for trc=1:num_trains
        dummy=spikes(trc,1:num_spikes);
        spikes(trc,1:num_spikes)=[dummy(dummy>0) dummy(dummy==0)];
    end
    d_para.interval_separators=500:500:d_para.tmax-500; % Edges of subsections
    d_para.interval_strings={'3 Cluster - ABBC';'4 Cluster - ABCD';'8 Cluster - ABCDEFGH';'Random Spiking'}; % Captions for subsections
    d_para.markers=d_para.interval_separators; % Edges of subsections
    d_para.comment_string='Clustering-2';

elseif dataset==71                                      % Clustering III (All)

    d_para.tmin=0; d_para.tmax=4000;
    d_para.dts=1;

    d_para.all_train_group_names={'G1';'G2';'G3';'G4'};
    d_para.all_train_group_sizes=[10 10 10 10];

    num_all_trains=40; num_spikes=16;
    noise=[0.1 0.15 0.2 0.25 0.2 0.15 0.1 0.1];

    spikes=zeros(num_all_trains,num_spikes);
    d_para.interval_str{1}='2 Cluster - AABB';
    for nc=1:num_spikes/8
        spikes(1:num_all_trains/2,nc)=(nc-0.5)/num_spikes*d_para.tmax+50*noise(1).*randn(1,num_all_trains/2)';
        spikes(num_all_trains/2+(1:num_all_trains/2),nc)=nc/num_spikes*d_para.tmax+50*noise(1).*randn(1,num_all_trains/2)';
    end

    d_para.interval_str{2}='2 Cluster - ABBA';
    for nc=num_spikes/8+(1:num_spikes/8)
        spikes(num_all_trains/4+(1:num_all_trains/2),nc)=(nc-0.5)/num_spikes*d_para.tmax+50*noise(2).*randn(1,num_all_trains/2)';
        spikes([1:num_all_trains/4 num_all_trains*3/4+(1:num_all_trains/4)],nc)=nc/num_spikes*d_para.tmax+50*noise(2).*randn(1,num_all_trains/2)';
    end

    d_para.interval_str{3}='2 Cluster - ABAB';
    for nc=num_spikes/4+(1:num_spikes/8)
        spikes([1:num_all_trains/4 num_all_trains/2+(1:num_all_trains/4)],nc)=(nc-0.5)/num_spikes*d_para.tmax+50*noise(3).*randn(1,num_all_trains/2)';
        spikes([num_all_trains/4+(1:num_all_trains/4) num_all_trains*3/4+(1:num_all_trains/4)],nc)=nc/num_spikes*d_para.tmax+50*noise(3).*randn(1,num_all_trains/2)';
    end

    d_para.interval_str{4}='2 Cluster - Random association';
    rand_st=randperm(num_all_trains);
    for nc=num_spikes*3/8+(1:num_spikes/8)
        spikes(rand_st(1:num_all_trains/2),nc)=(nc-0.5)/num_spikes*d_para.tmax+50*noise(4).*randn(1,num_all_trains/2)';
        spikes(rand_st(num_all_trains/2+(1:num_all_trains/2)),nc)=nc/num_spikes*d_para.tmax+50*noise(4).*randn(1,num_all_trains/2)';
    end

    d_para.interval_str{5}='3 Cluster - ABBC';
    for nc=num_spikes/2+(1:num_spikes/8)
        spikes(1:num_all_trains/4,nc)=(nc-0.25)/num_spikes*d_para.tmax+50*noise(5).*randn(1,num_all_trains/4)';
        spikes(num_all_trains/4+(1:num_all_trains/2),nc)=(nc-0.5)/num_spikes*d_para.tmax+50*noise(5).*randn(1,num_all_trains/2)';
        spikes(num_all_trains*3/4+(1:num_all_trains/4),nc)=nc/num_spikes*d_para.tmax+50*noise(5).*randn(1,num_all_trains/4)';
    end
    spikes(spikes>0)=spikes(spikes>0)-60;

    d_para.interval_str{6}='4 Cluster - ABCD';
    for nc=num_spikes*5/8+(1:num_spikes/8)
        spikes(1:num_all_trains/4,nc)=nc/num_spikes*d_para.tmax+50*noise(6).*randn(1,num_all_trains/4)'-30;
        spikes(num_all_trains/4+(1:num_all_trains/4),nc)=(nc-0.25)/num_spikes*d_para.tmax+50*noise(6).*randn(1,num_all_trains/4)'-30;
        spikes(num_all_trains/2+(1:num_all_trains/4),nc)=(nc-0.5)/num_spikes*d_para.tmax+50*noise(6).*randn(1,num_all_trains/4)'-30;
        spikes(num_all_trains*3/4+(1:num_all_trains/4),nc)=(nc-0.75)/num_spikes*d_para.tmax+50*noise(6).*randn(1,num_all_trains/4)'-30;
    end

    d_para.interval_str{7}='8 Cluster - ABCDEFGH';
    for nc=num_spikes*6/8+(1:num_spikes/8)
        spikes(1:num_all_trains/8,nc)=(nc-0.11)/num_spikes*d_para.tmax+50*noise(7)/2.*randn(1,num_all_trains/8)';
        spikes(num_all_trains/8+(1:num_all_trains/8),nc)=(nc-0.22)/num_spikes*d_para.tmax+50*noise(7)/2.*randn(1,num_all_trains/8)';
        spikes(num_all_trains/4+(1:num_all_trains/8),nc)=(nc-0.33)/num_spikes*d_para.tmax+50*noise(7)/2.*randn(1,num_all_trains/8)';
        spikes(num_all_trains*3/8+(1:num_all_trains/8),nc)=(nc-0.44)/num_spikes*d_para.tmax+50*noise(7)/2.*randn(1,num_all_trains/8)';
        spikes(num_all_trains/2+(1:num_all_trains/8),nc)=(nc-0.55)/num_spikes*d_para.tmax+50*noise(7)/2.*randn(1,num_all_trains/8)';
        spikes(num_all_trains*5/8+(1:num_all_trains/8),nc)=(nc-0.66)/num_spikes*d_para.tmax+50*noise(7)/2.*randn(1,num_all_trains/8)';
        spikes(num_all_trains*3/4+(1:num_all_trains/8),nc)=(nc-0.77)/num_spikes*d_para.tmax+50*noise(7)/2.*randn(1,num_all_trains/8)';
        spikes(num_all_trains*7/8+(1:num_all_trains/8),nc)=(nc-0.88)/num_spikes*d_para.tmax+50*noise(7)/2.*randn(1,num_all_trains/8)';
    end

    d_para.interval_str{8}='Random Spiking';
    for nc=num_spikes*7/8+(1:num_spikes/8)
        spikes(1:num_all_trains,nc)=nc/num_spikes*d_para.tmax-d_para.tmax/num_spikes.*rand(1,num_all_trains)';
    end

    d_para.select_train_mode=1;    % 1:all,2:selected groups,3:selected trains
    d_para.select_averages=[];
    d_para.trigger_averages=[];
    d_para.markers=500:500:d_para.tmax-500;
    d_para.interval_separators=500:500:d_para.tmax-500; % Edges of subsections
    d_para.comment_string='Clustering';

elseif dataset==81                                % Poisson (Divergence)

    num_all_spikes=200; num_trig_trac1_spikes=12;
    %num_all_spikes=10; num_trig_trac1_spikes=3;

    num_trains=20;
    trig_trac1=1; trig_tracs1=[4 8 11 16 19];
    d_para.tmin=0; d_para.tmax=100; %min(spikes(:,num_all_spikes))*1.0001;
    d_para.dts=0.001;

    spikes=zeros(num_trains,num_all_spikes);
    for trc=setdiff(1:num_trains,trig_trac1)
        dummy=f_poisson(num_all_spikes,1,0);
        dummy(dummy<d_para.dts)=d_para.dts;
        spikes(trc,1:num_all_spikes)=cumsum(dummy);
    end
    dummy=f_poisson(num_trig_trac1_spikes,1,0);
    dummy(dummy<d_para.dts)=d_para.dts;
    spikes(trig_trac1,1:num_trig_trac1_spikes)=cumsum(dummy);

    spikes(spikes<d_para.tmin | spikes>d_para.tmax)=0;

    num_spikes=zeros(1,num_trains);
    for trc=1:num_trains
        num_spikes(trc)=find(spikes(trc,1:num_all_spikes),1,'last');
    end

    spikes(trig_trac1,1:num_trig_trac1_spikes)=spikes(trig_trac1,1:num_trig_trac1_spikes)/max(spikes(trig_trac1,1:num_trig_trac1_spikes))*d_para.tmax*0.97;
    for trc=setdiff(1:num_trains,trig_trac1)
        spikes(trc,1:num_spikes(trc))=spikes(trc,1:num_spikes(trc))/max(spikes(trc,1:num_spikes(trc)))*d_para.tmax*0.995;
        spikes(trc,num_spikes(trc))=spikes(trc,num_spikes(trc))-rand*5;
        spikes(trc,1:num_spikes(trc))=sort(spikes(trc,1:num_spikes(trc)));
    end

    max_num_spikes=max(num_spikes);
    spikes=spikes(:,1:max_num_spikes);

    for trc=trig_tracs1
        indy=[];
        for spic=1:num_spikes(trig_trac1)
            rem_indy=setdiff(1:max_num_spikes,indy);
            [dummy,index]=min(abs(spikes(trc,rem_indy)-spikes(trig_trac1,spic)));
            indy=[indy rem_indy(index)];
            spikes(trc,rem_indy(index))=spikes(trig_trac1,spic)+0.05*rand;
        end
        spikes(trc,1:num_spikes(trc))=sort(spikes(trc,1:num_spikes(trc)));
    end
    spikes(spikes<d_para.tmin | spikes>d_para.tmax)=0;

    num_spikes2=zeros(1,num_trains);
    for trc=1:num_trains
        num_spikes2(trc)=find(spikes(trc,1:max_num_spikes),1,'last');
    end
    d_para.comment_string='Poisson-Driver';

elseif dataset==111                                % 2011 paper, Fig. 2a   (splay state vs. identical)

    num_all_trains=20; num_spikes=10;
    spikes=zeros(num_all_trains,num_spikes);
    spikes(1,1:num_spikes)=0:100:(num_spikes-1)*100;
    for trc=2:num_all_trains
        spikes(trc,1)=spikes(1,1);
        spikes(trc,2:5)=spikes(1,1:4)+(trc-1)*100/num_all_trains;
        spikes(trc,6:num_spikes)=spikes(1,5:num_spikes-1);
    end
    spikes=sort(spikes,2);

    d_para.dts=1; d_para.tmin=0; d_para.tmax=800;
    spikes=spikes(:,1:find(spikes(1,:),1,'last'));
    spikes(spikes>d_para.tmax)=0;
    d_para.comment_string='Splay';

elseif dataset==121                                       % Poisson Test - Scalefree

    %rate=[10 31 100 310 1000 3100 10000];
    rate=round(10.^(1.8:0.1:3.3));

    num_trains=length(rate);
    max_num_spikes=max(rate);
    spikes=zeros(num_trains,max_num_spikes);
    max_spikes=zeros(1,num_trains);
    for trc=1:num_trains
        spikes(trc,1:rate(trc))=cumsum(f_poisson(rate(trc),rate(trc),0));
        max_spikes(trc)=max(spikes(trc,1:rate(trc)));
    end

    d_para.tmin=0; d_para.tmax=min(max_spikes)*1.0001; d_para.dts=0.00001;
    spikes(spikes<d_para.tmin | spikes>d_para.tmax)=0;

    num_spikes=zeros(1,num_trains);
    for trc=1:num_trains
        num_spikes(trc)=find(spikes(trc,1:max_num_spikes),1,'last');
    end
    d_para.comment_string='Poisson-Scalefree';

elseif dataset==122                                       % Poisson Test - Scalefree 2

    rate=round(10.^[1.9*ones(1,12) 2:0.1:3.3]);

    num_trains=length(rate);
    max_num_spikes=max(rate);
    spikes=zeros(num_trains,max_num_spikes);
    max_spikes=zeros(1,num_trains);
    for trc=1:num_trains
        spikes(trc,1:rate(trc))=cumsum(f_poisson(rate(trc),rate(trc),0));
        max_spikes(trc)=max(spikes(trc,1:rate(trc)));
    end

    d_para.tmin=0; d_para.tmax=min(max_spikes)*1.0001; d_para.dts=0.00001;
    spikes(spikes<d_para.tmin | spikes>d_para.tmax)=0;

    num_spikes=zeros(1,num_trains);
    for trc=1:num_trains
        num_spikes(trc)=find(spikes(trc,1:max_num_spikes),1,'last');
    end

    [snum_spikes,sindex]=sort(num_spikes(1:12));
    dspikes1=spikes(sindex(11),1:num_spikes(sindex(11)));
    dspikes2=spikes(sindex(12),1:num_spikes(sindex(11)));
    dspikes2=dspikes2/max(dspikes2)*0.99*d_para.tmax;

    spikes(1,1:num_spikes(sindex(11)))=dspikes1(1:num_spikes(sindex(11)));
    for trc=2:12
        spikes(trc,:)=0;
        spikes(trc,1:num_spikes(sindex(11)))=(12-trc)*0.1*dspikes1(1:num_spikes(sindex(11)))+(trc-2)*0.1*dspikes2(1:num_spikes(sindex(11)));
    end

    num_spikes=zeros(1,num_trains);
    for trc=1:num_trains
        num_spikes(trc)=find(spikes(trc,:),1,'last');
    end

    d_para.comment_string='Poisson-Scalefree2';

elseif dataset==131                                       % Numerical calculation of expectation values for Poisson spike trains

    num_trains=100; num_spikes=1000;
    %num_trains=5; num_spikes=100;
    spikes=zeros(num_trains,num_spikes);
    for trc=1:num_trains
        spikes(trc,1:num_spikes)=cumsum(f_poisson(num_spikes,1,0));
    end

    d_para.tmin=0; d_para.tmax=min(spikes(:,num_spikes))*1.0001; d_para.dts=0.001;
    spikes(spikes<d_para.tmin | spikes>d_para.tmax)=0;
    d_para.comment_string='Poisson-ExpectationValue';

elseif dataset==141                                       % Reliable

    num_trains=10; num_spikes=9; dura=100;
    spikes=zeros(num_trains,num_spikes);
    for trc=1:num_trains
        spikes(trc,1:num_spikes)=((1:num_spikes)+0.2*rand(1,num_spikes))*dura;
    end
    d_para.tmin=0; d_para.tmax=(num_spikes+1)*dura; d_para.dts=1;
    spikes=round(spikes/d_para.dts)*d_para.dts;
    spikes(spikes<d_para.tmin | spikes>d_para.tmax)=0;
    d_para.comment_string='Reliable';

elseif dataset==142                                       % Random

    num_trains=10; num_spikes=9; dura=100;
    spikes=zeros(num_trains,num_spikes);
    for trc=1:num_trains
        spikes(trc,1:num_spikes)=((1:num_spikes)+0.2*rand(1,num_spikes))*dura;
    end
    d_para.tmin=0; d_para.tmax=(num_spikes+1)*dura; d_para.dts=1;
    spikes=round(spikes/d_para.dts)*d_para.dts;
    spikes(spikes<d_para.tmin | spikes>d_para.tmax)=0;
    all_spikes=reshape(spikes,1,num_trains*num_spikes);
    all_spikes2=all_spikes(randperm(num_trains*num_spikes));
    spikes=sort(reshape(all_spikes2,num_trains,num_spikes),2);
    d_para.comment_string='Random';

end
