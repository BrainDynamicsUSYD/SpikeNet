function [spikes,d_para]=f_get_data(dataset,d_para)

if dataset==1                                % Replace by your data, load files, set data parameters, ...

    num_all_trains=2; num_spikes=12;
    spikes=zeros(num_all_trains,num_spikes);
    spikes(1,1:num_spikes)=(100:100:1200); %+rand(1,num_spikes)*50;
    spikes(2,1:num_spikes-1)=[100:110:1200]; %+rand(1,num_spikes-1)*50;
    d_para.dts=1; d_para.tmin=0; d_para.tmax=1300;

    spikes=spikes(:,1:find(spikes(1,:),1,'last'));
    spikes(spikes>d_para.tmax)=0;

    d_para.select_averages=[d_para.tmin d_para.tmax];
    d_para.trigger_averages=[];
    %d_para.trigger_averages{1}=[250 750];
    d_para.comment_string='Bi-Frequency-Mismatch';

end
