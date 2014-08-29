function [ave,ave2,fave,fresult_str,yt,ytl]=f_measure_profiles(fave,fresult_str,fyt,fytl,fx,fy,fheadline,fsp_paras,fmaxval,fdata_type,s_para)

fsp_posi=fsp_paras(1);
fsp_size=fsp_paras(2);
fsp_start=fsp_paras(3);
%fsp_index=fsp_paras(4);

if s_para.window_mode==1
    dfcols=[s_para.colors([3 3])];    % individual traces (within window, outside/without window)
    dmfcols=[s_para.colors([1 2])];   % multiple traces (within window, outside/without window)
else
    dfcols=[s_para.colors([3 4])];
    dmfcols=[s_para.colors([1 2])];
end

for nmac=1:length(s_para.nma)

    if s_para.nma(nmac)==2 && size(fy,1)==1                                                       % moving average
        if s_para.causal==0
            if round(sum(fx)/s_para.dtm)*s_para.dtm==round((s_para.itmax-s_para.itmin)/s_para.dtm)*s_para.dtm
                for trac=1:size(fy,1)
                    fy(trac,:)=f_moving_average_weighted(fy(trac,:),fx,s_para.spike_mao);
                end
            else
                if size(fy,1)==1
                    fy=f_moving_average(fy,s_para.time_mao);
                else
                    fy=f_moving_average_para(fy,s_para.time_mao);
                end
            end
        else
            if round(sum(fx)/s_para.dtm)*s_para.dtm==round((s_para.itmax-s_para.itmin)/s_para.dtm)*s_para.dtm
                for trac=1:size(fy,1)
                    fy(trac,:)=f_moving_average_weighted_p(fy(trac,:),fx,s_para.spike_mao);
                end
            else
                if size(fy,1)==1
                    fy=f_moving_average_p(fy,s_para.time_mao);
                else
                    fy=f_moving_average_para_p(fy,s_para.time_mao);
                end
            end
        end
        %fheadline=[fheadline,'*'];
        s_para.line_width=s_para.line_width+1;
        if s_para.window_mode==1
            dfcols=[s_para.colors([5 5])];
        else
            dfcols=[s_para.colors([5 6])];
        end
    end

    if fdata_type==1                                        % using vectors of length "num_isi"
        ctfx=cumsum([0 fx]);
        pfx=s_para.itmin+unique([ctfx(1:end) ctfx(2:end-1)-s_para.dtm/100]);

        if size(fy,1)==1
            pfy=reshape([fy; fy],1,length(fy)*2);
            if s_para.window_mode==1
                
                ave=sum(fy.*repmat(fx,size(fy,1),1),2)/sum(fx);
                ave2=mean(fy);
            end
        end
    else                                                        % using vectors of length "len"
        pfx=fx;
        pfy=fy;
        if size(fy,1)==1 && s_para.window_mode==1
            ave=mean(fy,2);
            ave2=ave;
        end
    end

    if size(fy,1)==1                                                            % one single profile
        if s_para.window_mode>1
            first_winspike_ind=find(pfx>=s_para.wmin,1,'first');
            last_winspike_ind=find(pfx<=s_para.wmax,1,'last');
            win_pfx=pfx(first_winspike_ind:last_winspike_ind);
            win_pfy=pfy(first_winspike_ind:last_winspike_ind);
            if fdata_type<3
                if s_para.wmin<pfx(first_winspike_ind) && s_para.wmin>s_para.itmin     % interval to first spike
                    win_pfx=[s_para.wmin win_pfx];
                    win_pfy=[pfy(first_winspike_ind-1) win_pfy];
                end
                if s_para.wmax>pfx(last_winspike_ind) && s_para.wmax<s_para.itmax      % interval after last spike
                    win_pfx=[win_pfx s_para.wmax];
                    win_pfy=[win_pfy pfy(last_winspike_ind+1)];
                end
            end

            win_pfx_isi=diff(win_pfx);
            win_pfy_isi=win_pfy(1:end-1);

            if fdata_type==1                                        % using vectors of length "num_isi"
                ave=sum(win_pfy_isi.*repmat(win_pfx_isi,size(win_pfy_isi,1),1),2)/sum(win_pfx_isi);
                ave2=mean(win_pfy_isi);
            else                                                    % using vectors of length "len"
                ave=mean(win_pfy);
                ave2=ave;
            end
        end
    else
        ave=0;
        ave2=ave;
    end


    if mod(s_para.plot_mode,2)>0 && fsp_posi>0                                                              % plotting
        num_vals=2;
        if s_para.num_subplots<6
            intvals=[0.5 1];
        elseif s_para.num_subplots<10
            intvals=[0.4 0.8];
        else
            intvals=[0.3 0.6];
        end
        intlab=unique([0 f_lab(intvals*fmaxval,num_vals,1,1)]);
        %intlab=[0 0.5];
        yt=[fyt s_para.yl(2)-fsp_start+(0.05+intlab/fmaxval)/1.1*fsp_size];
        ytl=[fytl intlab];
        text(s_para.xl(1)-0.11*(s_para.xl(2)-s_para.xl(1)),s_para.yl(2)-fsp_start+0.6/1.1*fsp_size,fheadline,'Color','k','FontSize',s_para.font_size+2)
        line(s_para.xl,s_para.yl(2)-fsp_start+0.05/1.1*fsp_size*ones(1,2),'Color','k','LineStyle',':')
        line(s_para.xl,s_para.yl(2)-fsp_start+1.05/1.1*fsp_size*ones(1,2),'Color','k','LineStyle',':')

        if size(fy,1)==1                                                                % one single profile
            plot(pfx,s_para.yl(2)-fsp_start+(0.05+pfy/fmaxval)/1.1*fsp_size,[s_para.line_style,dfcols(2)],'LineWidth',s_para.line_width)
            %line([s_para.itmin s_para.itmax],s_para.yl(2)-fsp_start+(0.05+ave/fmaxval)/1.1*fsp_size*ones(1,2),'LineStyle','--','Color','g','LineWidth',s_para.line_width+1)
            %text(s_para.xl(2)-0.15*(s_para.xl(2)-s_para.xl(1)),s_para.yl(2)-fsp_start+0.93/1.1*fsp_size,num2str(ave,3),'Color','g','FontSize',s_para.font_size+1,'FontWeight','bold')
            if s_para.window_mode>1
                plot(win_pfx,s_para.yl(2)-fsp_start+(0.05+win_pfy/fmaxval)/1.1*fsp_size,[s_para.line_style,dfcols(1)],'LineWidth',s_para.line_width)
                line([s_para.wmin s_para.wmax],s_para.yl(2)-fsp_start+(0.05+ave/fmaxval)/1.1*fsp_size*ones(1,2),'LineStyle','--','Color',dfcols(1),'LineWidth',s_para.line_width)
            end
            %if nmac==length(s_para.nma)
            %    text(s_para.xl(2)-0.1*(s_para.xl(2)-s_para.xl(1)),s_para.yl(2)-fsp_start+0.5/1.1*fsp_size,num2str(ave,3),'Color',dfcols(1),'FontSize',s_para.font_size+1,'FontWeight','bold')
            %end

        elseif nmac<2                                                                   % ISIs and rates for individual trials

            if fdata_type==1
                if s_para.window_mode>1
                    first_winspike_ind=find(pfx>=s_para.wmin,1,'first');
                    last_winspike_ind=find(pfx<=s_para.wmax,1,'last');
                    win_pfx=pfx(first_winspike_ind:last_winspike_ind);
                    if fdata_type<3
                        if s_para.wmin<pfx(first_winspike_ind)
                            win_pfx=[s_para.wmin win_pfx];
                        end
                        if s_para.wmax>pfx(last_winspike_ind)
                            win_pfx=[win_pfx s_para.wmax];
                        end
                    end
                end
                for trac=1:size(fy,1)
                    fyy=fy(trac,:);
                    pfyy=reshape([fyy; fyy],1,length(fyy)*2);
                    plot(pfx,s_para.yl(2)-fsp_start+(0.05+pfyy/fmaxval)/1.1*fsp_size,[s_para.line_style,dmfcols(2)],'LineWidth',s_para.line_width)
                    if s_para.window_mode>1
                        win_pfyy=pfyy(first_winspike_ind:last_winspike_ind);
                        if s_para.wmin<pfx(first_winspike_ind) && first_winspike_ind>1
                            win_pfyy=[pfyy(first_winspike_ind-1) win_pfyy];
                        end
                        if s_para.wmax>pfx(last_winspike_ind) && length(pfyy)>last_winspike_ind
                            win_pfyy=[win_pfyy pfyy(last_winspike_ind+1)];
                        end
                        plot(win_pfx,s_para.yl(2)-fsp_start+(0.05+win_pfyy/fmaxval)/1.1*fsp_size,[s_para.line_style,dmfcols(1)],'LineWidth',s_para.line_width)
                    end
                end
            else
                if s_para.window_mode>1
                    first_winspike_ind=find(pfx>=s_para.wmin,1,'first');
                    last_winspike_ind=find(pfx<=s_para.wmax,1,'last');
                    win_pfx=pfx(first_winspike_ind:last_winspike_ind);
                    win_pfy=pfy(:,first_winspike_ind:last_winspike_ind);
                    if fdata_type<3
                        if s_para.wmin<pfx(first_winspike_ind)     % interval to first spike
                            win_pfx=[s_para.wmin win_pfx];
                            win_pfy=[pfy(first_winspike_ind-1) win_pfy];
                        end
                        if s_para.wmax>pfx(last_winspike_ind)      % interval after last spike
                            win_pfx=[win_pfx s_para.wmax];
                            win_pfy=[win_pfy pfy(last_winspike_ind+1)];
                        end
                    end
                end
                for trac=1:size(fy,1)
                    pfyy=fy(trac,:);
                    plot(pfx,s_para.yl(2)-fsp_start+(0.05+pfyy/fmaxval)/1.1*fsp_size,[s_para.line_style,dmfcols(2)],'LineWidth',s_para.line_width)
                    if s_para.window_mode>1
                        win_pfyy=win_pfy(trac,:);
                        plot(win_pfx,s_para.yl(2)-fsp_start+(0.05+win_pfyy/fmaxval)/1.1*fsp_size,[s_para.line_style,dmfcols(1)],'LineWidth',s_para.line_width)
                    end
                end
            end
        end
    else
        yt=[];
        ytl=[];
    end
end
if size(fy,1)==1
    fave=[fave ave];
    fresult_str=[fresult_str,'  ',fheadline,' = ',num2str(ave,7),','];
end

