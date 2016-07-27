if mod(f_para.plot_mode,32)>15 && ruc==1
    moviename=[f_para.moviespath,f_para.filename,d_para.comment_string,f_para.comment_string,'.avi'];
    tmov = avifile(moviename,'fps',f_para.mov_frames_per_second);   % ,'compression','none'
    if f_para.dendrograms && f_para.dendro_color_mode>0
        if num_trains>50
            dlw=1;
        elseif num_trains>10
            dlw=1.5;
        else
            dlw=2;
        end
    end
end



if num_frame_select+(num_select_averages+num_trigger_averages)*(1+(f_para.mov_num_average_frames-1)*(mod(f_para.plot_mode,32)>15))*(ruc==num_runs)>0

    if ruc==1 || (ruc==num_runs && num_frame_select==0)
        if num_all_subplots>8
            error('Please reduce the number of subplots to 8 !')
        end

        max_cols=[1 2 3 4 3 3 4 4];

        max_col=max_cols(num_all_subplots);

        rows=ceil((1:num_all_subplots)/max_col);
        cols=(1:num_all_subplots)-max_col*((1:num_all_subplots)>max_col);
        supos=zeros(num_all_subplots,4);
        for matc=1:num_all_subplots
            if num_all_subplots==1
                supos(matc,1:4)=[0.35  0.075  0.33   0.4];
            elseif num_all_subplots==2
                if f_para.publication
                    supos(matc,1:4)=[0.11+0.45*(cols(matc)-1)  0.075  0.35+0.0439*(matc==num_all_subplots)*(f_para.color_norm_mode<3)  0.4];
                else
                    supos(matc,1:4)=[0.125+0.4*(cols(matc)-1)  0.075  0.36  0.4];
                end
            elseif num_all_subplots==3
                supos(matc,1:4)=[0.075+(1/max(cols)-0.03)*(cols(matc)-1)  0.075  0.2+0.0438*(matc==num_all_subplots)*(f_para.color_norm_mode<3)  0.4];
            elseif num_all_subplots==4
                if max(cols)==4
                    supos(matc,1:4)=[0.075+(1/max(cols)-0.03)*(cols(matc)-1)  0.075  0.17+0.0439*((matc==num_mat_subplots)*...
                        (f_para.dendrograms==1)+(matc==num_all_subplots)*(f_para.dendrograms==0))*(f_para.color_norm_mode<3)  0.4];
                else
                    supos(matc,1:4)=[0.23+0.32*(cols(matc)-1)  0.075+0.28*(rows(matc)<max(rows))...
                        0.2+0.0439*(matc==num_all_subplots)*(f_para.color_norm_mode<3)   0.23];
                end
            elseif num_all_subplots==5
                supos(matc,1:4)=[0.07+(1/max(cols)-0.03)*(cols(matc)-1)+0.14*(rows(matc)-1)  0.05+0.29*(matc<=max(cols))...
                    1/max(cols)-0.1+0.06*(f_para.color_norm_mode<3)  0.22];
            elseif num_all_subplots==6
                supos(matc,1:4)=[0.065+(1/max(cols)-0.02-0.01*(f_para.color_norm_mode<3))*(cols(matc)-1)  0.05+0.3*(matc<=max(cols))...
                    1/max(cols)-0.1+0.06*(f_para.color_norm_mode<3)  0.21];
            elseif num_all_subplots==7
                supos(matc,1:4)=[0.065+(1/max(cols)-0.03)*(cols(matc)-1)+0.12*(rows(matc)-1)  0.05+0.29*(matc<=max(cols))...
                    1/max(cols)-0.1+0.06*(f_para.color_norm_mode<3)  0.22];
            elseif num_all_subplots==8
                supos(matc,1:4)=[0.06+(1/max(cols)-0.01-0.02*(f_para.color_norm_mode<3))*(cols(matc)-1)  0.05+0.29*(matc<=max(cols))...
                    1/max(cols)-0.09+(0.02+0.0439*((matc==num_all_subplots)*(f_para.dendrograms==0)+(matc==num_mat_subplots)*(f_para.dendrograms==1)))*(f_para.color_norm_mode<3)  0.22];
            end
        end
    end

    if ruc>1 && exist('liha','var')
        delete(liha)
        clear liha
    end

    for frc=1:num_frame_select+(num_select_averages+num_trigger_averages)*(1+(f_para.mov_num_average_frames-1)*(mod(f_para.plot_mode,32)>15))*(ruc==num_runs)

        fig=figure(4000+f_para.num_fig+frc*f_para.multi_figure*(mod(f_para.plot_mode,16)>7 && mod(f_para.plot_mode,32)<=15));
        set(gcf,'Name',[f_para.filename,'  ',d_para.comment_string,'  ',f_para.comment_string,'  ',f_para.title_string])
        set(gcf,'Position',f_para.pos_fig)
        if f_para.multi_figure && (mod(f_para.plot_mode,16)>7 && mod(f_para.plot_mode,32)<=15)
            set(fig,'DoubleBuffer','on');
            set(gca,'NextPlot','replace','Visible','off')
        end
        if ruc==1 && frc==1 clf; end

        supo1=[0.075 0.56+(num_all_subplots>max_col)*0.1 0.85 0.39-(num_all_subplots>max_col)*0.1];
        subplot('position',supo1)
        if (ruc==1 || (ruc==num_runs && num_frame_select==0)) && (frc==1 || (f_para.multi_figure && (mod(f_para.plot_mode,16)>7 && mod(f_para.plot_mode,32)<=15)))
            cla
            hold on
            xlim([s_para.itmin-0.02*itrange s_para.itmax+0.02*itrange])
            ylim([-0.1 1.2])
            xl=xlim; yl=ylim;
            box on
            for trac=1:num_trains
                for sc=1:num_pspikes(trac)
                    line(pspikes(trac,sc)*ones(1,2),0.05+(num_trains-1-(trac-1)+[0.05 0.95])/num_trains,'Color','k')
                end
            end
            spikesvals=[0.5 1];
            line(xl,0.05*ones(1,2),'Color','k','LineStyle',':')
            line(xl,1.05*ones(1,2),'Color','k','LineStyle',':')
            line(d_para.tmin*ones(1,2),yl,'Color','k','LineStyle',':')
            line(d_para.tmax*ones(1,2),yl,'Color','k','LineStyle',':')
            if num_select_train_groups>1
                for sgc=1:num_select_train_groups
                    if sgc<num_select_train_groups
                        line([s_para.itmin s_para.itmax],(1.05-cum_num_select_group_trains(sgc)/num_trains)*ones(1,2),'Color','k','LineStyle',':')
                    end
                    text(xl(1)-0.06*(xl(2)-xl(1)),1.05-select_group_center(sgc)/num_trains,select_group_names{sgc},'Color','k','FontSize',s_para.font_size)
                end
                set(gca,'YTick',fliplr(1.05-cum_num_select_group_trains/num_trains),'YTickLabel',fliplr(cum_num_select_group_trains))
            else
                ylabel('Spike trains','FontSize',s_para.font_size)
                if mod(num_trains,2)==0
                    set(gca,'YTick',0.05+[0 0.5],'YTickLabel',[num_trains num_trains/2])
                else
                    set(gca,'YTick',0.05+[0 1-(num_trains-1)/2/num_trains],'YTickLabel',[num_trains (num_trains-1)/2])
                end
            end
            if isfield(f_para,'xfact')
                if f_para.xfact>1
                    set(gca,'XTick',[0:f_para.xfact:s_para.itmax],'XTickLabel',0:fix(s_para.itmax/f_para.xfact))
                end
            end
            xlabel(['Time ',f_para.timeunit_string],'FontSize',s_para.font_size)
            if isfield(d_para,'separators') && ~isempty(d_para.separators)
                for sec=1:length(d_para.separators)
                    line([s_para.itmin s_para.itmax],(1.05-d_para.separators(sec)/num_trains)*ones(1,2),'Color','k','LineStyle','--','LineWidth',1)
                end
            end
            if isfield(d_para,'separators2') && ~isempty(d_para.separators2)
                for sec=1:length(d_para.separators2)
                    line([s_para.itmin s_para.itmax],(1.05-d_para.separators2(sec)/num_trains)*ones(1,2),'Color','k','LineStyle','-','LineWidth',2)
                end
            end
            if isfield(d_para,'markers') && ~isempty(d_para.markers)
                for mac=1:length(d_para.markers)
                    line(d_para.markers(mac)*ones(1,2),[0.05 1.05],'Color','k','LineStyle','--','LineWidth',1)
                end
            end
            if isfield(d_para,'markers2') && ~isempty(d_para.markers2)
                for mac=1:length(d_para.markers2)
                    line(d_para.markers2(mac)*ones(1,2),[0.05 1.05],'Color','r','LineStyle','-','LineWidth',2)
                end
            end
            if f_para.dendrograms==1 && f_para.dendro_color_mode>0
                cm=colormap;
                if f_para.dendro_color_mode==1 && num_select_train_groups>1
                    dcol_indy=round(1:63/(num_select_train_groups-1):64);
                    dcols=cm(dcol_indy(select_group_vect),:);
                else
                    dcol_indy=round(1:63/(num_trains-1):64);
                    dcols=cm(dcol_indy,:);
                end
                for trac=1:num_trains
                    %line(d_para.tmin*ones(1,2),1.05-(trac-1+[0 1])/num_trains,'Color',dcols(trac,:),'LineWidth',2)
                    %line(d_para.tmax*ones(1,2),1.05-(trac-1+[0 1])/num_trains,'Color',dcols(trac,:),'LineWidth',2)
                    ph1=patch([d_para.tmin d_para.tmin-0.01*(d_para.tmax-d_para.tmin)*ones(1,2) d_para.tmin],[1.05-(trac-1+[0 0 1 1])/num_trains],dcols(trac,:));
                    ph2=patch([d_para.tmax d_para.tmax+0.01*(d_para.tmax-d_para.tmin)*ones(1,2) d_para.tmax],[1.05-(trac-1+[0 0 1 1])/num_trains],dcols(trac,:));
                    set(ph1,'EdgeColor',dcols(trac,:),'FaceColor',dcols(trac,:))
                    set(ph2,'EdgeColor',dcols(trac,:),'FaceColor',dcols(trac,:))
                end
            end
        end


        if frc<=num_frame_select
            ms=frame_select(frc);
            mi=frc;
            if isfield(d_para,'interval_separators') && ~isempty(d_para.interval_separators) && isfield(d_para,'interval_strings') && ...
                    ~isempty(d_para.interval_strings) && length(d_para.interval_strings)==length(d_para.interval_separators)+1
                inty=find(ms>d_para.interval_separators,1,'last')+1;
                if isempty(inty)
                    inty=1;
                end
                if f_para.publication==0
                    title_str=d_para.interval_strings{inty};
                    title([num2str(ms),' (',num2str(len),')    ---    ',title_str],'Color','k','FontSize',s_para.font_size+3,'FontWeight','bold')
                end
            else
                if f_para.publication==0
                    title([num2str(ms),' (',num2str(len),')'],'Color','k','FontSize',s_para.font_size+2,'FontWeight','bold')
                end
            end
            if frc==1
                liha=zeros(1,3);
                liha(1)=line(frame_select(frc)*ones(1,2),[-0.1 0.05],'Color','g','LineStyle','-','LineWidth',2);
                liha(2)=line(frame_select(frc)*ones(1,2),[1.05 1.2],'Color','g','LineStyle','-','LineWidth',2);
                liha(3)=line(frame_select(frc)*ones(1,2),[0.05 1.05],'Color','g','LineStyle',':','LineWidth',1);
            else
                set(liha(1),'XData',frame_select(frc)*ones(1,2))
                set(liha(2),'XData',frame_select(frc)*ones(1,2))
                set(liha(3),'XData',frame_select(frc)*ones(1,2))
            end
        else
            if frc==num_frame_select+1 && num_frame_select>0 && exist('liha','var') && ~(f_para.multi_figure && (mod(f_para.plot_mode,16)>7 && mod(f_para.plot_mode,32)<=15))
                delete(liha)
                clear liha
            end
            if frc<=num_frame_select+num_select_averages*f_para.mov_num_average_frames
                savc=ceil((frc-num_frame_select)/f_para.mov_num_average_frames);
                savc2=mod((frc-num_frame_select-1),f_para.mov_num_average_frames)+1;
                mi=num_frame_select+savc;
                %[frc savc savc2 mi]
                if f_para.publication==0
                    title(['Mean Select ',num2str(savc)],'Color','k','FontSize',s_para.font_size+2,'FontWeight','bold')
                end
                if savc2==1
                    if savc>1
                        delete(seliha)
                        clear seliha*
                    end
                    seliha=zeros(2,num_sel_ave(savc));
                    for selc=1:num_sel_ave(savc)
                        seliha(1,selc)=line(d_para.select_averages{savc}(2*selc-1:2*selc),-0.025*ones(1,2),'Color','g','LineStyle','-','LineWidth',3);
                        seliha(2,selc)=line(d_para.select_averages{savc}(2*selc-1:2*selc),1.125*ones(1,2),'Color','g','LineStyle','-','LineWidth',3);
                    end
                end
            else
                tavc=ceil((frc-num_frame_select-num_select_averages*f_para.mov_num_average_frames)/f_para.mov_num_average_frames);
                tavc2=mod((frc-num_frame_select-num_select_averages*f_para.mov_num_average_frames-1),f_para.mov_num_average_frames)+1;
                mi=num_frame_select+num_select_averages+tavc;   % f_para.mov_num_average_frames
                %[frc tavc tavc2 mi]
                if frc==num_frame_select+num_select_averages*f_para.mov_num_average_frames+1 && num_select_averages>0 && tavc2==1 && exist('seliha','var') && ~(f_para.multi_figure && (mod(f_para.plot_mode,16)>7 && mod(f_para.plot_mode,32)<=15))
                    delete(seliha)
                    clear seliha
                end
                if f_para.publication==0
                    title(['Mean Trigger ',num2str(tavc)],'Color','k','FontSize',s_para.font_size+2,'FontWeight','bold')
                end
                if frc==num_frame_select+num_select_averages*f_para.mov_num_average_frames+1 || (f_para.multi_figure && (mod(f_para.plot_mode,16)>7 && mod(f_para.plot_mode,32)<=15))
                    trplha=zeros(1,2);
                    trplha(1)=plot(d_para.trigger_averages{tavc},-0.025*ones(1,length(d_para.trigger_averages{tavc})),'g^','LineWidth',2,'MarkerFaceColor','g');
                    trplha(2)=plot(d_para.trigger_averages{tavc},1.125*ones(1,length(d_para.trigger_averages{tavc})),'gv','LineWidth',2,'MarkerFaceColor','g');
                else
                    set(trplha(1),'XData',d_para.trigger_averages{tavc},'YData',-0.025*ones(1,length(d_para.trigger_averages{tavc})))
                    set(trplha(2),'XData',d_para.trigger_averages{tavc},'YData',1.125*ones(1,length(d_para.trigger_averages{tavc})))
                end
            end
        end

        for matc=1:num_mat_subplots

            subplot('position',supos(matc,1:4)); hold on; cla;
            if matc<=num_sel_measures
                plot_mat=flipud(shiftdim(mov_mat(matc,mi,:,:),2));
            else
                plot_mat=flipud(shiftdim(block_mov_mat(matc-num_sel_measures,mi,select_group_vect,select_group_vect),2));
            end
            imagesc(plot_mat);

            if ((ruc==1 || (ruc==num_runs && num_frame_select==0)) && (frc==1 || (num_frame_select==0 && frc==f_para.mov_num_average_frames))) || ...
                    (matc==num_mat_subplots || f_para.color_norm_mode==3)
                if f_para.color_norm_mode==1
                    set(gca,'CLim',[0 1])
                elseif f_para.color_norm_mode==2
                    set(gca,'CLim',[0 max(max(max(max(mov_mat))))])
                elseif f_para.color_norm_mode==3
                    if max(max(plot_mat))>0
                        set(gca,'CLim',[0 max(max(plot_mat))])
                    end
                end
                xlim([0.5 num_trains+0.5])
                ylim([0.5 num_trains+0.5])
                if matc<=num_sel_measures
                    if num_select_train_groups>1
                        set(gca,'XTick',0.5+sort(cum_num_select_group_trains),'XTickLabel',cum_num_select_group_trains)
                        set(gca,'YTick',sort(num_trains-cum_num_select_group_trains),'YTickLabel',fliplr(cum_num_select_group_trains))
                    else
                        set(gca,'YTick',num_trains+1-fliplr(get(gca,'XTick')),'YTickLabel',flipud(get(gca,'XTickLabel')))
                    end
                else
                    set(gca,'XTick',0.5+select_group_center,'XTickLabel',d_para.all_train_group_names)
                    set(gca,'YTick',num_trains-fliplr(select_group_center),'YTickLabel',flipud(d_para.all_train_group_names))
                end
                if length(get(gca,'YTick'))>num_trains
                    set(gca,'YTick',(1:num_trains),'YTickLabel',fliplr(1:num_trains));
                end
                axis square
                if num_all_subplots<4
                    title(mat_str{matc},'Color','k','FontSize',s_para.font_size+4,'FontWeight','bold')
                else
                    title(mat_str{matc},'Color','k','FontSize',s_para.font_size+2,'FontWeight','bold')
                end
                if rows(matc)==max(rows) || f_para.dendrograms
                    xlabel('Spike trains','FontSize',s_para.font_size)
                end
                if cols(matc)==1 || max(cols)<=2
                    ylabel('Spike trains','FontSize',s_para.font_size)
                end
                % [matc get(gca,'Position')]
            end
            xl=xlim; yl=ylim;
            if num_select_train_groups>1
                for sgc=1:num_select_train_groups-1
                    line(xl,(0.5+num_trains-cum_num_select_group_trains(sgc))*ones(1,2),'Color','w','LineStyle',':','LineWidth',1)
                    line((0.5+cum_num_select_group_trains(sgc))*ones(1,2),yl,'Color','w','LineStyle',':','LineWidth',1)
                end
            end
            if isfield(d_para,'separators2') && ~isempty(d_para.separators2)
                for sec=1:length(d_para.separators2)
                    line(xl,(0.5+num_trains-d_para.separators2(sec))*ones(1,2),'Color','k','LineStyle','-','LineWidth',1.5)
                    line((0.5+d_para.separators2(sec))*ones(1,2),yl,'Color','k','LineStyle','-','LineWidth',1.5)
                end
            end
            if matc==num_mat_subplots || f_para.color_norm_mode==3
                colorbar
            end

            if f_para.dendrograms==1
                fig=figure(4000+f_para.num_fig+frc*f_para.multi_figure*(mod(f_para.plot_mode,16)>7 && mod(f_para.plot_mode,32)<=15));
                subplot('position',supos(num_mat_subplots+matc,1:4));
                if matc<=num_sel_measures
                    dmat=shiftdim(mov_mat(matc,mi,:,:),2);
                    mm=dmat(logical(tril(ones(num_trains),-1)));
                else
                    dmat=shiftdim(block_mov_mat(matc-num_sel_measures,mi,:,:),2);
                    mm=dmat(logical(tril(ones(num_select_train_groups),-1)));
                end
                if ~any(any(isnan(dmat)))
                    sdm_linkage=linkage(mm');
                    dh=dendrogram(sdm_linkage,0);

                    if f_para.dendro_color_mode>0
                        cm=colormap;
                        if matc>num_sel_measures || (f_para.dendro_color_mode==1 && num_select_train_groups>1)
                            dcol_indy=round(1:63/(num_select_train_groups-1):64);
                            if matc<=num_sel_measures
                                dcols=cm(dcol_indy(select_group_vect),:);
                            else
                                dcols=cm(dcol_indy,:);
                            end
                        else
                            dcol_indy=round(1:63/(num_trains-1):64);
                            dcols=cm(dcol_indy,:);
                        end
                        yl=ylim;
                        xtl=get(gca,'XTickLabel');
                        xtln=zeros(1,num_trains);
                        for trac=1:size(dmat,1)
                            xtln(trac)=str2double(xtl(trac,:));
                        end

                        trac=0;
                        for lic=1:size(sdm_linkage,1)
                            for cc=1:2
                                if sdm_linkage(lic,cc)<=size(dmat,1)
                                    trac=trac+1;
                                    line(find(xtln==sdm_linkage(lic,cc))*ones(1,2),[yl(1) sdm_linkage(lic,3)],'Color',dcols(sdm_linkage(lic,cc),:),'LineWidth',2.5)
                                end
                            end
                        end
                    end
                end
                axis square
                %set(gca,'YTick',[])
                if num_all_subplots>3
                    set(gca,'YTick',[])
                end
                if matc<=num_sel_measures
                    if num_trains>30
                        set(gca,'XTickLabel',[])
                    elseif num_trains>15
                        set(gca,'FontSize',8)
                    end
                else
                    if num_select_train_groups>30
                        set(gca,'XTickLabel',[])
                    end
                end
                %hold on
                if num_all_subplots<4
                    %title(mat_str{matc},'Color','k','FontSize',s_para.font_size+2,'FontWeight','bold')
                elseif num_all_subplots~=8
                    title(mat_str{matc},'Color','k','FontSize',s_para.font_size,'FontWeight','bold')
                end
                if rows(num_mat_subplots+matc)==max(rows) || num_all_subplots<4
                    if matc<=num_sel_measures
                        xlabel('Spike trains','FontSize',s_para.font_size)
                    else
                        xl=xlim; yl=ylim;
                        text(xl(1)+0.13*(diff(xl)),yl(1)-0.17*(diff(yl)),'Spike train groups','Color','k','FontSize',s_para.font_size)
                        %xlabel('Spike train groups','FontSize',s_para.font_size)
                    end
                end
                axis square
            end
            %get(gca,'Position')
        end

        if mod(f_para.plot_mode,16)>7 && f_para.print_mode                                                                          % Create postscript file
            if f_para.publication==1
                set(gcf,'PaperOrientation','Portrait');set(gcf,'PaperType', 'A4');
                set(gcf,'PaperUnits','Normalized','PaperPosition',[0 0.25 1.0 0.7]);
            else
                set(gcf,'PaperOrientation','Landscape');set(gcf,'PaperType', 'A4');
                set(gcf,'PaperUnits','Normalized','PaperPosition', [0 0 1.0 1.0]);
            end
            if frc<=num_frame_select
                psname=[f_para.imagespath,f_para.filename,d_para.comment_string,f_para.comment_string,num2str(frc),'.ps']
                print(gcf,'-dpsc',psname);
            elseif frc<=num_frame_select+num_select_averages*f_para.mov_num_average_frames
                if savc2==1
                    psname=[f_para.imagespath,f_para.filename,d_para.comment_string,f_para.comment_string,'_SelAve',num2str(savc),'.ps']
                    print(gcf,'-dpsc',psname);
                end
            else
                if tavc2==1
                    psname=[f_para.imagespath,f_para.filename,d_para.comment_string,f_para.comment_string,'_TrigAve',num2str(tavc),'.ps']
                    print(gcf,'-dpsc',psname);
                end
            end
        end

        if mod(f_para.plot_mode,32)>15
            F = getframe(fig);
            tmov = addframe(tmov,F);
        end
    end
end

if mod(f_para.plot_mode,32)>15 && ruc==num_runs
    tmov = close(tmov);
    moviename
end
