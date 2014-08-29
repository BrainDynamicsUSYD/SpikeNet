% different cuts (subplots) for each measure (figure)

for matc=1:num_sel_measures
    figure(2000+matc);
    clf; hold on
    set(gcf,'Position',[5 39 1912 1072])
    
    subplot('position',[0.075 0.55+(num_frame_select>3)*0.15 0.85 0.4-(num_frame_select>3)*0.15])

    xlim([s_para.itmin-0.02*itrange s_para.itmax+0.02*itrange])
    ylim([0 1.1])
    xl=xlim; yl=ylim;
    box on
    for trac=1:num_trains
        for sc=1:num_pspikes(trac)
            line(pspikes(trac,sc)*ones(1,2),0.05+(num_trains-1-(trac-1)+[0.05 0.95])/num_trains,'Color','b')
        end
    end
    spikesvals=[0.5 1];
    line(xl,0.05*ones(1,2),'Color','k','LineStyle',':')
    line(xl,1.05*ones(1,2),'Color','k','LineStyle',':')
    line(d_para.tmin*ones(1,2),yl,'Color','k','LineStyle',':')
    line(d_para.tmax*ones(1,2),yl,'Color','k','LineStyle',':')
    if matc<=num_mat_subplots
        set(gcf,'Name',[mat_str{matc}])
        if f_para.publication==0
            title([mat_str{matc}],'Color','k','FontSize',s_para.font_size+2,'FontWeight','bold')
        end
    else
        set(gcf,'Name',['Dendrogram ',mat_str{matc-num_mat_subplots}])
        if f_para.publication==0
            title(['Dendrogram ',mat_str{matc-num_mat_subplots}],'Color','k','FontSize',s_para.font_size+2,'FontWeight','bold')
        end
    end
    for frc=1:num_frame_select
        ms=frame_select(frc);
        line(time(ms)*ones(1,2),[0.05 1.05],'Color','g','LineStyle','-','LineWidth',2)
    end
    %text(xl(1)-0.05*(xl(2)-xl(1)),yl(2)+0.07*(yl(2)-yl(1)),'A','Color','k','FontSize',s_para.font_size,'FontWeight','bold')
    xlabel(['Time ',f_para.timeunit_string],'FontSize',s_para.font_size,'FontWeight','bold')
    if isfield(f_para,'xfact')
        set(gca,'XTick',[0:f_para.xfact:s_para.itmax],'XTickLabel',0:fix(s_para.itmax/f_para.xfact))
    end
    ylabel('Spike trains','FontSize',s_para.font_size,'FontWeight','bold')
    set(gca,'YTick',0.05+[0 0.5],'YTickLabel',[num_trains num_trains/2])
    
    max_cols=4;
    rows=ceil((1:num_frame_select)/max_cols);
    cols=(1:num_frame_select)-max_cols*((1:num_frame_select)>max_cols);
    for frc=1:num_frame_select
        
        if num_frame_select==1
            supos=[0.35  0.075  0.33   0.36];
        elseif num_frame_select==2
            supos=[0.125+0.4*(cols(frc)-1)  0.075  0.36   0.4];
        elseif num_frame_select==3
            supos=[0.075+(1/max(cols)-0.03)*(cols(frc)-1)  0.075  0.2+0.0438*(frc==num_frame_select)   0.4];
        elseif num_frame_select==4
            if max(cols)==4
                supos=[0.075+(1/max(cols)-0.03)*(cols(frc)-1)  0.075  0.17+0.0438*(frc==num_frame_select)   0.4];
            else
                supos=[0.075+(1/max(cols)-0.03)*(cols(frc)-1)  0.075+0.3*(frc<=max(cols))  0.17   0.23];
            end
        elseif num_frame_select==5
            supos=[0.075+(1/max(cols)-0.03)*(cols(frc)-1)+0.15*(rows(frc)-1)  0.075+0.3*(frc<=max(cols))  1/max(cols)-0.04   0.23];
        elseif num_frame_select==6
            supos=[0.075+(1/max(cols)-0.03)*(cols(frc)-1)  0.075+0.3*(frc<=max(cols))  1/max(cols)-0.04   0.23];
        elseif num_frame_select==7
            supos=[0.075+(1/max(cols)-0.03)*(cols(frc)-1)+0.12*(rows(frc)-1)  0.075+0.3*(frc<=max(cols))  1/max(cols)-0.04   0.23];
        elseif num_frame_select==8
            supos=[0.075+(1/max(cols)-0.03)*(cols(frc)-1)  0.075+0.3*(frc<=max(cols))  1/max(cols)-0.04   0.23];
        end
        subplot('position',supos)
        
        if matc<=num_sel_measures
            imagesc(shiftdim(mov_mat(matc,frc,:,:),2));
            xlim([0.5 num_trains+0.5])
            ylim([0.5 num_trains+0.5])
            if rows(frc)==max(rows)
                xlabel('Spike trains','FontSize',s_para.font_size,'FontWeight','bold')
            end
            if cols(frc)==1 || num_frame_select<4
                %text(xl(1)-0.22*(xl(2)-xl(1)),yl(1)-0.2*(yl(2)-yl(1)),alphabet(frc+1),'Color','k','FontSize',s_para.font_size,'FontWeight','bold')
                ylabel('Spike trains','FontSize',s_para.font_size,'FontWeight','bold')
            end
        elseif matc<=num_mat_subplots
            imagesc(shiftdim(block_mov_mat(matc-num_sel_measures,frc,:,:),2));
            xlim([0.5 num_select_train_groups+0.5])
            ylim([0.5 num_select_train_groups+0.5])
            set(gca,'XTick',1:num_select_train_groups,'XTickLabel',d_para.all_train_group_names)
            set(gca,'YTick',1:num_select_train_groups,'YTickLabel',d_para.all_train_group_names)
            if rows(frc)==max(rows)
                xlabel('Spike train groups','FontSize',s_para.font_size,'FontWeight','bold')
            end
            if cols(frc)==1 || num_frame_select<4
                %text(xl(1)-0.22*(xl(2)-xl(1)),yl(1)-0.2*(yl(2)-yl(1)),alphabet(frc+1),'Color','k','FontSize',s_para.font_size,'FontWeight','bold')
                ylabel('Spike train groups','FontSize',s_para.font_size,'FontWeight','bold')
            end
        else
            if matc<=num_sel_measures
                dmat=shiftdim(mov_mat(matc-num_mat_subplots,frc,:,:),2);
                mm=dmat(logical(tril(ones(num_trains),-1)));
            else
                dmat=shiftdim(block_mov_mat(matc-num_mat_subplots-num_sel_measures,frc,:,:),2);
                mm=dmat(logical(tril(ones(num_select_train_groups),-1)));
            end
            sdm_linkage=linkage(mm');
            dh=dendrogram(sdm_linkage,0);
            axis square
            %hold on
            if f_para.publication
                title(mat_str_short{matc-num_mat_subplots},'Color','k','FontSize',s_para.font_size+2,'FontWeight','bold')
            else
                if num_sel_measures<4
                    title(mat_str_long{matc-num_mat_subplots},'Color','k','FontSize',s_para.font_size+2,'FontWeight','bold')
                else
                    title(mat_str_long{matc-num_mat_subplots},'Color','k','FontSize',s_para.font_size,'FontWeight','bold')
                end
            end
            if rows(matc-num_mat_subplots)==max(rows-num_mat_subplots) || num_mat_subplots<4
                if matc-num_mat_subplots<=num_sel_measures
                    xlabel('Spike trains','FontSize',s_para.font_size,'FontWeight','bold')
                else
                    xlabel('Spike train groups','FontSize',s_para.font_size,'FontWeight','bold')
                end
            end
            if num_frame_select>3
                set(gca,'YTick',[])
            end
            if matc<=num_sel_measures
                if num_trains>30
                    set(gca,'XTickLabel',[])
                end
            else
                if num_select_train_groups>30
                    set(gca,'XTickLabel',[])
                end
            end
        end
        if matc<=num_mat_subplots
            set(gca,'CLim',[0 1])
            xl=xlim;yl=ylim;
            axis square
            if f_para.publication==0
                title(['T = ',num2str(frame_select(frc))],'Color','k','FontSize',s_para.font_size+1,'FontWeight','bold')
            end
        end
        if frc==num_frame_select || f_para.color_norm_mode==3
            colorbar
        end

    end
    if f_para.print_mode                                                                         % Create postscript file
        if f_para.publication==1
            set(gcf,'PaperOrientation','Portrait');set(gcf,'PaperType', 'A4');
            set(gcf,'PaperUnits','Normalized','PaperPosition',[0 0.45 1.0 0.5]);
        else
            set(gcf,'PaperOrientation','Landscape');set(gcf,'PaperType', 'A4');
            set(gcf,'PaperUnits','Normalized','PaperPosition', [0 0 1.0 1.0]);
        end
        if matc<=num_sel_measures
            psname=[f_para.imagespath,f_para.filename,f_para.comment_string,'_Matrix_Cuts_',mat_str{matc},'.ps'];
        elseif matc<=num_mat_subplots
            psname=[f_para.imagespath,f_para.filename,f_para.comment_string,'_Matrix_Cuts_',mat_str{matc-num_sel_measures},'G.ps'];
        elseif matc<=num_mat_subplots+num_sel_measures
            psname=[f_para.imagespath,f_para.filename,f_para.comment_string,'_Matrix_Cuts_Dendro_',mat_str{matc-num_mat_subplots},'.ps'];
        else
            psname=[f_para.imagespath,f_para.filename,f_para.comment_string,'_Matrix_Cuts_Dendro_',mat_str{matc-num_mat_subplots-num_sel_measures},'G.ps'];
        end
        print(gcf,'-dpsc',psname);
    end

end

