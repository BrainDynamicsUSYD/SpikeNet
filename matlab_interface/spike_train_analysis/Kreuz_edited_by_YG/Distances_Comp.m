% different measures (rows) and cuts (columns)

num_mat_select=size(mov_mat,2);

num_matrices=num_sel_measures*num_mat_select;
max_cols=4;
rows=ceil((1:num_matrices)/max_cols);
cols=(1:num_matrices)-max_cols*((1:num_matrices)>max_cols);
supos=zeros(num_matrices,4);
for matc=1:num_matrices
    if num_matrices==1
        supos(matc,1:4)=[0.35  0.075  0.33   0.36];
    elseif num_matrices==2
        supos(matc,1:4)=[0.125+0.4*(cols(matc)-1)  0.075  0.36   0.4];
    elseif num_matrices==3
        supos(matc,1:4)=[0.075+(1/max(cols)-0.03)*(cols(matc)-1)  0.075  0.2+0.0438*(matc==num_matrices)   0.4];
    elseif num_matrices==4
        if max(cols)==4
            supos(matc,1:4)=[0.075+(1/max(cols)-0.03)*(cols(matc)-1)  0.075  0.17+0.0438*(matc==num_matrices)   0.4];
        else
            supos(matc,1:4)=[0.25+0.3*(cols(matc)-1)  0.075+0.3*(matc<=max(cols))  0.2   0.23];
        end
    elseif num_matrices==5
        supos(matc,1:4)=[0.075+(1/max(cols)-0.03)*(cols(matc)-1)+0.15*(rows(matc)-1)  0.075+0.3*(matc<=max(cols))  1/max(cols)-0.04   0.23];
    elseif num_matrices==6
        supos(matc,1:4)=[0.075+(1/max(cols)-0.03)*(cols(matc)-1)  0.075+0.3*(matc<=max(cols))  1/max(cols)-0.04   0.23];
    elseif num_matrices==7
        supos(matc,1:4)=[0.075+(1/max(cols)-0.03)*(cols(matc)-1)+0.12*(rows(matc)-1)  0.075+0.3*(matc<=max(cols))  1/max(cols)-0.04   0.23];
    elseif num_matrices==8
        supos(matc,1:4)=[0.12+(1/max(cols)-0.05)*(cols(matc)-1)  0.095+0.24*(matc<=max(cols))  1/max(cols)-0.07   0.17];
    end
end

figure(1000+f_para.num_fig);
clf; hold on
set(gcf,'Position',f_para.pos_fig)
set(gcf,'Name',[f_para.filename,'  ',d_para.comment_string,'  ',f_para.comment_string,'  ',f_para.title_string])

supo1=[0.12 0.48+(num_sel_measures-1)*0.1 0.8 0.47-(num_sel_measures-1)*0.1];

subplot('position',supo1)
%subplot('position',[0.075 0.55+(num_matrices>3)*0.15 0.85 0.4-(num_matrices>3)*0.15])
cla

xlim([s_para.itmin-0.02*itrange s_para.itmax+0.02*itrange])
ylim([-0.3 1.4])
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
%text(xl(1)-0.13*(xl(2)-xl(1)),yl(2)+0.07*(yl(2)-yl(1)),'A','Color','k','FontSize',s_para.font_size+5,'FontWeight','bold')
set(gca,'YTick',0.05+[0 0.25 0.5 0.75],'YTickLabel',[num_trains 3*num_trains/4 num_trains/2 num_trains/4])
text(xl(1)+0.44*(diff(xl)),yl(1)-0.12*(diff(yl)),['Time  ',f_para.timeunit_string],'Color','k','FontSize',s_para.font_size)
text(xl(1)-0.1*(xl(2)-xl(1)),yl(2)-0.43*(yl(2)-yl(1)),'Spike','Color','k','FontSize',s_para.font_size)
text(xl(1)-0.1*(xl(2)-xl(1)),yl(2)-0.59*(yl(2)-yl(1)),'trains','Color','k','FontSize',s_para.font_size)
set(gca,'FontSize',s_para.font_size-1)

for frc=1:num_frame_select
    line(f_para.mov_frames(frc)*ones(1,2),[-0.20 0.05],'Color','g','LineStyle','-','LineWidth',2)
    line(f_para.mov_frames(frc)*ones(1,2),[1.05 1.3],'Color','g','LineStyle','-','LineWidth',2)
    %line(f_para.mov_frames(frc)*ones(1,2),[0.05 1.05],'Color','g','LineStyle',':','LineWidth',2)
end
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
if num_select_averages>0
    for savc=1:num_select_averages
        for selc=1:num_sel_ave(savc)
            seliha1(selc)=line(d_para.select_averages{savc}(2*selc-1:2*selc),(-0.025-0.06*(savc-1))*ones(1,2),'Color','g','LineStyle','-','LineWidth',3);
            seliha2(selc)=line(d_para.select_averages{savc}(2*selc-1:2*selc),(1.125+0.06*(num_select_averages-savc))*ones(1,2),'Color','g','LineStyle','-','LineWidth',3);
            %seliha3(selc)=line(d_para.select_averages{savc}(2*selc-1)*ones(1,2),[-0.025 0.025]-0.05*(savc-1),'Color','g','LineStyle','-','LineWidth',3);
            %seliha4(selc)=line(d_para.select_averages{savc}(2*selc-1)*ones(1,2),[1.075 1.125]+0.05*(savc-1),'Color','g','LineStyle','-','LineWidth',3);
            %seliha5(selc)=line(d_para.select_averages{savc}(2*selc)*ones(1,2),[-0.025 0.025]-0.05*(savc-1),'Color','g','LineStyle','-','LineWidth',3);
            %seliha6(selc)=line(d_para.select_averages{savc}(2*selc)*ones(1,2),[1.075 1.125]+0.05*(savc-1),'Color','g','LineStyle','-','LineWidth',3);
        end
    end
end

max_cols=min([num_mat_select 10]);

rows=reshape(cumsum(ones(num_sel_measures,num_mat_select))',1,num_sel_measures*num_mat_select);
cols=reshape(cumsum(ones(num_mat_select,num_sel_measures)),1,num_sel_measures*num_mat_select);

for matc=1:num_sel_measures
    for frc=1:num_mat_select
        grc=(matc-1)*num_mat_select+frc;

        if grc<8
            subplot('position',supos(grc,:))
        else
            supos(grc,3)=supos(grc,3)*1.2434;
            subplot('position',supos(grc,:))
        end

        plot_mat=shiftdim(mov_mat(matc,frc,:,:),2);
        isc=imagesc(plot_mat);
        color_norm_mode=1;
        if color_norm_mode==1
            set(gca,'CLim',[0 1])
        elseif color_norm_mode==2
            set(gca,'CLim',[0 max(max(max(max(mov_mat))))])
        elseif color_norm_mode==3
            if max(max(plot_mat))>0
                set(gca,'CLim',[0 max(max(plot_mat))])
            end
        end
        xlim([0.5 num_trains+0.5])
        ylim([0.5 num_trains+0.5])
        xl=xlim; yl=ylim;
        if matc==num_sel_measures && frc==num_mat_select
            colorbar
        end
        axis square
        if f_para.publication==0 && rows(grc)==1 && frc<=num_frame_select
            title(['T = ',num2str(frame_select(frc))],'Color','k','FontSize',s_para.font_size+1,'FontWeight','bold')
        end
        %if matc==num_sel_measures
        xlabel('Spike trains','FontSize',s_para.font_size-1)
        %end
        if frc==1
            ylabel('Spike trains','FontSize',s_para.font_size-1)
            if num_sel_measures>1
                text(xl(1)-0.73*(xl(2)-xl(1)),yl(1)+0.46*(yl(2)-yl(1)),mat_str{matc},'Color','k','FontSize',s_para.font_size+3,'FontWeight','bold')
            end
            if matc==1
                %text(xl(1)-0.6*(xl(2)-xl(1)),yl(1)-0.2*(yl(2)-yl(1)),'B','Color','k','FontSize',s_para.font_size+5,'FontWeight','bold')
            end
        end
        if frc>1
            set(gca,'YTick',[])
        end

        set(gca,'FontSize',s_para.font_size-2)
        %[grc get(gca,'Position')]
    end

end

if f_para.print_mode                                                                         % Create postscript file
    if f_para.publication==1
        set(gcf,'PaperOrientation','Portrait');set(gcf,'PaperType', 'A4');
        set(gcf,'PaperUnits','Normalized','PaperPosition',[0 0.15 1.0 0.8]);
    else
        set(gcf,'PaperOrientation','Landscape');set(gcf,'PaperType', 'A4');
        set(gcf,'PaperUnits','Normalized','PaperPosition', [0 0 1.0 1.0]);
    end
    psname=[f_para.imagespath,f_para.filename,f_para.comment_string,'_Matrix_Cuts.ps'];
    print(gcf,'-dpsc',psname);
end


