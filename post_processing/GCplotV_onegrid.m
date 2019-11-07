function GCplotV_onegrid(Result_cell)

% record the video of membrane potential dynamics
disp('GCplotV...');

Result_num = length(Result_cell);

for r_num = 1:Result_num
    R = Result_cell{r_num};
%     R = Result_cell;
    Num_neu = R.N;
    Num_pop = R.Num_pop;
    dt = R.dt;
    gsize = sqrt(sum(Num_neu));%sqrt(Num_neu);
    V = cell(Num_pop,1);
    
    for pop=1:Num_pop
        file_name = [R.stamp(1:end-3), num2str(pop-1), '_neurosamp.mat'];
        R_samp = load(file_name);
        V{pop}=R_samp.V;        
    end    

    % set the figure invisible, so one can run it without GUI
    figure('NumberTitle', 'off','Name','membrane potential', 'visible', 'off');
% figure
    clf
   % hold on %  add this will be much slower

    caxis('manual')
    colormap('jet')

    vmin=-75.625;
    vmax=-40;
    caxis([vmin,vmax])
    video1 = VideoWriter([file_name,'_V.avi'],'Motion JPEG AVI');
    video1.FrameRate = 15;
    video1.Quality = 50;
    open(video1);
%     input_file = h5read([R.stamp(1:end-18),'.h5'],'/config/pops/pop0/file_current_input/fname'); 
%     start_step = h5read(char(input_file),'/start_step')-200; 
%     end_step = h5read(char(input_file),'/end_step')+200;
%     start_step = 1; 
     [~,end_step] = size(V{1});
    MemP = zeros(gsize);
    MemP(2:2:gsize,2:2:gsize) = 1;
    ind_E = MemP == 0;
    ind_I = MemP == 1;
    for i=2:10:end_step
        %i
        MemP(ind_E) = V{1}(:,i);        
        MemP(ind_I) = V{2}(:,i);                    
        imagesc(MemP,[vmin,vmax]);            
        title(['Vm, ',sprintf('Time: %g ms',(-4002+i)/10)]);           
        axis equal
        axis tight
        
        drawnow
        vframe = getframe(gcf);
        writeVideo(video1,vframe);
    %         pause(0.5)
    %         waitforbuttonpress
    end

     close(video1);
     clear R
end 
end
