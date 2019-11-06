function GCplotSpike_onegrid(Result_cell)

% record the video of membrane potential dynamics
disp('GCplotSpike...');
Result_num = length(Result_cell);
for r_num = 1:Result_num
    R = Result_cell{r_num};
    Num_neu = R.N;
%     Num_pop = R.Num_pop;
    spike_hist = R.Analysis.spike_hist_combined;
    file_name = R.stamp(1:end-3);
    gridsize = sqrt(sum(Num_neu));
    
    % set the figure invisible, so one can run it without GUI
    figure('NumberTitle', 'off','Name','membrane potential', 'visible', 'off');
    clf
%     hold off
    video1 = VideoWriter([file_name,'_S.avi'],'Motion JPEG AVI');
    video1.FrameRate = 15;
    video1.Quality = 50;
    open(video1);
    %     input_file = h5read([R.stamp(1:end-18),'.h5'],'/config/pops/pop0/file_current_input/fname');
    %     start_step = h5read(input_file,'start_step')-200;
    %     end_step = h5read(input_file,'start_step')+200;
    start_step = 10000;
    end_step = 22000;
    step = 10;
    for i=start_step:step:end_step
        grid = reshape(spike_hist(:,i),gridsize,gridsize);
        [y,x] = find(grid);
        scatter(x,abs(gridsize-y),'k.');
        axis equal
        axis tight
        axis([1,gridsize,1,gridsize]);
        title(sprintf('Time: %f ms', (start_step + (i-start_step)*step)*R.dt))
        drawnow
        vframe = getframe(gcf);
        writeVideo(video1,vframe);
        %         pause(0.5)
        %         waitforbuttonpress
    end
    
    close(video1);
    
end
end

