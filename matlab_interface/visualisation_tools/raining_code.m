% This code generates digital rain as in Matrix and embeds a message
% Yifan GU, 30 Jan 2017
% yigu8115@gmail.com 
clc;
%clear;
close all;

% random seed
seed = 20170218;
rng(seed);

% generate video file
gen_video = 0;
if gen_video == 1
    vidObj = VideoWriter('test_raining_code.avi');
    open(vidObj);
end
	
% define the frame size and number of time steps
fig_width = 1000;
fig_height = 600;
h = 50;
w = 100;
t = 1000;

% define the character library
slCharacterEncoding('UTF-8')
% c_lib1 = 'SJTUSJTUSJTUSJTUSJTU';
% c_lib2 contains Chinese characters that needs to be defined in command window
c_lib = [c_lib1,c_lib2];
c_length = length(c_lib);
c_map = randi([1, c_length], 1, h*w); % the map between characters and h-w positions; a vector (use vec2ind)
c_change_rate = 0.005; % the charactor-hw map random update rate

% define the message
msg = 'MAY YOUR LOVE LAST FOREVER';
msg = fliplr(msg); % remember to flip the video file horizontally later on
msg_len_tmp = length(msg);
msg_pos = round((w-msg_len_tmp)/2)+(1:msg_len_tmp);
msg_del = isspace(msg); % delete space from msg so that rains can fall through
msg(msg_del) = [];
msg_pos(msg_del) = [];
msg_len = length(msg);
msg_t = 300; % when message start to fall
other_rain_stop_t = 365; % when other new rains stop to fall
msg_drawn = zeros(1, msg_len);
msg_h = floor(h/2-1);
msg_fade_steps = unifrnd(3,25,1,msg_len);
msg_rain_stop = false(1, msg_len);
msg_color = ones(msg_len, 3);

% define the rain
matrix_green = [1/3 0.9 1]; % in HSV color space
head_len = 3; % rain head length
tail_len_range = [15 22]; % random rain tail length range
life_range = [26 40]; % random rain length range
new_rate = 0.025; % rate of generating new rain
v_range = [0.3 0.6]; % random falling speed range, should be less than 1
if life_range(1) < head_len + max(tail_len_range) % check definition consistency
	error('life_range(1) < head_len + tail_len!')
end
			
% initialize the rain state vectors
rain_pos = [];
rain_pos_rev = [];
h_pos = [];
w_pos = [];
rain_tail_len = [];
v_pos = []; % different speed of the tick
v_tick = []; % when the tick runs over an integer, the rain falls by 1 step

% define the color-hw map
hsv_v_map_rand = unifrnd(0.5, 1, 1, w*h); % hsv v-value map

% intialize the figure
figure(1);
set(gcf,'color','k','Position', [100, 100, 100+fig_width, 100+fig_height])%, 'visible','off')
axis off;
axis([0 w+1 0 h+1]);
set(gca,'YDir','reverse')
set(gca,'color','k')
hold on;
set(gcf,'Units','normal')
set(gca,'Position',[0 0 1 1])

% start the simulation
tic
for i = 1:t
    % generate new rains
    if i < msg_t
        new_rain = setdiff(1:w, unique(w_pos(h_pos < h/3)));
    elseif i < other_rain_stop_t
        other = setdiff(1:w, msg_pos);
        new_rain = setdiff(msg_pos(~msg_rain_stop), unique(w_pos(h_pos < h/3)));
        new_rain = [new_rain other];
    else
        new_rain = setdiff(msg_pos(~msg_rain_stop), unique(w_pos(h_pos < h/3)));
    end
    new_rain(rand(1,length(new_rain)) > new_rate) = [];
    
    % add new rains
    w_pos = [w_pos new_rain]; %#ok<*AGROW>
    h_pos = [h_pos ones(1,length(new_rain))];
    v_pos = [v_pos unifrnd(v_range(1), v_range(2), 1,length(new_rain))];
    v_tick = [v_tick zeros(1, length(new_rain))];
    new_rain_length = randi(life_range, 1,length(new_rain));
    rain_pos = [rain_pos new_rain_length]; % counting from N to 1 (from head to tail)
    rain_pos_rev = [rain_pos_rev ones(1,length(new_rain))];% counting from 1 to N (from head to tail)
    rain_tail_len = [rain_tail_len randi(tail_len_range, 1,length(new_rain))]; % random tail length
    
    
    % falling down
    f_fall =  h_pos == 1 & (rain_pos > 1) & ((floor(v_tick+v_pos)-floor(v_tick)) == 1);
    
    w_pos_new = {w_pos, w_pos(f_fall)};
    rain_pos_new = {rain_pos, rain_pos(f_fall)-1};
    h_pos_new = h_pos;
    h_pos_new((floor(v_tick+v_pos)-floor(v_tick)) == 1) = h_pos_new((floor(v_tick+v_pos)-floor(v_tick)) == 1) + 1;
    h_pos_new = {h_pos_new, ones(1, sum(f_fall))};
    v_pos_new = {v_pos, v_pos(f_fall)};
    v_tick_new = {v_tick, v_tick(f_fall)};
    
    rain_pos_rev_new = {rain_pos_rev, rain_pos_rev(f_fall)+1};
    rain_tain_len_new = {rain_tail_len, rain_tail_len(f_fall)};
    
    w_pos = cell2mat(w_pos_new);
    rain_pos = cell2mat(rain_pos_new);
    h_pos = cell2mat(h_pos_new);
    v_pos = cell2mat(v_pos_new);
    v_tick = cell2mat(v_tick_new)+ v_pos;
    rain_pos_rev = cell2mat(rain_pos_rev_new);
    rain_tail_len = cell2mat(rain_tain_len_new);
    
    % delete those falling out of screen
    w_pos(h_pos > h) = [];
    rain_pos(h_pos > h) = [];
    v_pos(h_pos > h) = [];
    v_tick(h_pos > h) = [];
    rain_pos_rev(h_pos > h) = [];
    rain_tail_len(h_pos > h) = [];
    h_pos(h_pos > h) = [];
    
    % update character-to-hw map
    c_change = rand(1,h*w) < c_change_rate;
    c_map( c_change ) = randi([1, c_length], 1, sum(c_change));
    c_pos = c_lib(c_map(sub2ind([h w], h_pos, w_pos)));
    
    % generate the brightness-to-hw map based on head and tail effects
    hsv_v_pos = hsv_v_map_rand(sub2ind([h w], h_pos, w_pos));
    % head effect
    head_effect = linspace(0, 0.5, head_len);
    hsv_s_pos_head = ones(size(rain_pos_rev));
    is_rain_head = rain_pos_rev <= head_len;
    hsv_s_pos_head(is_rain_head) =  head_effect(rain_pos_rev(is_rain_head));
    % tail effect
    hsv_v_pos_tail = ones(size(rain_pos_rev));
    is_rain_tail = rain_pos <= rain_tail_len;
    hsv_v_pos_tail(is_rain_tail) = (rain_pos(is_rain_tail)-1)./rain_tail_len(is_rain_tail);

    % draw
    cla;
    msg_j_del = [];
    for j = 1:length(w_pos)
        % embed the message
        if i >= msg_t && h_pos(j) == msg_h && sum(w_pos(j) == msg_pos) == 1 % hitting massage positions
            msg_pos_ind = find(w_pos(j) == msg_pos);
            if msg_drawn(msg_pos_ind) == 1 % if the message has been drawn
                msg_j_del = [msg_j_del, j];
                continue; % skipping the draw loop
            elseif rain_pos_rev(j) == 1 % if not drawn yet
                msg_drawn(msg_pos_ind) = 1;
                msg_rain_stop(msg_pos_ind) = true;
                msg_j_del = [msg_j_del, j];
                continue; % skipping the draw loop
            end
        end
        
        % draw the rains
        color_tmp = matrix_green;
        if rain_pos_rev(j) <= head_len % head effect
            color_tmp(2) = hsv_s_pos_head(j);
            color_tmp = hsv2rgb(color_tmp);
            text( w_pos(j), h_pos(j), c_pos(j), 'color', color_tmp , ...
                'FontWeight', 'bold', 'HorizontalAlignment','center','VerticalAlignment','top',...
                'FontName','Courier')
        else
            color_tmp(3) = color_tmp(3)*hsv_v_pos(j)*hsv_v_pos_tail(j);
            color_tmp = hsv2rgb(color_tmp);
            text( w_pos(j), h_pos(j), c_pos(j), 'color', color_tmp , ...
                'FontWeight', 'bold', 'HorizontalAlignment','center','VerticalAlignment','top',...
                'FontName','Courier')
        end
    end
    
    % delete those falling beyond the message
    w_pos(msg_j_del) = [];
    rain_pos(msg_j_del) = [];
    v_pos(msg_j_del) = [];
    v_tick(msg_j_del) = [];
    rain_pos_rev(msg_j_del) = [];
    rain_tail_len(msg_j_del) = [];
    h_pos(msg_j_del) = [];

    % draw the message
    for m = 1:length(msg_drawn)
        if msg_drawn(m) == 1
            text( msg_pos(m), msg_h, msg(m), 'color',  msg_color(m,:), ...
                'FontWeight', 'bold', 'HorizontalAlignment','center','VerticalAlignment','top',...
                'FontName','Courier')
        end
    end
    
    % the messaage fades
    if isempty(h_pos) && nnz(msg_drawn-1) == 0
        for m = 1:length(msg_drawn)
            msg_color(m,:) = msg_color(m,:).*exp(-1/msg_fade_steps(m));
        end
    end
    
    % save to the video file
    drawnow;
    if gen_video == 1
        writeVideo(vidObj, getframe(gcf));
    end
    
    % show progress
    i/t %#ok<NOPTS>
end
toc; % show the total time taken

% close the video file
if gen_video == 1
    close(vidObj);
end





