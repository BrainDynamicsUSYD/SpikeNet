function main_example_cc_half_module_size(varargin)
% Do it!!!
% Find it!!!
% Hunt it down!!!


% varargin is for PBS arrary job
if nargin == 0
    clc;clear all;close all;
    cd /import/yossarian1/yifan/Project1/
    addpath(genpath(cd));
    cd tmp_data
end % Basic parameters
N = [4000; 1000]; %!!!
dt = 0.1;
sec = round(10^3/dt); % 1*(10^3/dt) = 1 sec
step_tot = 200*sec; % use 10 second!

% Loop number for PBS array job
Num_pop = length(N);
loop_num = 0;
discard_transient = 500; % ms
EE_factor = 0.6;
II_factor = 0.8;
rr = 0.7; % this is different from 0.6!!
% Mnum = 16; %!!!
kk = 1; %2:5; % use 2 to roughly compensate synaptic saturation
phi_E = 1;
phi_I = 1.05;

Mnum = 16; % don't delete it! It will be useful later!



for lesion_0 = 0.2
    
for lesion_1 = 1.1      %+(-0.05:0.05:0.05) %1.1:0.1:1.4 % range [0-1]
    for lesion_2 = 1    %+(-0.05:0.05:0.05)
        for lesion_3 = 1    %+(-0.05:0.05:0.05)
            for lesion_4 = 0.6  %+(-0.05:0.05:0.05)

                for I_ext_strength = repmat(1.1:0.1:1.4, 1, 10) %1.2*ones(1,40)
                    
                    
                    loop_num = loop_num + 1;
                    
                    % For PBS array job
                    if nargin ~= 0
                        PBS_ARRAYID = varargin{1};
                        if loop_num ~=  PBS_ARRAYID
                            continue;
                        end
                    end
                    
                    % seed the matlab rand function! The seed is global.
                    % Be very careful about that!!!!!!!!!!!
                    rand_seed = loop_num*10^5+eval(datestr(now,'SSFFF'));
                    rng(rand_seed,'twister');
                    
                    % Creat ygin file
                    % using loop_num in filename to ensure unique naming!
                    % Otherwise overwriting may occur when using PBS.
                    name = [ sprintf('%03g-', loop_num), datestr(now,'yyyymmddHHMM-SSFFF')];
                    
                    fprintf('Data file name is: /n%s/n', strcat(name,'.ygin') ); % write the file name to stdout and use "grep ygin" to extract it
                    FID = fopen([name,'.ygin'], 'w'); % creat file
                    FID_syn = fopen([name,'.ygin_syn'], 'w'); % creat file
                    
                    % write basic parameters
                    writeBasicPara(FID, dt, step_tot, N)
                    % write pop para
                    writePopPara(FID, 1,  'tau_ref', 2);
                    writePopPara(FID, 2,  'tau_ref', 2);
                    % write synapse para
                    writeSynPara(FID, 'tau_decay_GABA', 3);
                    
                    %%%%%%% write runaway killer
                    runaway_steps = round(50/dt); % 50 ms
                    runaway_mean_num_ref = 0.2;
                    writeRunawayKiller(FID, runaway_steps, runaway_mean_num_ref);
                    %%%%%%%%%%%%%%%%%%%%%%%
                    
                    
                    
                    Kmat = [2.4*EE_factor  1.4;
                        4.5  5.7*II_factor]*kk*10^-3; % miuSiemens
                    
                    Kmat(1,:) = Kmat(1,:)*phi_E;
                    Kmat(2,:) = Kmat(2,:)*phi_I;
                    
                    Pmat = [0.2 0.5;
                        0.5 0.5];
                    
                    TYPEmat = [1 2];
                    
                    
                    
                    % External current
                    writeExtCurrentSettings(FID, 1, I_ext_strength, 0)
                    writeExtCurrentSettings(FID, 2, I_ext_strength, 0)
                    
                    
                           
                    %%%%%%% data sampling
                    writeNeuronSampling(FID, 1, [1,1,1,1,0,0,1],[100:500:4000]);
                    writeNeuronSampling(FID, 2, [1,1,1,1,0,0,1],[100;600]);
                    
                    
                    %%%%%%% random initial condition settings (int pop_ind, double p_fire)
                    p_fire = 0.00*ones(size(N)); % between [0,1], 0.05
                    writeInitV(FID, p_fire);
                    
                    %%%%%%%%%%%%%%%%%%% Chemical Connections %%%%%%%%%%%%%%%%%%%%%%%
                    % type(1:AMAP, 2:GABAa, 3:NMDA)
                    
                    %
                    for i_pre = 1:2
                        for j_post = 1:2
                            
                            if i_pre == 1 && j_post == 1
                                
                                %%%%%%%% Hierarchical Connection and Lesion%%%%%%%%%%%%%%%%%
                                % Get the guts out of
                                % it so I can do lesion
                                % study
                                P0 = Pmat(i_pre,j_post);
                                
                               
                                Msize8 = Mnum_2_Msize( 8, N(i_pre) );
                                [P8, CL] = inter_module_Pmatrix(Msize8, P0, rr);
                                P8(CL==4) = P8(CL==4)*lesion_4; % highest level connection
                                P8(CL==3) = P8(CL==3)*lesion_3;  %*lesion_left;  % second-highest level connection
                                P8(CL==2) = P8(CL==2)*lesion_2; % third-highest level connection
                                P8(CL==1) = P8(CL==1)*lesion_1;
                                
                                P16 = zeros(16, 16);
                                for ip = 1:16
                                    for jp = 1:16
                                        P16(ip,jp) = P8( floor((ip-1)/2)+ 1, floor((jp-1)/2) + 1 );
                                    end
                                end
                                
                                for ip = 1:16
                                    for jp = 1:16
                                        i8 = floor((ip-1)/2) + 1;
                                        j8 = floor((jp-1)/2) + 1;
                                        if i8 == j8
                                            if ip == jp
                                                P16(ip,jp) = P16(ip,jp)*(1+lesion_0);
                                            else
                                                P16(ip,jp) = P16(ip,jp)*2 - P16(ip,jp)*(1+lesion_0);
                                            end
                                        end
                                    end
                                end
                                
                                
                                Msize16 = Mnum_2_Msize( 16, N(i_pre) );
                                
                                A11 = P_2_A(P16,Msize16);
                                
                                % Display
                               
                                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                                [I, J, ~] = find(A11);
                                % % save A11
                                % save([name,'A11.mat'], 'A11');
                            else
                                [I, J, ~] = find(MyRandomGraphGenerator('E_R_pre_post','N_pre',N(i_pre),'N_post',N(j_post),'p',Pmat(i_pre,j_post)));
                            end
                            
                            
                            K = ones(size(I))*Kmat(i_pre,j_post);
                            D = rand(size(I))*1;
                            writeChemicalConnection(FID_syn, TYPEmat(i_pre),  i_pre,j_post,   I,J,K,D); % (FID, type, i_pre, j_post, I, J, K, D)
                            clear I J K D;
                        end
                    end
                    
                    
                    
                    % Explanatory (ExplVar) and response variables (RespVar) for cross-simulation data gathering and post-processing
                    % Record explanatory variables, also called "controlled variables"
                    
                    writeExplVar(FID, 'discard_transient', discard_transient, ...
                        'loop_num', loop_num, ...
                        'r', rr, ...
                        'Mnum', Mnum, ...
                        'I_ext_strength', I_ext_strength, ...
                        'lesion_4',lesion_4,'lesion_3',lesion_3,'lesion_2',lesion_2,'lesion_1',lesion_1, 'lesion_0', lesion_0,...
                        'phi_E', phi_E, 'phi_I', phi_I);
                    
                    
                    % Adding comments in raster plot
                    comment1 = 'p=[0.2 0.5 0.5 0.5], k = [2.4*EE_fatcor 1.4;4.5 5.7*II_factor]*k*10^-3, tau_decay_GABA=3';
                    comment2 = datestr(now,'dd-mmm-yyyy-HH:MM');
                    writeExplVar(FID, 'comment1', comment1, 'comment2', comment2);
                    
                    
                    % append this file self into .ygin for future reference
                    appendThisMatlabFile(FID)
                    
                    
                end
            end
            
        end
    end
end

end

end


% This function must be here!
function appendThisMatlabFile(FID)
breaker = ['>',repmat('#',1,80)];
fprintf(FID, '%s\n', breaker);
fprintf(FID, '%s\n', '> MATLAB script generating this file: ');
fprintf(FID, '%s\n', breaker);
Fself = fopen([mfilename('fullpath'),'.m'],'r');
while ~feof(Fself)
    tline = fgetl(Fself);
    fprintf(FID, '%s\n', tline);
end
fprintf(FID, '%s\n', breaker);
fprintf(FID, '%s\n', breaker);
fprintf(FID, '%s\n', breaker);
end

