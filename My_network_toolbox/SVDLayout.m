function [x, y, z] = SVDLayout(A,varargin)
% A is adjacency matrix, can be weighted and directed. But it will be
% symmetrised.
%
% Based on Bernard's code.

%% Settings
% default settings
ShowPlots = 1;
ShowModularity = 1;
ShowDegree = 1;
ShowConnection = 0;
Method = 'Laps'; % SDE: Spectral distance embedding
ShowNode = 1;
DefaultNodeSize = 10;
DefaultNodeColor = [0 1 0];
maxMarkerSize = 20;
Modularity_threshold = 0.3; % Only above this value modularity will be shown
MaxConnection = 4e5; % Connection more than this will be slow to plot
NiceEdgeColor = [35 163 200]/255; % Nice blue-ish: [35 163 200]/255

% read settings
for i = 1:length(varargin)/2
    temp = varargin{2*i};
    eval([varargin{2*i-1}, '= temp;']);
end

%% Solve layout problem
N = size(A,1);% number of nodes
[uxA,wA,uvA] = svd(A);
svd_error = sum(sum((uxA*wA*uvA'-A).^2)) % check the accuracy of svd() solution!
        
%degree matrix
Deg = diag(sum(A));   %  sum each Col,  then place on diagonals of a matrix
if min(diag(Deg)) == 0
    disp('Note that there is zero(s) in Deg and thus NaN/Inf in degree-normalised symmetric Laplacian!');
end
        
switch Method
    case 'Laps'
        % Calcualte Lap, Lapn, Laps
        Iden = eye(N);  % Identity matrix
        %Laplacian matrix
        Lap = Deg - A;
        % % simple degree-normalized graph Laplacian - wrt row averageas - asymm.
        % % Lapn = inv(Deg)*A + Iden;  % nb. this is working from the "unsigned" def'n of L = D+A
        % Lapn = Deg\Lap; %inv(Deg)*Lap;  % use L=D-A (as above) - the is L-rw
        
        % symmetrise this - via deg-i, pre&post multiply :
        Degh = Iden;  Deginv = eye(N,N);% start w identity matrix
        for i=1:N
            Degh(i,i) = 1/sqrt(Deg(i,i));
            Deginv(i,i) = 1/(Deg(i,i));
        end
        Laps = (Degh*Lap)*Degh;  % this is Lap-sym
        
        % SVD with Laps
        [ux,w,uv] = svd(Laps);
        svd_error = sum(sum((ux*w*uv'-Laps).^2)) % check the accuracy of svd() solution!
        rank_laps = rank(Laps)
        NumZeroEigVal = N-rank_laps % number of zero eigenvalues, corresponding to number of components in the graph, usually one
        %[ux,w] = sortem(ux,w); % not necessary since svd outputs decreasing order
        
        % generate the 3D layout coords
        x = diag(Degh).*ux(:,N-NumZeroEigVal);  %  take smallest EigVals (here the later ones) drop all the zeros
        y = diag(Degh).*ux(:,N-NumZeroEigVal-1);  %  & normalise by  sqrt(deg-i) - not deg-i
        z = diag(Degh).*ux(:,N-NumZeroEigVal-2);
        
    case 'SDE'
        % SDE ( Ali Civril'06)
        %  uses Paths
        Paths = zeros(N); Pred = zeros(N);
        h_wb = waitbar(0,'Please wait for breath(A,i)...');
        for i=1:N
            [Paths(i,:), Pred(:,i)] = breadth(A,i);
            Paths(i,i) = 0;  % need to eliminate souce (spurious entry??)
            waitbar(i/N);
        end
        delete(h_wb);
        % and need CleanPaths to elim. Inf's
        onevec = ones(1,N);
        L = Paths.^2; % array of (d-ij)^2, graph dist's.
        gamma = eye(N)- onevec'*onevec/N;  % form projection matrix, to get centered coords (mean 0)
        M = -gamma*L*gamma/2;  %  normalised dist matrix (of graph path lengths)
        % SVD
        [ux,w,uv] = svd(M);
        svd_error = sum(sum((ux*w*uv'-M).^2)) % check the accuracy of svd() solution!
        % for g-unweighted:
        sl1 = sqrt(w(1,1)); sl2=sqrt(w(2,2)); sl3=sqrt(w(3,3)); % sqrt(Lamda-1) etc
        x = sl1*ux(:,1);  %  take largest EigVals
        y = sl2*ux(:,2);
        z = sl3*ux(:,3);
        clear onevec;      % saves space
        
end % switch method
Method

%% Plots
if ShowPlots == 1
    switch Method
        case 'Laps'
            % show Laps
            figure('Name','Spy(Ajacency Matrix)','NumberTitle','off');
            spy(Laps);
            title('Spy(), SymLap','FontWeight','bold','FontSize',14 )
    end
    
    % Show eigen spectrum
    figure('Name','Eigen Spectrum','NumberTitle','off');
    subplot(2,1,1);
    plot(1:length(diag(wA)),diag(wA),'x');legend('Eigen spectrum of adjcency matrix');
    subplot(2,1,2);
    plot(1:length(diag(w)),diag(w),'x');legend('Eigen spectrum of layout problem');
            
    % Show 3D layout
    % Prepare Modularyity
    if ShowModularity == 1
        [Ci, Q ] = modularity_newman_dir(A);
        % [Ci, Q ] = modularity_louvain_dir(A);
        % [Ci, Q ] = modularity_finetune_dir(A,Ci0);
        Modularity = Q
        if Q >= Modularity_threshold % If significant enough
            [~,A_re] = reorder_mod(A, Ci);
            figure('Name','Reordered Ajacency Matrix','NumberTitle','off');
            spy(A_re);
            ColorLabel = {'r.','b.','w.','g.','y.','c.','m.','r*','b*','w*','g*','y*','c*','m*','ro','bo','wo','go','yo','co','mo'};
            Num_of_modules = max(Ci)
            if length(ColorLabel) < Num_of_modules
                ShowModularity = 0;
            end
        else
            ShowModularity = 0;
        end
    end
    
    % Prepare large figure
    figwidth = 768;
    figheight = 768;
    figposition = [100, 100, figwidth, figheight];
    figure('Name','SVD Layout','NumberTitle','off','position',figposition, 'units','pixels');
    hold on;
    set(gcf,'Renderer','OpenGL');  % to exploit grahpics cards (using built-in OpenGL in graphic card, i.e., firmware), making things magnitudes faster
    view(59,32);
    set(gca, 'Color',[0,0,0]); % background black
    
    % Show Connections
    if nnz(A) > MaxConnection
        ShowConnection = 0;
        disp('Too many edges! Visualisation result will be compromised if plot all.');
    end
    if ShowConnection == 1
        % Plot the connections (if too many, then not working!)
        % Strategy: only show the important ones! (pairs with nodes of highest degree/betweenness centrality, etc. Usually no more than 10k pairs)
        % Use "patch" to gain more control. alpha = 0.03
        [ii,jj] = find(A);
        patchline([x(ii)'; x(jj)'],[y(ii)'; y(jj)'], [z(ii)'; z(jj)'], 'EdgeColor', NiceEdgeColor, 'LineSmoothing','on', 'EdgeAlpha', 0.02, 'LineWidth',0.1);
        % Line smoothing is anti-aliasing technique
    end
    
    % Prepare Degree
    DegVec = diag(Deg);
    
    % Plot the nodes
    if ShowNode == 1
        if ShowDegree == 1 && ShowModularity == 1
            for ii = 1:N
                plot3(x(ii),y(ii),z(ii),ColorLabel{Ci(ii)},'MarkerSize',maxMarkerSize*DegVec(ii)/max(DegVec))
            end
        elseif ShowDegree == 1 && ShowModularity == 0
            for ii = 1:N
                plot3(x(ii),y(ii),z(ii),'.','Color',DefaultNodeColor,'MarkerSize',maxMarkerSize*DegVec(ii)/max(DegVec))
            end
        elseif ShowDegree == 0 && ShowModularity == 1
            for ii = 1:N
                plot3(x(ii),y(ii),z(ii),ColorLabel{Ci(ii)},'MarkerSize',DefaultNodeSize)
            end
        elseif ShowDegree == 0 && ShowModularity == 0
            for ii = 1:N
                plot3(x(ii),y(ii),z(ii),'.','Color',DefaultNodeColor,'MarkerSize',DefaultNodeSize)
            end
        end
    end
    
    xlabel('EigVec X','FontWeight','bold','FontSize',14)
    ylabel('EigVec Y','FontWeight','bold','FontSize',14)
    zlabel('EigVec Z','FontWeight','bold','FontSize',14)
end % ShowPlots


end % SVDLayout







%% Supplementary functions
