function Vini = get_initialV(Result_cell,varargin)

% Get initial membrane potential
% please make sure the main function sampled the initial V
disp('Get initial membrane potential...');

Result_num = length(Result_cell);

if nargin == 1
    T = 1;
else
    T = varargin{1};
end

for r_num = 1:Result_num
    R = Result_cell{r_num};
%     R = Result_cell;    
    Num_pop = R.Num_pop;    
    Vini = cell(Num_pop,1);
    
    for pop=1:Num_pop
        file_name = [R.stamp(1:end-3), num2str(pop-1), '_neurosamp.mat'];
        R_samp = load(file_name);
        Vini{pop}=R_samp.V(:,T);        
    end    
    V=Vini{1};
    save([R.stamp(1:4),'_init_V_1.mat'],'V')
    V=Vini{2};
    save([R.stamp(1:4),'_init_V_2.mat'],'V')
end 
end
