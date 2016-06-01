function [syn_cell] = read_ygin_syn(filename, pop_pre, pop_post)
% refer to the function writeChemicalConnection

syn_cell = {};

FID = fopen(filename,'r');

while ~feof(FID)
    tline = fgetl(FID);
    if isempty(tline)
        continue;
    elseif strcmp(tline(1), '>')

        if strfind(tline,'INIT006') % fprintf(FID, '%s\n', '> INIT006');
            tline = fgetl(FID);
            scan_temp = textscan(tline,'%f','Delimiter',',');
            syn.type = scan_temp{1}(1)+1; % be careful here!
            syn.pop_pre = scan_temp{1}(2)+1; % be careful here!
            syn.pop_post = scan_temp{1}(3)+1; % be careful here!
            
            if syn.pop_pre == pop_pre && syn.pop_post == pop_post
                fprintf('\t Reading connectivity from pop %d to pop %d...\n', syn.pop_pre, syn.pop_post);
                
                tline = fgetl(FID);
                scan_temp = textscan(tline,'%f','Delimiter',',');
                syn.I = scan_temp{1} + 1;
                
                tline = fgetl(FID);
                scan_temp = textscan(tline,'%f','Delimiter',',');
                syn.J = scan_temp{1} + 1;
                
                tline = fgetl(FID);
                scan_temp = textscan(tline,'%f','Delimiter',',');
                syn.K = scan_temp{1};
                
                tline = fgetl(FID);
                scan_temp = textscan(tline,'%f','Delimiter',',');
                syn.D = scan_temp{1};
                
                if isempty(syn_cell)
                    syn_cell = {syn};
                else
                    syn_cell = [syn_cell, syn];
                end
                
            end
            

        end
        
    end
end
           
fclose(FID);
               

end