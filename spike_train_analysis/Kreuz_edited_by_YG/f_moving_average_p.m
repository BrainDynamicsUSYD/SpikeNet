function z=f_moving_average_p(y,o)    % only past (as far as possible at the beginning), order of moving average o = window length

if size(y,length(size(y)))==1
    y=y'; 
end
num_data=length(y);
if num_data<o
    o=num_data;
end

if o>1
    z=zeros(1,num_data);
    
    for cc=1:o-1
        z(cc)=sum(y(1:cc))/cc;
    end

    dummy=zeros(o,num_data-o+1);
    for rc=1:o
        dummy(rc,1:num_data-o+1)=y(rc-1+(1:num_data-o+1));
    end
    z(o:num_data)=mean(dummy);
else
    z=y; 
end
