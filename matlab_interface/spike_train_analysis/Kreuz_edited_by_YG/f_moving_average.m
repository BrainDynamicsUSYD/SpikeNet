function z=f_moving_average(y,o)    % past+future, inwards averaging at the edges (asymmetric), order of moving average o = window length

if size(y,length(size(y)))==1
    y=y'; 
end
num_data=length(y);
if num_data<o
    o=num_data;
end
h=fix(o/2);

if o>1
    z=zeros(1,num_data);
    
    for cc=1:h
        z(cc)=sum(y(1:cc+h))/(cc+h);
    end

    dummy=zeros(o,num_data-2*h);
    for rc=1:o
        dummy(rc,1:num_data-2*h)=y(rc-1+(1:num_data-2*h));
    end
    z(h+1:num_data-h)=mean(dummy);
    
    for cc=num_data-h+1:num_data
        z(cc)=sum(y(cc-h:num_data))/(h+1+(num_data-cc));
    end 
else
    z=y; 
end
