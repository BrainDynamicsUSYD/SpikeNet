function z=f_moving_average_weighted_p(y,x,o)    % past+future, inwards averaging at the edges (asymmetric), order of moving average o = window length

if size(y,length(size(y)))==1 
    y=y';
end
if size(x,length(size(x)))==1 
    x=x';
end

num_data=length(y);
if num_data<o
    o=num_data;
end

if o>1
    z=zeros(1,num_data);
    
    for cc=1:o-1
        z(cc)=sum(y(1:cc).*x(1:cc))/sum(x(1:cc));
    end
    
    y2=zeros(o,num_data-o+1);
    x2=zeros(o,num_data-o+1);
    for rc=1:o
        y2(rc,1:num_data-o+1)=y(rc-1+(1:num_data-o+1));
        x2(rc,1:num_data-o+1)=x(rc-1+(1:num_data-o+1));
    end
    z(o:num_data)=sum(y2.*x2)./sum(x2);
else
    z=y; 
end