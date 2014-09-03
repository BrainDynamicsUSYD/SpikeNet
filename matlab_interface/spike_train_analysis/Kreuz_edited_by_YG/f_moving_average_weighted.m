function z=f_moving_average_weighted(y,x,o)    % past+future, inwards averaging at the edges (asymmetric), order of moving average o = window length

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
h=fix(o/2);

if o>1
    z=zeros(1,num_data);
    
    for cc=1:h
        z(cc)=sum(y(1:cc+h).*x(1:cc+h))/sum(x(1:cc+h));
    end

    y2=zeros(o,num_data-2*h);
    x2=zeros(o,num_data-2*h);
    for rc=1:o
        y2(rc,1:num_data-2*h)=y(rc-1+(1:num_data-2*h));
        x2(rc,1:num_data-2*h)=x(rc-1+(1:num_data-2*h));
    end
    z(h+1:num_data-h)=sum(y2.*x2)./sum(x2,1);
    
    for cc=num_data-h+1:num_data
        z(cc)=sum(y(cc-h:num_data).*x(cc-h:num_data))/sum(x(cc-h:num_data));
    end
else
    z=y; 
end
