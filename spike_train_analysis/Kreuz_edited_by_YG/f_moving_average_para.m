function z=f_moving_average_para(y,o)    % past+future, inwards averaging at the edges (asymmetric), calculates many rows in parallel, order of moving average o = window length

if size(y,length(size(y)))==1
    y=y';
end
num_rows=size(y,1);
num_data=size(y,2);
if num_data<o
    o=num_data;
end
h=fix(o/2);

if o>1
    z=zeros(num_rows,num_data);
    
    for cc=1:h+1
        z(1:num_rows,cc)=sum(y(1:num_rows,1:cc+h),2)/(cc+h);
    end
    
    for cc=h+2:num_data-h
        z(1:num_rows,cc)=z(1:num_rows,cc-1)+( y(1:num_rows,cc+h) - y(1:num_rows,cc-h-1) )/o;
        %yy=z(1:num_rows,cc)
    end
    
    for cc=num_data-h+1:num_data
        z(1:num_rows,cc)=sum(y(1:num_rows,cc-h:num_data),2)/(h+1+(num_data-cc));
    end
else
    z=y;
end

