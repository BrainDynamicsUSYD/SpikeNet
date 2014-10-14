function z=f_moving_average_weighted_para(y,x,o)    % past+future, inwards averaging at the edges (asymmetric), calculates many rows in parallel, order of moving average o = window length

if size(y,length(size(y)))==1 
    y=y';
end
if size(x,length(size(x)))==1 
    x=x';
end
num_rows=size(y,1);
num_data=size(y,2);
if num_data<o
    o=num_data;
end
h=fix(o/2);

if o>1
    z=zeros(num_rows,num_data);
    
    for cc=1:h
        z(1:num_rows,cc)=sum(y(1:num_rows,1:cc+h).*repmat(x(1:cc+h),num_rows,1),2)/sum(x(1:cc+h));
    end

    y2=zeros(num_rows,o,num_data-2*h);
    x2=zeros(o,num_data-2*h);
    for rc=1:o
        y2(1:num_rows,rc,1:num_data-2*h)=y(1:num_rows,rc-1+(1:num_data-2*h));
        x2(rc,1:num_data-2*h)=x(rc-1+(1:num_data-2*h));
    end
    z(1:num_rows,h+1:num_data-h)=shiftdim(permute(sum(y2.*permute(repmat(x2,[1,1,num_rows]),[3 1 2]),2),[2 1 3]),1)./repmat(sum(x2,1),num_rows,1);
    
    for cc=num_data-h+1:num_data
        z(1:num_rows,cc)=sum(y(1:num_rows,cc-h:num_data).*repmat(x(cc-h:num_data),num_rows,1),2)/sum(x(cc-h:num_data));
    end
else
    z=y; 
end
