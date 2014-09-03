function z=f_moving_average_weighted_para_p(y,x,o)    % past+future, inwards averaging at the edges (asymmetric), calculates many rows in parallel, order of moving average o = window length

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

if o>1
    z=zeros(num_rows,num_data);
    
    for cc=1:o-1
        z(1:num_rows,cc)=sum(y(1:num_rows,1:cc).*repmat(x(1:cc),num_rows,1),2)/sum(x(1:cc));
    end
    
    y2=zeros(num_rows,o,num_data-o+1);
    x2=zeros(o,num_data-o+1);
    for rc=1:o
        y2(1:num_rows,rc,1:num_data-o+1)=y(1:num_rows,rc-1+(1:num_data-o+1));
        x2(rc,1:num_data-o+1)=x(rc-1+(1:num_data-o+1));
    end
    z(1:num_rows,o:num_data)=shiftdim(permute(sum(y2.*permute(repmat(x2,[1,1,num_rows]),[3 1 2]),2),[2 1 3]),1)./repmat(sum(x2,1),num_rows,1);
else
    z=y; 
end