function DM=GridDM(InGrid, OutGrid, pbc)
% Calculates the distance between all pairs of points from one grid to
% another.
% Grid is a vector of [ no rows, no columns, grid step size]
NoIn=InGrid(1)*InGrid(2);
NoOut=OutGrid(2)*OutGrid(2);

DM=zeros(NoIn,NoOut);

for i=0:NoIn-1
    % Convert 1D index to 2D spat coords
    yi=floor(i/InGrid(1));
    xi=i-yi*InGrid(1);
    
    % Scale for grid step size
    yi=yi*InGrid(3);
    xi=xi*InGrid(3);
    
    for j=0:NoOut-1
        % Convert 1D index to 2D spat coords
        yj=floor(j/OutGrid(1));
        xj=j-yj*OutGrid(1);

        % Scale for grid step size
        yj=yj*OutGrid(3);
        xj=xj*OutGrid(3);
        
        delx=abs(xi-xj);
        dely=abs(yi-yj);
        
        if pbc
            %If grid has periodic boundaries then check to see which
            % distance is shortest (i.e via the periodic boundary or not
            delx=min([delx,abs(InGrid(2)*InGrid(3)-delx)]);
            dely=min([dely,abs(InGrid(1)*InGrid(3)-dely)]);
        end
        
        DM(i+1,j+1)=sqrt(delx^2+dely^2);
        
    end 
end
end
