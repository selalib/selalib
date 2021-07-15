function eval = evaluate_splines_grid_points( dofs, degree )

n = size(dofs);
eval = zeros(n(1),n(2));
for i=1:n(1)
    
    if ( degree == 1 )
        eval(i,2:n(2)) = dofs(i,1:n(2)-1);
        eval(i,1) = dofs(i,n(2));
    elseif (degree == 2) 
        eval(i,1) = 0.5 * ( dofs(i,n(2)-1) + dofs(i,n(2)) );
        eval(i,2) = 0.5 * ( dofs(i,1) + dofs(i,n(2)) );
        for j=3:n(2)
            eval(i,j) = 0.5*( dofs(i,j-2)+dofs(i,j-1) );
        end
    elseif (degree == 3) 
        eval(i,1) = (dofs(i,n(2)-2)+dofs(i,n(2)))/6+2/3*dofs(i,n(2)-1);
        eval(i,2) = (dofs(i,n(2)-1)+dofs(i,1))/6+2/3*dofs(i,n(2));
        eval(i,3) = (dofs(i,n(2))+dofs(i,2))/6+2/3*dofs(i,1);
        for j=4:n(2)
            eval(i,j) = (dofs(i,j-3)+dofs(i,j-1))/6+2/3*dofs(i,j-2);
        end
        
    else
        disp('Not implemented.');
    end
end