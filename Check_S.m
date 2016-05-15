function [shear_check] = Check_S(x, y, z, dz, q)

%%Initialize

delta_x_sea = zeros(1,length(x));
delta_y_sea = zeros(1,length(x));

%% Use shear_check.delta_x to check against Sx, and _y to Sy.
for k = 1:5
    for j = 1:dz:length(z) 
        for i = 1:length(x)-1
           delta_x_sea(i,j,k) = q.sea(i,j,k)*(x(i+1)-x(i));
           delta_y_sea(i,j,k) = q.sea(i,j,k)*(y(i+1)-y(i));
        end
        shear_check.x_sea(j,k) = sum(delta_x_sea(:,j,k));
        shear_check.y_sea(j,k) = sum(delta_y_sea(:,j,k));
    end
end



