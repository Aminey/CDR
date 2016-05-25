function [shear_check] = Check_S(x, y, z, dz, q, q01, q02, Sx, Sy, spar)

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
        
        delta_x_sea(length(x),j,k) = q.sea(length(x),j,k) * (x(1) - x(length(x)));      % should be zero, cuz x(1) = x(end)
        delta_y_sea(length(x),j,k) = q.sea(length(x),j,k) * (y(1) - y(length(x)));    
        
        shear_check.x_sea(j,k) = sum(delta_x_sea(:,j,k));
        shear_check.y_sea(j,k) = sum(delta_y_sea(:,j,k)) + (q.sea_spar(j,k))*spar.h(1) + (q02.sea(j,k)*spar.h(2));
    end
end

disp('Sx.sea(1,1) = ');
disp(Sx.sea(1,1));
disp('shear_check.x_sea(1,1) = ');
disp(vpa(shear_check.x_sea(1,1)));
disp('Sx percent error = ');
disp(100 * abs((Sx.sea(1,1)-shear_check.x_sea(1,1))/Sx.sea(1,1)));

disp('Sy.sea(1,1) = ');
disp(Sy.sea(1,1));
disp('shear_check.y_sea(1,1) = ');
disp(vpa(shear_check.y_sea(1,1)));
disp('Sy percent error = ');
disp(vpa(100 * abs((Sy.sea(1,1)-shear_check.y_sea(1,1))/Sy.sea(1,1))));
