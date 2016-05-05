%% Basic Shear Flow qb

x = zeros(1,length(airfoil.booms));
y = zeros(1,length(airfoil.booms));

for i = 1:length(airfoil.booms)
    x(i) = airfoil.booms(i).x_coordinate;                %%
    y(i) = airfoil.booms(i).y_coordinate;
end

qb_sea = zeros(length(x),length(z),5);
qb_alt = zeros(length(x),length(z),5);

for k = 1:5
    for j = 1:length(z)
        for i = 1:length(x)-1   %% for n number of points, there will be n flows, but the first is cut. 
                                 % So moving counterclockwise, you will have point 1, flow 1, point 2, flow 2, ...
            qb_sea(i,j,k) = (((Sx.sea(j,k)*structure.inertias(1) - Sy.sea(j,k)*stucture.inertias(3))/...
                                         (structure.inertias(1)*structure.inertias(2) - structure.inertias(3)))...
                                         *BOOM_A_SEA(i,j,k)*x(i)) - (((Sy.sea(j,k)*structure.inertias(2) -...
                                         Sx.sea(j,k)*structure.inertias(3))/(structure.inertias(1)*structure.inertias(2)...
                                         -structure.inertias(3)))*BOOM_A_SEA(i,j,k)*y(i));
            qb_alt(i,j,k) = (((Sx.alt(j,k)*structure.inertias(1) - Sy.alt(j,k)*stucture.inertias(3))/...
                                         (structure.inertias(1)*structure.inertias(2) - structure.inertias(3)))...
                                         *BOOM_A_ALT(i,j,k)*x(i)) - (((Sy.alt(j,k)*structure.inertias(2) -...
                                         Sx.alt(j,k)*structure.inertias(3))/(structure.inertias(1)*structure.inertias(2)...
                                         -structure.inertias(3)))*BOOM_A_ALT(i,j,k)*y(i));
       
        end
    end
end


