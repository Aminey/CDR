%% Initialize
x = zeros(1,length(airfoil.booms));
y = zeros(1,length(airfoil.booms));
delta_A = zeros(1,length(airfoil.booms));
delta_Term_2_alt = zeros(length(x),length(z),5);
delta_Term_2_sea = zeros(length(x),length(z),5);
Term_2.alt = zeros(length(z),5);
Term_2.sea = zeros(length(z),5);

for i = 1:length(airfoil.booms)
    x(i) = airfoil.booms(i).x_coordinate;                %%
    y(i) = airfoil.booms(i).y_coordinate;
end

qb_sea = zeros(length(x),length(z),5);
qb_alt = zeros(length(x),length(z),5);

%% Find incremental and total area
for i = 1:length(x)-1
    delta_A(i) = abs(x(i+1)*y(i) - y(i+1)*x(i))/2;
end

A_total = sum(delta_A);

%% BASIC SHEAR FLOW
for k = 1:5
    for j = 1:length(z)
            qb_sea(1,j,k) = (((Sx.sea(j,k)*structure.inertias(1) - Sy.sea(j,k)*stucture.inertias(3))/...
                                         (structure.inertias(1)*structure.inertias(2) - structure.inertias(3)))...
                                         *Booms.sea(1,j,k)*x(i)) - (((Sy.sea(j,k)*structure.inertias(2) -...
                                         Sx.sea(j,k)*structure.inertias(3))/(structure.inertias(1)*structure.inertias(2)...
                                         -structure.inertias(3)))*Booms.sea(1,j,k)*y(i));
            qb_alt(1,j,k) = (((Sx.alt(j,k)*structure.inertias(1) - Sy.alt(j,k)*stucture.inertias(3))/...
                                         (structure.inertias(1)*structure.inertias(2) - structure.inertias(3)))...
                                         *Booms.alt(1,j,k)*x(i)) - (((Sy.alt(j,k)*structure.inertias(2) -...
                                         Sx.alt(j,k)*structure.inertias(3))/(structure.inertias(1)*structure.inertias(2)...
                                         -structure.inertias(3)))*Booms.alt(1,j,k)*y(i));
            delta_Term_2_alt(1,j,k) = 2*delta_A(1)*qb_alt(1,j,k);
            delta_Term_2_sea(1,j,k) = 2*delta_A(1)*qb_sea(1,j,k);
        for i = 2:length(x)-1   %% for n number of points, there will be n flows, but the first is cut. 
                                 % So moving counterclockwise, you will
                                 % have point 1, flow 1, point 2, flow 2,
                                 % etc. since an i-1 term is utilized, the
                                 % i=1 term is initialized, and the loop
                                 % starts at i = 2.
            qb_sea(i,j,k) = qb_sea(i-1,j,k) + (((Sx.sea(j,k)*structure.inertias(1) - Sy.sea(j,k)*stucture.inertias(3))/...
                                         (structure.inertias(1)*structure.inertias(2) - structure.inertias(3)))...
                                         *Booms.sea(i,j,k)*x(i)) - (((Sy.sea(j,k)*structure.inertias(2) -...
                                         Sx.sea(j,k)*structure.inertias(3))/(structure.inertias(1)*structure.inertias(2)...
                                         -structure.inertias(3)))*Booms.sea(i,j,k)*y(i));
            qb_alt(i,j,k) = qb_alt(i-1,j,k) + (((Sx.alt(j,k)*structure.inertias(1) - Sy.alt(j,k)*stucture.inertias(3))/...
                                         (structure.inertias(1)*structure.inertias(2) - structure.inertias(3)))...
                                         *Booms.alt(i,j,k)*x(i)) - (((Sy.alt(j,k)*structure.inertias(2) -...
                                         Sx.alt(j,k)*structure.inertias(3))/(structure.inertias(1)*structure.inertias(2)...
                                         -structure.inertias(3)))*Booms.alt(i,j,k)*y(i));
            delta_Term_2_alt(i,j,k) = 2*delta_A(i)*qb_alt(i,j,k);
            delta_Term_2_sea(i,j,k) = 2*delta_A(i)*qb_sea(i,j,k);           
        end
        Term_2.alt(j,k) = sum(delta_Term_2_alt(:,j,k));
        Term_2.sea(j,k) = sum(delta_Term_2_sea(:,j,k));
       
    end
end


