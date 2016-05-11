function [Term_2, A1, A2, A_total, qb] = Shear_Flow_Basic(x, y, z, Ixx, Iyy, Ixy, Booms, Sx, Sy, spar)

%% Initialize

delta_A = zeros(1,length(x));
delta_A1 = zeros(1,length(x));
delta_Term_2_alt = zeros(length(x),length(z),5);
delta_Term_2_sea = zeros(length(x),length(z),5);
Term_2.alt = zeros(length(z),5);
Term_2.sea = zeros(length(z),5);

qb.sea = zeros(length(x),length(z),5);
qb.alt = zeros(length(x),length(z),5);

x_temp = zeros(1,length(x));

%% Find incremental, cell, and total area

for i = 1:length(x)-1
    delta_A(i) = abs(x(i+1)*y(i) - y(i+1)*x(i))/2;
end
    
for i = 1:length(x)
    x_temp(i) = x(i) - spar.x(1); % temporary coordinate shift with x=0 on front spar.
                                  % need this to calculate area of each cell.
end

for i = spar.i_CCW(2):spar.i_CCW(3)-1
    delta_A1(i) = abs(x_temp(i+1)*y(i) - y(i+1)*x_temp(i))/2;
end

A_total = sum(delta_A);
A1 = sum(delta_A1);
A2 = A_total - A1;
%% BASIC SHEAR FLOW
for k = 1:5
    for j = 1:length(z)
        
            qb.sea(1,j,k) = (((Sx.sea(j,k)*Ixx - Sy.sea(j,k)*Ixy)/(Ixx*Iyy - Ixy))*Booms.sea(1,j,k)*x(i)) -...
                                         (((Sy.sea(j,k)*Iyy -Sx.sea(j,k)*Ixy)/(Ixx*Iyy-Ixy^2))*Booms.sea(1,j,k)*y(i));
            qb.alt(1,j,k) = (((Sx.alt(j,k)*Ixx - Sy.alt(j,k)*Ixy)/(Ixx*Iyy - Ixy))*Booms.alt(1,j,k)*x(i)) -...
                                         (((Sy.alt(j,k)*Iyy -Sx.alt(j,k)*Ixy)/(Ixx*Iyy-Ixy^2))*Booms.alt(1,j,k)*y(i));
                                     
            delta_Term_2_alt(1,j,k) = 2*delta_A(1)*qb.alt(1,j,k);
            delta_Term_2_sea(1,j,k) = 2*delta_A(1)*qb.sea(1,j,k);
            
        for i = 2:length(x)-1   %% for n number of points, there will be n flows, but the first is cut. 
                                 % So moving counterclockwise, you will
                                 % have point 1, flow 1, point 2, flow 2,
                                 % etc. since an i-1 term is utilized, the
                                 % i=1 term is initialized, and the loop
                                 % starts at i = 2.
            qb.sea(i,j,k) = qb.sea(i-1,j,k) + (((Sx.sea(j,k)*Ixx - Sy.sea(j,k)*Ixy)/(Ixx*Iyy - Ixy))*Booms.sea(i,j,k)*x(i)) -...
                                         (((Sy.sea(j,k)*Iyy -Sx.sea(j,k)*Ixy)/(Ixx*Iyy-Ixy^2))*Booms.sea(i,j,k)*y(i));
            qb.alt(i,j,k) = qb.alt(i-1,j,k) + (((Sx.alt(j,k)*Ixx - Sy.alt(j,k)*Ixy)/(Ixx*Iyy - Ixy))*Booms.alt(i,j,k)*x(i)) -...
                                         (((Sy.alt(j,k)*Iyy -Sx.alt(j,k)*Ixy)/(Ixx*Iyy-Ixy^2))*Booms.alt(i,j,k)*y(i));
            delta_Term_2_alt(i,j,k) = 2*delta_A(i)*qb.alt(i,j,k);
            delta_Term_2_sea(i,j,k) = 2*delta_A(i)*qb.sea(i,j,k);           
        end
        Term_2.alt(j,k) = sum(delta_Term_2_alt(:,j,k));
        Term_2.sea(j,k) = sum(delta_Term_2_sea(:,j,k));
       
    end
end


