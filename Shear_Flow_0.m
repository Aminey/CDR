function [q, q01, q02, tau] = Shear_Flow_0(x, y, z, Ixx, Iyy, Ixy, Booms, Sx, Sy, A1, A2, qb, structure, Term_2, M_0, x_quarterchord, skin, str, caps, spar, dz);
%% Initialize

% dtheta_dz_1
delta_term_a1 = zeros(length(x),length(z),5);
delta_term_a2.sea = zeros(length(x),length(z),5);
delta_term_a2.alt = zeros(length(x),length(z),5);
delta_term_a3.sea = zeros(length(x),length(z),5);
delta_term_a3.alt = zeros(length(x),length(z),5);

% dtheta_dz 2
delta_term_b1 = zeros(length(x),length(z),5);
delta_term_b2.sea = zeros(length(x),length(z),5);
delta_term_b2.alt = zeros(length(x),length(z),5);
delta_term_b3.sea = zeros(length(x),length(z),5);
delta_term_b3.alt = zeros(length(x),length(z),5);

sum_term_a1 = zeros(length(z),5);
sum_term_a2.sea = zeros(length(z),5);
sum_term_a2.alt = zeros(length(z),5);
sum_term_a3.sea = zeros(length(z),5);
sum_term_a3.alt = zeros(length(z),5);

sum_term_b1 = zeros(length(z),5);
sum_term_b2.sea = zeros(length(z),5);
sum_term_b2.alt = zeros(length(z),5);
sum_term_b3.sea = zeros(length(z),5);
sum_term_b3.alt = zeros(length(z),5);

dtheta_dz_1.sea = sym(zeros(length(z),5));
dtheta_dz_1.alt = sym(zeros(length(z),5));

dtheta_dz_2.sea = sym(zeros(length(z),5));
dtheta_dz_2.alt = sym(zeros(length(z),5));

eq1.sea = sym(zeros(length(z),5));
eq1.alt = sym(zeros(length(z),5));

eq2.sea = sym(zeros(length(z),5));
eq2.alt = sym(zeros(length(z),5));

q01.sea = sym(zeros(length(z),5));
q02.sea = sym(zeros(length(z),5));
q01.alt = sym(zeros(length(z),5));
q02.alt = sym(zeros(length(z),5));

q.sea = zeros(length(x),length(z),5);
q.alt = zeros(length(x),length(z),5);
q.sea_spar = zeros(length(z),5);
q.alt_spar = zeros(length(z),5);

tau.sea = zeros(length(x),length(z),5);
tau.alt = zeros(length(x),length(z),5);
tau.sea_spar = zeros(length(x),length(z),5);
tau.alt_spar = zeros(length(x),length(z),5);

for k = 1:1
for j = 1:1
    
% q01_alt = sym(zeros(length(z),5));
% q01_sea = sym(zeros(length(z),5));
% q02_alt = sym(zeros(length(z),5));
% q02_sea = sym(zeros(length(z),5));
syms q01_alt
syms q01_sea
syms q02_alt
syms q02_sea
    
%% EQ 1
eq1.sea(j,k) = M_0.sea(k) + Sy.sea(j,k)*(x_quarterchord) - Sx.sea(j,k)*(0) - (2*A1*q01_sea + 2*A2*q02_sea + Term_2.sea(j,k));
eq1.alt(j,k) = M_0.alt(k) + Sy.alt(j,k)*(x_quarterchord) - Sx.alt(j,k)*(0) - (2*A1*q01_alt + 2*A2*q02_alt + Term_2.alt(j,k));
%% EQ 2

% dtheta_dz_1
m = 1;
for i = spar.i_CCW(2):spar.i_CCW(3)-1
    delta_term_a1(m,j,k) = (((x(i+1)-x(i))^2 + (y(i+1) - y(i))^2)^0.5)/skin.t;
    delta_term_a2.sea(m,j,k) = Booms.sea(i,j,k) * x(i) * (((x(i+1)-x(i))^2 + (y(i+1) - y(i))^2)^0.5)/skin.t;
    delta_term_a2.alt(m,j,k) = Booms.alt(i,j,k) * x(i) * (((x(i+1)-x(i))^2 + (y(i+1) - y(i))^2)^0.5)/skin.t;
    delta_term_a3.sea(m,j,k) = Booms.sea(i,j,k) * y(i) * (((x(i+1)-x(i))^2 + (y(i+1) - y(i))^2)^0.5)/skin.t;
    delta_term_a3.alt(m,j,k) = Booms.alt(i,j,k) * y(i) * (((x(i+1)-x(i))^2 + (y(i+1) - y(i))^2)^0.5)/skin.t;
    m = m+1;
end
sum_term_a1(j,k) = sum(delta_term_a1(:,j,k));
sum_term_a2.sea(j,k) = sum(delta_term_a2.sea(:,j,k));
sum_term_a2.alt(j,k) = sum(delta_term_a2.alt(:,j,k));
sum_term_a3.sea(j,k) = sum(delta_term_a3.sea(:,j,k));
sum_term_a3.alt(j,k) = sum(delta_term_a3.alt(:,j,k));

dtheta_dz_1.sea(j,k) = (1/(2*A1*structure.G)) * (q01_sea*(sum_term_a1(j,k)) + (q01_sea-q02_sea)*(caps.y_low_c(1) - caps.y_upp_c(1))/spar.t + ((Sy.sea(j,k)*Ixy - Sx.sea(j,k)*Ixx)/(Ixx*Iyy-Ixy^2))*(sum_term_a2.sea(j,k)) + ((Sx.sea(j,k)*Ixy-Sy.sea(j,k)*Iyy)/(Ixx*Iyy-Ixy^2))*(sum_term_a3.sea(j,k)));
dtheta_dz_1.alt(j,k) = (1/(2*A1*structure.G)) * (q01_alt*(sum_term_a1(j,k)) + (q01_alt-q02_alt)*(caps.y_low_c(1) - caps.y_upp_c(1))/spar.t + ((Sy.alt(j,k)*Ixy - Sx.alt(j,k)*Ixx)/(Ixx*Iyy-Ixy^2))*(sum_term_a2.alt(j,k)) + ((Sx.alt(j,k)*Ixy-Sy.alt(j,k)*Iyy)/(Ixx*Iyy-Ixy^2))*(sum_term_a3.alt(j,k)));

counter = 1;
for i = 1:spar.i_CCW(2)-1 %top rear skin
    delta_term_b1(i,j,k) = (((x(i+1)-x(i))^2 + (y(i+1) - y(i))^2)^0.5)/skin.t;
    delta_term_b2.sea(i,j,k) = Booms.sea(i,j,k) * x(i) * (((x(i+1)-x(i))^2 + (y(i+1) - y(i))^2)^0.5)/skin.t;
    delta_term_b2.alt(i,j,k) = Booms.alt(i,j,k) * x(i) * (((x(i+1)-x(i))^2 + (y(i+1) - y(i))^2)^0.5)/skin.t;
    delta_term_b3.sea(i,j,k) = Booms.sea(i,j,k) * y(i) * (((x(i+1)-x(i))^2 + (y(i+1) - y(i))^2)^0.5)/skin.t;
    delta_term_b3.alt(i,j,k) = Booms.alt(i,j,k) * y(i) * (((x(i+1)-x(i))^2 + (y(i+1) - y(i))^2)^0.5)/skin.t;
    counter = counter + 1;
end
for i = spar.i_CCW(3):length(x)-1 %bottom rear skin
    delta_term_b1(counter,j,k) = (((x(i+1)-x(i))^2 + (y(i+1) - y(i))^2)^0.5)/skin.t;
    delta_term_b2.sea(counter,j,k) = Booms.sea(i,j,k) * x(i) * (((x(i+1)-x(i))^2 + (y(i+1) - y(i))^2)^0.5)/skin.t;
    delta_term_b2.alt(counter,j,k) = Booms.alt(i,j,k) * x(i) * (((x(i+1)-x(i))^2 + (y(i+1) - y(i))^2)^0.5)/skin.t;
    delta_term_b3.sea(counter,j,k) = Booms.sea(i,j,k) * y(i) * (((x(i+1)-x(i))^2 + (y(i+1) - y(i))^2)^0.5)/skin.t;
    delta_term_b3.alt(counter,j,k) = Booms.alt(i,j,k) * y(i) * (((x(i+1)-x(i))^2 + (y(i+1) - y(i))^2)^0.5)/skin.t;
    counter = counter + 1;
end
    
sum_term_b1(j,k) = sum(delta_term_b1(:,j,k));
sum_term_b2.sea(j,k) = sum(delta_term_b2.sea(:,j,k));
sum_term_b2.alt(j,k) = sum(delta_term_b2.alt(:,j,k));
sum_term_b3.sea(j,k) = sum(delta_term_b3.sea(:,j,k));
sum_term_b3.alt(j,k) = sum(delta_term_b3.alt(:,j,k));

dtheta_dz_2.sea(j,k) = q02_sea*(sum_term_b1(j,k)) + (q02_sea-q01_sea)*(caps.y_upp_c(1)-caps.y_low_c(1))/spar.t + q02_sea*(caps.y_low_c(1) - caps.y_upp_c(1))/(spar.t) + ((Sy.sea(j,k)*Ixy - Sx.sea(j,k)*Ixx)/(Ixx*Iyy - Ixy^2))*(sum_term_b2.sea(j,k)) + ((Sx.sea(j,k)*Ixy - Sy.sea(j,k)*Iyy)/(Ixx*Iyy-Ixy^2))*(sum_term_b3.sea(j,k));
dtheta_dz_2.alt(j,k) = q02_alt*(sum_term_b1(j,k)) + (q02_alt-q01_alt)*(caps.y_upp_c(1)-caps.y_low_c(1))/spar.t + q02_alt*(caps.y_low_c(1) - caps.y_upp_c(1))/(spar.t) + ((Sy.alt(j,k)*Ixy - Sx.alt(j,k)*Ixx)/(Ixx*Iyy - Ixy^2))*(sum_term_b2.alt(j,k)) + ((Sx.alt(j,k)*Ixy - Sy.sea(j,k)*Iyy)/(Ixx*Iyy-Ixy^2))*(sum_term_b3.alt(j,k));

% Equating dtheta_dz_1 and _2
eq2.sea(j,k) = dtheta_dz_2.sea(j,k) - dtheta_dz_1.sea(j,k);
eq2.alt(j,k) = dtheta_dz_2.alt(j,k) - dtheta_dz_1.alt(j,k);

%% Solve
[q01.sea(j,k),q02.sea(j,k)] = solve([eq1.sea(j,k)==0,eq2.sea(j,k)==0], [q01_sea,q02_sea]);
[q01.alt(j,k),q02.alt(j,k)] = solve([eq1.alt(j,k)==0,eq2.alt(j,k)==0], [q01_alt,q02_alt]);

%% Shear Flow and Shear Stress Results
for i = 1:spar.i_CCW(2)-1
    q.sea(i,j,k) = qb.sea(i,j,k) + q02.sea(j,k);
    q.alt(i,j,k) = qb.alt(i,j,k) + q02.alt(j,k);
    tau.sea(i,j,k) = q.sea(i,j,k) / skin.t;
    tau.alt(i,j,k) = q.alt(i,j,k) / skin.t;
end
for i = spar.i_CCW(2):spar.i_CCW(3)-1
    q.sea(i,j,k) = qb.sea(i,j,k) + q01.sea(j,k);
    q.alt(i,j,k) = qb.alt(i,j,k) + q01.alt(j,k);
    tau.sea(i,j,k) = q.sea(i,j,k) / skin.t;
    tau.alt(i,j,k) = q.alt(i,j,k) / skin.t;
end
for i = spar.i_CCW(3):length(x)
    q.sea(i,j,k) = qb.sea(i,j,k) + q02.sea(j,k);
    q.alt(i,j,k) = qb.alt(i,j,k) + q02.alt(j,k);
    tau.sea(i,j,k) = q.sea(i,j,k) / skin.t;
    tau.alt(i,j,k) = q.alt(i,j,k) / skin.t;
end
    q.sea_spar(j,k) = q01.sea(j,k) - q02.sea(j,k);
    q.alt_spar(j,k) = q01.alt(j,k) - q02.alt(j,k);
    tau.sea_spar(j,k) = q.sea_spar(j,k) / spar.t;
    tau.alt_spar(j,k) = q.alt_spar(j,k) / spar.t;
    
disp('q01 = ');
disp(vpa(q01.sea(1,1)));
    
disp('q02 =');
disp(vpa(q02.sea(1,1)));

disp('dtheta_dz.sea(1,1) = ');
disp(vpa((1/(2*A1*structure.G)) * (q01.sea(1,1)*(sum_term_a1(1,1)) + (q01.sea(1,1)-q02.sea(1,1))*(caps.y_low_c(1) - caps.y_upp_c(1))/spar.t + ((Sy.sea(1,1)*Ixy - Sx.sea(1,1)*Ixx)/(Ixx*Iyy-Ixy^2))*(sum_term_a2.sea(1,1)) + ((Sx.sea(1,1)*Ixy-Sy.sea(1,1)*Iyy)/(Ixx*Iyy-Ixy^2))*(sum_term_a3.sea(1,1)))));

% figure;
% plot(z,q01.sea(:,:));
% xlabel('z')
% ylabel('q01.sea');

figure;
plot(1:length(x),tau.sea(:,1,1)/100000);
xlabel('x');
ylabel('Shear Stress (MPa)');


end
end