[Term_1] = Shear_Flow_0(x, y, z, Booms, Sx, Sy, A1, A2, qb, G, Term_2, M_0, x_quarterchord, skin, str, caps, spar);

% check what thickness should be

%% Initialize
syms q01.alt;
syms q01.sea;
syms q02.alt;
syms q02.sea;

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

qtot.sea = zeros(1,length(x)+1);           % the +1 is for the spar
qtot.alt = zeros(1,length(x)+1);

sum_term_a1(j,k) = zeros(length(z),5);
sum_term_a2.sea(j,k) = zeros(length(z),5);
sum_term_a2.alt(j,k) = zeros(length(z),5);
sum_term_a3.sea(j,k) = zeros(length(z),5);
sum_term_a3.alt(j,k) = zeros(length(z),5);

sum_term_b1(j,k) = zeros(length(z),5);
sum_term_b2.sea(j,k) = zeros(length(z),5);
sum_term_b2.alt(j,k) = zeros(length(z),5);
sum_term_b3.sea(j,k) = zeros(length(z),5);
sum_term_b3.alt(j,k) = zeros(length(z),5);

dtheta_dz_1.sea(j,k) = zeros(length(z),5);
dtheta_dz_1.alt(j,k) = zeros(length(z),5);

dtheta_dz_2.sea(j,k) = zeros(length(z),5);
dtheta_dz_2.alt(j,k) = zeros(length(z),5);

eq1.sea(j,k) = zeros(length(z),5);
eq1.alt(j,k) = zeros(length(z),5);

eq2.sea(j,k) = zeros(length(z),5);
eq2.alt(j,k) = zeros(length(z),5);

q.sea(i,j,k) = zeros(length(x),length(z),5);
q.alt(i,j,k) = zeros(length(x),length(z),5);
q.sea_spar(j,k) = zeros(length(z),5);
q.alt_spar(j,k) = zeros(length(z),5);

for k = 1:5
for j = 1:length(z)
    
%% EQ 1
eq1.sea(j,k) = M_0.sea(k) + Sy.sea(j,k)*(x_quarterchord) - Sx.sea(j,k)*(0) - (2*A1*q01.sea(j,k) + 2*A2*q02.sea(j,k) + Term_2.sea(j,k));
eq1.alt(j,k) = M_0.alt(k) + Sy.alt(j,k)*(x_quarterchord) - Sx.alt(j,k)*(0) - (2*A1*q01.alt(j,k) + 2*A2*q02.alt(j,k) + Term_2.alt(j,k));

%% EQ 2

% dtheta_dz_1
m = 1;
for i = spar.i_CCW(2):spar.i_CCW(3)-1
    delta_term_a1(m,j,k) = (((x(i+1)-x(i))^2 + (y(i+1) - y(i))^2)^0.5)/skin.t;
    delta_term_a2.sea(m,j,k) = Booms.sea(i,j,k) * x(i) * (((x(i+1)-x(i))^2 + (y(i+1) - y(i))^2)^0.5)/thickness(i);
    delta_term_a2.alt(m,j,k) = Booms.alt(i,j,k) * x(i) * (((x(i+1)-x(i))^2 + (y(i+1) - y(i))^2)^0.5)/thickness_1;
    delta_term_a3.sea(m,j,k) = Booms.sea(i,j,k) * y(i) * (((x(i+1)-x(i))^2 + (y(i+1) - y(i))^2)^0.5)/thickness(i);
    delta_term_a3.alt(m,j,k) = Booms.alt(i,j,k) * y(i) * (((x(i+1)-x(i))^2 + (y(i+1) - y(i))^2)^0.5)/thickness_1;
    m = m+1;
end
sum_term_a1(j,k) = sum(delta_term_a1(:,j,k));
sum_term_a2.sea(j,k) = sum(delta_term_a2.sea(:,j,k));
sum_term_a2.alt(j,k) = sum(delta_term_a2.alt(:,j,k));
sum_term_a3.sea(j,k) = sum(delta_term_a3.sea(:,j,k));
sum_term_a3.alt(j,k) = sum(delta_term_a3.alt(:,j,k));

dtheta_dz_1.sea(j,k) = (1/(2*A1*G)) * (q01.sea(j,k)*(sum_term_a1(j,k)) + (q01.sea(j,k)-q02.sea(j,k))*(caps.y_low - y_top)/spar.t + ((Sy.sea(j,k)*Ixy - Sx.sea(j,k)*Ixx)/(Ixx*Iyy-Ixy^2))*(sum_term_a2.sea(j,k)) + ((Sx.sea(j,k)*Ixy-Sy.sea(j,k)*Iyy)/(Ixx*Iyy-Ixy^2))*(sum_term_a3.sea(j,k)));
dtheta_dz_1.alt(j,k) = (1/(2*A1*G)) * (q01.sea(j,k)*(sum_term_a1(j,k)) + (q01.alt(j,k)-q02.alt(j,k))*(caps.y_low - y_top)/spar.t + ((Sy.alt(j,k)*Ixy - Sx.alt(j,k)*Ixx)/(Ixx*Iyy-Ixy^2))*(sum_term_a2.alt(j,k)) + ((Sx.alt(j,k)*Ixy-Sy.alt(j,k)*Iyy)/(Ixx*Iyy-Ixy^2))*(sum_term_a3.alt(j,k)));

counter = 1;
for i = 1:spar.i_CCW(2)-1 %top rear skin
    delta_term_b1(i,j,k) = (((x(i+1)-x(i))^2 + (y(i+1) - y(i))^2)^0.5)/skin.t;
    delta_term_b2.sea(i,j,k) = Booms.sea(i,j,k) * x(i) * (((x(i+1)-x(i))^2 + (y(i+1) - y(i))^2)^0.5)/thickness(i);
    delta_term_b2.alt(i,j,k) = Booms.alt(i,j,k) * x(i) * (((x(i+1)-x(i))^2 + (y(i+1) - y(i))^2)^0.5)/thickness_1;
    delta_term_b3.sea(i,j,k) = Booms.sea(i,j,k) * y(i) * (((x(i+1)-x(i))^2 + (y(i+1) - y(i))^2)^0.5)/thickness(i);
    delta_term_b3.alt(i,j,k) = Booms.alt(i,j,k) * y(i) * (((x(i+1)-x(i))^2 + (y(i+1) - y(i))^2)^0.5)/thickness_1;
    counter = counter + 1;
end
for i = spar.i_CCW(3):length(x) %bottom rear skin
    delta_term_b1(counter,j,k) = (((x(i+1)-x(i))^2 + (y(i+1) - y(i))^2)^0.5)/skin.t;
    delta_term_b2.sea(counter,j,k) = Booms.sea(i,j,k) * x(i) * (((x(i+1)-x(i))^2 + (y(i+1) - y(i))^2)^0.5)/thickness(i);
    delta_term_b2.alt(counter,j,k) = Booms.alt(i,j,k) * x(i) * (((x(i+1)-x(i))^2 + (y(i+1) - y(i))^2)^0.5)/thickness_1;
    delta_term_b3.sea(counter,j,k) = Booms.sea(i,j,k) * y(i) * (((x(i+1)-x(i))^2 + (y(i+1) - y(i))^2)^0.5)/thickness(i);
    delta_term_b3.alt(counter,j,k) = Booms.alt(i,j,k) * y(i) * (((x(i+1)-x(i))^2 + (y(i+1) - y(i))^2)^0.5)/thickness_1;
    counter = counter + 1;
end
    
sum_term_b1(j,k) = sum(delta_term_b1(:,j,k));
sum_term_b2.sea(j,k) = sum(delta_term_b2.sea(:,j,k));
sum_term_b2.alt(j,k) = sum(delta_term_b2.alt(:,j,k));
sum_term_b3.sea(j,k) = sum(delta_term_b3.sea(:,j,k));
sum_term_b3.alt(j,k) = sum(delta_term_b3.alt(:,j,k));

dtheta_dz_2.sea(j,k) = q02.sea(j,k)*(sum_term_b1(j,k)) + (q02.sea(j,k)-q01.sea(j,k))*(y_top-t_bottom)/spar.t + q02.sea(j,k)*(y_bot-y_top)/(skin.thickness + spar.thickness) + ((Sy.sea(j,k)*Ixy - Sx.sea(j,k)*Ixx)/(Ixx*Iyy - Ixy^2))*(sum_term_b2.sea(j,k)) + ((Sx.sea(j,k)*Ixy - Sy.sea(j,k)*Iyy)/(Ixx*Iyy-Ixy^2))*(sum_term_b3.sea(j,k));
dtheta_dz_2.alt(j,k) = q02.alt(j,k)*(sum_term_b1(j,k)) + (q02.alt(j,k)-q01.alt(j,k))*(y_top-t_bottom)/spar.t + q02.alt(j,k)*(y_bot-y_top)/(skin.thickness + spar.thickness) + ((Sy.sea(j,k)*Ixy - Sx.sea(j,k)*Ixx)/(Ixx*Iyy - Ixy^2))*(sum_term_b2.alt(j,k)) + ((Sx.sea(j,k)*Ixy - Sy.sea(j,k)*Iyy)/(Ixx*Iyy-Ixy^2))*(sum_term_b3.alt(j,k));

% Equating dtheta_dz_1 and _2
eq2.sea(j,k) = dtheta_dz_2.sea(j,k) - dtheta_dz_1.sea(j,k);
eq2.alt(j,k) = dtheta_dz_2.alt(j,k) - dtheta_dz_1.alt(j,k);

%% Solve
[q01.sea(j,k),q02.sea(j,k)] = solve(eq1.sea(j,k)==0,eq2.sea(j,k)==0);
[q01.alt(j,k),q02.alt(j,k)] = solve(eq1.alt(j,k)==0,eq2.alt(j,k)==0);

for i = 1:spar.i_CCW(2)-1
    q.sea(i,j,k) = qb.sea(i,j,k) + q02.sea(j,k);
    q.alt(i,j,k) = qb.alt(i,j,k) + q02.alt(j,k);
end
for i = spar.i_CCW(2):spar.i_CCW(3)-1
    q.sea(i,j,k) = qb.sea(i,j,k) + q01.sea(j,k);
    q.alt(i,j,k) = qb.alt(i,j,k) + q01.alt(j,k);
end
for i = spar.i_CCW(3):length(x)
    q.sea(i,j,k) = qb.sea(i,j,k) + q01.sea(j,k);
    q.alt(i,j,k) = qb.alt(i,j,k) + q01.alt(j,k);
end
    q.sea_spar(j,k) = q01.sea(j,k) - q02.sea(j,k);
    q.alt_spar(j,k) = q01.alt(j,k) - q02.alt(j,k);

% % shear stress tau
% % please calculate the shear stress from the shear flow
end
end