[Term_1] = Shear_Flow_0(x, y, z, Booms, Sx, Sy, A1, A2, qb, Term_2, M_0, x_quarterchord);

%% Initialize

syms q01.alt;
syms q01.sea;
syms q02.alt;
syms q02.sea;


for k = 1:5
for j = 1:length(z)
    
    
%% EQ 1

eq1.sea = M_0.sea + Sy.sea(j,k)*(x_quarterchord) - Sx.sea(j,k)*(0) - (2*A1*q01.sea(j,k) + 2*A2*q02.sea(j,k) + Term_2.sea);
eq1.alt = M_0.alt + Sy.alt(j,k)*(x_quarterchord) - Sx.alt(j,k)*(0) - (2*A1*q01.alt(j,k) + 2*A2*q02.alt(j,k) + Term_2.alt);

%% EQ 2
% dtheta_dz_1
delta_term_a1 = zeros(1,length(x));
delta_term_a2.sea = zeros(1,length(x));
delta_term_a2.alt = zeros(1,length(x));
delta_term_a3.sea = zeros(1,length(x));
delta_term_a3.alt = zeros(1,length(x));

for i = spar.i_CCW(2):spar.i_CCW(3)-1
    m = 1;
    delta_term_a1(m) = (((x(i+1)-x(i))^2 + (y(i+1) - y(i))^2)^0.5)/skin.t;
    delta_term_a2.sea(m) = Booms.sea * x(i) * (((x(i+1)-x(i))^2 + (y(i+1) - y(i))^2)^0.5)/thickness(i);
    delta_term_a2.alt(m) = Booms.alt * x(i) * (((x(i+1)-x(i))^2 + (y(i+1) - y(i))^2)^0.5)/thickness_1;
    delta_term_a3.sea(m) = Booms.sea * y(i) * (((x(i+1)-x(i))^2 + (y(i+1) - y(i))^2)^0.5)/thickness(i);
    delta_term_a3.alt(m) = Booms.alt * y(i) * (((x(i+1)-x(i))^2 + (y(i+1) - y(i))^2)^0.5)/thickness_1;
    m = m+1;
end
sum_term_a1 = sum(delta_term_a1(:));
sum_term_a2.sea = sum(delta_term_a2.sea(:));
sum_term_a2.alt = sum(delta_term_a2.alt(:));
sum_term_a3.sea = sum(delta_term_a3.sea(:));
sum_term_a3.alt = sum(delta_term_a3.alt(:));

dtheta_dz_1.sea = (1/(2*A1*G)) * (q01.sea*(sum_term_a1) + (q01.sea-q02.sea)*(y_bot - y_top)/spar.t + ((Sy.sea*Ixy - Sx.sea*Ixx)/(Ixx*Iyy-Ixy^2))*(sum_term_a2.sea) + ((Sx.sea*Ixy-Sy.sea*Iyy)/(Ixx*Iyy-Ixy^2))*(sum_term_a3.sea));
dtheta_dz_1.alt = (1/(2*A1*G)) * (q01.sea*(sum_term_a1) + (q01.alt-q02.alt)*(y_bot - y_top)/spar.t + ((Sy.alt*Ixy - Sx.alt*Ixx)/(Ixx*Iyy-Ixy^2))*(sum_term_a2.alt) + ((Sx.alt*Ixy-Sy.alt*Iyy)/(Ixx*Iyy-Ixy^2))*(sum_term_a3.alt));

% dtheta_dz 2
delta_term_b1 = zeros(1,length(x));
delta_term_b2.sea = zeros(1,length(x));
delta_term_b2.alt = zeros(1,length(x));
delta_term_b3.sea = zeros(1,length(x));
delta_term_b3.alt = zeros(1,length(x));

for i = 1:spar.i_CCW(2)-1 %top rear skin
    counter = 1;
    delta_term_b1(i) = (((x(i+1)-x(i))^2 + (y(i+1) - y(i))^2)^0.5)/skin.t;
    delta_term_b2.sea(i) = Booms.sea * x(i) * (((x(i+1)-x(i))^2 + (y(i+1) - y(i))^2)^0.5)/thickness(i);
    delta_term_b2.alt(i) = Booms.alt * x(i) * (((x(i+1)-x(i))^2 + (y(i+1) - y(i))^2)^0.5)/thickness_1;
    delta_term_b3.sea(i) = Booms.sea * y(i) * (((x(i+1)-x(i))^2 + (y(i+1) - y(i))^2)^0.5)/thickness(i);
    delta_term_b3.alt(i) = Booms.alt * y(i) * (((x(i+1)-x(i))^2 + (y(i+1) - y(i))^2)^0.5)/thickness_1;
    counter = counter + 1;
end
for i = spar.i_CCW(3):length(x) %bottom rear skin
    delta_term_b1(counter) = (((x(i+1)-x(i))^2 + (y(i+1) - y(i))^2)^0.5)/skin.t;
    delta_term_b2.sea(counter) = Booms.sea * x(i) * (((x(i+1)-x(i))^2 + (y(i+1) - y(i))^2)^0.5)/thickness(i);
    delta_term_b2.alt(counter) = Booms.alt * x(i) * (((x(i+1)-x(i))^2 + (y(i+1) - y(i))^2)^0.5)/thickness_1;
    delta_term_b3.sea(counter) = Booms.sea * y(i) * (((x(i+1)-x(i))^2 + (y(i+1) - y(i))^2)^0.5)/thickness(i);
    delta_term_b3.alt(counter) = Booms.alt * y(i) * (((x(i+1)-x(i))^2 + (y(i+1) - y(i))^2)^0.5)/thickness_1;
    counter = counter + 1;
end
    
sum_term_b1 = sum(delta_term_b1(:));
sum_term_b2.sea = sum(delta_term_b2.sea(:));
sum_term_b2.alt = sum(delta_term_b2.alt(:));
sum_term_b3.sea = sum(delta_term_b3.sea(:));
sum_term_b3.alt = sum(delta_term_b3.alt(:));

dtheta_dz_2.sea = q02.sea*(sum_term_b1) + (q02.sea-q01.sea)*(y_top-t_bottom)/spar.t + q02.sea*(y_bot-y_top)/(skin.thickness + spar.thickness) + ((Sy.sea*Ixy - Sx.sea*Ixx)/(Ixx*Iyy - Ixy^2))*(sum_term_b2.sea) + ((Sx.sea*Ixy - Sy.sea*Iyy)/(Ixx*Iyy-Ixy^2))*(sum_term_b3.sea);
dtheta_dz_2.alt = q02.alt*(sum_term_b1) + (q02.alt-q01.alt)*(y_top-t_bottom)/spar.t + q02.alt*(y_bot-y_top)/(skin.thickness + spar.thickness) + ((Sy.sea*Ixy - Sx.sea*Ixx)/(Ixx*Iyy - Ixy^2))*(sum_term_b2.alt) + ((Sx.sea*Ixy - Sy.sea*Iyy)/(Ixx*Iyy-Ixy^2))*(sum_term_b3.alt);

% Equating dtheta_dz_1 and _2
eq2.sea = dtheta_dz_2.sea - dtheta_dz_1.sea;
eq2.alt = dtheta_dz_2.alt - dtheta_dz_1.alt;

%% Solve
[q01.sea,q02.sea] = solve(eq1.sea==0,eq2.sea==0);
[q01.alt,q02.alt] = solve(eq1.alt==0,eq2.alt==0);

qtot.sea = zeros(1,length(x)+1);           % the +1 is for the spar
qtot.alt = zeros(1,length(x)+1);

for i = 1:spar.i_CCW(2)-1
    q.sea(i,j,k) = qb.sea(i,j,k) + q02.sea;
    q.alt(i,j,k) = qb.alt(i,j,k) + q02.alt;
end
for i = spar.i_CCW(2):spar.i_CCW(3)-1
    q.sea(i,j,k) = qb.sea(i,j,k) + q01.sea;
    q.alt(i,j,k) = qb.alt(i,j,k) + q01.alt;
end
for i = spar.i_CCW(3):length(x)
    q.sea(i,j,k) = qb.sea(i,j,k) + q01.sea;
    q.alt(i,j,k) = qb.alt(i,j,k) + q01.alt;
end
    q.sea.spar(j,k) = q01.sea - q02.sea;
    q.alt.spar(j,k) = q01.alt - q02.alt;


% % shear stress tau
% % please calculate the shear stress from the shear flow
end
end