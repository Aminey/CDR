function [] = Shear_Flow_0()
%% Initialize



%% Find area of each section
i_A1 = find(ismember(x,spars.x(1)));        % i don't think this is right.... need to fiddle with it.
A1 = sum(delta_A(i_A1(1):i_A1(2)));
A2 = A_total - A1;
%%
syms q01.sea(j,k) q02.sea(j,k)
syms q01.alt(j,k) q02.alt(j,k)
%% EQ 1
for k = 1:5
    for j = 1:length(z)
        eq1.sea = M_0.sea + Sy.sea(j,k)*(x_quarterchord) - Sx.sea(j,k)*(0) - (2*A1*q01.sea(j,k) + 2*A2*q02.sea(j,k) + Term_2.sea);
        eq1.alt = M_0.alt + Sy.alt(j,k)*(x_quarterchord) - Sx.alt(j,k)*(0) - (2*A1*q01.alt(j,k) + 2*A2*q02.alt(j,k) + Term_2.alt);
    end
end
%% EQ 2
% dtheta_dz_1
delta_term_a1 = zeros(1,length(x)); %% this is wrong. doesnt go from 1 to length(x). goes over front section
delta_term_a2.sea = zeros(1,length(x));
delta_term_a2.alt = zeros(1,length(x));
delta_term_a3.sea = zeros(1,length(x));
delta_term_a3.alt = zeros(1,length(x));

for i = 1:length(x)-1 % "Sum CCW over front section" also wrong...
    delta_term_a1(i) = (((x(i+1)-x(i))^2 + (y(i+1) - y(i))^2)^0.5)/skin.thickness;
    delta_term_a2.sea(i) = Booms.sea * x(i) * (((x(i+1)-x(i))^2 + (y(i+1) - y(i))^2)^0.5)/thickness(i);
    delta_term_a2.alt(i) = Booms.alt * x(i) * (((x(i+1)-x(i))^2 + (y(i+1) - y(i))^2)^0.5)/thickness_1;
    delta_term_a3.sea(i) = Booms.sea * y(i) * (((x(i+1)-x(i))^2 + (y(i+1) - y(i))^2)^0.5)/thickness(i);
    delta_term_a3.alt(i) = Booms.alt * y(i) * (((x(i+1)-x(i))^2 + (y(i+1) - y(i))^2)^0.5)/thickness_1;
end
sum_term_a1 = sum(delta_term_a1(:));
sum_term_a2.sea = sum(delta_term_a2.sea(:));
sum_term_a2.alt = sum(delta_term_a2.alt(:));
sum_term_a3.sea = sum(delta_term_a3.sea(:));
sum_term_a3.alt = sum(delta_term_a3.alt(:));

dtheta_dz_1.sea = (1/(2*A1*G)) * (q01.sea*(sum_term_a1) + (q01.sea-q02.sea)*(y_bot - y_top)/spar.thickness + ((Sy.sea*Ixy - Sx.sea*Ixx)/(Ixx*Iyy-Ixy^2))*(sum_term_a2.sea) + ((Sx.sea*Ixy-Sy.sea*Iyy)/(Ixx*Iyy-Ixy^2))*(sum_term_a3.sea));
dtheta_dz_1.alt = (1/(2*A1*G)) * (q01.sea*(sum_term_a1) + (q01.alt-q02.alt)*(y_bot - y_top)/spar.thickness + ((Sy.alt*Ixy - Sx.alt*Ixx)/(Ixx*Iyy-Ixy^2))*(sum_term_a2.alt) + ((Sx.alt*Ixy-Sy.alt*Iyy)/(Ixx*Iyy-Ixy^2))*(sum_term_a3.alt));

% dtheta_dz 2
delta_term_b1 = zeros(1,length(x)); %% this is wrong. doesnt go from 1 to length(x). goes over top and bottom
delta_term_b2.sea = zeros(1,length(x));
delta_term_b2.alt = zeros(1,length(x));
delta_term_b3.sea = zeros(1,length(x));
delta_term_b3.alt = zeros(1,length(x));

for i = 1:length(x)-1 % "Sum CCW over front section" also wrong...
    delta_term_b1(i) = (((x(i+1)-x(i))^2 + (y(i+1) - y(i))^2)^0.5)/skin.thickness;
    delta_term_b2.sea(i) = Booms.sea * x(i) * (((x(i+1)-x(i))^2 + (y(i+1) - y(i))^2)^0.5)/thickness(i);
    delta_term_b2.alt(i) = Booms.alt * x(i) * (((x(i+1)-x(i))^2 + (y(i+1) - y(i))^2)^0.5)/thickness_1;
    delta_term_b3.sea(i) = Booms.sea * y(i) * (((x(i+1)-x(i))^2 + (y(i+1) - y(i))^2)^0.5)/thickness(i);
    delta_term_b3.alt(i) = Booms.alt * y(i) * (((x(i+1)-x(i))^2 + (y(i+1) - y(i))^2)^0.5)/thickness_1;
end
sum_term_b1 = sum(delta_term_b1(:));
sum_term_b2.sea = sum(delta_term_b2.sea(:));
sum_term_b2.alt = sum(delta_term_b2.alt(:));
sum_term_b3.sea = sum(delta_term_b3.sea(:));
sum_term_b3.alt = sum(delta_term_b3.alt(:));

dtheta_dz_2.sea = q02*(sum_term_b1) + (q02-q01)*(y_top-t_bottom)/spar.thickness + q02*(y_bot-y_top)/(skin.thickness + spar.thickness) + ((Sy.sea*Ixy - Sx.sea*Ixx)/(Ixx*Iyy - Ixy^2))*(sum_term_b2.sea) + ((Sx.sea*Ixy - Sy.sea*Iyy)/(Ixx*Iyy-Ixy^2))*(sum_term_b3.sea);
dtheta_dz_2.alt = q02*(sum_term_b1) + (q02-q01)*(y_top-t_bottom)/spar.thickness + q02*(y_bot-y_top)/(skin.thickness + spar.thickness) + ((Sy.sea*Ixy - Sx.sea*Ixx)/(Ixx*Iyy - Ixy^2))*(sum_term_b2.alt) + ((Sx.sea*Ixy - Sy.sea*Iyy)/(Ixx*Iyy-Ixy^2))*(sum_term_b3.alt);

% Equating dtheta_dz_1 and _2
eq2.sea = dtheta_dz_2.sea - dtheta_dz_1.sea;
eq2.alt = dtheta_dz_2.alt - dtheta_dz_1.alt;

%% Solve
[q01.sea,q02.sea] = solve(eq1.sea==0,eq2.sea==0);
[q01.alt,q02.alt] = solve(eq1.alt==0,eq2.alt==0);

% % shear flow along airfoil contour from top right corner to top right corner CCW 
% q = zeros(1,nq);   
% for i = 1:i_A1(1)-1
%     q(i) = qb2(i) + q02;
% end
% 
% for i = 1 : length(qb1)
%     j = i + i_A1(1) - 1;
%     q(j) = qb1(i) + q01;
% end
% 
% for i = 1 : nq - i_A1(2) + 1
%     j = i + i_A1(2) - 1;
%     k = i + i_A1(1) - 1;
%     q(j) = qb2(k) + q02;
% end
% % please calculate the shear flow in the central spar here
% % verification:
% % please verify your results in 4 ways (see instructions)
% % shear stress tau
% % please calculate the shear stress from the shear flow