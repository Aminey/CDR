% %% shear flow
% % form cell 1: boom area, distance between booms and coordinates of each boom
% B = [fliplr(BU),BL(2:end)];
% x_boom = [fliplr(x_boomU),x_boomL(2:end)];
% y_boom = [fliplr(yU(i_BU(:))),yL(i_BL(2:end))];
% L_boom = [fliplr(L_boomU),L_boomL,h_spar(2)];
% 
% %% calculate area of triangle formed by two nodes on the airfoil profile and point(x_spar(1),0)
% nq = length(L_boom);
% A = zeros(1,nq);
% for i = 1:nq-1
%     A(i) = abs(x_boom(i)*(y_boom(i+1)-0) + x_boom(i+1)*(0-y_boom(i)) + x_spar(1)*(y_boom(i)-y_boom(i+1)))/2;
% end
% A(end) = abs(x_boom(end)*(y_boom(1)-0) + x_boom(1)*(0-y_boom(end)) + x_spar(1)*(y_boom(end)-y_boom(1)))/2;
% Asum = sum(A);


% % separate cell 1 and cell 2 
% qb1 = qb(i_A1(1):i_A1(2)-1);
% L_boom1 = L_boom(i_A1(1):i_A1(2)-1);
% L1sum = sum(L_boom1);
% 
% qb2 = [qb(1:i_A1(1)-1),qb(i_A1(2):end)];
% L_boom2 = [L_boom(1:i_A1(1)-1),L_boom(i_A1(2):end)];

% % modify distance between booms at the rear spar to compensate the change of the thickness

% L_boom2(end) = L_boom2(end)*t_skin/t_spar; 
% L2sum = sum(L_boom2);
% 
% % equations
% syms q01 q02
% % eq1: equating moments of applied shear and pitch moment to moments of internal shear flow
% %
% % please set up equation 1 here by using eq1=...
% %
% 
% % eq2: set the angle of twist of each cell to be the same
% %
% % please set up equation 2 here by using eq2=...
% %
%   
% [q01,q02] = solve(eq1==0,eq2==0);
% 
% 
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

%% Initialize
x = zeros(1,length(airfoil.booms));
y = zeros(1,length(airfoil.booms));

%% Find area of each section
i_A1 = find(ismember(x,spars.x(1)));        % i don't think this is right.... need to fiddle with it.
A1 = sum(delta_A(i_A1(1):i_A1(2)));
A2 = A_total - A1;
%%
syms q01.sea q02.sea
syms q01.alt q02.alt
%% EQ 1
for k = 1:5
    for j = 1:length(z)
        M_0.sea + Sy.sea(j,k)*(x_quarterchord) - Sx.sea(j,k)*(0) = 2*A1*q01.sea(j,k) + 2*A2*q02.sea(j,k) + Term_2.sea;
        M_0.alt + Sy.alt(j,k)*(x_quarterchord) - Sx.alt(j,k)*(0) = 2*A1*q01.alt(j,k) + 2*A2*q02.alt(j,k) + Term_2.alt;
    end
end
%% EQ 2
% dtheta_dz_1
delta_term_a1(i) = zeros(1,length(x)); %% this is wrong. doesnt go from 1 to length(x). goes over front section
delta_term_a2.sea(i) = zeros(1,length(x));
delta_term_a2.alt(i) = zeros(1,length(x));
delta_term_a3.sea(i) = zeros(1,length(x));
delta_term_a3.alt(i) = zeros(1,length(x));

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

dtheta_dz.sea = (1/(2*A1*G)) * (q01.sea*(sum_term_a1) + (q01.sea-q02.sea)*(y_bot - y_top)/spar.thickness + ((Sy.sea*Ixy - Sx.sea*Ixx)/(Ixx*Iyy-Ixy^2))*(sum_term_a2.sea) + ((Sx.sea*Ixy-Sy.sea*Iyy)/(Ixx*Iyy-Ixy^2))*(sum_term_a3.sea));
dtheta_dz.alt = (1/(2*A1*G)) * (q01.sea*(sum_term_a1) + (q01.alt-q02.alt)*(y_bot - y_top)/spar.thickness + ((Sy.alt*Ixy - Sx.alt*Ixx)/(Ixx*Iyy-Ixy^2))*(sum_term_a2.alt) + ((Sx.alt*Ixy-Sy.alt*Iyy)/(Ixx*Iyy-Ixy^2))*(sum_term_a3.alt));

% dtheta_dz 2
delta_term_b1(i) = zeros(1,length(x)); %% this is wrong. doesnt go from 1 to length(x). goes over top and bottom
delta_term_b2.sea(i) = zeros(1,length(x));
delta_term_b2.alt(i) = zeros(1,length(x));
delta_term_b3.sea(i) = zeros(1,length(x));
delta_term_b3.alt(i) = zeros(1,length(x));

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

dtheta_dz.sea = q02*(sum_term_b1) + (q02-q01)*(y_top-t_bottom)/spar.thickness + q02*(y_bot-y_top)/(skin.thickness + spar.thickness) + ((Sy.sea*Ixy - Sx.sea*Ixx)/(Ixx*Iyy - Ixy^2))*(sum_term_b2.sea) + ((Sx.sea*Ixy - Sy.sea*Iyy)/(Ixx*Iyy-Ixy^2))*(sum_term_b3.sea);
dtheta_dz.alt = q02*(sum_term_b1) + (q02-q01)*(y_top-t_bottom)/spar.thickness + q02*(y_bot-y_top)/(skin.thickness + spar.thickness) + ((Sy.sea*Ixy - Sx.sea*Ixx)/(Ixx*Iyy - Ixy^2))*(sum_term_b2.alt) + ((Sx.sea*Ixy - Sy.sea*Iyy)/(Ixx*Iyy-Ixy^2))*(sum_term_b3.alt);
