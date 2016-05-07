% %% shear flow
% % form cell 1: boom area, distance between booms and coordinates of each boom
% B = [fliplr(BU),BL(2:end)];
% x_boom = [fliplr(x_boomU),x_boomL(2:end)];
% y_boom = [fliplr(yU(i_BU(:))),yL(i_BL(2:end))];
% L_boom = [fliplr(L_boomU),L_boomL,h_spar(2)];
% 

% 
% % calculate area of triangle formed by two nodes on the airfoil profile and point(x_spar(1),0)
% nq = length(L_boom);
% A = zeros(1,nq);
% for i = 1:nq-1
%     A(i) = abs(x_boom(i)*(y_boom(i+1)-0) + x_boom(i+1)*(0-y_boom(i)) + x_spar(1)*(y_boom(i)-y_boom(i+1)))/2;
% end
% A(end) = abs(x_boom(end)*(y_boom(1)-0) + x_boom(1)*(0-y_boom(end)) + x_spar(1)*(y_boom(end)-y_boom(1)))/2;
% Asum = sum(A);
% 
% % calculate the area of each cell
% i_A1 = find(ismember(x_boom,x_spar(1)));
% A1 = A(i_A1(1):i_A1(2)-1);
% A1sum = sum(A1);
% A2sum = Asum - A1sum;
% 

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
x = zeros(1,length(airfoil.booms));
y = zeros(1,length(airfoil.booms));
delta_A = zeros(1,length(airfoil.booms));
delta_Term_2_alt = zeros(length(x),length(z),5);
delta_Term_2_sea = zeros(length(x),length(z),5);
Term_2_alt = zeros(length(z),5);
Term_2_sea = zeros(length(z),5);

for i = 1:length(airfoil.booms)
    x(i) = airfoil.booms(i).x_coordinate;                %%
    y(i) = airfoil.booms(i).y_coordinate;
end


for k = 1:5
    for j = 1:length(z)
        for i = 1:length(x)-1 % for x points there will be x flows, but the first is cut
            delta_A(i) = abs(x(i+1)*y(i) - y(i+1)*x(i))/2;
            delta_Term_2_alt(i,j,k) = 2*delta_A(i)*qb_alt(i,j,k);
            delta_Term_2_sea(i,j,k) = 2*delta_A(i)*qb_sea(i,j,k);           
        end
        Term_2_alt(j,k) = sum(delta_Term_2_alt);
        Term_2_sea(j,k) = sum(delta_Term_2_sea);
        
    end
end