% q2 - q1 = (Sx Ixx - Sy Ixy / IxxIyy - Ixy2) * Br * xr - (Sy Iyy - SxIxy / IxxIyy - Ixy2) * Br *yr

qb_sea = zeros(length(x)*2,length(z),5);
qb_alt = zeros(length(x)*2,length(z),5);

for k = 1:5
    for j = 1:length(z)
        for i = length(x):-1:1
            qb_sea(length(x)-i+1,j,k) = (((Sx.sea(j,k)*structure.inertias(1) - Sy.sea(j,k)*stucture.inertias(3))/...
                                         (structure.inertias(1)*structure.inertias(2) - structure.inertias(3)))...
                                         *BOOM_A_SEA_UPPER(i,j,k)*x(i)) - (((Sy.sea(j,k)*structure.inertias(2) -...
                                         Sx.sea(j,k)*structure.inertias(3))/(structure.inertias(1)*structure.inertias(2)...
                                         -structure.inertias(3)))*BOOM_A_SEA_UPPER(i,j,k)*yU(i));
        for l = 1:length(x)-1
            qb_sea(length(x)+l,j,k) = (((Sx.sea(j,k)*structure.inertias(1) - Sy.sea(j,k)*stucture.inertias(3))/...
                                         (structure.inertias(1)*structure.inertias(2) - structure.inertias(3)))...
                                         *BOOM_A_SEA_UPPER(i,j,k)*x(i)) - (((Sy.sea(j,k)*structure.inertias(2) -...
                                         Sx.sea(j,k)*structure.inertias(3))/(structure.inertias(1)*structure.inertias(2)...
                                         -structure.inertias(3)))*BOOM_A_SEA_UPPER(i,j,k)*yU(i));
        end
        end
    end
end


% %% shear flow
% % form cell 1: boom area, distance between booms and coordinates of each boom
% B = [fliplr(BU),BL(2:end)];
% x_boom = [fliplr(x_boomU),x_boomL(2:end)];
% y_boom = [fliplr(yU(i_BU(:))),yL(i_BL(2:end))];
% L_boom = [fliplr(L_boomU),L_boomL,h_spar(2)];
% 
% figure
% plot(x_boom,y_boom,'ro','markersize',6)
% hold on
% plot(x,yU,'k',x,yL,'k','linewidth',1.5)
% plot([x(end),x(end)],[yU(end),yL(end)],'b',[x_spar(1),x_spar(1)],[yU(i_spar(1)),yL(i_spar(1))],'b','linewidth',2)
% ylim([-0.3 0.3]) 
% xlabel('x (m)')
% ylabel('y (m)')
% title('Boom Distribution')
% grid on
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
% % calcualte qb at the root of the wing 
% %
% %
% % please finish this part
% %
% %
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
% 
% %
% % please calculate the shear flow in the central spar here
% %
% 
% % verification:
% %
% % please verify your results in 4 ways (see instructions)
% %
% 
% 
% % shear stress tau
% %
% % please calculate the shear stress from the shear flow

