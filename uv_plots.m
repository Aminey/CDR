function [u, v] = uv_plots(z, Ixx, Iyy, Ixy, structure, My, Mx)

% Ixx = structure.inertias(1);
% Iyy = structure.inertias(2);
% Ixy = structure.inertias(3);

%% Preallocation
% u's and v's
% How to read: _sea = dA, .sea = integral
du_sea = zeros(length(z),5);    du.sea = zeros(length(z),5);
u_sea = zeros(length(z),5);     u.sea = zeros(length(z),5);
du_alt = zeros(length(z),5);    du.alt = zeros(length(z),5);
u_alt = zeros(length(z),5);     u.alt = zeros(length(z),5);
dv_sea = zeros(length(z),5);    dv.sea = zeros(length(z),5);
v_sea = zeros(length(z),5);     v.sea = zeros(length(z),5);
dv_alt = zeros(length(z),5);    dv.alt = zeros(length(z),5);
v_alt = zeros(length(z),5);     v.alt = zeros(length(z),5);


%% Displacements
K = 1/(structure.E*(Ixx*Iyy - Ixy^2));
ddu_sea = -K.*(-Ixy.*Mx.sea + Ixx.*My.sea);
ddu_alt = -K.*(-Ixy.*Mx.alt + Ixx.*My.alt);
ddv_sea = -K.*(Iyy.*Mx.sea - Ixy.*My.sea);
ddv_alt = -K.*(Iyy.*Mx.alt - Ixy.*My.alt);

% duu = cumtrapz(z,ddu_sea);
% uu = cumtrapz(z,duu);
% figure;
% plot(z,uu);
% legend('PHAA','PLAA','NHAA','Maximum Downward Gust','NLAA Gust','Location','Best');
% xlabel('Spanwise Length (m)');
% ylabel('Displacement u, at sea level (m)');


% dz = z(2) - z(1);
% 
% for i = 1:5
%     for j = 1:length(z)-1
%         du_sea(j+1,i) = (z(j+1) - z(j))*(ddu_sea(j,i));
%         du.sea(j+1,i) = du.sea(j,i) + du_sea(j+1,i);
%         dv_sea(j+1,i) = (z(j+1) - z(j))*(ddv_sea(j,i));
%         dv.sea(j+1,i) = dv.sea(j,i) + dv_sea(j+1,i);
%         
%         du_alt(j+1,i) = (z(j+1) - z(j))*(ddu_alt(j,i));
%         du.alt(j+1,i) = du.alt(j,i) + du_alt(j+1,i);  
%         dv_alt(j+1,i) = (z(j+1) - z(j))*(ddv_alt(j,i));
%         dv.alt(j+1,i) = dv.alt(j,i) + dv_alt(j+1,i);  
%         
%     end
% end
% 
% for i = 1:5
%     for j = 1:length(z)-1
%         u_sea(j,i) = (z(j+1) - z(j))*(du_sea(j,i));
%         u.sea(j+1,i) = u.sea(j,i) + u_sea(j,i);
%         u_alt(j,i) = (z(j+1) - z(j))*(du_alt(j,i));
%         u.alt(j+1,i) = u.alt(j,i) + u_alt(j,i);        
%         
%         v_sea(j,i) = (z(j+1) - z(j))*(dv_sea(j,i));
%         v.sea(j+1,i) = v.sea(j,i) + v_sea(j,i);
%         v_alt(j,i) = (z(j+1) - z(j))*(dv_alt(j,i));
%         v.alt(j+1,i) = v.alt(j,i) + v_alt(j,i);        
%                 
%     end
% end

du.sea = cumtrapz(z,ddu_sea);
u.sea = cumtrapz(z,du.sea);

dv.sea = cumtrapz(z,ddv_sea);
v.sea = cumtrapz(z,dv.sea);

du.alt = cumtrapz(z,ddu_alt);
u.alt = cumtrapz(z,du.alt);

dv.alt = cumtrapz(z,ddv_alt);
v.alt = cumtrapz(z,dv.alt);



%% Plot u and v, at Sea level and altitude
figure;
plot(z,u.sea);
legend('PHAA','PLAA','NHAA','Maximum Downward Gust','NLAA Gust','Location','Best');
xlabel('Spanwise Length (m)');
ylabel('Displacement u, at sea level (m)');

figure;
plot(z,u.alt);
legend('PHAA','PLAA','NHAA','Maximum Downward Gust','NLAA Gust','Location','Best');
xlabel('Spanwise Length (m)');
ylabel('Displacement u, at altitude (m)');

figure;   
plot(z,v.sea);
legend('PHAA','PLAA','NHAA','Maximum Downward Gust','NLAA Gust','Location','Best');
xlabel('Spanwise Length (m)');
ylabel('Displacement v, at sea level (m)');

figure;
plot(z,v.alt);
legend('PHAA','PLAA','NHAA','Maximum Downward Gust','NLAA Gust','Location','Best');
xlabel('Spanwise Length (m)');
ylabel('Displacement v, at altitude (m)');
