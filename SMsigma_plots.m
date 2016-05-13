function [Sx, Sy, Mx, My, sigma_z, wx_sea] = SMsigma_plots(alpha, x, y, z, Ixx, Iyy, Ixy, L_distribution, D_distribution)

%% Preallocate Matrices
wx_sea = zeros(length(z),5);
wx_alt = zeros(length(z),5);
wy_sea = zeros(length(z),5);
wy_alt = zeros(length(z),5);

Sx_sea = zeros(length(z),5);
Sy_sea = zeros(length(z),5);
Sx_alt = zeros(length(z),5);
Sy_alt = zeros(length(z),5);

% Sum Matrix of Sx and Sy
Sx.sea = zeros(length(z),5);
Sy.sea = zeros(length(z),5);
Sx.alt = zeros(length(z),5);
Sy.alt = zeros(length(z),5);

Mx_sea = zeros(length(z),5);
My_sea = zeros(length(z),5);
Mx_alt = zeros(length(z),5);
My_alt = zeros(length(z),5);

% Sum Matrix of Mx and My
Mx.sea = zeros(length(z),5);
My.sea = zeros(length(z),5);
Mx.alt = zeros(length(z),5);
My.alt = zeros(length(z),5);

sigma_z.alt = zeros(length(x),length(z),5);
sigma_z.sea = zeros(length(x),length(z),5);
sigma_plot.alt = zeros(length(x),5);
sigma_plot.sea = zeros(length(x),5);

%transform forces
for i = 1:5
    %axial
    wx_sea(1:length(z),i) = -L_distribution.sea(1:length(z),i)*sind(alpha.sea(i))+D_distribution.sea(1:length(z),i)*cosd(alpha.sea(i));
    %normal
    wy_sea(1:length(z),i) = L_distribution.sea(1:length(z),i)*cosd(alpha.sea(i))+D_distribution.sea(1:length(z),i)*sind(alpha.sea(i));
    
    %axial
    wx_alt(1:length(z),i) = -L_distribution.alt(1:length(z),i)*sind(alpha.alt(i))+D_distribution.alt(1:length(z),i)*cosd(alpha.alt(i));
    %normal
    wy_alt(1:length(z),i) = L_distribution.alt(1:length(z),i)*cosd(alpha.alt(i))+D_distribution.alt(1:length(z),i)*sind(alpha.alt(i));
end

%% Plots
figure; 
plot(z,wx_sea);
legend('PHAA','PLAA','NHAA','Maximum Downward Gust','NLAA Gust','Location','Best');
xlabel('Spanwise Length (m)');
ylabel('W_x Force at Sea Level (N/m)');

figure; 
plot(z,wy_sea);
legend('PHAA','PLAA','NHAA','Maximum Downward Gust','NLAA Gust','Location','Best');
xlabel('Spanwise Length (m)');
ylabel('W_y Force at Sea Level (N/m)');

figure; 
plot(z,wx_alt);
legend('PHAA','PLAA','NHAA','Maximum Downward Gust','NLAA Gust','Location','Best');
xlabel('Spanwise Length (m)');
ylabel('W_x Force at Service Ceiling (N/m)');

figure; 
plot(z,wy_alt);
legend('PHAA','PLAA','NHAA','Maximum Downward Gust','NLAA Gust','Location','Best');
xlabel('Spanwise Length (m)');
ylabel('W_y Force at Service Ceiling (N/m)');

for i = 1:5;
    for j = length(z):-1:2;
        Sx_sea(j,i) = -(z(j)-z(j-1)) * (wx_sea(j,i) + wx_sea(j-1,i))/2;
        Sx.sea(j-1,i) = Sx.sea(j,i) + Sx_sea(j,i);
        
        Sx_alt(j,i) = -(z(j)-z(j-1)) * (wx_alt(j,i) + wx_alt(j-1,i))/2;
        Sx.alt(j-1,i) = Sx.alt(j,i) + Sx_alt(j,i);
        
        Sy_sea(j,i) = -(z(j)-z(j-1)) * (wy_sea(j,i) + wy_sea(j-1,i))/2;
        Sy.sea(j-1,i) = Sy.sea(j,i) + Sy_sea(j,i);
        
        Sy_alt(j,i) = -(z(j)-z(j-1)) * (wy_alt(j,i) + wy_alt(j-1,i))/2;
        Sy.alt(j-1,i) = Sy.alt(j,i) + Sy_alt(j,i);
        %
        %
        Mx_sea(j,i) = (z(j)-z(j-1)) * (Sy_sea(j,i) + Sy_sea(j-1,i))/2;
        Mx.sea(j-1,i) = Mx.sea(j,i) + Mx_sea(j,i);
        
        Mx_alt(j,i) = (z(j)-z(j-1)) * (Sy_alt(j,i) + Sy_alt(j-1,i))/2;
        Mx.alt(j-1,i) = Mx.alt(j,i) + Mx_alt(j,i);
        
        My_sea(j,i) = (z(j)-z(j-1)) * (Sx_sea(j,i) + Sx_sea(j-1,i))/2;
        My.sea(j-1,i) = My.sea(j,i) + My_sea(j,i);
        
        My_alt(j,i) = (z(j)-z(j-1)) * (Sx_alt(j,i) + Sx_alt(j-1,i))/2;
        My.alt(j-1,i) = My.alt(j,i) + My_alt(j,i);
    end
end

figure;
plot(z,Sx.sea);
legend('PHAA','PLAA','NHAA','Maximum Downward Gust','NLAA Gust','Location','Best');
xlabel('Spanwise Length (m)');
ylabel('Force Sx at Sea Level (N)');

figure;
plot(z,Sy.sea);
legend('PHAA','PLAA','NHAA','Maximum Downward Gust','NLAA Gust','Location','Best');
xlabel('Spanwise Length (m)');
ylabel('Force Sy at Sea Level (N)');

figure;
plot(z,Mx.sea);
legend('PHAA','PLAA','NHAA','Maximum Downward Gust','NLAA Gust','Location','Best');
xlabel('Spanwise Length (m)');
ylabel('Moment Mx at Sea Level (N*m)');

figure;
plot(z,My.sea);
legend('PHAA','PLAA','NHAA','Maximum Downward Gust','NLAA Gust','Location','Best');
xlabel('Spanwise Length (m)');
ylabel('Moment My at Sea Level (N*m)');

figure;
plot(z,Sx.alt);
legend('PHAA','PLAA','NHAA','Maximum Downward Gust','NLAA Gust','Location','Best');
xlabel('Spanwise Length (m)');
ylabel('Force Sx at Altitude (N)');

figure;
plot(z,Sy.alt);
legend('PHAA','PLAA','NHAA','Maximum Downward Gust','NLAA Gust','Location','Best');
xlabel('Spanwise Length (m)');
ylabel('Force Sy at Altitude (N)');

figure;
plot(z,Mx.alt);
legend('PHAA','PLAA','NHAA','Maximum Downward Gust','NLAA Gust','Location','Best');
xlabel('Spanwise Length (m)');
ylabel('Moment Mx at Altitude (N*m)');

figure;
plot(z,My.alt);
legend('PHAA','PLAA','NHAA','Maximum Downward Gust','NLAA Gust','Location','Best');
xlabel('Spanwise Length (m)');
ylabel('Moment My at Altitude (N*m)');


%% Sigma Z's

for k = 1:5
    for j = 1:length(z) 
        for i = 1:length(x)
            sigma_z.alt(i,j,k) = Mx.alt(j,k)*(Iyy*y(i) - Ixy*x(i))/(Ixx*Iyy - Ixy^2) + ...
                                 My.alt(j,k)*(Ixx*x(i) - Ixy*y(i))/(Ixx*Iyy - Ixy^2);
            sigma_z.sea(i,j,k) = Mx.sea(j,k)*(Iyy*y(i) - Ixy*x(i))/(Ixx*Iyy - Ixy^2) + ...
                                 My.sea(j,k)*(Ixx*x(i) - Ixy*y(i))/(Ixx*Iyy - Ixy^2);
        end
    end    
end

    for i = 1:length(x)
        sigma_plot.alt(i) = sigma_z.alt(i,1,1);
        sigma_plot.sea(i) = sigma_z.sea(i,1,1);
    end


figure;
plot(x,sigma_plot.alt);
legend('PHAA','PLAA','NHAA','Maximum Downward Gust','NLAA Gust','Location','Best');
xlabel('Chordwise Length (m)');
ylabel('Sigma_zz at altitude (N/m^2');

figure;
plot(x,sigma_plot.sea);
legend('PHAA','PLAA','NHAA','Maximum Downward Gust','NLAA Gust','Location','Best');
xlabel('Chordwise Length (m)');
ylabel('Sigma_zz at sea level (N/m^2');

disp('SMsigma_plots complete');


