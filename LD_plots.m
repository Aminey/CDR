function [L_distribution, D_distribution] = LD_plots()



%airfoil shape
% airfoil_x = xlsread('NACA 2415.xlsx','A1:A99');
% airfoil_y = xlsread('NACA 2415.xlsx','B1:B99');

%constants
g = 9.8; %m/s^2                 gravity
rho_sea = 1.225; %kg/m^3        density of air at sea level
rho_alt = 0.78205; %kg/m^3      density of air at 14,600 ft
e = 0.79; %                     oswald efficiency

%wing geometry
b = 10.82; %m                   span
c = 1.5; %m                     chord
S = b*c; %m^2                   surface area
m = 1100; %kg                   aircraft maximum gross mass
W = 1100*9.8; %N                  aircraft maximum gross weight
AR = b^2/S; %                   aspect ratio

%aerodynamic data
CD0_sea = 0.00562;
CD0_alt = 0.0058;

C_la_sea = 6.537; %1/rad
C_l0_sea = 0.2411;
C_la_alt = 6.526; %1/rad
C_l0_alt = 0.2411;

C_La_sea = C_la_sea/(1+C_la_sea/(AR*pi*e)); %1/rad
C_L0_sea = (C_La_sea/C_la_sea)*C_l0_sea;
C_La_alt = C_la_alt/(1+C_la_alt/(AR*pi*e)); %1/rad
C_L0_alt = (C_La_alt/C_la_alt)*C_l0_alt;

CD0_sea = 0.00586;
CD0_alt = 0.00608;

%critical flight condition loading
%[PHAA PLAA NHAA Downward_Gust NLAA]
n_sea = [4.4 4.4 -1.76 -1.82 -1.117];
v_sea = [59.7 95.8 39.9 63.9 95.83];
n_alt = [4.4 4.4 -1.76 -2.103 -1.33];
v_alt = [75.9 95.8 51.3 63.9 95.83];

%define wing semispan
z = 0:0.01:b/2;

% Preallocate L and D matrices
L_distribution.sea = zeros(length(z),5);
L_distribution.alt = zeros(length(z),5);
D_distribution.sea = zeros(length(z),5);
D_distribution.alt = zeros(length(z),5);

for i = 1:5
    %altitude sea level
    L_sea(i) = n_sea(i)*W;
    CL_sea(i) = L_sea(i)/(0.5*rho_sea*v_sea(i)^2*S);
    alpha_sea(i) = ((CL_sea(i)-C_L0_sea)/C_La_sea)*180/pi; %degrees
    
    %lift distribution
    L_distribution.sea(1:length(z),i) = (L_sea(i)/b + 4*L_sea(i)/(pi*b)*sqrt(1-(2*z/b).^2))/2;
    
    %altitude 14600 ft
    L_alt(i) = n_alt(i)*W;
    CL_alt(i) = L_alt(i)/(0.5*rho_alt*v_alt(i)^2*S);
    alpha_alt(i) = ((CL_alt(i)-C_L0_alt)/C_La_alt)*180/pi; %degrees
    
    %lift distribution
    L_distribution.alt(1:length(z),i) = (L_alt(i)/b + 4*L_alt(i)/(pi*b)*sqrt(1-(2*z/b).^2))/2;
end

% CL_seat(4) = -1.76*W/(0.5*rho_sea*v_sea(i)^2*S);
% CL_seat(5) = -1*W/(0.5*rho_sea*v_sea(i)^2*S);
% alpha_sea(4) = ((CL_seat(4)-C_L0_sea)/C_La_sea)*180/pi - atand(15.24/63.9); %degrees
% alpha_sea(5) = ((CL_seat(5)-C_L0_sea)/C_La_sea)*180/pi - atand(15.24/63.9); %degrees
% alpha_sea(4) = alpha_sea(4) - atand(15.24/63.9);
% alpha_sea(5) = alpha_sea(5) - atand(7.62/95.8);
% 
% CL_altt(4) = -1.76*W/(0.5*rho_alt*v_alt(i)^2*S);
% CL_altt(5) = -1*W/(0.5*rho_alt*v_alt(i)^2*S);
% alpha_alt(4) = ((CL_altt(4)-C_L0_alt)/C_La_alt)*180/pi - atand(15.24/63.9); %degrees
% alpha_alt(5) = ((CL_altt(5)-C_L0_alt)/C_La_alt)*180/pi - atand(15.24/63.9); %degrees
% alpha_alt(4) = alpha_alt(4) - atand(15.24/63.9);
% alpha_alt(5) = alpha_alt(5) - atand(7.62/95.8);

%% Plots figure 1 to 8
figure; 
plot(z,L_distribution.sea);
legend('PHAA','PLAA','NHAA','Maximum Downward Gust','NLAA Gust','Location','Best');
xlabel('Spanwise Length (m)');
ylabel('Lift Force at Sea Level (N/m)');

figure; 
plot(z,L_distribution.alt);
legend('PHAA','PLAA','NHAA','Maximum Downward Gust','NLAA Gust','Location','Best');
xlabel('Spanwise Length (m)');
ylabel('Lift Force at Altitude (N/m)'); 

%generate drag distributions
for i = 1:5
    %altitude sea level
    for j = 1:length(z)
        if z(j) < 0.9*b/2
            D_distribution.sea(j,i) = (0.5*rho_sea*v_sea(i)^2*S*(CD0_sea+CL_sea(i)^2/(pi*AR*e)))/2;
        else
            D_distribution.sea(j,i) = (1.2*0.5*rho_sea*v_sea(i)^2*S*(CD0_sea+CL_sea(i)^2/(pi*AR*e)))/2;
        end
    end

    %altitude 14600 ft
    for j = 1:length(z)
        if z(j) < 0.9*b/2
            D_distribution.alt(j,i) = (0.5*rho_alt*v_alt(i)^2*S*(CD0_alt+CL_alt(i)^2/(pi*AR*e)))/2;
        else
            D_distribution.alt(j,i) = (1.2*0.5*rho_alt*v_alt(i)^2*S*(CD0_alt+CL_alt(i)^2/(pi*AR*e)))/2;
        end
    end
end

%sea level
figure; 
plot(z,D_distribution.sea);
legend('PHAA','PLAA','NHAA','Maximum Downward Gust','NLAA Gust','Location','Best');
xlabel('Spanwise Length (m)');
ylabel('Drag Force at Sea Level (N/m)'); 

%service ceiling
figure; 
plot(z,D_distribution.alt);
legend('PHAA','PLAA','NHAA','Maximum Downward Gust','NLAA Gust','Location','Best');
xlabel('Spanwise Length (m)');
ylabel('Drag Force at Altitude (N/m)'); 

disp('LD_plots complete');

