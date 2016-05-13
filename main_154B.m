%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     UCLA MAE 154B Main Script    %
%           Cessna 177B            %
%       Elhafsi, Fung, Kam         %
%       Spring Quarter 2016        %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all;
close all;
clc;

%airfoil shape
% airfoil_x = xlsread('NACA 2415.xlsx','A1:A99');
% airfoil_y = xlsread('NACA 2415.xlsx','B1:B99');

%% Constants
g = 9.8;                    % m/s^2                  gravity
rho_sea = 1.225;            % kg/m^3        density of air at sea level
rho_alt = 0.78205;          % kg/m^3      density of air at 14,600 ft
mu_sea = 1.789E-5;          % Pa*s
mu_alt = 2.4905E-5;         % Pa*s
theta = 30;                 % deg
alpha = 10;                 % deg


%% Wing Geometry
b = 10.82;                  % m        span
c = 1.5;                    % m        chord
S = b*c;                    % m^2      surface area
m = 1100;                   % kg       aircraft maximum gross mass
W = m*9.8;               % N        aircraft maximum gross weight
AR = b^2/S;                 %          aspect ratio
e = 0.79;                   %          oswald efficiency
z = 0:0.01:b/2;             % m        semi span vector
z_intervals = 5;
dz = round(length(z)/z_intervals);

%% Aerodynamic Data
CD0_sea = 0.00562;
CD0_alt = 0.0058;
CD0_sea = 0.00586;
CD0_alt = 0.00608;

C_la_sea = 6.537;           % 1/rad
C_l0_sea = 0.2411;
C_la_alt = 6.526;           % 1/rad
C_l0_alt = 0.2411;

C_La_sea = C_la_sea/(1+C_la_sea/(AR*pi*e));     %1/rad
C_L0_sea = (C_La_sea/C_la_sea)*C_l0_sea;
C_La_alt = C_la_alt/(1+C_la_alt/(AR*pi*e));     %1/rad
C_L0_alt = (C_La_alt/C_la_alt)*C_l0_alt;

Cm.sea = [-.021 -.054 -.055 -.055 -.055];   % first two terms are legit. the rest are bs
Cm.alt = [-.055 -.055 -.055 -.055 -.055];   % all bs

%critical flight condition loading
%[PHAA PLAA NHAA Downward_Gust NLAA]
n_sea = [4.4 4.4 -1.76 -1.82 -1.117];
v_sea = [59.7 95.8 39.9 63.9 95.83];
n_alt = [4.4 4.4 -1.76 -2.103 -1.33];
v_alt = [75.9 95.8 51.3 63.9 95.83];

Re.sea = rho_sea*v_sea*c/mu_sea;
Re.alt = rho_alt*v_alt*c/mu_sea;

M_0.sea = 0.5*rho_sea*v_sea.^2*S.*Cm.sea*c;                      %replace
M_0.alt = 0.5*rho_alt*v_alt.^2*S.*Cm.alt*c;   

i = 1; %select condition

% CL_seat(4) = -1.76*W/(0.5*rho_sea*v_sea(i)^2*S);
% CL_seat(5) = -1*W/(0.5*rho_sea*v_sea(i)^2*S);
% alpha.sea(4) = ((CL_seat(4)-C_L0_sea)/C_La_sea)*180/pi - atand(15.24/63.9); %degrees
% alpha.sea(5) = ((CL_seat(5)-C_L0_sea)/C_La_sea)*180/pi - atand(15.24/63.9); %degrees
% alpha.sea(4) = alpha.sea(4) - atand(15.24/63.9);
% alpha.sea(5) = alpha.sea(5) - atand(7.62/95.8);
% 
% CL_altt(4) = -1.76*W/(0.5*rho_alt*v_alt(i)^2*S);
% CL_altt(5) = -1*W/(0.5*rho_alt*v_alt(i)^2*S);
% alpha.alt(4) = ((CL_altt(4)-C_L0_alt)/C_La_alt)*180/pi - atand(15.24/63.9); %degrees
% alpha.alt(5) = ((CL_altt(5)-C_L0_alt)/C_La_alt)*180/pi - atand(15.24/63.9); %degrees
% alpha.alt(4) = alpha.alt(4) - atand(15.24/63.9);
% alpha.alt(5) = alpha.alt(5) - atand(7.62/95.8);

alpha.sea = [13.90527,4.10748,-16.457656,-31.4285542,-21.488200]; %degrees
alpha.alt = [13.42583,7.63933,-15.72596,-32.8453936,-22.294760]; %degrees

%% Structural Elements and Parameters
%box beam parameters
structure.kt = 0.001016;              % m   skin thickness
structure.sh = 0.08;                  % m
structure.st = 0.0025;                % m   spar thickness
structure.bl = 0.012;                 % m   bracket height
structure.bt = 0.0025;                % m   bracket thickness`
structure.E  = 70E9;                  % Pa  Young's Modulus
structure.G  = 28E9;                  % Pa
%wing = build_wing(c,structure,theta); % Build wing
%[wing_section_centroid, structure, component_moments] = calculate_geometry(wing,c,structure,theta);

%% Functions
[x, y, Ixx, Iyy, Ixy, skin, spar, str, caps, x_quarterchord] = build_airfoil();
[L_distribution, D_distribution] = LD_plots();
[Sx, Sy, Mx, My, sigma_z, wx_sea] = SMsigma_plots(alpha, x, y, z, Ixx, Iyy, Ixy, L_distribution, D_distribution);
[u, v] = uv_plots(z,Ixx, Iyy, Ixy, structure, My, Mx);
[Booms] = Booms(x, y, z, sigma_z, skin, str, caps, spar, dz);
[Term_2, A1, A2, A_total, qb] = Shear_Flow_Basic(x, y, z, Ixx, Iyy, Ixy, Booms, Sx, Sy, spar, dz);
[q, tau] = Shear_Flow_0(x, y, z, Ixx, Iyy, Ixy, Booms, Sx, Sy, A1, A2, qb, structure, Term_2, M_0, x_quarterchord, skin, str, caps, spar, dz); %% need M_0

