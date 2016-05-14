% function [c,Ixx,Iyy,Ixy,x,yU,yL,x_strU,x_strL,x_skinU,x_skinL,L_skinU,L_skinL,x_spar,h_spar,i_spar,dx] = airfoil_section(A_cap,A_str,t_spar,t_skin,x_spar,x_strU,x_strL)
%% airfoil section profile
% NACA 2415
close all; clc; clear all;

m = 0.02;               % 
p = 0.4;                % max camber
t = 0.15;               % max thickness 
c = 1.5;                % chord length
nx = 500;               % number of increments
dx = c/nx;  
x = 0:dx:c;             % even spacing

yc = zeros(1,nx+1);
yt = zeros(1,nx+1);
yU0 = zeros(1,nx+1);
yL0 = zeros(1,nx+1);
theta = zeros(1,nx+1);
xb = p*c;
i_xb = xb/dx + 1;  

% NACA equation source: http://www.aerospaceweb.org/question/airfoils/q0041.shtml    
for i = 1:nx+1
    yt(i) = 5*t*c*(0.2969*sqrt(x(i)/c) - 0.1260*(x(i)/c) - 0.3516*(x(i)/c)^2 + 0.2843*(x(i)/c)^3 - 0.1015*(x(i)/c)^4);
    if i <= i_xb
        yc(i) = m*x(i)/p^2*(2*p - x(i)/c);
        theta(i) = atan(2*m/p^2*(p - x(i)/c));
    else
        yc(i) = m*(c - x(i))/(1-p)^2*(1 + x(i)/c - 2*p);
        theta(i) = atan(2*m/(1-p)^2*(p-x(i)/c));
    end
    yU0(i) = yc(i) + yt(i)*cos(theta(i));
    yL0(i) = yc(i) - yt(i)*cos(theta(i));
end

% The last 20% of the chord length of the airfoil was neglected under the assumption that
% this section contained flaps and ailerons, and would therefore not support aerodynamic loads.
i_xend = round(0.8*c/dx)+1;
x = x(1:i_xend);
yU = yU0(1:i_xend);
yL = yL0(1:i_xend);

%% Spars
spar.x = [x(20), x(30)];                     % choose x locations
spar.t = 0.5;                                % thickness of spar
spar.x = [spar.x,x(end)];                    % x location of spar
spar.n = length(spar.x);                     % number of spars
spar.i = round(spar.x./dx)+1;                % index of spars  
spar.h = yU(spar.i) - yL(spar.i);            % height of spar
spar.Cy = (yU(spar.i) + yL(spar.i))/2;       % y centroid of spar  
spar.A = spar.t*spar.h;                      % area of spar


%% Stringers
str.x_upp = [x(5), x(10), x(15), x(20), x(25),x(30)];     % x locations of upper stringers
str.n_upp = length(str.x_upp);                            % number of upper stringers
str.i_upp = round(str.x_upp./dx)+1;                       % indices of upper stringers
str.A_upp = 0.5;                                          % area of each upper stringer

str.x_low = [x(5), x(10), x(15), x(20), x(25),x(30)];     % x locations of lower stringers
str.n_low = length(str.x_low);                            % number of lower stringers
str.i_low = round(str.x_low./dx)+1;                       % indices of lower stringers
str.A_low = 0.5;                                          % area of each lower stringer


%% Nodes and Skins
% Nodes include spar caps and stringers
node.x_upp = [str.x_upp,spar.x];            % x location of upper stringers and spars
node.x_upp = sort(node.x_upp);              % sort ascending

% skin should be broken into smaller elements for higher accuracy of calculation
% break one skin element into two by adding one more node in between

skin.t = 5;                                  % thickness of skin
skin.x_upp = zeros(1,2*length(node.x_upp)-1);
for i = 1:length(node.x_upp)-1
    skin.x_upp(2*i-1) = node.x_upp(i);
    skin.x_upp(2*i) = (node.x_upp(i) + node.x_upp(i+1))/2;
end
skin.x_upp(end) = node.x_upp(end);

skin.i_upp = round(skin.x_upp/dx)+1;       % indices of skins
skin.n_upp = length(skin.x_upp)-1;         % number of upper skins
skin.l_upp = zeros(1,skin.n_upp);          % length of upper skins
skin.A_upp = zeros(1,skin.n_upp);          % area of upper skins
skin.Cx_upp = zeros(1,skin.n_upp);         % x centroid of upper skins
skin.Cy_upp = zeros(1,skin.n_upp);         % y centroid of upper skins

for i = 1:skin.n_upp
    skin.l_upp(i) = sqrt((yU(skin.i_upp(i+1)) - yU(skin.i_upp(i)))^2 + (skin.x_upp(i+1) - skin.x_upp(i))^2);  
    skin.A_upp(i) = skin.t*skin.l_upp(i);
    skin.Cx_upp(i) = (skin.x_upp(i+1) + skin.x_upp(i))/2;
    skin.Cy_upp(i) = (yU(skin.i_upp(i+1)) + yU(skin.i_upp(i)))/2;
end

% lower part
node.x_low = [str.x_low,spar.x];
node.x_low = sort(node.x_low);

skin.x_low = zeros(1,2*length(node.x_low)-1);
for i = 1:length(node.x_low)-1
    skin.x_low(2*i-1) = node.x_low(i);
    skin.x_low(2*i) = (node.x_low(i) + node.x_low(i+1))/2;
end
skin.x_low(end) = node.x_low(end);

skin.i_low = round(skin.x_low/dx)+1;
skin.n_low = length(skin.x_low)-1;
skin.l_low = zeros(1,skin.n_low);
skin.A_low = zeros(1,skin.n_low);
skin.Cx_low = zeros(1,skin.n_low);
skinL.Cy_low = zeros(1,skin.n_low);

for i = 1:skin.n_low
    skin.l_low(i) = sqrt((yL(skin.i_low(i+1)) - yL(skin.i_low(i)))^2 + (skin.x_low(i+1) - skin.x_low(i))^2);
    skin.A_low(i) = skin.t*skin.l_low(i);
    skin.Cx_low(i) = (skin.x_low(i+1) + skin.x_low(i))/2;
    skin.Cy_low(i) = (yL(skin.i_low(i+1)) + yL(skin.i_low(i)))/2;
end

%% centroid of the wing section
% initial value
Cx_sum = sum(str.x_upp) + sum(str.x_low) + sum(spar.x) + sum(skin.Cx_low) + sum(skin.Cx_upp) ;
Cy_sum = sum(yU(str.i_upp)) + sum(yL(str.i_low)) + sum(spar.Cy) + sum(skin.Cy_low) + sum(skin.Cy_upp);
A_sum = str.A_upp*str.n_upp + str.A_low*str.n_low + sum(spar.A) + sum(skin.A_low) + sum(skin.A_upp);

% spars
% for i = 1:n
%     Cx_sum = Cx_sum + x(i)*area(i);
%     Cy_sum = Cy_sum + Cy(i)*area(i);
%     A_sum = A_sum + area(i);
% end

%
%
% please finish the part for stringers, spar caps and skins
%
%

Cx = Cx_sum/A_sum;
Cy = Cy_sum/A_sum;

figure
plot(x,yU,'k',x,yL,'k','Linewidth',2);
ylim([-0.3 0.3])
hold on
xlim([0,c]);
plot(str.x_upp,yU(str.i_upp),'or',str.x_low,yL(str.i_low),'or','markersize',5);
plot([x(1),x(1)],[yU(i(1)),yL(i(1))],'b',[x(end),x(end)],[yU(end),yL(end)],'b','Linewidth',3);
% plot(x,yU(i),'sg',x,yL(i),'sg','markersize',7);
% plot([x(end),x(end)],[yU(end),yL(end)]);
scatter(Cx,Cy,'m*')
ylabel('y (m)')
xlabel('x (m)')
grid on

%% Area moments of inertia
% initial value
Ixx = 0;
Iyy = 0;
Ixy = 0;

% spars
% for i = 1:n
%     Ixx = Ixx + spar.t*height(i)^3/12 + area(i)*(Cy(i)-Cy)^2;
%     Iyy = Iyy + spar.t^3*height(i)/12 + area(i)*(x(i)-Cx)^2;
%     Ixy = Ixy + area(i)*(Cy(i)-Cy)*(x(i)-Cx);
% end

%
%
% please add contrbution of spar caps, stringers and skin panels. You can use the previous code as a reference.
%
%

%% coordiante transformation: new origin at centroid
%
%
% You need to update x, yU, yL, x_skinU, x_skinL, x_spar, x_strU, x_strL
%
%





