function [x_c, y_c, Ixx, Iyy, Ixy, skin, spar, str, caps, x_quarterchord] = build_airfoil()


% clear all;
%% 2415 Airfoil 
m = 0.02;               % 
p = 0.4;                % max camber
t = 0.15;               % max thickness 
c = 1.5;                % chord length
nx = 500;               % number of increments
dx = c/nx;  
x = 0:dx:c;             % even spacing

%% Structural elements (arbitrary inputs for now)
skin.t = 0.008;                               % thickness of skin
spar.t = 0.01;                               % thickness of spar
str.A_upp = 5E-6;                            % area of each upper stringer
str.A_low = 5E-6;                            % area of each lower stringer
caps.A = 1E-5;                                % area of spar caps
spar.x = 0.3*1.5;                                                      % choose x location of spar
str.x_upp = [x(1), x(50), x(150), x(200), x(250),x(300), x(350)];     % choose x locations of upper stringers
str.x_low = [x(1), x(50), x(150), x(200), x(250),x(300), x(350)];     % choose x locations of lower stringers

%% Create Airfoil from equation
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
spar.x = [spar.x,x(end)];                    % x location of spar
spar.n = length(spar.x);                     % number of spars
spar.i = round(spar.x./dx)+1;                % index of spars  
spar.h = yU(spar.i) - yL(spar.i);            % height of spar
spar.Cy = (yU(spar.i) + yL(spar.i))/2;       % y centroid of spar  
spar.A = spar.t*spar.h;                      % area of spar
spar.y_upp = yU(spar.i);
spar.y_low = yL(spar.i);

%% Stringers
str.n_upp = length(str.x_upp);                            % number of upper stringers
str.i_upp = round(str.x_upp./dx)+1;                       % indices of upper stringers
str.n_low = length(str.x_low);                            % number of lower stringers
str.i_low = round(str.x_low./dx)+1;                       % indices of lower stringers
str.y_upp = yU(str.i_upp);
str.y_low = yL(str.i_low);

%% Nodes (based on TA's definition: node = skin and spar cap)
node.x_upp = [str.x_upp,spar.x];            % x location of upper stringers and spars
node.x_upp = sort(node.x_upp);              % sort ascending
% skin should be broken into smaller elements for higher accuracy of calculation
% break one skin element into two by adding one more node in between
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

% lower nodes
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
skin.Cy_low = zeros(1,skin.n_low);

for i = 1:skin.n_low
    skin.l_low(i) = sqrt((yL(skin.i_low(i+1)) - yL(skin.i_low(i)))^2 + (skin.x_low(i+1) - skin.x_low(i))^2);
    skin.A_low(i) = skin.t*skin.l_low(i);
    skin.Cx_low(i) = (skin.x_low(i+1) + skin.x_low(i))/2;
    skin.Cy_low(i) = (yL(skin.i_low(i+1)) + yL(skin.i_low(i)))/2;
end

%% Spar Caps
caps.x = spar.x;
caps.n = 2*length(caps.x);
caps.y_upp = yU(spar.i);
caps.y_low = yL(spar.i);
caps.i = round(caps.x/dx)+1;
% Centroids of caps are at the x and y locations: caps.Cy = caps.C,

%% Centroid of the overall wing section
% initial value
Cx_sum = 2*(sum(caps.x.*caps.A)) + sum(str.x_upp.*str.A_upp) + sum(str.x_low.*str.A_low) + sum(spar.x.*spar.A) ...
         + sum(skin.Cx_low.*skin.A_low) + sum(skin.Cx_upp.*skin.A_upp);
Cy_sum = sum(caps.y_low.*caps.A) + sum(caps.y_upp.*caps.A) + sum(yU(str.i_upp).*str.A_upp) + sum(yL(str.i_low).*str.A_low) ... 
         + sum(spar.Cy.*spar.A) + sum(skin.Cy_low.*skin.A_low) + sum(skin.Cy_upp.*skin.A_upp);
A_sum = caps.A*caps.n + str.A_upp*str.n_upp + str.A_low*str.n_low + sum(spar.A) + sum(skin.A_low) + sum(skin.A_upp);

Cx = Cx_sum/A_sum;
Cy = Cy_sum/A_sum;

%% Leading Edge Origin Plot
% figure
% plot(x,yU,'k',x,yL,'k','Linewidth',2);
% ylim([-0.3 0.3])
% hold on
% xlim([0,c]);
% plot(str.x_upp,yU(str.i_upp),'or',str.x_low,yL(str.i_low),'or','markersize',5);
% % plot([x(1),x(1)],[yU(i(1)),yL(i(1))],'b',[x(end),x(end)],[yU(end),yL(end)],'b','Linewidth',3);
% rectangle('Position',[x(spar.i(1))-spar.t/2, yL(spar.i(1)), spar.t, yU(spar.i(1))-yL(spar.i(1))]);
% rectangle('Position',[x(spar.i(2))-spar.t/2, yL(spar.i(2)), spar.t, yU(spar.i(2))-yL(spar.i(2))]);
% plot(spar.x,yU(spar.i),'sg',spar.x,yL(spar.i),'sg','markersize',7);
% % plot([x(end),x(end)],[yU(end),yL(end)]);
% scatter(Cx,Cy,'m*')
% ylabel('y (m)')
% xlabel('x (m)')
% grid on

%% Area moments of inertia
for i = 1:length(spar.x)
    spar.Ixx(i) = (spar.t*spar.h(i)^3)/12 + spar.A(i)*(spar.Cy(i) - Cy)^2;
    spar.Iyy(i) = (spar.h(i)*spar.t^3)/12 + spar.A(i)*(spar.x(i) - Cx)^2;
    spar.Ixy(i) = spar.A(i)*(spar.Cy(i) - Cy)*(spar.x(i) - Cx);
end

for i = 1:skin.n_upp
    skin.Ixx_upp(i) = skin.A_upp(i)*(skin.Cy_upp(i) - Cy)^2;
    skin.Iyy_upp(i) = skin.A_upp(i)*(skin.Cx_upp(i) - Cx)^2;
    skin.Ixy_upp(i) = skin.A_upp(i)*(skin.Cy_upp(i) - Cy)*(skin.Cx_upp(i) - Cx);
end

for i = 1:skin.n_low
    skin.Ixx_low(i) = skin.A_low(i)*(skin.Cy_low(i) - Cy)^2;
    skin.Iyy_low(i) = skin.A_low(i)*(skin.Cx_low(i) - Cx)^2;
    skin.Ixy_low(i) = skin.A_low(i)*(skin.Cy_low(i) - Cy)*(skin.Cx_low(i) - Cx);
end

for i = 1:str.n_upp
    str.Ixx_upp(i) = str.A_upp*(yU(str.i_upp(i)) - Cy)^2;
    str.Iyy_upp(i) = str.A_upp*(str.x_upp(i) - Cx)^2;
    str.Ixy_upp(i) = str.A_upp*(yU(str.i_upp(i)) - Cy)*(str.x_upp(i) - Cx);
end

for i = 1:str.n_low
    str.Ixx_low(i) = str.A_low*(yL(str.i_low(i)) - Cy)^2;
    str.Iyy_low(i) = str.A_low*(str.x_low(i) - Cx)^2;
    str.Ixy_low(i) = str.A_low*(yL(str.i_low(i)) - Cy)*(str.x_low(i) - Cx);
end

for i = 1:length(caps.x)
    caps.Ixx_upp(i) = caps.A*(caps.y_upp(i) - Cy)^2;
    caps.Iyy_upp(i) = caps.A*(caps.x(i) - Cx)^2;
    caps.Ixy_upp(i) = caps.A*(caps.y_upp(i) - Cy)*(caps.x(i) - Cx);
end

for i = 1:length(caps.x)
    caps.Ixx_low(i) = caps.A*(caps.y_low(i) - Cy)^2;
    caps.Iyy_low(i) = caps.A*(caps.x(i) - Cx)^2;
    caps.Ixy_low(i) = caps.A*(caps.y_low(i) - Cy)*(caps.x(i) - Cx);
end

Ixx = sum(spar.Ixx)+sum(skin.Ixx_upp)+sum(skin.Ixx_low)+sum(str.Ixx_upp)+sum(str.Ixx_low)+sum(caps.Ixx_low);
Iyy = sum(spar.Iyy)+sum(skin.Iyy_upp)+sum(skin.Iyy_low)+sum(str.Iyy_upp)+sum(str.Iyy_low)+sum(caps.Iyy_low);
Ixy = sum(spar.Ixy)+sum(skin.Ixy_upp)+sum(skin.Ixy_low)+sum(str.Ixy_upp)+sum(str.Ixy_low)+sum(caps.Ixy_low);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% BELOW CALCULATIONS ARE COORDINATE TRANSFORMATIONS FOR JARED's CODE %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Shift coordinate system to centroid-centered frame
% x's
spar.x_c = spar.x - Cx;
skin.x_upp_c = skin.Cx_upp - Cx;
skin.x_low_c = skin.Cx_low - Cx;
str.x_upp_c = str.x_upp - Cx;
str.x_low_c = str.x_low - Cx;
caps.x_upp_c = caps.x - Cx;
caps.x_low_c = caps.x - Cx;

% y's (Cy's are y's for point-components)
skin.y_upp_c = skin.Cy_upp - Cy;
skin.y_low_c = skin.Cy_low - Cy;
str.y_upp_c = str.y_upp - Cy;
str.y_low_c = str.y_low - Cy;
caps.y_upp_c = caps.y_upp - Cy;
caps.y_low_c = caps.y_low - Cy;

%% Make CCW coordinate system for shear flow
x_CCW = zeros(1,length(x));
yU_CCW = zeros(1,length(x));
yL_CCW = zeros(1,length(x));

for i = 1:length(x)
    x_CCW(i) = x(end - i + 1);
    yU_CCW(i) = yU(end - i + 1);
    yL_CCW(i) = yL(end - i + 1);
end
y_CCW = [yU_CCW yL_CCW(end-1:-1:1)];         
x_CCW = [x_CCW x(2:end)];

spar.i_CCW = length(x) - spar.i + 1;
spar.i_CCW = sort([spar.i_CCW, (length(x) + spar.i - 1)]); 
skin.i_CCW = length(x) - skin.i_upp + 1;
skin.i_CCW = sort([skin.i_CCW, (length(x) + skin.i_low - 1)]);
str.i_CCW  = length(x) - str.i_upp + 1;
str.i_CCW  = sort([str.i_CCW, (length(x) + str.i_low - 1)]);
caps.i_CCW = spar.i_CCW;

x_c = x_CCW - Cx;
y_c = y_CCW - Cy;

x_quarterchord = c/4 - Cx;


%% Centroid-centered plot
figure
hold on
plot(x-Cx,yU-Cy,'k',x-Cx,yL-Cy,'k','Linewidth',2);
ylim([-0.3 0.3])
xlim([1.1*(-Cx),1.1*(x(end)-Cx)]);
p.skin_upp = plot(skin.x_low_c, skin.y_low_c, 'g.', 'markersize',8);
p.skin_low = plot(skin.x_upp_c, skin.y_upp_c, 'g.','markersize',8);
p.str_upp  = plot(str.x_upp_c, str.y_upp_c, 'or', 'markersize',5);
p.str_low  = plot(str.x_low_c, str.y_low_c, 'or', 'markersize',5);
p.spar1    = plot(spar.x_c(1), 0.5*(spar.y_upp(1) + spar.y_low(1)) - Cy, 'bx', 'markersize', 8);
p.spar2    = plot(spar.x_c(2), 0.5*(spar.y_upp(2) + spar.y_low(2)) - Cy, 'bx', 'markersize', 8);
p.caps_upp = plot(caps.x_upp_c, caps.y_upp_c, 'sm', 'markersize', 10);
p.caps_low = plot(caps.x_low_c, caps.y_low_c, 'sm', 'markersize', 10);
rectangle('Position',[x(spar.i(1))-spar.t/2-Cx, yL(spar.i(1))-Cy, spar.t, yU(spar.i(1))-yL(spar.i(1))]);
rectangle('Position',[x(spar.i(2))-spar.t/2-Cx, yL(spar.i(2))-Cy, spar.t, yU(spar.i(2))-yL(spar.i(2))]);
p.centroid = plot(0,0,'r*', 'markersize', 10);
legend([p.centroid, p.skin_upp, p.str_upp, p.spar1, p.caps_upp],'Overall Centroid', 'Skin', 'Stringer', 'Spar', 'Spar Cap');
ylabel('y (m)')
xlabel('x (m)')
grid on

disp('build_airfoil complete');