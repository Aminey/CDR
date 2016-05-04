% function [airfoil] = build_airfoil()

    %read NACA 2415 coordinates from excel file
    x_coordinates = 1.5*xlsread('NACA 2415.xlsx','A1:A199');
    y_coordinates = 1.5*xlsread('NACA 2415.xlsx','B1:B199');
    
    airfoil.booms.x_coordinates = x_coordinates;
    airfoil.booms.y_coordinates = y_coordinates;
    
%     figure;
%     plot(x_coordinates,y_coordinates);
%     axis equal;
    
    %determine number of points
    n_points = length(x_coordinates);

    %
    A_cap = ;
    A_str = ;
    t_spar = ;
    t_skin = ;
    % locations of spars, spar caps and stringers (nose at the origin of the coordinate)
    x_spar0 = ;                       % front spar (2 cell beam)
    x_strU0 = [];                     % upper surface
    x_strL0 = [];                     % lower surface
    % new coordinate with origin at the centroid is used for the output below
    [c,Ixx,Iyy,Ixy,x,yU,yL,x_strU,x_strL,x_boomU,x_boomL,L_boomU,L_boomL,x_spar,h_spar,i_spar,dx] = airfoil_section(A_cap,A_str,t_spar,t_skin,x_spar0,x_strU0,x_strL0);


    
    
% end