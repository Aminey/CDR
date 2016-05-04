function [airfoil] = build_airfoil()

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

    
    
end