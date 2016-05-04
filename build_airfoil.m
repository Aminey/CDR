function [airfoil] = build_airfoil()

    %read NACA 2415 coordinates from excel file
    x_coordinates = 1.5*xlsread('NACA 2415.xlsx','A1:A199');
    y_coordinates = 1.5*xlsread('NACA 2415.xlsx','B1:B199');
    
    for i = 1:length(x_coordinates)
        airfoil.booms(i).x_coordinate = x_coordinates(i);
        airfoil.booms(i).y_coordinate = y_coordinates(i);
        airfoil.booms(i).isStringer = false;
        airfoil.booms(i).isSpar = false;
        airfoil.booms(i).isSparCap = false;
    end
    
    %set stringer, spar, spar cap locations
    %...
    %    
    
end