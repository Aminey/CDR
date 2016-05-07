% boom area = (skin thickness * distance between booms / 6) * (2 + stress2/stress1)

% need all x-values used to plot airfoil. need stringer area, spar cap area,
% and stringer, spar, and spar cap x-locations
x = zeros(1,length(airfoil.booms));
y = zeros(1,length(airfoil.booms));

for i = 1:length(airfoil.booms)
    x(i) = airfoil.booms(i).x_coordinate;                %%
    y(i) = airfoil.booms(i).y_coordinate;
end

BOOM_A_ALT = zeros(length(x),length(z),5);
BOOM_A_SEA = zeros(length(x),length(z),5);

for k = 1:5
    for j = 1:length(z)    
        for i = 2:length(x)-1
            dist_rear = ((x(i)+x(i-1))^2 + (y(i)+y(i-1))^2)^0.5;
            dist_fore = ((x(i)+x(i+1))^2 + (y(i)+y(i+1))^2)^0.5;
            if airfoil.booms(i).isStringer
                area_term = stringers.area;
            elseif airfoil.booms(i).isCap
                area_term = caps.area;
            else
                area_term = 0;
            end
            BOOM_A_ALT(i,j,k) = (kt * dist_rear/6) * (2 + (sigma_z.alt(i-1,j,k)/sigma_z.alt(i,j,k))) + ...
                          (kt * dist_fore/6) * (2 + (sigma_z.alt(i+1,j,k)/sigma_z.alt(i,j,k))) + area_term;
            BOOM_A_SEA(i,j,k) = (kt * dist_rear/6) * (2 + (sigma_z.sea(i-1,j,k)/sigma_z.sea(i,j,k))) + ...
                          (kt * dist_fore/6) * (2 + (sigma_z.sea(i+1,j,k)/sigma_z.sea(i,j,k))) + area_term;
        end
        
%% First Coordinate
        BOOM_A_ALT(1,j,k) = kt * (((x(end)-x(1))^2 + (y(end)-y(1))^2)^0.5)/6 *...
                              (2 + sigma_z.alt(end,j,k)/sigma_z.alt(1,j,k))+...
                              kt * (((x(2)-x(1))^2 + (y(2)-y(1))^2)^0.5)/6*...
                              (2 + sigma_z.alt(2,j,k)/sigma_z.alt(1,j,k));
        BOOM_A_SEA(1,j,k) = kt * (((x(end)-x(1))^2 + (y(end)-y(1))^2)^0.5)/6 *...
                              (2 + sigma_z.sea(end,j,k)/sigma_z.sea(1,j,k))+...
                              kt * (((x(2)-x(1))^2 + (y(2)-y(1))^2)^0.5)/6*...
                              (2 + sigma_z.sea(2,j,k)/sigma_z.sea(1,j,k));
       

%% Last Coordinate
        BOOM_A_ALT(end,j,k) = kt * (((x(end)-x(end-1))^2 + (y(end)-y(end-1))^2)^0.5)/6 *...
                              (2 + sigma_z.alt(end-1,j,k)/sigma_z.alt(end,j,k))+...
                              kt * (((x(end)-x(1))^2 + (y(end)-y(1))^2)^0.5)/6*...
                              (2 + sigma_z.alt(1,j,k)/sigma_z.alt(end,j,k));
        BOOM_A_SEA(end,j,k) = kt * (((x(end)-x(end-1))^2 + (y(end)-y(end-1))^2)^0.5)/6 *...
                              (2 + sigma_z.sea(end-1,j,k)/sigma_z.sea(end,j,k))+...
                              kt * (((x(end)-x(1))^2 + (y(end)-y(1))^2)^0.5)/6*...
                              (2 + sigma_z.sea(1,j,k)/sigma_z.sea(end,j,k));
       
       

%% Stringers on first and last coordinates?
        if airfoil.booms(1).isStringer
           BOOM_A_ALT(1,j,k) = BOOM_A_ALT(1,j,k) + stringers.area; 
           BOOM_A_SEA(1,j,k) = BOOM_A_SEA(1,j,k) + stringers.area;
        end
        if airfoil.booms(end).isStringer
           BOOM_A_ALT(end,j,k) = BOOM_A_ALT(end,j,k) + stringers.area; 
           BOOM_A_SEA(end,j,k) = BOOM_A_SEA(end,j,k) + stringers.area;
        end
       


%% Add the contribution of spar's path for shear flow
    m = 1;
    spar_holder = zeros(1,4);  
    
    for i = 1:length(airfoil.booms)
        if airfoil.booms(i).isSpars
           spar_holder(m) = i;
        end
    end
    
    for m = 1:4   % two spars will have four connection points. going counterclockwise, connection 1 will be the top of the rear 
                  % spar and 4 will be the bottom. 2 and 3 will be the top
                  % and bottom of the forward spar. thus, index (5-m) is
                  % used to relate spar connections 1 and 4 and connections 2 and 3
         BOOM_A_SEA(spar_holder(m),j,k) = BOOM_A_SEA(spar_holder(m),j,k) + kt * spars.height(1)/6 * ...
                                                        (2 + sigma_z.sea(spar_holder(5-m),j,k)/sigma_z.sea(spar_holder(m),j,k));



         BOOM_A_ALT(spar_holder(m),j,k) = BOOM_A_ALT(spar_holder(m),j,k) + kt * spars.height(1)/6 * ...
                                                        (2 + sigma_z.alt(spar_holder(5-m),j,k)/sigma_z.alt(spar_holder(5-m),j,k));

    end
    end
end


