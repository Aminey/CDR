% boom area = (skin thickness * distance between booms / 6) * (2 + stress2/stress1)

% need all x-values used to plot airfoil. need stringer area, spar cap area,
% and stringer, spar, and spar cap x-locations
x = zeros(1,length(x_CCW));
y = zeros(1,length(y_CCW));

for i = 1:length(x_CCW)
    x(i) = x_CCW(i);                %%
    y(i) = y_CCW(i);
end

Booms.alt = zeros(length(x),length(z),5);
Booms.sea = zeros(length(x),length(z),5);

for k = 1:5
    for j = 1:length(z)    
        for i = 2:length(x)-1
            dist_rear = ((x(i)+x(i-1))^2 + (y(i)+y(i-1))^2)^0.5;
            dist_fore = ((x(i)+x(i+1))^2 + (y(i)+y(i+1))^2)^0.5;
            if airfoil.booms(i).isStringer
                area_term = str.A;
            elseif airfoil.booms(i).isCap
                area_term = caps.A;
            else
                area_term = 0;
            end
            Booms.alt(i,j,k) = (kt * dist_rear/6) * (2 + (sigma_z.alt(i-1,j,k)/sigma_z.alt(i,j,k))) + ...
                          (kt * dist_fore/6) * (2 + (sigma_z.alt(i+1,j,k)/sigma_z.alt(i,j,k))) + area_term;
            Booms.sea(i,j,k) = (kt * dist_rear/6) * (2 + (sigma_z.sea(i-1,j,k)/sigma_z.sea(i,j,k))) + ...
                          (kt * dist_fore/6) * (2 + (sigma_z.sea(i+1,j,k)/sigma_z.sea(i,j,k))) + area_term;
        end
        
%% First Coordinate
        Booms.alt(1,j,k) = kt * (((x(end)-x(1))^2 + (y(end)-y(1))^2)^0.5)/6 *...
                              (2 + sigma_z.alt(end,j,k)/sigma_z.alt(1,j,k))+...
                              kt * (((x(2)-x(1))^2 + (y(2)-y(1))^2)^0.5)/6*...
                              (2 + sigma_z.alt(2,j,k)/sigma_z.alt(1,j,k));
        Booms.sea(1,j,k) = kt * (((x(end)-x(1))^2 + (y(end)-y(1))^2)^0.5)/6 *...
                              (2 + sigma_z.sea(end,j,k)/sigma_z.sea(1,j,k))+...
                              kt * (((x(2)-x(1))^2 + (y(2)-y(1))^2)^0.5)/6*...
                              (2 + sigma_z.sea(2,j,k)/sigma_z.sea(1,j,k));
       

%% Last Coordinate
        Booms.alt(end,j,k) = kt * (((x(end)-x(end-1))^2 + (y(end)-y(end-1))^2)^0.5)/6 *...
                              (2 + sigma_z.alt(end-1,j,k)/sigma_z.alt(end,j,k))+...
                              kt * (((x(end)-x(1))^2 + (y(end)-y(1))^2)^0.5)/6*...
                              (2 + sigma_z.alt(1,j,k)/sigma_z.alt(end,j,k));
        Booms.sea(end,j,k) = kt * (((x(end)-x(end-1))^2 + (y(end)-y(end-1))^2)^0.5)/6 *...
                              (2 + sigma_z.sea(end-1,j,k)/sigma_z.sea(end,j,k))+...
                              kt * (((x(end)-x(1))^2 + (y(end)-y(1))^2)^0.5)/6*...
                              (2 + sigma_z.sea(1,j,k)/sigma_z.sea(end,j,k));
       
       

%% Stringers on first and last coordinates?
        if airfoil.booms(1).isStringer
           Booms.alt(1,j,k) = Booms.alt(1,j,k) + str.A; 
           Booms.sea(1,j,k) = Booms.sea(1,j,k) + str.A;
        end
        if airfoil.booms(end).isStringer
           Booms.alt(end,j,k) = Booms.alt(end,j,k) + str.A; 
           Booms.sea(end,j,k) = Booms.sea(end,j,k) + str.A;
        end
       


%% Add the contribution of spar's path for shear flow
    m = 1;
    spar_holder = zeros(1,4);  
    
    for i = 1:length(x_CCW)
        if airfoil.booms(i).isSpars
           spar_holder(m) = i;
           m = m+1;
        end
    end
    
    for m = 1:4   % two spars will have four connection points. going counterclockwise, connection 1 will be the top of the rear 
                  % spar and 4 will be the bottom. 2 and 3 will be the top
                  % and bottom of the forward spar. thus, index (5-m) is
                  % used to relate spar connections 1 and 4 and connections 2 and 3
         Booms.sea(spar_holder(m),j,k) = Booms.sea(spar_holder(m),j,k) + kt * spars.height(1)/6 * ...
                                                        (2 + sigma_z.sea(spar_holder(5-m),j,k)/sigma_z.sea(spar_holder(m),j,k));



         Booms.alt(spar_holder(m),j,k) = Booms.alt(spar_holder(m),j,k) + kt * spars.height(1)/6 * ...
                                                        (2 + sigma_z.alt(spar_holder(5-m),j,k)/sigma_z.alt(spar_holder(5-m),j,k));

    end
    end
end


