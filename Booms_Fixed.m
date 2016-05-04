% boom area = (skin thickness * distance between booms / 6) * (2 + stress2/stress1)

% need all x-values used to plot airfoil. need stringer area, spar cap area,
% and stringer, spar, and spar cap x-locations

x = airfoil.booms.x_coordinates;                %%
y = airfoil.booms.y_coordinates;


BOOM_A_ALT = zeros(length(x),length(z),5);
BOOM_A_SEA = zeros(length(x),length(z),5);

for k = 1:5
    for j = 1:length(z)    
        for i = 2:length(x)-1
            dist_rear = ((x(i)+x(i-1))^2 + (y(i)+y(i-1))^2)^0.5;
            dist_fore = ((x(i)+x(i+1))^2 + (y(i)+y(i+1))^2)^0.5;
            if ismember(x(i),stringers.x) && ismember(y(i),stringers.y) == 1;
                area_term = stringers.area;
            elseif ismember(x(i),caps.x) == 1;
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
        if ismember(x(1),stringers.x)==1
           BOOM_A_ALT(1,j,k) = BOOM_A_ALT(1,j,k) + stringers.area; 
           BOOM_A_SEA(1,j,k) = BOOM_A_SEA(1,j,k) + stringers.area;
        end
        if ismember(x(end),stringers.x)==1
           BOOM_A_ALT(end,j,k) = BOOM_A_ALT(end,j,k) + stringers.area; 
           BOOM_A_SEA(end,j,k) = BOOM_A_SEA(end,j,k) + stringers.area;
        end
       


%% Add the contribution of spars
%% i think this is wrong :(
        if_cap = ismember(x,caps.x);
        i_Bspar = find(if_cap);

        for i = 1:2
            BOOM_A_SEA(i_Bspar(i),j,k) = BOOM_A_SEA(i_Bspar(i),j,k) + kt * spars.height(i)/6 * (2 + sigma_z.sea(i_BsparL(i),j,k)/sigma_z.sea_upper(i_BsparU(i)));
            BOOM_A_ALT(i_Bspar(i),j,k) = BOOM_A_ALT(i_Bspar(i),j,j) + kt * spars.height(i)/6 * (2 + sigma_z.alt(i_BsparU(i),j,k)/sigma_z.alt_lower(i_BsparL(i)));
        end

    end
end


