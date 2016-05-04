% boom area = (skin thickness * distance between booms / 6) * (2 + stress2/stress1)

% need all x-values used to plot airfoil. need stringer area, spar cap area,
% and stringer, spar, and spar cap x-locations

BOOM_A_ALT_U = zeros(length(x),length(z),5);
BOOM_A_SEA_U = zeros(length(x),length(z),5);
BOOM_A_ALT_L = zeros(length(x),length(z),5);
BOOM_A_SEA_L = zeros(length(x),length(z),5);


for k = 1:5
    for j = 1:length(z)
        
        %%Upper Surface        
        for i = 2:length(x)-1
            dist_rear = ((x(i)+x(i-1))^2 + (yU(i)+yU(i-1))^2)^0.5;
            dist_fore = ((x(i)+x(i+1))^2 + (yU(i)+yU(i+1))^2)^0.5;
            if ismember(x(i),x_stringersU) == 1;
                area_term = stringer_area;
            elseif ismember(x(i),x_sparsU) == 1;
                area_term = spar_area;
            else
                area_term = 0;
            end
            BOOM_A_ALT_U(i,j,k) = (kt * dist_rear/6) * (2 + (sigma_z.alt_upper(i-1,j,k)/sigma_z.alt_upper(i,j,k))) + ...
                          (kt * dist_fore/6) * (2 + (sigma_z.alt_upper(i+1,j,k)/sigma_z.alt_upper(i,j,k))) + area_term;
            BOOM_A_SEA_U(i,j,k) = (kt * dist_rear/6) * (2 + (sigma_z.sea_upper(i-1,j,k)/sigma_z.sea_upper(i,j,k))) + ...
                          (kt * dist_fore/6) * (2 + (sigma_z.sea_upper(i+1,j,k)/sigma_z.sea_upper(i,j,k))) + area_term;
        end

        BOOM_A_ALT_U(end,j,k) = kt * (((x(end)+x(end-1))^2 + (yU(end)+yU(end-1))^2)^0.5) / 6 * (2 + sigma_z.alt_upper(end-1,j,k)/sigma_z.alt_upper(end,j,k));
        BOOM_A_SEA_U(end,j,k) = kt * (((x(end)+x(end-1))^2 + (yU(end)+yU(end-1))^2)^0.5) / 6 * (2 + sigma_z.sea_upper(end-1,j,k)/sigma_z.sea_upper(end,j,k));

        %% Lower Surface

        for i = 2:length(x)-1
            dist_rear = ((x(i)+x(i-1))^2 + (yU(i)+yU(i-1))^2)^0.5;
            dist_fore = ((x(i)+x(i+1))^2 + (yU(i)+yU(i+1))^2)^0.5;
            if ismember(x(i),x_stringersL) == 1;
                area_term = stringer_area;
            elseif ismember(x(i),x_sparsL) == 1;
                area_term = spar_area;
            else
                area_term = 0;
            end
            BOOM_A_ALT_L(i,j,k) = (kt * dist_rear/6) * (2 + (sigma_z.alt_lower(i-1,j,k)/sigma_z.alt_lower(i,j,k))) + ...
                          (kt * dist_fore/6) * (2 + (sigma_z.alt_lower(i+1,j,k)/sigma_z.alt_lower(i,j,k))) + area_term;
            BOOM_A_SEA_L(i,j,k) = (kt * dist_rear/6) * (2 + (sigma_z.sea_lower(i-1,j,k)/sigma_z.sea_lower(i,j,k))) + ...
                          (kt * dist_fore/6) * (2 + (sigma_z.sea_lower(i+1,j,k)/sigma_z.sea_lower(i,j,k))) + area_term;
        end

        BOOM_A_ALT_L(end,j,k) = kt * (((x(end)+x(end-1))^2 + (yL(end)+yL(end-1))^2)^0.5) / 6 * (2 + sigma_z.alt_lower(end-1,j,k)/sigma_z.alt_lower(end,j,k));
        BOOM_A_SEA_L(end,j,k) = kt * (((x(end)+x(end-1))^2 + (yL(end)+yL(end-1))^2)^0.5) / 6 * (2 + sigma_z.sea_lower(end-1,j,k)/sigma_z.sea_lower(end,j,k));

        %% Single Shared Boom between upper and lower surfaces
        BOOM_A_ALT_U(1,j,k) = kt * (((x(2)+x(1))^2 + (yU(2)+yU(1))^2)^0.5)/6 * (2 + sigma_z.alt_upper(2,j,k)/sigma_z.alt_upper(1,j,k))...
                            + kt * (((x(2)+x(1))^2 + (yL(2)+yL(1))^2)^0.5)/6 * (2 +sigma_z.alt_lower(2,j,k)/sigma_z.alt_lower(1,j,k));
        BOOM_A_SEA_U(1,j,k) = kt * (((x(2)+x(1))^2 + (yU(2)+yU(1))^2)^0.5)/6 * (2 + sigma_z.sea_upper(2,j,k)/sigma_z.sea_upper(1,j,k))...
                            + kt * (((x(2)+x(1))^2 + (yL(2)+yL(1))^2)^0.5)/6 * (2 +sigma_z.sea_lower(2,j,k)/sigma_z.sea_lower(1,j,k));                

        if ismember(x(1),x_stringersU)==1
           BOOM_A_ALT_U(1,j,k) = BOOM_A_ALT_U(1,j,k) + stringer_area; 
        end
        BOOM_A_ALT_L(1,j,k) = BOOM_A_ALT_U(1,j,k);


        if ismember(x(1),x_stringersU)==1
           BOOM_A_SEA_U(1,j,k) = BOOM_A_SEA_U(1,j,k) + stringer_area; 
        end
        BOOM_A_SEA_L(1,j,k) = BOOM_A_SEA_U(1,j,k);


        % add the contribution of spars

        if_capU = ismember(x,x_sparsU);
        if_capL = ismember(x,x_sparsL);
        i_BsparU = find(if_capU);
        i_BsparL = find(if_capL);

        for i = 1:2
            BOOM_A_SEA_U(i_BsparU(i),j,k) = BOOM_A_SEA_U(i_BsparU(i),j,k) + kt * h_spars(i)/6 * (2 + sigma_z.sea_lower(i_BsparL(i),j,k)/sigma_z.sea_upper(i_BsparU(i)));
            BOOM_A_SEA_L(i_BsparL(i),j,k) = BOOM_A_SEA_L(i_BsparL(i),j,k) + kt * h_spars(i)/6 * (2 + sigma_z.sea_upper(i_BsparU(i),j,k)/sigma_z.sea_lower(i_BsparL(i)));
            BOOM_A_ALT_U(i_BsparU(i),j,k) = BOOM_A_ALT_U(i_BsparU(i),j,k) + kt * h_spars(i)/6 * (2 + sigma_z.alt_lower(i_BsparL(i),j,k)/sigma_z.alt_upper(i_BsparU(i)));
            BOOM_A_SEA_L(i_BsparL(i),j,k) = BOOM_A_ALT_L(i_BsparL(i),j,j) + kt * h_spars(i)/6 * (2 + sigma_z.alt_upper(i_BsparU(i),j,k)/sigma_z.alt_lower(i_BsparL(i)));
        end

    end
end


