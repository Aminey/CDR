function [yield] = Von_Mises(sigma_z, tau, x, dz, z)


%% Initialize
yield_stress = 0; % change
yield.alt = zeros(length(x), length(z)/dz, 5);
yield.sea = zeros(length(x), length(z)/dz, 5);


%% Check for Yielding at all points, along span, for each condition
for k = 1:5
    for j = 1:dz:length(z)
        for i = 1:length(x)
            if yield_stress <= ((( 2 * (sigma_z.alt(i,j,k))^2 + 6 * (tau.alt(i,j,k)^2))/2)^0.5)
                yield.alt(i,j,k) = 1;
            end
            if yield_stress <= ((( 2 * (sigma_z.sea(i,j,k))^2 + 6 * (tau.sea(i,j,k)^2))/2)^0.5)
                yield.sea(i,j,k) = 1;
            end
        end
        if yield_stress <= ((( 2 * (%sigma_z.sea_spar(j,k))^2 + 6 * (tau.sea_spar(j,k)^2))/2)^0.5)
            yield.sea_spar(j,k) = 1;
        end
        if yield_stress <= ((( 2 * (%sigma_z.alt_spar(j,k))^2 + 6 * (tau.alt_spar(j,k)^2))/2)^0.5)
            yield.alt_spar(j,k) = 1;
        end
        % need rear spar here           if yield_stress <= ((( 2 * (sigma_z.sea_spar(j,k))^2 + 6 * ((q02.sea(j,k)/spar.t)^2))/2)^0.5)
           % yield.sea_spar(j,k) = 1;
        end
    end
end


%% Plots
plot(1:length(x), yield.alt(:,1,:));
plot(1:length(x), yield.alt(:,1,:));

% need to plot yielding on spars as well

if ismember(1, yield);
    disp('Structure Yields')
else
    disp('Structure Does Not Yield')
end