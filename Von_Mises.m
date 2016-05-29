function [yield] = Von_Mises(sigma_z, tau, x, dz, z)


%% Initialize
yield_stress = 324*10^6; % Pa
yield.alt = zeros(length(x), round(length(z)/dz), 5);
yield.sea = zeros(length(x), round(length(z)/dz), 5);


%% Check for Yielding at all points, along span, for each condition
for k = 1:5
    for j = 1:dz:length(z)
        for i = 1:length(x)
            if yield_stress <= 1.5 * ((( 2 * (sigma_z.alt(i,j,k))^2 + 6 * (tau.alt(i,j,k)^2))/2)^0.5)
                yield.alt(i,j,k) = 1;
            end
            if yield_stress <= 1.5 * ((( 2 * (sigma_z.sea(i,j,k))^2 + 6 * (tau.sea(i,j,k)^2))/2)^0.5)
                yield.sea(i,j,k) = 1;
            end
        end       
    end
end


%% Plots

if ismember(1, yield.sea);
    disp('Structure Yields')
    plot(1:length(x), yield.sea(:,1,1));
elseif ismember(1, yield.alt);
    disp('Structure Yields')
    plot(1:length(x), yield.alt(:,1,1));
else
    disp('Structure Does Not Yield')
end