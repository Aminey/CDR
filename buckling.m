function [rib_spacing] = buckling(sigma_z,str,structure)

    %critical flight condition loading
    %[PHAA PLAA NHAA Downward_Gust NLAA]
    
    % i is index ccw
    % j is location along span
    % k is condition
    
    %load in variables
    stringer_index = str.i_CCW;
    stringer_x = [fliplr(str.x_upp) str.x_low]; %m
    stringer_y = [fliplr(str.y_upp) str.y_low]; %m
    stringer_Ixx = [fliplr(str.Ixx_upp) str.Ixx_low]; %m^4
    stringer_Iyy = [fliplr(str.Iyy_upp) str.Iyy_low]; %m^4
    stringer_Ixy = [fliplr(str.Ixy_upp) str.Ixy_low]; %m^4    
    E = structure.E; %Pa
    
    %find max stress
    for i = 1:length(str.I_CCW)
        for j = 1:5        
            %rows correspond to max bending stress in a given stringer, columns 
            %correspond to critical flight conditions as following convention above
            stringer_stress(i,j) = sigma_z(i,0,j);
        end
    end     
    P_crit = max(max(stringer_stress));
    
    
    
end