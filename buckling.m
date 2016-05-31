function [rib_spacing] = buckling(sigma_z,str,structure)

    %critical flight condition loading
    %[PHAA PLAA NHAA Downward_Gust NLAA]
    
    % i is index ccw
    % j is location along span
    % k is condition
    
    %load in variables
    stringer_area = str.A_upp; %m^2
    stringer_l = sqrt(str.A_upp); %m approximate stringers as boxes, so length is sqrt(area)
    stringer_I = (stringer_l^4)/12; %m^4
    E = structure.E; %Pa
    
    %find max stress
    for i = 1:length(str.i_CCW)
        for j = 1:5        
            %rows correspond to max bending stress in a given stringer, columns 
            %correspond to critical flight conditions as following convention above
            stress_sea = sigma_z.sea(i,1,j);
            stress_alt = sigma_z.alt(i,1,j);            
            stringer_stress(i,j) = max([stress_sea stress_alt]);
        end
    end     
    max_stress = max(max(stringer_stress));
    
    l = sqrt(pi^2*E*stringer_I/(1.5*max_stress*stringer_area*0.5^2)); %m    
    rib_spacing = l; %m
    
    disp([num2str(5.41/rib_spacing) ' ribs need to be placed ' num2str(rib_spacing) ' meters apart to prevent buckling']);
    disp('Buckling complete');
    disp('You receive AIDS!');
    
end