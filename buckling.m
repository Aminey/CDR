function [rib_spacing] = buckling(sigma_z,str,structure,skin,z)

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
    skin_t = skin.t;
    
    z_current = 0; %m 
    z_rib = [0];
    z_index = [];
    
    while z_rib < 5.41
        %find z index corresponding to spanwise location
        [val index] = min(abs(z_current-z));
        z_index = index;
        %find max stress
        for i = 1:length(str.i_CCW)
            for j = 1:5        
                %rows correspond to max bending stress in a given stringer, columns 
                %correspond to critical flight conditions as following convention above
                stress_sea = sigma_z.sea(i,z_index,j);
                stress_alt = sigma_z.alt(i,z_index,j);            
                stringer_stress(i,j) = max([stress_sea stress_alt]);
            end
        end     
        max_stress = max(max(stringer_stress));

        l = sqrt(pi^2*E*stringer_I/(1.5*max_stress*stringer_area*0.5^2)); %m    
        z_current = z_current+l;
        z_rib = [z_rib z_current];
    end
    
    rib_spacing = z_rib;
    disp(['Ribs need to be placed at ' num2str(z_rib) ' to prevent buckling']);
    disp('Stringer buckling complete');
    
%     %%stringer buckling
%     %find max stress
%     for i = 1:length(str.i_CCW)
%         for j = 1:5        
%             %rows correspond to max bending stress in a given stringer, columns 
%             %correspond to critical flight conditions as following convention above
%             stress_sea = sigma_z.sea(i,1,j);
%             stress_alt = sigma_z.alt(i,1,j);            
%             stringer_stress(i,j) = min([stress_sea stress_alt]);
%         end
%     end     
%     max_stress = min(min(stringer_stress));
%     
%     l = sqrt(pi^2*E*stringer_I/(1.5*-1*max_stress*stringer_area*0.5^2)); %m    
%     rib_spacing = l; %m
%     
%     disp([num2str(5.41/rib_spacing) ' ribs need to be placed ' num2str(rib_spacing) ' meters apart to prevent buckling']);
    
    %%skin buckling (upper side only)
    stringer_distances = [];
    for i = 1:length(str.i_upp)-1
        d = sqrt((str.x_upp(i+1)-str.x_upp(i))^2-(str.y_upp(i+1)-str.y_upp(i))^2);
        stringer_distances = [stringer_distances d];
    end 
    
    skin_stress_sea = [];
    skin_stress_alt = [];
    
    for j = 1:5        
        %rows correspond to max bending stress in a given stringer, columns 
        %correspond to critical flight conditions as following convention above
        skin_stress_sea = [skin_stress_sea min(sigma_z.sea(:,1,j))];
        skin_stress_alt = [skin_stress_alt min(sigma_z.alt(:,1,j))];            
    end
    
    max_skin_stress = min([skin_stress_sea skin_stress_alt]);    
    k = 9;
    v = 0.33;
    
    b = sqrt((skin_t^2*k*pi^2*E)/(12*(1-v^2)*1.5*-1*max_skin_stress));
    
    if max(stringer_distances) <= b
        disp(['Skin buckling satisfied. Stringer spacing needs to be no more than ' ...
            'more than ' num2str(b) ' m apart. The largest stringer spacing'...
            ' is ' num2str(max(stringer_distances)) ' m.']);
    else
        disp(['Skin buckling not satisfied. Stringer spacing needs to be no ' ...
            'more than ' num2str(b) ' m apart. The largest stringer spacing'...
            ' is ' num2str(max(stringer_distances)) ' m.']);
    end
    
    
    disp('Buckling complete');

end