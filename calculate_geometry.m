function [wing_section_centroid,structure,...
         component_inertias] = calculate_geometry(wing,c,structure,theta)
    %calculates and plots centroids of components and wing section, then
    %calculates component and section moments of area
    
    % Initialized for better readability in the loops-
    bl = structure.bl;
    bt = structure.bt;
    sh = structure.sh;
    st = structure.st;
    kt = structure.kt;
    
    component = 1;
    cx = 0;
    cy = 0;
    total_area = 0;
    X = (-c*cosd(theta)+sqrt(c^2*cosd(theta)^2+3*c^2))/2;
    beta = -theta;
    %skins' cross section centroids
    for i = 1:4
       xcentroid(component) = 0.25*(wing.skin(i).x(1)+wing.skin(i).x(2)+...
           wing.skin(i).x(3)+wing.skin(i).x(4));
       ycentroid(component) = 0.25*(wing.skin(i).y(1)+wing.skin(i).y(2)+...
           wing.skin(i).y(3)+wing.skin(i).y(4));
       area(component) = 0.5*abs(wing.skin(i).x(1)*wing.skin(i).y(3)-...
           wing.skin(i).y(1)*wing.skin(i).x(3)+wing.skin(i).x(3)*...
           wing.skin(i).y(4)-wing.skin(i).y(3)*wing.skin(i).x(4)+...
           wing.skin(i).x(4)*wing.skin(i).y(2)-wing.skin(i).y(4)*...
           wing.skin(i).x(2)+wing.skin(i).x(2)*wing.skin(i).y(1)-...
           wing.skin(i).y(2)*wing.skin(i).x(1));
       cx = cx + xcentroid(component)*area(component);
       cy = cy + ycentroid(component)*area(component);
       total_area = total_area + area(component);
       component = component+1;
    end
    
    %spars' cross section centroids
    for i = 1:3
       xcentroid(component) = 0.25*(wing.spar(i).x(1)+wing.spar(i).x(2)+...
           wing.spar(i).x(3)+wing.spar(i).x(4));
       ycentroid(component) = 0.25*(wing.spar(i).y(1)+wing.spar(i).y(2)+...
           wing.spar(i).y(3)+wing.spar(i).y(4));
       area(component) = 0.5*abs(wing.spar(i).x(1)*wing.spar(i).y(3)-...
           wing.spar(i).y(1)*wing.spar(i).x(3)+wing.spar(i).x(3)*...
           wing.spar(i).y(4)-wing.spar(i).y(3)*wing.spar(i).x(4)+...
           wing.spar(i).x(4)*wing.spar(i).y(2)-wing.spar(i).y(4)*...
           wing.spar(i).x(2)+wing.spar(i).x(2)*wing.spar(i).y(1)-...
           wing.spar(i).y(2)*wing.spar(i).x(1));
       cx = cx + xcentroid(component)*area(component);
       cy = cy + ycentroid(component)*area(component);
       total_area = total_area + area(component);
       component = component+1;
    end  
    
    %brackets' cross section centroids
    for i = 1:16
       xcentroid(component) = 0.25*(wing.bracket(i).x(1)+wing.bracket(i).x(2)+...
           wing.bracket(i).x(3)+wing.bracket(i).x(4));
       ycentroid(component) = 0.25*(wing.bracket(i).y(1)+...
           wing.bracket(i).y(2)+wing.bracket(i).y(3)+wing.bracket(i).y(4));
       area(component) = 0.5*abs(wing.bracket(i).x(1)*wing.bracket(i).y(3)-...
           wing.bracket(i).y(1)*wing.bracket(i).x(3)+wing.bracket(i).x(3)*...
           wing.bracket(i).y(4)-wing.bracket(i).y(3)*wing.bracket(i).x(4)+...
           wing.bracket(i).x(4)*wing.bracket(i).y(2)-wing.bracket(i).y(4)*...
           wing.bracket(i).x(2)+wing.bracket(i).x(2)*wing.bracket(i).y(1)-...
           wing.bracket(i).y(2)*wing.bracket(i).x(1));
       cx = cx + xcentroid(component)*area(component);
       cy = cy + ycentroid(component)*area(component);
       total_area = total_area + area(component);
       component = component+1;
    end
    component = component-1;
    
    %calculate overall centroid
    cx = cx/total_area;
    cy = cy/total_area;
    
    %plot centroids
    hold on    
    for i = 1:component
        h1 = plot(xcentroid(i),ycentroid(i),'r*');
    end
    h2 = plot(cx,cy,'b*');
    legend([h1 h2],{'Component Centroids','Wing Section Centroid'});    
    hold off
    
    %calculate area moments of inertia
    component = 1;

    %skins
    for i = 1:4
       if (i == 1) || (i == 3)
           Ixx(component) = 1/12*((c/2)*kt^3)+...
               area(component)*(ycentroid(component)-cy)^2;
           Iyy(component) = 1/12*((c/2)^3*kt)+...
               area(component)*(xcentroid(component)-cx)^2;
           Ixy(component) = (c/2*kt)*(xcentroid(component)-cx)*(ycentroid(component)-cy);
       elseif (i == 2) || (i == 4)
           Ixx(component) = 1/12*(X^3*kt)*sind(beta)^2+...
               area(component)*(ycentroid(component)-cy)^2;
           Iyy(component) = 1/12*(X^3*kt)*cosd(beta)^2+...
               area(component)*(xcentroid(component)-cx)^2;
           Ixy(component) = 1/24*(X^3*kt)*sind(2*beta)+area(component)*...
               ((xcentroid(component)-cx)^2+(ycentroid(component)-cy)^2);
       end
       component = component+1;
    end
    
    %spars
    for i = 1:3
       if (i == 1) || (i == 2)
           Ixx(component) = 1/12*(sh^3*st)+...
               area(component)*(ycentroid(component)-cy)^2;
           Iyy(component) = 1/12*(sh*st^3)+...
               area(component)*(xcentroid(component)-cx)^2;
           Ixy(component) = area(component)*(xcentroid(component)-cx)*(ycentroid(component)-cy);
       elseif (i == 3)
           Ixx(component) = 1/12*(st^3*sh)*sind(beta)^2+...
               area(component)*(ycentroid(component)-cy)^2;
           Iyy(component) = 1/12*(st^3*sh)*cosd(beta)^2+...
               area(component)*(xcentroid(component)-cx)^2;
           Ixy(component) = 1/24*(st^3*sh)*sind(2*beta)+area(component)*...
               ((xcentroid(component)-cx)^2+(ycentroid(component)-cy)^2);
       end
       component = component+1;
    end
    
    %brackets
    for i = 1:16
       if (i == 1) || (i == 4) || (i == 9) || (i == 12)
           Ixx(component) = 1/12*(bl^3*bt)+...
               area(component)*(ycentroid(component)-cy)^2;
           Iyy(component) = 1/12*(bl*bt^3)+...
               area(component)*(xcentroid(component)-cx)^2;
           Ixy(component) = area(component)*(xcentroid(component)-cx)*(ycentroid(component)-cy);
       elseif (i == 2) || (i == 3) || (i == 10) || (i == 11)
           Ixx(component) = 1/12*(bl*bt^3)+...
               area(component)*(ycentroid(component)-cy)^2;
           Iyy(component) = 1/12*(bl^3*bt)+...
               area(component)*(xcentroid(component)-cx)^2;
           Ixy(component) = area(component)*(xcentroid(component)-cx)*(ycentroid(component)-cy);
       elseif (i == 6) || (i == 7) || (i == 14) || (i == 15)
           Ixx(component) = 1/12*(bl^3*bt)*sind(beta)^2+...
               area(component)*(ycentroid(component)-cy)^2;
           Iyy(component) = 1/12*(bl^3*bt)*cosd(beta)^2+...
               area(component)*(xcentroid(component)-cx)^2;
           Ixy(component) = 1/24*(bl^3*bt)*sind(2*beta)+area(component)*...
               ((xcentroid(component)-cx)^2+(ycentroid(component)-cy)^2);
       elseif (i == 5) || (i == 8) || (i == 13) || (i == 16)
           Ixx(component) = 1/12*(bl*bt^3)*sind(beta)^2+...
               area(component)*(ycentroid(component)-cy)^2;
           Iyy(component) = 1/12*(bl*bt^3)*cosd(beta)^2+...
               area(component)*(xcentroid(component)-cx)^2;
           Ixy(component) = 1/24*(bl*bt^3)*sind(2*beta)+area(component)*...
               ((xcentroid(component)-cx)^2+(ycentroid(component)-cy)^2);
       end
       component = component+1;
    end
    component = component-1;

    %calculate overall area moments about overall centroid
    Ixx_wing = sum(Ixx);
    Iyy_wing = sum(Iyy);
    Ixy_wing = sum(Ixy);
    
    %function output
    wing_section_centroid = [cx,cy]; % centroid coordinate units m
    structure.component_centroids = [xcentroid', ycentroid'];
    structure.inertias = [Ixx_wing',Iyy_wing',Ixy_wing']; % moment units m^4
    component_inertias = [Ixx',Iyy',Ixy'];
    
end