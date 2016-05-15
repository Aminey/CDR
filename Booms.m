function [Booms] = Booms(x,y,z,sigma_z,skin, str, caps, spar, dz)

Booms.alt = zeros(length(x),length(z),5);
Booms.sea = zeros(length(x),length(z),5);
str.A = str.A_upp;

for k = 1:5
    for j = 1:dz:length(z)    
        for i = 2:length(x)-1           % first and last coordinates will be calculated separately
            dist_rear = ((x(i)+x(i-1))^2 + (y(i)+y(i-1))^2)^0.5;
            dist_fore = ((x(i)+x(i+1))^2 + (y(i)+y(i+1))^2)^0.5;
            if ismember(i, str.i_CCW)
                area_term = str.A;
            elseif ismember(i, spar.i_CCW)
                area_term = 2*caps.A;
            else
                area_term = 0;
            end
            Booms.alt(i,j,k) = (skin.t * dist_rear/6) * (2 + (sigma_z.alt(i-1,j,k)/sigma_z.alt(i,j,k))) + ...
                          (skin.t * dist_fore/6) * (2 + (sigma_z.alt(i+1,j,k)/sigma_z.alt(i,j,k))) + area_term;
            Booms.sea(i,j,k) = (skin.t * dist_rear/6) * (2 + (sigma_z.sea(i-1,j,k)/sigma_z.sea(i,j,k))) + ...
                          (skin.t * dist_fore/6) * (2 + (sigma_z.sea(i+1,j,k)/sigma_z.sea(i,j,k))) + area_term;
        end
        
%% First Coordinate
        Booms.alt(1,j,k) = spar.t * spar.h(2)/6 *...
                              (2 + sigma_z.alt(end,j,k)/sigma_z.alt(1,j,k))+...
                              skin.t * (((x(2)-x(1))^2 + (y(2)-y(1))^2)^0.5)/6*...
                              (2 + sigma_z.alt(2,j,k)/sigma_z.alt(1,j,k)) + caps.A;
        Booms.sea(1,j,k) = spar.t * spar.h(2)/6 *...
                              (2 + sigma_z.sea(end,j,k)/sigma_z.sea(1,j,k))+...
                              skin.t * (((x(2)-x(1))^2 + (y(2)-y(1))^2)^0.5)/6*...
                              (2 + sigma_z.sea(2,j,k)/sigma_z.sea(1,j,k)) + caps.A;
       

%% Last Coordinate
        Booms.alt(end,j,k) = skin.t * (((x(end)-x(end-1))^2 + (y(end)-y(end-1))^2)^0.5)/6 *...
                              (2 + sigma_z.alt(end-1,j,k)/sigma_z.alt(end,j,k))+...
                              spar.t * spar.h(2)/6*...
                              (2 + sigma_z.alt(1,j,k)/sigma_z.alt(end,j,k)) + caps.A;
        Booms.sea(end,j,k) = skin.t * (((x(end)-x(end-1))^2 + (y(end)-y(end-1))^2)^0.5)/6 *...
                              (2 + sigma_z.sea(end-1,j,k)/sigma_z.sea(end,j,k))+...
                              spar.t * spar.h(2)/6*...
                              (2 + sigma_z.sea(1,j,k)/sigma_z.sea(end,j,k)) + caps.A;
       
       

%% Stringers on first and last coordinates?
        if ismember(1,str.i_CCW)
           Booms.alt(1,j,k) = Booms.alt(1,j,k) + str.A; 
           Booms.sea(1,j,k) = Booms.sea(1,j,k) + str.A;
        end
        if ismember(length(x),str.i_CCW)
           Booms.alt(end,j,k) = Booms.alt(end,j,k) + str.A; 
           Booms.sea(end,j,k) = Booms.sea(end,j,k) + str.A;
        end
       


%% Add the contribution of spar's path for shear flow
%     m = 1;
%     spar_holder = zeros(1,4);  % will hold the indices of the spar connection locations
%     
%     for i = 1:length(x)
%         if ismember(i,spar.i_CCW) 
%            spar_holder(m) = i;
%            m = m+1;
%         end
%     end
%     
%     for m = 1:4   % two spars will have four connection points. going counterclockwise, connection 1 will be the top of the rear 
%                   % spar and 4 will be the bottom. 2 and 3 will be the top
%                   % and bottom of the forward spar. thus, index (5-m) is
%                   % used to relate spar connections 1 and 4 and connections 2 and 3
%          Booms.sea(spar_holder(m),j,k) = Booms.sea(spar_holder(m),j,k) + skin.t * spar.h(1)/6 * ...
%                                                         (2 + sigma_z.sea(spar_holder(5-m),j,k)/sigma_z.sea(spar_holder(m),j,k));
%          Booms.alt(spar_holder(m),j,k) = Booms.alt(spar_holder(m),j,k) + skin.t * spar.h(1)/6 * ...
%                                                         (2 + sigma_z.alt(spar_holder(5-m),j,k)/sigma_z.alt(spar_holder(5-m),j,k));
% 
%     end
    for i = 2:3
       Booms.sea(spar.i_CCW(i),j,k) = Booms.sea(spar.i_CCW(i),j,k) + spar.t * spar.h(1)/6 * ...
                                                ((2 + sigma_z.sea(spar.i_CCW(5-i),j,k))/sigma_z.sea(spar.i_CCW(i),j,k));
       Booms.alt(spar.i_CCW(i),j,k) = Booms.alt(spar.i_CCW(i),j,k) + spar.t * spar.h(1)/6 * ...
                                                ((2 + sigma_z.alt(spar.i_CCW(5-i),j,k))/sigma_z.alt(spar.i_CCW(i),j,k));

    end
    end
end

figure;
plot(x,Booms.alt(:,1,1));
xlabel('x');
ylabel('Booms.alt');

figure;
plot(x,sigma_z.alt(:,1,1).*Booms.alt(:,1,1));
xlabel('x');
ylabel('sigma_z.alt*Booms.alt');

disp('sum sigma*boom = ');
sum(sigma_z.alt(:,1,1).*Booms.alt(:,1,1))

disp('Booms complete');