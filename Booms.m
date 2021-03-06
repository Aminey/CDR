function [Booms] = Booms(x,y,z,sigma_z,skin, str, caps, spar, dz)

Booms.alt = zeros(length(x),length(z),5);
Booms.sea = zeros(length(x),length(z),5);
str.A = str.A_upp;
dist_rear = zeros(1, length(x)-2);
dist_fore = zeros(1, length(x)-2);

for k = 1:5
    for j = 1:dz:length(z)    
        for i = 2:length(x)-1           % first and last coordinates will be calculated separately
            dist_rear(i) = ((x(i)-x(i-1))^2 + (y(i)-y(i-1))^2)^0.5;
            dist_fore(i) = ((x(i)-x(i+1))^2 + (y(i)-y(i+1))^2)^0.5;
            if ismember(i, str.i_CCW)
                area_term = str.A;
            elseif ismember(i, spar.i_CCW)
                area_term = caps.A;
            else
                area_term = 0;
            end
            Booms.alt(i,j,k) = (skin.t * dist_rear(i)/6) * (2 + (sigma_z.alt(i-1,j,k)/sigma_z.alt(i,j,k))) + ...
                          (skin.t * dist_fore(i)/6) * (2 + (sigma_z.alt(i+1,j,k)/sigma_z.alt(i,j,k))) + area_term;
            Booms.sea(i,j,k) = (skin.t * dist_rear(i)/6) * (2 + (sigma_z.sea(i-1,j,k)/sigma_z.sea(i,j,k))) + ...
                          (skin.t * dist_fore(i)/6) * (2 + (sigma_z.sea(i+1,j,k)/sigma_z.sea(i,j,k))) + area_term;
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
%         Booms.alt(1,j,k) = skin.t * (((x(2)-x(1))^2 + (y(2)-y(1))^2)^0.5)/6*...
%                               (2 + sigma_z.alt(2,j,k)/sigma_z.alt(1,j,k)) + caps.A;
%         Booms.sea(1,j,k) = skin.t * (((x(2)-x(1))^2 + (y(2)-y(1))^2)^0.5)/6*...
%                               (2 + sigma_z.sea(2,j,k)/sigma_z.sea(1,j,k)) + caps.A;

%% Last Coordinate
        Booms.alt(end,j,k) = skin.t * (((x(end)-x(end-1))^2 + (y(end)-y(end-1))^2)^0.5)/6 *...
                              (2 + sigma_z.alt(end-1,j,k)/sigma_z.alt(end,j,k))+...
                              spar.t * spar.h(2)/6*...
                              (2 + sigma_z.alt(1,j,k)/sigma_z.alt(end,j,k)) + caps.A;
        Booms.sea(end,j,k) = skin.t * (((x(end)-x(end-1))^2 + (y(end)-y(end-1))^2)^0.5)/6 *...
                              (2 + sigma_z.sea(end-1,j,k)/sigma_z.sea(end,j,k))+...
                              spar.t * spar.h(2)/6*...
                              (2 + sigma_z.sea(1,j,k)/sigma_z.sea(end,j,k)) + caps.A;
%         Booms.alt(end,j,k) = skin.t * (((x(end)-x(end-1))^2 + (y(end)-y(end-1))^2)^0.5)/6 *...
%                               (2 + sigma_z.alt(end-1,j,k)/sigma_z.alt(end,j,k))+...
%                               caps.A;
%         Booms.sea(end,j,k) = skin.t * (((x(end)-x(end-1))^2 + (y(end)-y(end-1))^2)^0.5)/6 *...
%                               (2 + sigma_z.sea(end-1,j,k)/sigma_z.sea(end,j,k))+...
%                               caps.A;
   
%% Add the contribution of spar's path for shear flow

    for i = 2:3
       Booms.sea(spar.i_CCW(i),j,k) = Booms.sea(spar.i_CCW(i),j,k) + spar.t * spar.h(1)/6 * ...
                                                (2 + (sigma_z.sea(spar.i_CCW(5-i),j,k)/sigma_z.sea(spar.i_CCW(i),j,k)));
       Booms.alt(spar.i_CCW(i),j,k) = Booms.alt(spar.i_CCW(i),j,k) + spar.t * spar.h(1)/6 * ...
                                                (2 + (sigma_z.alt(spar.i_CCW(5-i),j,k)/sigma_z.alt(spar.i_CCW(i),j,k)));

    end
    end
end

figure;
plot(1:length(x),Booms.alt(:,1,1));
xlabel('index');
ylabel('Boom areas (m^2)');

figure;
plot(x,sigma_z.alt(:,1,1).*Booms.alt(:,1,1));
xlabel('x');
ylabel('sigma_z.alt*Booms.alt');

disp('sum sigma*boom = ');
disp(sum(sigma_z.alt(:,1,1).*Booms.alt(:,1,1)));

disp('rms sigma*boom = ');
disp(rms(sigma_z.alt(:,1,1).*Booms.alt(:,1,1)));

disp('                              percent of final to ave = ');
disp( 100*(abs(sum(sigma_z.alt(:,1,1).*Booms.alt(:,1,1)))) / rms(sigma_z.alt(:,1,1).*Booms.alt(:,1,1)));

disp('Booms complete');