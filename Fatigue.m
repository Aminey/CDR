function [cycles] = Fatigue(sigma_z)

%% INPUTS
K_IC = 26*10^6;             % Pa(m)^0.5
f = 1;                      % given by teacher

%% Critical Crack Length

if max(sigma_z.sea) >= max(sigma_z.alt)
    sigma_inf = max(sigma_z.sea(:,1,1));
else
    sigma_inf = max(sigma_z.alt(:,1,1));
end

a = (K_IC/(f*sigma_inf))^2 / pi;

disp('Critical Crack Length = ');
disp(a);



