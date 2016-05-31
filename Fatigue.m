function [cycles] = Fatigue(x, sigma_z)

%% INPUTS
K_IC = 26*10^6;             % Pa(m)^0.5
f = 1;                      % given by teacher
c = 1.6*10^-11;             % given by teacher
m = 3.59;                   % given by teacher

%% Critical Crack Length

% loop will find critical crack length "a" at each point along the airfoil
for i = 1:length(x)
    a.sea(i) = (K_IC/(f*sigma_z.sea(:,1,1)))^2/pi;
    a.alt(i) = (K_IC/(f*sigma_z.alt(:,1,1)))^2/pi;
end

%% Cycles to Failure

cycles = (1/(c*((sigma_z.sea(:,1,1)-sigma_z.sea(:,1,4))*pi^0.5)^m)) * ()



