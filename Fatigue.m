function [cycles, a] = Fatigue(x, sigma_z)

%% INPUTS
K_IC = 26*10^6;             % Pa(m)^0.5
f = 1;                      % given by teacher
c = 1.6*10^-11;             % given by teacher
m = 3.59;                   % given by teacher
frac = 0.1;                 % initial crack length fraction, given by teacher

%% Initialize
a.sea = zeros(1,length(x));
a.alt = zeros(1,length(x));
a.sea_frac = zeros(1,length(x));
a.alt_frac = zeros(1,length(x));
cycles.sea = zeros(1,length(x));
cycles.alt = zeros(1,length(x));

%% Critical Crack Length

% loop will find critical crack length "a" at each point along the airfoil
for i = 1:length(x)
    a.sea(i) = (K_IC/(f*sigma_z.sea(i,1,1)))^2/pi;
    a.alt(i) = (K_IC/(f*sigma_z.alt(i,1,1)))^2/pi;
end

a.sea_frac(:) = a.sea(:) * frac;
a.alt_frac(:) = a.alt(:) * frac;

%% Cycles to Failure

for i = 1:length(x)
    cycles.sea(i) = (1/(c*(abs(sigma_z.sea(i,1,1)-sigma_z.sea(i,1,4))*pi^0.5)^m)) * (((a.sea(i)^(1-m/2))/(1-m/2))-((a.sea_frac(i)^(1-m/2))/(1-m/2)));
    cycles.alt(i) = (1/(c*(abs(sigma_z.alt(i,1,1)-sigma_z.alt(i,1,4))*pi^0.5)^m)) * (((a.alt(i)^(1-m/2))/(1-m/2))-((a.alt_frac(i)^(1-m/2))/(1-m/2)));
end

if max(cycles.sea) >= max(cycles.alt)
    cycles.fail = max(cycles.sea);
else
    cycles.fail = max(cycles.alt);
end

disp('Cycles to Failure =')
disp(cycles.fail)