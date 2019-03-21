k = 1.9;
A = 8.5;
rho = 1.225;
max_power_coeff = 0.496;
blade_diameter = 150;
cutin_speed = 4;
cutout_speed = 25;
area = 0.5 * pi * 75^2;
%% PART 1
%A direct drive variable-speed generator with maximum power capacity of
%6MW, used with the turbine. Efficiency is assumed to always be 0.9.

%% QUESTION 1
%Obtain characteristic curve and determine the rated speed of the WTG.
p_a = 0.5 * rho * area * w_s^3;



%% QUESTION 2


%X = A[-ln(1-u)^(1/k)]
%8760 hours in a year
%make dist of 1000 values between 0 and 1, pick randomly for each hour and
%save to variable called vw using above formula for X
%plot a 23 bin histogram and determine the average wind speed
%save vw to a .mat file