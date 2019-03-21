%% Initial Set up of parameters

clear variables;
clc;

A = 8.5;
d = 150;
area = pi*((d/2)^2);
cp = 0.496;
efficiency = 0.9;
filename = 'Harshil_Sumaria.mat';
k = 1.9;
steps = [0:0.05:30]; %can change 0.1 to control the step size and smoothness of plots
n = length (steps);
no_bins = 23;
pout_6 = 6000000;
pmax_6 = pout_6/efficiency;
pout_12 = 12000000;
pmax_12 = pout_12 / efficiency;
rho = 1.225;


 %% PART 1 Question 1
figure
% do your job here
drawnow;
set(get(handle(gcf),'JavaFrame'),'Maximized',1);
subplot(2,3,1);
%Characteristic curve and rated speed for 6MW WTG
[power_6, ratedspeed_6] = characteristic_curve ( steps, rho, cp, pmax_6, efficiency, area); 
ratedspeed_6
plot (steps, power_6);
axis([0 steps(n) 0 pmax_6*5]);
hold on;

%Total power available
[power_no_limit] = power_betz (steps, rho,1,area); 

plot (steps, power_no_limit);

%Total power available with Betz Limit
[power_betz_limit] = power_betz(steps, rho, cp, area); 
plot (steps, power_betz_limit);
title('Characteristic Curves (6MW WTG)');
xlabel('Wind speed (m/s)');
ylabel('Power (MW)');
legend('6MW WTG','Available Power','With Betz Limit')
hold off;

%% Part 1 Question 2
sample_size = 365*24;        %Number of hours in a year

%vw = wblrnd(A, k, [sample_size, 1]);     %random weibull data for wind speeds
%save(filename,'vw')                        %save windspeeds to file for repeatability
load(filename,'vw')                     %load windspeeds from file for repeatability

%figure
subplot(2,3,2);
hist(vw, no_bins);     %histogram of wind speeds
histval = hist(vw, no_bins);       %save quantity in each bin for later
title ('Histogram of Wind Speeds over One Year');
xlabel('Wind speed (m/s)');
ylabel('Number of hours');
avgwindspeed = mean(vw)
medwindspeed = median(vw)%save average speed


%% Part 1 Question 3

%get the middle of each histogram bin
mids = steps(n)/no_bins;
for p = 1:no_bins
    wind_speeds(p) = (mids*(p-1) + mids*(p))/2;
end

%Available energy for ideal, betz limit and real WTG
[sum_energy_ideal, energy_ideal] = e_for_extraction (rho, 1,  area, no_bins, histval, wind_speeds);
[~, energy_betz] = e_for_extraction (rho, cp,  area, no_bins, histval, wind_speeds);
[sum_energy_real, energy_real] = real_e_for_extraction (rho, cp, pmax_6, efficiency, area,no_bins,  histval, wind_speeds, ratedspeed_6);

%Plot available energy on same graph for three cases
subplot(2,3,3);
bar(wind_speeds, energy_ideal,'red');
%figure
hold on;
bar(wind_speeds, energy_betz,'green');
%figure
hold on;
bar(wind_speeds, energy_real,'blue');
title('Energy Available for Extraction (6MW WTG)');
xlabel('Wind speed (m/s)');
ylabel('Energy (Wh)');
legend('Total Available Energy','With Betz Limit','Energy Extracted')
legend('Location','northwest')

%% Part 1 Question 4
%Percentage change between available energy and energy extracted by 6MW WTG
perc_energy_produced_actual = 100*(sum_energy_real/sum_energy_ideal)

%Energy produced if the WTG were to run at average speed all year
energy_prod_avg = 0.5*rho*area*(avgwindspeed^3)*365*24*cp*efficiency;

%Percentage change between available energy and energy extracted by 6MW WTG running at average speed all year
perc_energy_produced_average = 100*(energy_prod_avg/sum_energy_ideal)

%% Part 2 Question 1

subplot(2,3,4);
[power_12, ratedspeed_12] = characteristic_curve (steps, rho, cp, pmax_12, efficiency, area); %Generating Characteristic Curve with rated speed and power being returned
ratedspeed_12
plot (steps, power_12);
axis([0 steps(n) 0 (pmax_12)*5]);

hold on;
plot (steps, power_6); %plot the s6MW characteristic curve for reference
plot (steps, power_no_limit); %plot the no limit curve
plot (steps, power_betz_limit);%plot the betz limit curve
title('Characteristic Curves (12MW WTG)');
xlabel('Wind speed (m/s)');
ylabel('Power Output');
legend('12MW WTG', '6MW WTG','Available Power','With Betz Limit')
legend('Location','northwest')
hold off;

%% Part 2 Question 2

%bar graph as in part 1 question 3
[sum_energy_real_12, energy_real_12] = real_e_for_extraction ( rho, cp, pmax_12, efficiency, area, no_bins, histval, wind_speeds, ratedspeed_12);
subplot(2,3,5);
bar(wind_speeds, energy_real_12,'blue');
hold on;
bar(wind_speeds, energy_real,'green');
title('Energy Availability Comparison');
ylabel('Energy (Wh)');
xlabel('Wind speed (m/s)');
legend('12MW WTG', '6MW WTG')
hold off;

%% Part 2 Question 3
%Percentage change between available energy and energy extracted by 12MW WTG
perc_energy_produced_actual_12 = 100*(sum_energy_real_12/sum_energy_ideal)


%% Part 2 Question 4

%percentage difference between yield of 6MW and 12MW WTG
yield_change = 100*(sum_energy_real_12 - sum_energy_real)/sum_energy_real 

%% FUNCTION DEFINITIONS

%generate the characteristic curve
function [power, rated_speed] = characteristic_curve (steps, rho, Cp, power_max, efficiency, a)
n = length (steps);
rated_speed = 0;
flag = 0;
for ii = 1:n
    power(ii) = 0.5*rho*a*(steps(ii)^3)*Cp*efficiency;
            
    if steps(ii)>= 25 || steps(ii) < 4
        power(ii)=0;
    
    elseif power(ii) > efficiency*power_max
        if flag == 0
            rated_speed = steps(ii);
        end
        power(ii) = power_max*efficiency;
        flag = 1;
    end

end
end


%generate the bar charts of available energy
function [sum_energy, energy] = e_for_extraction ( rho, cp,  a, no_bins,  no_hours, wind_speeds)

for ii = 1:no_bins
    energy(ii) = 0.5*rho*a*(wind_speeds(ii)^3)*no_hours(ii)*cp;
end
sum_energy = sum(energy);

end


%generate the bar charts of available energy for the WTG
function [sum_energy, energy] = real_e_for_extraction ( rho, cp, power_max, efficiency, a, no_bins, no_hours, wind_speeds, rated_speed)

for ii = 1:no_bins
    energy(ii) = 0.5*rho*a*(wind_speeds(ii)^3)*no_hours(ii)*cp*efficiency;
    if wind_speeds(ii)>=25 || wind_speeds(ii) < 4       %cut-out and cut in condition
        energy(ii)=0;
    elseif wind_speeds(ii)>=rated_speed  %rated maximum speed
       energy(ii) = power_max*efficiency*no_hours(ii);
    end
end
sum_energy = sum(energy);

end


%plot ideal and betz limit curves
function [Pwbetz] = power_betz( steps, rho, cp,  a)

n = length (steps);
for ii = 1:n
    Pwbetz(ii) = 0.5*rho*a*(steps(ii)^3)*cp;

end
end