%{
Author: 
Assignment:
Creation Date: 
Inputs:
Outputs: 
Purpose: 
%}

clear;
clc;
close all;

%% Load data

% this is going to load with the names 
% here im using variblenamerule and setting it to preserve 
%this makes it read colum headers i think 
data_ISS_A = readtable('Data_ISS_A.csv', 'VariableNamingRule','preserve');
data_ISS_B = readtable('Data_ISS_B.csv','VariableNamingRule','preserve');

%i think this giving varible names 
data_ISS_A.Properties.VariableNames = {'Time_s','X_km', 'Y_km'};
data_ISS_B.Properties.VariableNames = {'Time_s','X_km', 'Y_km'};

%adjust?? 
%im not sure if this is right 

%this is taking the column into sperate varibles 
time_a = data_ISS_A.Time_s;
x_a = data_ISS_A.X_km;
y_a = data_ISS_A.Y_km;

time_b = data_ISS_B.Time_s;
x_b = data_ISS_B.X_km;
y_b = data_ISS_B.Y_km;

%% linear fit 
%least square for x and y 
% makes y = mx +b 

%linear model for x and y over time so we can get the velocity for both 
%lin fit for A 
c_x_a = polyfit(time_a, x_a, 1);
c_y_a = polyfit(time_a, y_a, 1);
u_a = c_x_a(1); 
x_a0 = c_x_a(2);
v_a = c_y_a(1);
y_a0 = c_y_a(2);

%lin fit for B 
c_x_b = polyfit(time_b, x_b, 1);
c_y_b = polyfit(time_b, y_b, 1);
u_b = c_x_b(1); 
x_b0 = c_x_b(2);
v_b = c_y_b(1);
y_b0 = c_y_b(2);

%% time of the closest approach 
% this is using the formula that was given to cal time when they are close 

% this will be getting the numerator and demoninator first 
n = (x_b0-x_a0) * (u_b-u_a) + (y_b0-y_a0) * (v_b- v_a);
d = (u_b-u_a)^2 + (v_b- v_a)^2; 
T_ca = -n / d; 

%% min distance at T of the closet approach 

%T_ca goes back to the postion equation 
%positon at T_Ca 
x_a_Tca = x_a0 + u_a * T_ca;
y_a_Tca = y_a0 + v_a * T_ca;
x_b_Tca = x_b0 + u_b * T_ca;
y_b_Tca = y_b0 + v_b * T_ca;

% this is getting min distance 
D_min = sqrt((x_b_Tca - x_a_Tca)^2 + (y_b_Tca-y_a_Tca)^2);

%% Error prop 


%% descion making and warning 

% determining if we should use prevetive manuvers 
% use D_min thresholds 

%im going to use a if statments for warning 
if D_min < 1.8 
    warning_code = 'Red - Action must be taken';
elseif D_min < 28.2 
    warning_code = 'Yello - Plans are devloped';
else 
    warning_code = 'Green - All clear';
end


%displaying 
fprintf('T_ca: %.2f second\n', T_ca);
fprintf('D_min: %.2f Km\n', D_min);
fprintf( 'warning code: %s\n', warning_code);

%% Plotting 

% postiton and trajectories to show closets approach 

figure; 
hold on; 
plot(x_a, y_a);
plot(x_b, y_b);
plot([x_a0,x_a_Tca], [y_a0, y_a_Tca]);
plot([x_b0,x_b_Tca], [y_b0, y_b_Tca]);
plot(x_a_Tca, y_a_Tca);
grid on; 
hold off; 
