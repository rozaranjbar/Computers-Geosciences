% Test route generator

% Cleaning

clear all
close all
clc

% Start

segments = [1 3 5 7 9];

N = 5;

current_segment = segments(3);

% Distances (totally invented)
          %     1    3    5    7    9
distances = [  0.0  2.3  3.2  7.0 15.0;  % 1
               2.3  0.0  0.8  4.2 12.8;  % 2
               3.2  0.8  0.0  2.6 10.1;  % 3
               7.0  4.2  2.6  0.0  8.6;  % 4
              15.0 12.8 10.1  8.6  0.0]; % 5
          
distance_limit = 4;

% Battery

energy_consumption = 4; % 4 units of SOC (battery) for each unit of distance

initial_battery = 28;

% Function route generator

ST = clock;
[Routes,distance_step,distance_traveled,wasted_energy,energy_step] = limited_route_generator2(current_segment,segments,N,distances,initial_battery,energy_consumption,distance_limit,battery_rechage);

FT = clock;

TimeRoutes = etime(FT,ST);