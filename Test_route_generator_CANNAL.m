% Test route generator (REAL CANNAL)

% Cleaning

clear all
close all
clc

% Start

load('amplied_cannal.mat','A','segments')

water_reaches = [];
 
for i = 1:size(A,1)
     
    if A(i,i) == 1

        water_reaches = [water_reaches, i];

    end

end

segments = water_reaches;

current_segment = segments(3); % Initial position of the robot

N = 5; % Length of the route (time, allocation horizon)

% Distances (totally invented)

distances_reaches = [100 1200 400 800 2000 1700 1600 1700];

distances = zeros(length(distances_reaches));

for i = 1:size(distances,1)
    
    for j = 1:size(distances,2)
        
        if distances(i,j) == 0
        
            distances(i,j) = 0;

                for k = i:j-1

                    distances(i,j) = distances(i,j) + distances_reaches(k);

                end

            distances(j,i) = distances(i,j);
        
        end
        
    end
    
end

robot_velocity = 1; % in m/s

sampletime = 0.5; % hours
          
distance_limit = robot_velocity*(3600)*sampletime; % distance limit in meters

% Battery

energy_consumption = 0.004; % units of SOC (battery) for each unit of distance

initial_battery = 28; % Initial battery of the robot

battery_rechage = 0.3; % units of SOC (battery) for each sample time without moving 

% Function route generator

ST = clock;

[Routes,distance_step,distance_traveled,wasted_energy,energy_step] = limited_route_generator2(current_segment,segments,N,distances,initial_battery,energy_consumption,distance_limit,battery_rechage);

FT = clock;

TimeRoutes = etime(FT,ST);

% Compact data. Remember the first distance_traveled is the traveling from
% the current segment to the first water reach of the route.

CompactData = [Routes,distance_traveled,distance_step,wasted_energy,energy_step];
CompactDataDistance = [Routes,distance_traveled,distance_step];
CompactDataEnergy = [Routes,wasted_energy,energy_step];