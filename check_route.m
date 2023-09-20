function [route_OK,distance_step,distance_traveled,wasted_energy,energy_step] = check_route(New_route,distances,initial_battery,energy_consumption,distance_limit,segments,battery_rechage)

    route_OK = 1;
    
    N = length(New_route);
    
    distance_traveled = 0;
    distance_step = zeros(1,N);
    
    energy_step = zeros(1,N);
    energy_step(1) = initial_battery;
    
    for i = 1:length(New_route)-1
        
        from = find(segments==New_route(i));
        to = find(segments==New_route(i+1));
        
        distance_step(i) = distances(from,to);
        
        distance_traveled = distance_traveled + distance_step(i);
        
        energy_step(i+1) = energy_step(i) + battery_rechage - distance_step(i)*energy_consumption;
        
        if energy_step(i+1) > 100
            
            energy_step(i+1) = 100;
            
        end

        if distances(from,to) > distance_limit || energy_step(i+1) < 0

            route_OK = 0;

        end
    
    end
    
    wasted_energy = distance_traveled*energy_consumption;
    
    if initial_battery - distance_traveled*energy_consumption < 0
        
        route_OK = 0;
        
    end

end