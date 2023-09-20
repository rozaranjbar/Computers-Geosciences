function [Routes,distance_step,distance_traveled,wasted_energy,energy_step] = limited_route_generator2(current_segment, segments,N,distances,initial_battery,energy_consumption,distance_limit,battery_rechage)
% N = Np - Nred

    Routes = [];
    Number_of_possible_routes = length(segments)^N;
    
    if N == 1

        Routes = segments;
        
    else 
        
        Initial_Routes = [current_segment*ones(length(segments),1) segments'];
        
        Possible_Routes{1} = [];
        Possible_distance_step{1} = [];
        Possible_distance_traveled{1} = [];
        Possible_wasted_energy{1} = [];
        Possible_energy_step{1} = [];
        
        for i = 1:length(Initial_Routes)
        
            [route_OK,distance_step,distance_traveled,wasted_energy,energy_step] = check_route(Initial_Routes(i,:),distances,initial_battery,energy_consumption,distance_limit,segments,battery_rechage);
            
            if route_OK
                
                Possible_Routes{1} = [Possible_Routes{1}; Initial_Routes(i,:)];
                Possible_distance_step{1} = [Possible_distance_step{1}; distance_step];
                Possible_distance_traveled{1} = [Possible_distance_traveled{1}; distance_traveled];
                Possible_wasted_energy{1} = [Possible_wasted_energy{1}; wasted_energy];
                Possible_energy_step{1} = [Possible_energy_step{1}; energy_step];
                
            end
            
        end
        
        for k = 2:N
            
            Possible_Routes{k} = [];
            Possible_distance_step{k} = [];
            Possible_distance_traveled{k} = [];
            Possible_wasted_energy{k} = [];
            Possible_energy_step{k} = [];
            
            for i = 1:size(Possible_Routes{k-1},1)
                
                for j = 1:length(segments)
                    
                    New_route = [Possible_Routes{k-1}(i,:) segments(j)];
                    
                    [route_OK,distance_step,distance_traveled,wasted_energy,energy_step] = check_route(New_route,distances,initial_battery,energy_consumption,distance_limit,segments,battery_rechage);
                    
                    if route_OK
                
                        Possible_Routes{k} = [Possible_Routes{k}; New_route];
                        Possible_distance_step{k} = [Possible_distance_step{k}; distance_step];
                        Possible_distance_traveled{k} = [Possible_distance_traveled{k}; distance_traveled];
                        Possible_wasted_energy{k} = [Possible_wasted_energy{k}; wasted_energy];
                        Possible_energy_step{k} = [Possible_energy_step{k}; energy_step];
                    
                    end
                
                end

            end
        
        end
        
    end
    
    Routes = Possible_Routes{k}(:,2:N+1);
    distance_step = Possible_distance_step{k}(:,1:N);
    distance_traveled = Possible_distance_traveled{k};
    wasted_energy = Possible_wasted_energy{k};
    energy_step = Possible_energy_step{k}(:,2:N+1);
    
end