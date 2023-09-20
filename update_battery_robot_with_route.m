function new_battery = update_battery_robot_with_route(segments,previous_segment,current_segment,distances,current_battery,battery_rechage,energy_consumption)

    from = find(segments==previous_segment);
    to = find(segments==current_segment);

    distance_traveled = distances(from,to);
    
    new_battery = current_battery + battery_rechage - distance_traveled*energy_consumption;
       
end