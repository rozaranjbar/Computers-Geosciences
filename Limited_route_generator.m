
function [Routes] = Limited_route_generator(segments,N,distances,energy_consumption)
% N = Np - Nred

    Number_of_possible_routes = length(segments)^N;
    Routes = zeros(Number_of_possible_routes,N);

    if N == 1

        Routes = segments;

    elseif N == 2

        for ik = 1:length(segments):length(segments)^N

            inds = ik:ik+length(segments)-1;
            Routes(inds,1) = segments(ceil(ik/length(segments)));
            Routes(inds,2) = segments(:);

        end
        
    else
        
        row = 1;
        
        for ik = 1:length(segments)
            
            tmp = Limited_route_generator(segments,N-1,distances,energy_consumption);
            Routes(row:row+size(tmp,1)-1,end-N+1:end) = [repmat(segments(ik),size(tmp,1),1) tmp];
            row = row  + size(tmp,1);
            
        end
        
        Routes(row+1:Number_of_possible_routes,:) = [];
        
    end
    
end