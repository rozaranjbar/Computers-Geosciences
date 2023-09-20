function [completed_Routes] = complete_Routes(Routes,Np,segments)

    %% Complete routes
    
    length_route = size(Routes,2);
    
    N_elems = size(segments,2);
    
    completed_Routes = zeros(size(Routes, 1),Np);
    
    for i = 1:size(Routes, 1)
        
        aux = Routes(i,:);
        
        last_segment = Routes(i,length_route);
        
        start_pos_segments = find(segments==last_segment);
        
        cont = 1;
        
        while size(aux,2) < Np
            
            if mod(cont,2) == 1 && start_pos_segments ~= length_route
            
                aux = [aux segments(start_pos_segments+1:N_elems)];
                
            elseif mod(cont,2) == 0
                
                aux = [aux segments(N_elems-1:-1:1)];
                
            end
            
            cont = cont + 1;            
            
        end        
        
        completed_Routes(i, :) = aux(1:Np);
        
    end

end