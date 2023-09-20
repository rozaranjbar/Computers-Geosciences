function predicted_sigma = predict_sigma_V2(initial_sigma,Route_predicted,Np,A,D,sigma_disturbances)   
    
    predicted_sigma(:,1) = initial_sigma;
    
    for k = 1:Np-1
        
        for i = 1:length(initial_sigma)
                
            if Route_predicted(k) == i

                predicted_sigma(i,k+1) = 0;

            else

                predicted_sigma(i,k+1) = sum( (A(i,:).^2).*predicted_sigma(:,k)' ) + sum( (D(i,:).^2).*sigma_disturbances' );

            end 
                
        end
        
    end
    
end