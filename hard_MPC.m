function [U,J] =  hard_MPC(H, F, AU, bU)

    % PEPE : opciones para que quadprog imprima menos
    ops = optimset('quadprog');
    ops.Display = 'off';
    
    U = quadprog(H,F,AU,bU,[],[],[],[],[],ops);
    
    if isempty(find(U,1))
    
        J = Inf;
        U = zeros(size(H,2),1);
    
    else
        
        J = 1/2*U'*H*U+F*U;
        
    end

end
   