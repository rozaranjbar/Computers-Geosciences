function [U,J] =  soft_MPC(H, F, AU, bU, soft_weight)

    ns = length(bU);
    Hs = soft_weight*eye(ns);
    H = blkdiag(H, Hs);
    F(end+ns) = 0;
    AU = [AU eye(ns)];

    % PEPE : opciones para que quadprog imprima menos
    ops = optimset('quadprog');
    ops.Display = 'off';
    
    U = quadprog(H,F,AU,bU,[],[],[],[],[],ops);
    
    if isempty(find(U,1))
    
        J = Inf;
        U = zeros(size(H,1),1);
    
    else
        
        J = U'*H*U+F*U;
        
    end

end
   