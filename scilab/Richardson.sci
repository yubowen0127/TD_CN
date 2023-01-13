
function[x_sol,r_k,niter,info,rk_vec] = myRichardson(A,b,nmaxiter,epsilon,alpha,x_0)
    // vérification sur A : il faut que tous les a_ii soient non nuls
    if ~and(diag(A)) then
        error('Presence de zéro sur la diagonale de A')
    end
    // Décomposition de A = D-E-F
    D = diag(diag(A));
    E = -tril(A) + D;
    F = -triu(A) + D;
    
    // initialisation
    x_sol = x_0;
    niter = 0; // aucune itération
    info  = 0;  // pas de convergence
    r_k = b - A*x_sol; // erreur à l'iétaration k;
    
    // résolution 
    for k = 1:nmaxiter
        x_sol =x_sol + alpha*(b - A*x_sol); 
        r_k = b - A*x_sol;
        rk_vec(k) = norm(r_k,2); // Pour garder l'historique des normes de r_k
        if  (norm(r_k,2) < epsilon)
            info = 1;
            niter = k;
            break;
         end
    end
endfunction 
 
