//info pour signifier si la méthode converge ou non
//niter : nombre d'iteration à la convergence
// x solution à la convergence
// nmaxietr nombre maximal d'iétéartion que la méthode converge ou pas
// epsilonn: seuil de tolérance

function[x_sol,r_k,niter,info,rk_vec] = myJacobi(A,b,nmaxiter,epsilon,x_0)
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
    r_k = b - A*b; // erreur à l'iétaration k;
    // résolution 
    for k = 1:nmaxiter
        x_sol = inv(D)*((E+F)*x_sol + b); 
        r_k = b - A*x_sol;
        rk_vec(k) = norm(r_k,2); // Pour garder l'historique des normes de r_k
        if  (norm(r_k) < epsilon)
            info = 1; //convergence
            niter = k;
            break;
         end
    end
endfunction
