% Proseminar Numerik WS24/25 | Kjell Machalowsky
% Diese Funktion nutzt die Incomplete-Block-Cholesky Zerlegung zur Vorkonditionierung.
%
% INPUTS
%  - A:      sparse Systemmatrix
%  - b:      Vektor auf den der Vorkonditionierer angewandt werden soll
%  - approx: Function-Handle für die Approximation der Inversen von
%            tridiagonalen Matrizen innerhalb des Verfahrens
%
% OUTPUTS
%  - x:      Ergebnis der Anwendung des Vorkonditionierers
%  - kappa:  Konditionszahl der konditionierten Matrix

function [x, kappa] = ibc_cond(A,b,approx)

    % argument validation
    [m1,m2] = size(A);
    [m3,m4] = size(b);
    assert(m1==m2, "A muss quadratisch sein.")
    assert(sqrt(m1)==round(sqrt(m1)), "A hat die falsche Dimension");
    assert((m4==1) && (m1==m3), "b hat die falsche Dimension.")

    % Initialisierung
    n = sqrt(m1);
    D = A(1:n,1:n);
    Delta_1 = D;
    T = approx(Delta_1);
    L_1 = chol(Delta_1)';
    L_1_inv = chol(T);
    
    K = spalloc(m1,m1, nnz(A));
    K(1:n,1:n) = L_1;
    K(n+1:2*n,1:n) = -L_1_inv';

    % Iterationsschleife
    for k = 2:n-1
        Delta_k = sparse(D - T);
        L_k = chol(Delta_k)';
        T = approx(Delta_k);
        L_k_inv = chol(T);
        
        K((k-1)*n+1:k*n,(k-1)*n+1:k*n) = L_k;
        K(k*n+1:(k+1)*n,(k-1)*n+1:k*n) = -L_k_inv';
    end

    Delta_n = D - T;
    L_n = chol(Delta_n);
    K(m1-n+1:end,m1-n+1:end) = L_n;
    K = sparse(K);
    
    % Konditionierung
    K_inv = K^-1;
    y = K\b;
    x = K'\y;
    kappa = condest(K_inv'*K_inv*A,2);


end