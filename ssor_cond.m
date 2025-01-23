% Proseminar Numerik WS24/25 | Kjell Machalowsky
% Funktion zur Anwednung der SSOR-Vorkonditionierung.
%
% INPUTS
%  - A:     sparse Systemmatrix
%  - b:     Vektor auf den der Vorkonditionierer angewandt werden soll
%  - omega: Relaxationsparameter
%
% OUTPUTS
%  - x:     Ergebnis der Anwendung des Vorkonditionierers
%  - kappa: Konditionszahl der konditionierten Matrix

function [x,kappa] = ssor_cond(A,b,omega)

    % argument validation
    [m1,m2] = size(A);
    [m3,m4] = size(b);
    assert(m1==m2, "A muss quadratisch sein.")
    assert((m4==1) && (m1==m3), "b hat die falsche Dimension.")
    assert((omega<2) && (omega>0), "Relaxationsparamter nicht im Intervall (0,2).")
    
    % Initalisierung
    D = spdiags(diag(A),0, size(A,1),size(A,2));
    D_ = spdiags(sqrt(diag(D)),0, size(D,1),size(D,2));
    D_inv = spdiags(1./sqrt(diag(D)),0, size(D,1),size(D,2));
    L = tril(A)-D;

    % Berechnung (Dreiecksform von K fuer LR-Zerlegung nutzen)
    % K*K' * x = b
    K = (D_ + omega*L*D_inv);
    K_inv = K^-1;
    y = K\b;
    x = K'\y;
    kappa = condest(K_inv'*K_inv,2);

end