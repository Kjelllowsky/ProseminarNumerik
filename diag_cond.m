% Proseminar Numerik WS24/25 | Kjell Machalowsky
% Diese Funktion nutzt eine inverse Diagonalmatrix zur Vorkonditionierung
% des gegebenen Gleichungssystems.
%
% INPUTS
%  - A:     Systemmatrix
%  - b:     Vektor auf den der Vorkonditionierer angewandt werden soll
%
% OUTPUTS
%  - x:     Ergebnis der Anwendung des Vorkonditionierers
%  - kappa: Konditionszahl der konditionierten Matrix

function [x, kappa] = diag_cond(A,b)

    % argument validation
    [m1,m2] = size(A);
    [m3,m4] = size(b);
    assert(m1==m2, "A muss quadratisch sein.")
    assert((m4==1) && (m1==m3), "b hat die falsche Dimension.")

    % Berechnung (Cx=b <=> x = C^-1*b)
    K_inv = sparse(diag(sqrt(1./diag(A))));
    kappa = condest(K_inv'*K_inv*A,2);
    x = K_inv'*K_inv*b;
  
end