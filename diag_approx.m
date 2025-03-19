% Proseminar Numerik WS24/25 | Kjell Machalowsky
% Funktion zur Approximation einer Tridiagonal-Inversene über die Inverse
% der Diagonalmatrix.
%
% INPUTS
%  - T: sym. tridiag. Matrix deren Inverse approximiert werden soll
%
% OUTPUTS
%  - T_inv: approximierte Inverse von T
%
function T_inv = diag_approx(T)
    [n,m] = size(T);
    % argument validation
    assert(n==m,'The given sytem matrix must be quadratic!')
    
    % Approximation
    T_inv = sparse(diag(1./diag(T)));
end