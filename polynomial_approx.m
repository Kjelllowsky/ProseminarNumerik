% Proseminar Numerik WS24/25 | Kjell Machalowsky
% Funktion zur Approximation einer Tridiagonal-Inversen über einen
% Polynom-Ansatz.
%
% INPUTS
%  - T:     tridiag. Matrix deren Inverse bestimmt werden soll
%  - alpha: 1. Koeffizient des Polynom-Ansatzes
%  - beta:  2. Koeffizient des Polynom-Ansatzes
%
% OUTPUTS
%  - T_inv: approximierte Inverse von T
%
function T_inv = polynomial_approx(T,alpha,beta)
    [n,m] = size(T);
    % input validation
    assert(n==m,'The given sytem matrix must be quadratic!')

    % Approximation
    D_T = diag(diag(T));
    D_T_inv = sparse(diag(1./diag(T)));
    T_inv = alpha*D_T_inv + beta*D_T_inv*(T-D_T)*D_T_inv;
    T_inv = sparse(T_inv);

end