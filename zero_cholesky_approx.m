% Proseminar Numerik WS24/25 | Kjell Machalowsky
% Incomplete-Cholesky Zerlegung mit Neuordnung und Nullsetzen
%
% INPUTS
%  - T:        tridiag. Matrix deren Inverse bestimmt werden soll
%  - zero_tol: Toleranz für das Nullsetzen der Matrixeinträge
%
% OUTPUTS
% - T_inv: approximierte Inverse von T

function T_inv = zero_cholesky_approx(T, zero_tol)
    [n,m] = size(T);
    % argument validation
    assert(n==m,'The given sytem matrix must be quadratic!')
    assert(isequal(T,T'), 'The given matrix must be sym. and pos. definite.')

    L = incl_chol(T,zero_tol);
    L_inv = L^(-1);

    T_inv = L_inv'*L_inv;
end