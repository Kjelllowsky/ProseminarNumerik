% Function to approximate the inverse of a tridiagonal matrix by
% its diagonal inverse.
%
% INPUTS
%  - T: sym. tridiag. matrix, whose inverse is to be approximated
%
% OUTPUTS
%  - T_inv: approximate inverse of T
%
function T_inv = diag_approx(T)
    [n,m] = size(T);
    % argument validation
    assert(n==m,'The given sytem matrix must be quadratic!')
    
    % Approximation
    T_inv = sparse(diag(1./diag(T)));
end