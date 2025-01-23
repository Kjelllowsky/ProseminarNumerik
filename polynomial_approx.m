% Function to approximate the inverse of a tridiagonal matrix by a
% polynomial form.
%
% INPUTS
%  - T: tridiag. matrix, whose inverse is to be approximated
%  - alpha: 1. coefficient of the polynome
%  - beta: 2. coefficient of the polynome
%
% OUTPUTS
%  - T_inv: approximate inverse of T
%
function T_inv = polynomial_approx(T,alpha,beta)
    [n,m] = size(T);
    % input validation
    assert(n==m,'The given sytem matrix must be quadratic!')

    % approximation
    D_T = diag(diag(T));
    D_T_inv = sparse(diag(1./diag(T)));
    T_inv = alpha*D_T_inv + beta*D_T_inv*(T-D_T)*D_T_inv;
    T_inv = sparse(T_inv);

end