% Function to approximate the inverse of a tridiagonal matrix by a
% cholesky decomposition.
%
% INPUTS
%  - T: tridiag. matrix, whose inverse is to be approximated
%
% OUTPUTS
% - T_inv: approximate inverse of T
%
function T_inv = cholesky_approx(T)
    [n,m] = size(T);
    % argument validation
    assert(n==m,'The given sytem matrix must be quadratic!')
    
    % initialization
    L = sparse(chol(T)');
    gamma = diag(L);
    a = diag(T);
    b = -diag(T,1);
    T_inv = sparse(size(T));
    gamma_i_sq = gamma(1)^2;
    gamma_ip1_sq = gamma(2)^2;

    % approximation
    T_inv(1,1) = 1/(gamma_i_sq) + (a(2)-gamma_ip1_sq)/(gamma_i_sq*gamma_ip1_sq);
    T_inv(1,2) = b(1)/(gamma_i_sq*gamma_ip1_sq);
    for i=2:n-1
        gamma_i_sq = gamma_ip1_sq;
        gamma_ip1_sq = gamma(i+1)^2;
        T_inv(i,i-1) = T_inv(i-1,i);
        T_inv(i,i) = 1/(gamma_i_sq) + (a(i)-gamma_ip1_sq)/(gamma_i_sq*gamma_ip1_sq);
        T_inv(i,i+1) = b(i)/(gamma_i_sq*gamma_ip1_sq);
    end
    T_inv(n,n-1) = T_inv(n-1,n);
    T_inv(n,n) = 1/(gamma_i_sq) + (a(n)-gamma_i_sq)/(gamma_i_sq^2);

end