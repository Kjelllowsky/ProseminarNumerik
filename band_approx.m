% Function to approximate the inverse of a tridiagonal matrix by
% its sparsity.
%
% INPUTS
%  - T: sym. tridiag. matrix, whose inverse is to be approximated
%  - p: number of sub-diagonals
%
% OUTPUTS
%  - T_inv: approximate inverse of T
%
function T_inv = band_approx(T,p)
    [n,m] = size(T);
    % argument validation
    assert(n==m,'The given sytem matrix must be quadratic!')

    % Initialisierung
    a = diag(T);
    b = -diag(T,1);
    u = 1;
    u = [u, a(1)/b(1)];

    % approximation
    for i=3:n
        u(i) = (a(i-1)*u(i-1)-b(i-2)*u(i-2))/b(i-1);
    end
    v = zeros(size(u));
    v(end) = 1/(-b(n-1)*u(n-1)+a(n)*u(n));
    for i=n-1:-1:2
        v(i) = (1+b(i)*u(i)*v(i+1))/(a(i)*u(i)-b(i-1)*u(i-1));
    end
    v(1) = (1+b(1)*u(1)*v(2))/(a(1)*u(1));
    T_inv = sparse(size(T));
    for i=1:n
        T_inv(i,i) = u(i)*v(i);
        if i+p<=n
            mark = p;
        else
            mark = n-i;
        end
        for j=1:mark
            T_inv(i,i+j) = u(i)*v(j+i);
            T_inv(i+j,i) = u(i)*v(j+i);
        end
    end
    
end