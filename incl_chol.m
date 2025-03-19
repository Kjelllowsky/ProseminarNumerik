% Proseminar Numerik WS24/25 | Kjell Machalowsky
% Incomplete-Cholesky Zerlegung mit Neuordnung un Nullsetzen
%
% INPUTS
%  - A:        Matrix deren ICZ bestimmt werden soll
%  - zero_tol: Toleranz für das Nullsetzen der Matrixeinträge
%
% OUTPUTS
% - L: ICZ Faktor-Matrix

function L = incl_chol(A, zero_tol)
    [n,m] = size(A);
    % argument validation
    assert(n==m,'The given sytem matrix must be quadratic!')
    assert(isequal(A,A'), 'The given matrix must be sym. and pos. definite.')
    
    % Neuordnung von A durchführen
    p = amd(A);
    A_reordered = A(p,p);
    L = zeros(n,n);

    % Incomplete-Cholesky mit Nullsetzen
    for k=1:n
        % Diagonale
        L(k,k) = sqrt(A_reordered(k,k) - sum(L(k,1:k-1).^2));
        % Nebendiagonale
        for j=k+1:n
            if A_reordered(j,k) ~=0
                L(j,k) = (A_reordered(j,k) - sum(L(j,1:k-1).*L(k,1:k-1))) / L(k,k);
                if abs(L(j,k)) < zero_tol
                    L(j,k) = 0;
                end
            end
        end
    end

    L(p,:) = L;
    L = sparse(L);

end