% Proseminar Numerik WS24/25 | Kjell Machalowsky
% Funktion zum Erstellen der Systemmatrix in verschiedenen Größen
%
% INPUT:
% - n  Dimension von A
%
% OUTPUT:
% - A  Matrix A
function A = create_matrix_A(n)
    assert(sqrt(n) == round(sqrt(n)));
    m = sqrt(n);
    B = eye(m)*4 + diag(ones(1,m-1)*-1,-1) + diag(ones(1,m-1)*-1,1);
    A = kron(eye(m),B) + diag(ones(1,n-m)*-1,-m) + diag(ones(1,n-m)*-1,m);,
    A = sparse(A);
end