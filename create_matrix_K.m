% Proseminar Numerik WS24/25 | Kjell Machalowsky
% Funktion zum Erstellen einer nicht sym. pos. def. Systemmatrix
%
% INPUT:
% - n  Dimension von A
%
% OUTPUT:
% - K  Matrix K
function K = create_matrix_K(n)
    % erstelle Hilfsmatrizen
    mat1 = eye(n)*3;
    vec2 = ones(n-1,1)*(-4);
    mat2 = diag(vec2,-1);
    vec3 = ones(n-2,1);
    mat3 = diag(vec3,-2);

    % erstelle K matrix
    K = mat1+mat2+mat3;
end
