% Proseminar Numerik (Bachelor)
% Wintersemester 2024/25
% Kjell Machalowsky

% Dieses Skript vergleicht die Effektivität der Diagonal-Vorkonditionierung
% für verschiedene Strukturen an Matrizen.

%% Initialisierung
% Matrizen erstellen
m = 9;
A = create_matrix_A(m)*eye(m);
B = [1 3  1 9 0  1 0 1  9;
     4 10 7 7 5  3 5 4  7;
     4 1  1 7 8  9 3 4  6;
     7 7  8 7 6  8 4 2  5;
     0 7  8 5 1  7 3 7  3;
     6 7  2 6 4  9 8 6  5;
     3 4  0 3 3  2 4 10 1;
     4 3  0 5 8  8 5 5  3;
     8 9  2 1 10 2 9 6  6];
C = eye(m)*4 + diag(ones(1,m-1)*-0.25,-1) + diag(ones(1,m-1)*-0.25,1);

% ursprüngliche Konditionszahlen berechnen
kappa_A = cond(A,2);
kappa_B = cond(B,2);
kappa_C = cond(K,2);

%% Diagonal-Konditionierung
K_inv_A = diag(sqrt(1./diag(A)));
K_inv_B = diag(sqrt(1./diag(B)));
K_inv_C = diag(sqrt(1./diag(C)));

tilde_A = K_inv_A*A*K_inv_A';
tilde_B = K_inv_B*B*K_inv_B';
tilde_C = K_inv_C*C*K_inv_C';

% neue Konditionszahlen berechnen
kappa_A_tilde = cond(tilde_A,2);
kappa_B_tilde = cond(tilde_B,2);
kappa_C_tilde = cond(tilde_C,2);
