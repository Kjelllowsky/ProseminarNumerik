% Proseminar Numerik WS24/25 | Kjell Machalowsky
% Funktion zur Anwendung des PCG-Verfahrens.
% 
% INPUTS
%  - A:         Systemmatrix
%  - b:         Rechte-Seite-Vektor
%  - x:         Startvektor
%  - tol:       Toleranz für das Abbruchkriterium ||A*x_n - b|| < tol * ||A*x0 - b||
%  - max_iter:  maximale Anzahl an Iterationen
%  - precon:    Function-Handle auf den Vorkonditionierer. Der Vorkonditionierer
%               soll den Funktionskopf y = precon(A,d) besitzen:
%               Das Bedeutet: Der Vorkonditionierer bekommt die Systemmatrix
%               A übergeben und gibt die Anwendung C^(-1)*d zurück
%
% OUTPUTS
%  - x:         Lösungsvektor
%  - res:       Vektor der Residuen in der euklidischen Norm zu jeder Iteration
%               angefangen mit dem Startresiduum ||b - A * x_0||
%               Länge ist gegeben durch die maximale Anzahl an Iterationen plus 1
%  - kappa:     Vektor mit Konditionszahlen der konditionierten Matrizen

function [x,res,kappa] = my_pcg(A,b,x,tol,maxiter,precon)

    % argument validation
    [m1,m2] = size(A);
    [m3,m4] = size(b);
    [m5,m6] = size(x);
    assert(m1==m2, "Die Systemmatrix muss quadratisch sein.")
    assert((m4==1) && (m6==1) && (m1==m3) && (m3==m5), "Startwert oder Störvektor haben falsche Dimension.")
    assert((maxiter>0) && (tol > 0), "Max. Iterationszahl und Toleranz müssen größer null sein.")

    % Initialisierung
    iter = 0;
    g = A*x-b;
    [rho, kappa] = precon(A,g);
    d = -rho;
    A_d = A*d;
    g0_norm = norm(g);
    grho_SP2 = g' * rho;
    res = g0_norm;
    tol_tilde = tol*g0_norm;
    
    % PCG-Verfahren anwenden
    while (iter < maxiter)
        alpha = grho_SP2 / (d' * A_d);
        x = x + alpha*d;
        g = g + alpha*A_d;
        [rho, ~] = precon(A,g);
        grho_SP2_next = g' * rho;
        beta = grho_SP2_next/grho_SP2;
        grho_SP2 = grho_SP2_next;
        d = -rho + beta*d;
        A_d = A*d;
        res = [res, norm(g)];
        iter = iter+1;

        % teste Abbruchkriterium
        if res(end) <= tol_tilde
            break
        end
    end
end