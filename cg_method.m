% Proseminar Numerik WS24/25 | Kjell Machalowsky
% Funktion zur Anwendung des Conjugierte-Gradienten Verfahrens.
% 
% INPUTS
%  - A:       Systemmatrix
%  - b:       Rechte-Seite-Vektor
%  - x0:      Startwert
%  - tol:     Toleranz für das Abbruchkriterium
%             ||Ax_k -b||_2<tol*||Ax_0-b||_2
%  - maxiter: maximale Anzahl an Iterationen
%
% OUTPUTS
%  - x:       Lösungsvektor
%  - res:     Vektor der Länge #iter Normen der Residuen
%       ||Ax_k - b||_2

function [x, res] = cg_method(A,b,x0,tol,maxiter)

    % argument validation
    [m1,m2] = size(A);
    [m3,m4] = size(b);
    [m5,m6] = size(x0);
    assert(m1==m2, "A muss quadratisch sein.")
    assert((m4==1) && (m6==1) && (m1==m3) && (m3==m5), "Startwert oder Störvektor haben falsche Dimension.")
    assert((maxiter>0) && (tol > 0), "Max. Iterationszahl und Toleranz müssen größer null sein.")

    % Initialisierung
    x = x0;
    d = b - A*x;
    g = -d;
    norm_g = norm(g,2);
    k = 0;
    res = norm(b,2);
    tol_tilde = tol*res;

    % CG-Verfahren anwenden
    while k<maxiter
        alpha = norm_g^2 / (d'*(A*d));
        x = x + alpha*d;
        g_kp1 = g + alpha*A*d;
        norm_gkp1 = norm(g_kp1,2);
        beta = norm_gkp1^2 / norm_g^2;
        d = -g_kp1 + beta*d;

        g = g_kp1;
        norm_g = norm_gkp1;
        res = [res, norm(A*x-b)];
        k = k+1;

        % teste Abbruchkriterium
        if res(end) <= tol_tilde
            break
        end
    end
    
end