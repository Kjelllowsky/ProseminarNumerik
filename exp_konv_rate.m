% Proseminar Numerik WS24/25 | Kjell Machalowsky
% Diese Funktion berechnet die experimentelle Konvergenzrate für ein Verfahren
% mit gegeben Fehlervektor.
%
% INPUTS
%  - res:    Vektor der Residuen für gewähltes Verfahren
%
% OUTPUTS
%  - q:      experimentelle Konvergenzrate

function q = exp_konv_rate(res)
    n = length(res);

    % Nomen der Diferenzen berchnen
    d1_vec = zeros(n-3,1);
    d2_vec = zeros(n-3,1);
    d3_vec = zeros(n-3,1);
    for k=4:n
        d1_vec(k-3) = norm(res(k-1)-res(k),2);
        d2_vec(k-3) = norm(res(k-2)-res(k),2);
        d3_vec(k-3) = norm(res(k-3)-res(k),2);
    end

    % durchschnittliche Normen der Differenzen berchnen
    d1 = mean(d1_vec);
    d2 = mean(d2_vec);
    d3 = mean(d3_vec);

    % experimentelle Konvergenzrate berchnen
    q = log(d2/d1) / log(d3/d2);

end