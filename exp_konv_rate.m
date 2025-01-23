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
    % Normen der Differenzen berchnen
    d1 = norm(res(end-1)-res(end),2);
    d2 = norm(res(end-2)-res(end),2);
    d3 = norm(res(end-3)-res(end),2);

    % experimentelle Konvergenzrate berchnen
    q = log(d2/d1) / log(d3/d2);

end