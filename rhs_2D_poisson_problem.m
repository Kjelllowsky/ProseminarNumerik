% Proseminar Numerik WS24/25 | Kjell Machalowsky
% Funktion zum Erstellen der rechten Seite für das zu lösende LGS
%
% INPUTS
%  - n:  Dimension von A
%  - f:  Function-Handle der Störfunktion
%        f = 0 erzeugt das homogene Problem
%
% OUTPUTS
%  - b:  Vektor b (Rechte Seite 

function b = rhs_2D_poisson_problem(n,f)
    % argument validation
    assert(n~=0, "A muss eine Matrix sein.")

    % Initialisierung
    b = zeros(n,1);
    m = sqrt(n);
    h = 1/m;
    x_grid = linspace(0, 1, m+2);
    y_grid = linspace(0, 1, m+2);

    % Grundstruktur aufbauen
    for i=1:m
        for j=1:m
            index = i + (j-1)*m;
            b(index) = h^2 * f(x_grid(i+1),y_grid(j+1));
        end
    end

    % Randbedingungen
    RB_x0 = @(x) sin(2*pi*x);
    RB_x1 = @(x) sin(2*pi*x);
    RB_0y = @(x) sin(2*pi*x);
    RB_1y = @(x) sin(2*pi*x);
    for i=1:m
        % untere RB (x,0)
        b(i) = b(i) + RB_x0(x_grid(i+1));
        % obere RB (x,1)
        b(i+(m-1)*m) = b(i+(m-1)*m) + RB_x1(x_grid(i+1));
    end
    for j=1:m
        % linke RB (0,y)
        b((j-1)*m + 1) = b((j-1)*m + 1) + RB_0y(y_grid(i+1));
        % rechte RB (1,y)
        b(j*m) = b(j*m) + RB_1y(y_grid(i+1));
    end

end
