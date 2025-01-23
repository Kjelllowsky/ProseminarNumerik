% Proseminar Numerik (Bachelor)
% Wintersemester 2024/25
% Kjell Machalowsky

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         ACHTUNG          %
% Ausführung kann u. U.    %
% mehrere Stunden brauchen %
% (für Dim. bis 16384)     %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Dieses Skript vergleicht verschiedene Ansätze der Vorkonditionierung im 
% PCG-Verfahren für lineare Gleichungssysteme. 
% 1 | Basis (CG-Verfahren)
% 2 | Diagonale Konditionierung
% 3 | SSOR-Verfahren
% 4 | incomplete-cholesky Zerlegung (Block-vorkonditionierung)
%   | a) Diagonal-Approximation der tridiagonalen Inversen
%   | b) Band-Approximation der tridiagonalen Inversen
%   | c) Cholesky-Block-Zerlegung der tridiagonalen Inversen
%   | d) Polynom-Approximation der tridiagonalen Inversen

%% Initialisierung
tol = 1e-8;
maxiter = 1e5;
% getestete Diskretisierungen: 2,4,8,12,16,24,32,64,96,128
dim_vec = [16,64,144,256,576,1024];
% dim_vec = [16,64,144,256,576];
% Störfunktion der Poisson-Gleichung
rhs_fun = @(x,y) x^2+y^2;
% Werte für Plots
residual_mat = zeros(7,length(dim_vec));
kappa_mat = zeros(7,length(dim_vec));
kappa_mat_sca = zeros(7,length(dim_vec));
iter_mat = zeros(7,length(dim_vec));
konv_mat = zeros(7,length(dim_vec));


%% Test durchführen
for i=1:length(dim_vec)
    dim = dim_vec(i);
    A = create_matrix_A(dim);
    b = rhs_2D_poisson_problem(dim,rhs_fun);
    x0 = zeros(dim,1);
    kappa_A = condest(A,2);
    disp("##############################")
    disp("Test für Dimension "+dim)
    disp("##############################")

    % 1 | Basisfall (CG-Verfahren)
    tic;
    [x_CG, res_CG] = cg_method(A,b,x0,tol,maxiter);
    toc
    residual_mat(1,i) = res_CG(end);
    kappa_mat(1,i) = kappa_A;
    kappa_mat_sca(1,i) = 1;
    iter_mat(1,i) = length(res_CG);
    konv_mat(1,i) = exp_konv_rate(res_CG);
    disp("Fall 1 fertig")

    % 2 | Diagonale Konditionierung
    tic;
    [x_diag, res_diag, kappa_diag] = my_pcg(A,b,x0,tol,maxiter,@(C,g) diag_cond(C,g));
    toc;
    residual_mat(2,i) = res_diag(end);
    kappa_mat(2,i) = kappa_diag;
    kappa_mat_sca(2,i) = kappa_diag/kappa_A;
    iter_mat(2,i) = length(res_diag);
    konv_mat(2,i) = exp_konv_rate(res_diag);
    disp("Fall 2 fertig")

    % 3 | SSOR-Verfahren
    tic;
    [x_ssor, res_ssor, kappa_ssor] = my_pcg(A,b,x0,tol,maxiter,@(C,g) ssor_cond(C,g,1.6));
    toc;
    residual_mat(3,i) = res_ssor(end);
    kappa_mat(3,i) = kappa_ssor;
    kappa_mat_sca(3,i) = kappa_ssor/kappa_A;
    iter_mat(3,i) = length(res_ssor);
    konv_mat(3,i) = exp_konv_rate(res_ssor);
    disp("Fall 3 fertig")

    % 4 | IBC a) Diagonal-Approximation
    tic;
    [x_IBCa, res_IBCa, kappa_IBCa] = my_pcg(A,b,x0,tol,maxiter,@(C,g) ibc_cond(C,g,@(T) diag_approx(T)));
    toc;
    residual_mat(4,i) = res_IBCa(end);
    kappa_mat(4,i) = kappa_IBCa;
    kappa_mat_sca(4,i) = kappa_IBCa/kappa_A;
    iter_mat(4,i) = length(res_IBCa);
    konv_mat(4,i) = exp_konv_rate(res_IBCa);
    disp("Fall 4a) fertig")

    % 4 | IBC b) Band-Approximation
    tic
    [x_IBCb, res_IBCb, kappa_IBCb] = my_pcg(A,b,x0,tol,maxiter,@(C,g) ibc_cond(C,g,@(T) band_approx(T,1)));
    toc
    residual_mat(5,i) = res_IBCb(end);
    kappa_mat(5,i) = kappa_IBCb;
    kappa_mat_sca(5,i) = kappa_IBCb/kappa_A;
    iter_mat(5,i) = length(res_IBCb);
    konv_mat(5,i) = exp_konv_rate(res_IBCb);
    disp("Fall 4b) fertig")

    % 4 | IBC c) Cholesky-Block-Zerlegung
    tic;
    [x_IBCc, res_IBCc, kappa_IBCc] = my_pcg(A,b,x0,tol,maxiter,@(C,g) ibc_cond(C,g,@(T) cholesky_approx(T)));
    toc;
    residual_mat(6,i) = res_IBCc(end);
    kappa_mat(6,i) = kappa_IBCc;
    kappa_mat_sca(6,i) = kappa_IBCc/kappa_A;
    iter_mat(6,i) = length(res_IBCc);
    konv_mat(6,i) = exp_konv_rate(res_IBCc);
    disp("Fall 4c) fertig")

    % 4 | IBC d) Polynom-Approximation
    tic;
    [x_IBCd, res_IBCd, kappa_IBCd] = my_pcg(A,b,x0,tol,maxiter,@(C,g) ibc_cond(C,g,@(T) polynomial_approx(T,1.14,-1.14)));
    toc;
    residual_mat(7,i) = res_IBCd(end);
    kappa_mat(7,i) = kappa_IBCd;
    kappa_mat_sca(7,i) = kappa_IBCd/kappa_A;
    iter_mat(7,i) = length(res_IBCd);
    konv_mat(7,i) = exp_konv_rate(res_IBCd);
    disp("Fall 4d) fertig")
end

% Ergebnisse speichern
save("residual_mat.mat", "residual_mat")
save("kappa_mat.mat", "kappa_mat")
save("kappa_mat_sca.mat", "kappa_mat_sca")
save("iter_mat.mat", "iter_mat")
save("konv_mat.mat", "konv_mat")


%% Plots erzeugen

% optional: Matrizen laden
residual_mat = load("residual_mat.mat");
kappa_mat = load("kappa_mat.mat");
kappa_mat_sca = load("kappa_mat_sca.mat");
iter_mat = load("iter_mat.mat");
konv_mat = load("konv_mat.mat");

% plot 1: Anzahl Iterationen über Dimension für verschiedene Verfahren
fig1 = figure;
loglog(dim_vec,iter_mat(1,:), 'k-',Marker='o', DisplayName='CG-Verfahren')
hold on
loglog(dim_vec,iter_mat(2,:), 'r-',Marker='+', DisplayName='DIAG')
loglog(dim_vec,iter_mat(3,:), 'b-',Marker='*', DisplayName='SSOR')
loglog(dim_vec,iter_mat(4,:), 'g-',Marker='x', DisplayName='IBC-DIAG')
loglog(dim_vec,iter_mat(5,:), 'm-',Marker='square', DisplayName='IBC-BAND')
loglog(dim_vec,iter_mat(6,:), 'c-',Marker='diamond', DisplayName='IBC-CHOL')
loglog(dim_vec,iter_mat(7,:), '-',Color='#D95319', Marker='^', DisplayName='IBC-POLY')
hold off
title('Benötigte Iterationen nach Systemgröße')
xlabel('Systemgröße [-]')
ylabel('Anzahl Iterationen [-]')
xlim([0,dim_vec(end)*1.05])
ylim([0,2*max(max(iter_mat))])
legend(Location="northwest")


% plot 2: Verlauf des Residuums über Iterationszahl für versch. Verfahren
% plot 2a): CG-Verfahren
% fig2a = figure;
% semilogy


% plot 3: Konvergenzrate über Systemgröße für verschiedene Verfahren
fig3 = figure;
semilogx(dim_vec,konv_mat(1,:), 'k-',Marker='o', DisplayName='CG-Verfahren')
hold on
semilogx(dim_vec,konv_mat(2,:), 'r-',Marker='+', DisplayName='DIAG')
semilogx(dim_vec,konv_mat(3,:), 'b-',Marker='*', DisplayName='SSOR')
semilogx(dim_vec,konv_mat(4,:), 'g-',Marker='x', DisplayName='IBC-DIAG')
semilogx(dim_vec,konv_mat(5,:), 'm-',Marker='square', DisplayName='IBC-BAND')
semilogx(dim_vec,konv_mat(6,:), 'c-',Marker='diamond', DisplayName='IBC-CHOL')
semilogx(dim_vec,konv_mat(7,:), '-',Color='#D95319', Marker='^', DisplayName='IBC-POLY')
hold off
title('Konvergenzrate der Verfahren nach Systemgröße')
xlabel('Systemgröße [-]')
ylabel('Konvergenzrate [-]')
xlim([0,dim_vec(end)*1.05])
ylim([0,1.05*max(max(konv_mat))])
legend(Location="northeast")


% plot 4: Quotient der Kond.-zahlen über Systemgröße für versch. Verfahren
fig4 = figure;
loglog(dim_vec,kappa_mat_sca(1,:), 'k-',Marker='o', DisplayName='CG-Verfahren')
hold on
loglog(dim_vec,kappa_mat_sca(2,:), 'r-',Marker='+', DisplayName='DIAG')
loglog(dim_vec,kappa_mat_sca(3,:), 'b-',Marker='*', DisplayName='SSOR')
loglog(dim_vec,kappa_mat_sca(4,:), 'g-',Marker='x', DisplayName='IBC-DIAG')
loglog(dim_vec,kappa_mat_sca(5,:), 'm-',Marker='square', DisplayName='IBC-BAND')
loglog(dim_vec,kappa_mat_sca(6,:), 'c-',Marker='diamond', DisplayName='IBC-CHOL')
loglog(dim_vec,kappa_mat_sca(7,:), '-',Color='#D95319', Marker='^', DisplayName='IBC-POLY')
hold off
title('Quotient der Konditionszahl nach Systemgröße')
xlabel('Systemgröße [-]')
ylabel('$\displaystyle \frac{\kappa(\tilde{A})}{\kappa(A)}$', 'Interpreter','latex')
xlim([0,dim_vec(end)*1.05])
ylim([0,2*max(max(kappa_mat_sca))])
legend(Location="northwest")


