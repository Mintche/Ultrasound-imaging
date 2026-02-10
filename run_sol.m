run("mesh_out.m");       % Charge Coor, Tri, etc.
run("solution_out.m");   % Charge Ure, Uim

% --- Paramètres ---
h = 0.6; 
k0 = 30;

% Reconstruction du nombre complexe
U_total = Ure + 1i*Uim;

% --- Construction du champ incident ---
% Sur un maillage non-structuré, X est la première colonne de Coor
X = Coor(:,1); 
Phi_inc = (1/sqrt(h)) * exp(1i * k0 * X); 

% --- Calcul du champ diffracté (Scattered field) ---
U_scat = U_total - Phi_inc;

% --- Affichage ---
figure('Position', [100, 100, 1200, 400]);

% 1. Affichage du Champ Total (pour comparer)
subplot(1,2,1);
trisurf(Tri(:,1:3), Coor(:,1), Coor(:,2), real(U_total), 'EdgeColor', 'none');
view(2); axis equal; shading interp; colorbar; colormap jet;
title('Champ Total (Re)');
xlabel('x'); ylabel('y');

% 2. Affichage du Champ Diffracté (C'est ce qui nous intéresse pour LSM)
subplot(1,2,2);
trisurf(Tri(:,1:3), Coor(:,1), Coor(:,2), real(U_scat), 'EdgeColor', 'none');
view(2); axis equal; shading interp; colorbar; colormap jet;
title('Champ Diffracté u^s (Re)');
xlabel('x'); ylabel('y');