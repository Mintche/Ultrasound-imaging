run("mesh_out.m");   % charge Coor, Tri, RefTri, EdgeRef, IsDef

figure; hold on;

Z = zeros(size(Coor,1),1);

% Triangles "milieu" (non défaut)
trisurf(Tri(IsDef==0,1:3), Coor(:,1), Coor(:,2), Z, ...
        'FaceColor', [0.8 0.8 0.8], 'EdgeColor', [0.3 0.3 0.3]);

% Triangles "défaut"
trisurf(Tri(IsDef==1,1:3), Coor(:,1), Coor(:,2), Z, ...
        'FaceColor', [1 0.2 0.2], 'EdgeColor', [0.3 0.3 0.3]);

view(2); axis equal; box on;
title("Défaut en rouge");