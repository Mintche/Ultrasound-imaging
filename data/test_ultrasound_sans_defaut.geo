// test_ultrasound_nodefect.geo
// Domaine rectangulaire SANS défaut
// Physical tags:
//  - Surface MILIEU : 1
//  - SIGMA_PLUS_L (gauche) : 11
//  - SIGMA_MINUS_L (droite): 12
//  - BORD_HAUT : 13
//  - BORD_BAS  : 14

Mesh.MshFileVersion = 2.2;

SetFactory("OpenCASCADE");

// ----------------------
// Paramètres géométriques
// ----------------------
Lx = 1.0;    // largeur du rectangle
Ly = 0.6;    // hauteur du rectangle

// ----------------------
// Paramètres de maillage
// ----------------------
h_bulk = 0.04;   // taille dans le milieu (uniforme ici)

Mesh.CharacteristicLengthMin = h_bulk;
Mesh.CharacteristicLengthMax = h_bulk;

// ----------------------
// Géométrie
// ----------------------
Rectangle(1) = {0, 0, 0, Lx, Ly};   // surface unique

// ----------------------
// Physical groups (surfaces)
// ----------------------
Physical Surface(1) = {1}; // MILIEU

// ----------------------
// Physical groups (bords)
// ----------------------
// On récupère les courbes du contour du rectangle via BoundingBox.
// Gmsh ne garantit pas l'ordre des courbes, donc on les classe par position.

eps = 1e-6;

// Courbes sur x=0 (gauche)
cLeft[]  = Curve In BoundingBox {-eps, -eps, -1, eps, Ly+eps, 1};
// Courbes sur x=Lx (droite)
cRight[] = Curve In BoundingBox {Lx-eps, -eps, -1, Lx+eps, Ly+eps, 1};
// Courbes sur y=Ly (haut)
cTop[]   = Curve In BoundingBox {-eps, Ly-eps, -1, Lx+eps, Ly+eps, 1};
// Courbes sur y=0 (bas)
cBot[]   = Curve In BoundingBox {-eps, -eps, -1, Lx+eps, eps, 1};

Physical Curve(11) = {cLeft[]};   // SIGMA_PLUS_L
Physical Curve(12) = {cRight[]};  // SIGMA_MINUS_L
Physical Curve(13) = {cTop[]};    // BORD_HAUT
Physical Curve(14) = {cBot[]};    // BORD_BAS

// Optionnel : forcer triangles
Mesh.RecombineAll = 0;
