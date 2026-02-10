// test_ultrasound_funnel.geo
// Tube en forme d'entonnoir (trapèze), sans défaut
// Physical tags:
//  - Surface MILIEU : 1
//  - ENTREE (gauche, grand côté) : 11
//  - SORTIE (droite, petit côté) : 12
//  - PAROI_HAUT : 13
//  - PAROI_BAS  : 14

Mesh.MshFileVersion = 2.2;
SetFactory("OpenCASCADE");

// ----------------------
// Paramètres géométriques
// ----------------------
Lx = 1.0;      // longueur
Ly_in  = 0.6;  // hauteur à l'entrée (x=0)
Ly_out = 0.3;  // hauteur à la sortie (x=Lx)

// Pour centrer verticalement l'entonnoir autour de y = Ly_in/2 (optionnel)
y0_in  = 0.0;                    // bas à l'entrée
y1_in  = y0_in + Ly_in;          // haut à l'entrée
// On centre la sortie sur le même axe médian :
yc     = 0.5*(y0_in + y1_in);
y0_out = yc - 0.5*Ly_out;        // bas à la sortie
y1_out = yc + 0.5*Ly_out;        // haut à la sortie

// ----------------------
// Paramètres de maillage
// ----------------------
h_bulk = 0.04;
Mesh.CharacteristicLengthMin = h_bulk;
Mesh.CharacteristicLengthMax = h_bulk;

// ----------------------
// Géométrie (trapèze)
// ----------------------
Point(1) = {0,   y0_in,  0, h_bulk};   // entrée bas
Point(2) = {0,   y1_in,  0, h_bulk};   // entrée haut
Point(3) = {Lx,  y1_out, 0, h_bulk};   // sortie haut
Point(4) = {Lx,  y0_out, 0, h_bulk};   // sortie bas

Line(1) = {1,2};  // entrée (gauche)
Line(2) = {2,3};  // paroi haute
Line(3) = {3,4};  // sortie (droite)
Line(4) = {4,1};  // paroi basse

Curve Loop(1) = {1,2,3,4};
Plane Surface(1) = {1};

// ----------------------
// Physical groups
// ----------------------
Physical Surface(1) = {1}; // MILIEU

Physical Curve(11) = {1}; // ENTREE (SIGMA_PLUS_L)
Physical Curve(12) = {3}; // SORTIE (SIGMA_MINUS_L)
Physical Curve(13) = {2}; // PAROI_HAUT
Physical Curve(14) = {4}; // PAROI_BAS

// Optionnel : forcer triangles
Mesh.RecombineAll = 0;
