// test_ultrasound_small_defect_topwall.geo
// Tube rectangle + petit défaut circulaire tangent à la paroi supérieure
// Physical tags:
//  - Surface MILIEU : 1
//  - Surface DEFAUT : 2
//  - SIGMA_PLUS_L (gauche) : 11
//  - SIGMA_MINUS_L (droite): 12
//  - BORD_HAUT : 13
//  - BORD_BAS  : 14

Mesh.MshFileVersion = 2.2;

SetFactory("OpenCASCADE");

// ----------------------
// Paramètres géométriques
// ----------------------
Lx = 1.0;    // largeur du tube
Ly = 0.6;    // hauteur du tube

// Petit défaut collé à la paroi supérieure (tangent à y=Ly)
r  = 0.03;   // rayon du défaut (petit)
cx = 0.55;   // position en x du défaut (à ajuster)
cy = Ly - r; // pour être tangent à la paroi supérieure

// ----------------------
// Paramètres de maillage
// ----------------------
h_bulk   = 0.04;   // taille dans le milieu
h_defaut = 0.012;  // taille près du défaut

Mesh.CharacteristicLengthMin = h_defaut;
Mesh.CharacteristicLengthMax = h_bulk;

// ----------------------
// Géométrie
// ----------------------
Rectangle(1) = {0, 0, 0, Lx, Ly};          // Surface 1 (provisoire)
Disk(2)      = {cx, cy, 0, r, r};          // Surface 2 (provisoire)

// Découpe du rectangle par le disque (le disque devient la surface défaut)
out[] = BooleanFragments{ Surface{1}; Delete; }{ Surface{2}; Delete; };

// Récupération des surfaces résultantes
sDef[] = Surface In BoundingBox {cx-r-1e-6, cy-r-1e-6, -1, cx+r+1e-6, cy+r+1e-6, 1};
sAll[] = Surface In BoundingBox {0-1e-6, 0-1e-6, -1, Lx+1e-6, Ly+1e-6, 1};

sMil[] = sAll[];
For i In {0:#sDef[]-1}
  sMil[] -= {sDef[i]};
EndFor

// ----------------------
// Champs de taille (raffinement près du défaut)
// ----------------------
Field[1] = Distance;
Field[1].SurfacesList = {sDef[]};
Field[1].NumPointsPerCurve = 100;

Field[2] = Threshold;
Field[2].InField = 1;
Field[2].SizeMin = h_defaut;
Field[2].SizeMax = h_bulk;
Field[2].DistMin = r;
Field[2].DistMax = 3*r;

Background Field = 2;

// ----------------------
// Physical groups (surfaces)
// ----------------------
Physical Surface(1) = {sMil[]}; // MILIEU
Physical Surface(2) = {sDef[]}; // DEFAUT

// ----------------------
// Physical groups (bords)
// ----------------------
// Bords du rectangle (même logique que ton script)
eps = 1e-6;

cLeft[]  = Curve In BoundingBox {-eps, -eps, -1, eps, Ly+eps, 1};
cRight[] = Curve In BoundingBox {Lx-eps, -eps, -1, Lx+eps, Ly+eps, 1};
cTop[]   = Curve In BoundingBox {-eps, Ly-eps, -1, Lx+eps, Ly+eps, 1};
cBot[]   = Curve In BoundingBox {-eps, -eps, -1, Lx+eps, eps, 1};

Physical Curve(11) = {cLeft[]};   // SIGMA_PLUS_L
Physical Curve(12) = {cRight[]};  // SIGMA_MINUS_L
Physical Curve(13) = {cTop[]};    // BORD_HAUT
Physical Curve(14) = {cBot[]};    // BORD_BAS

Mesh.RecombineAll = 0;
