# Imagerie Ultrasonore par Méthode d'Échantillonnage Linéaire (LSM) et Éléments Finis (FEM)

Ce projet implémente une chaîne numérique complète en C++ pour la détection et la reconstruction de défauts (imagerie ultrasonore). Il couple un solveur direct par Méthode des Éléments Finis (FEM) pour simuler la propagation des ondes (équation d'Helmholtz) avec la Méthode d'Échantillonnage Linéaire (LSM) pour imager les défauts internes d'un milieu.

## Fonctionnalités Principales

* **Maillage & Géométrie** : 
  * Parseur natif pour les fichiers maillage Gmsh (`.msh` v2.2 ASCII).
  * Enrichissement automatique des éléments P1 vers des éléments triangulaires quadratiques (P2).
  * Optimisation de la numérotation des nœuds via l'algorithme **Reverse Cuthill-McKee (RCM)** pour réduire la largeur de bande des matrices.
* **Solveur FEM (Éléments Finis)** :
  * Assemblage des matrices de rigidité et de masse (quadrature de Gauss).
  * Stockage matriciel optimisé en format **Profil / Skyline** (`ProfileMatrix`).
  * Solveur direct intégré utilisant une **factorisation LDL*** pour les matrices symétriques complexes.
* **Imagerie LSM** :
  * Calcul de la méthode d'échantillonnage linéaire multi-fréquentielle.
  * Ajout de bruit synthétique paramétrable pour évaluer la robustesse de l'algorithme.
* **Interopérabilité** : Export des résultats et maillages vers Matlab/Octave et format texte pour Python.

## Prérequis

Pour compiler et exécuter ce projet, vous aurez besoin de :
* Un compilateur C++ supportant le standard **C++17** (GCC, Clang, MSVC).
* **CMake** (version 3.10 ou supérieure) ou **Make**.
* *Optionnel :* [Gmsh](https://gmsh.info/) pour créer/modifier les maillages (`.geo` vers `.msh`).
* *Optionnel :* Matlab, Octave ou Python pour visualiser les champs d'ondes et l'image reconstruite.

## Compilation

Le projet peut être compilé avec le Makefile fourni ou via CMake (recommandé).

### Option 1 : Avec CMake (Recommandé)
```bash
mkdir build
cd build
cmake ..
make -j4
```

### Option 2 : Avec GNU Make
Directement à la racine du projet :
```bash
make
```

## Utilisation

L'exécutable généré (`us_imaging.x`) prend 3 arguments en ligne de commande :

```bash
./us_imaging.x <chemin_vers_maillage.msh> <pourcentage_de_bruit> <nombre_de_frequences>
```

**Exemple d'exécution :**
```bash
./us_imaging.x data/test_ultrasound_defaut_centre.msh 0.05 3
```
*Cet exemple lance l'imagerie sur le maillage `test_ultrasound_defaut_centre.msh`, ajoute 5% de bruit au signal, et moyenne les résultats sur 3 fréquences.*

### Sorties générées
Lors de son exécution, le programme génère par défaut :
* `mesh_out.m` : Un script Matlab contenant les données du maillage pour vérification.
* `image_lsm.txt` : La grille de données contenant la valeur de l'indicateur LSM pour chaque point testé de la zone d'imagerie.

## Génération de données FEM pour un PINN

Le générateur résout un même problème pour plusieurs fréquences et modes incidents. Il exporte les mesures sur les deux ports, le champ sur tous les degrés de liberté P2 et, en option, une interpolation P2 sur une grille régulière.

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build --target generate_pinn_data.x -j4

./build/generate_pinn_data.x \
  --mesh data/test_us_barrehalf_centree.msh \
  --output-dir pinn_data \
  --dataset barrehalf_contrast20percent \
  --c0 340 --speed-ratio 0.8 \
  --frequencies 500,600,700 \
  --modes 0,1 \
  --grid 201x61
```

`--speed-ratio` désigne précisément `c_defaut / c0`; une valeur de `0.8` donne donc `c_defaut = 272 m/s`. Les tags par défaut sont 2 pour le défaut, 11 pour le port gauche et 12 pour le port droit. Ils peuvent être changés avec `--defect-tag`, `--left-tag` et `--right-tag`.

Pour chaque mode, les fichiers de bord suivent directement la convention de `pinn_waveguide_multi_mode.py` :

* `pinn_boundary_left_<dataset>_mode<N>.csv`
* `pinn_boundary_right_<dataset>_mode<N>.csv`

Les autres sorties sont :

* `fem_field_*` : solution complexe aux degrés de liberté P2 ;
* `fem_grid_*` : solution interpolée par les six fonctions de forme P2, avec célérité et tag de région ;
* `fem_evaluation_grid_*` : grille fixe, six indices nodaux et six poids de Lagrange P2 par point ;
* `fem_mesh_nodes_*` et `fem_mesh_elements_*` : coordonnées et connectivité P2, indices à partir de zéro ;
* `fem_metadata_*` : paramètres et dimensions du jeu de données.

Le couple `fem_field_*` + `fem_evaluation_grid_*` est la représentation de référence. Pour un point `p` de la grille, la reconstruction est exactement

```text
U_h[p] = phi0[p] * U[node0[p]] + ... + phi5[p] * U[node5[p]].
```

Le fichier `fem_grid_*` contient déjà ce calcul pour un usage immédiat. Il peut être régénéré en Python à partir des coefficients nodaux :

```bash
python3 tools/reconstruct_fem_p2.py \
  --grid-map pinn_data/fem_evaluation_grid_barrehalf_contrast20percent.csv \
  --field pinn_data/fem_field_barrehalf_contrast20percent_mode0.csv \
  --frequency 600 --mode 0 \
  --output reconstructed_fem_mode0_600Hz.csv
```

### Comparaison PINN–FEM

Évaluer le PINN aux coordonnées `x_norm,y_norm` de `fem_evaluation_grid_*`, puis exporter les prédictions avec les coordonnées physiques `x,y` et les colonnes `Re_U,Im_U`. Les colonnes `f` et `mode` sont recommandées pour un fichier multi-cas. Ainsi, FEM et PINN sont comparés point par point sur une grille strictement identique.

```bash
python3 tools/compare_pinn_fem.py \
  --fem pinn_data/fem_grid_barrehalf_contrast20percent_mode0.csv \
  --pinn predictions/pinn_mode0.csv \
  --frequency 600 --mode 0 \
  --output comparison_mode0_600Hz.png \
  --metrics-json comparison_mode0_600Hz.json
```

Par défaut, le script refuse des grilles différentes afin d'éviter d'ajouter une erreur d'interpolation à la comparaison. L'option `--interpolate-pinn` permet explicitement une interpolation linéaire si elle est souhaitée. Comme le PINN courant normalise chaque champ par `U_norm`, il faut soit réexporter `U_pred * U_norm`, soit passer cette valeur avec `--pinn-scale`. `--normalization max` permet aussi une comparaison de forme indépendante de l'amplitude, mais ne mesure plus l'erreur physique absolue.

## Architecture du Projet

* `src/` : Fichiers sources C++ (`main.cpp`, algorithmes FEM, maillage, LSM).
* `include/` : Fichiers d'en-tête contenant les définitions de classes et l'algèbre linéaire (`math.hpp`).
* `data/` : Dossier contenant les géométries `.geo` et maillages `.msh` de test.
* `test/` : Scripts de validation et de tests de sous-modules.
