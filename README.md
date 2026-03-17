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

## Architecture du Projet

* `src/` : Fichiers sources C++ (`main.cpp`, algorithmes FEM, maillage, LSM).
* `include/` : Fichiers d'en-tête contenant les définitions de classes et l'algèbre linéaire (`math.hpp`).
* `data/` : Dossier contenant les géométries `.geo` et maillages `.msh` de test.
* `test/` : Scripts de validation et de tests de sous-modules.
