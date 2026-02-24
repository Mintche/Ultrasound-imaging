# 1. Version minimale de CMake requise
cmake_minimum_required(VERSION 3.10)

# 2. Nom du projet et langage
project(UltrasoundImaging CXX)

# 3. Spécifier le standard C++ (équivalent de -std=c++17)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)

# 4. Indiquer où se trouvent les fichiers d'en-tête (.hpp)
# (équivalent de -Iinclude)
include_directories(include)

# 5. Lister manuellement tes fichiers sources
# C'est une meilleure pratique que d'utiliser GLOB pour que CMake
# détecte immédiatement l'ajout d'un nouveau fichier.
set(SOURCES
    src/main.cpp
    src/fem.cpp
    src/mesh.cpp
    src/linear_sampling.cpp
)

# 6. Déclarer l'exécutable final à générer
add_executable(us_imaging.x ${SOURCES})

# 7. Ajouter les flags de compilation (équivalent de -Wall -Wextra -O3)
target_compile_options(us_imaging.x PRIVATE -Wall -Wextra -O3)