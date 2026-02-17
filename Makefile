CC = g++
CXXFLAGS = -Iinclude -Wall -std=c++17 -O3 -fopenmp

SRC = src/main.cpp
OBJ = us_imaging.x

# Cible par défaut
all: $(OBJ)

# Compilation
$(OBJ): $(SRC) include/fem.hpp include/math.hpp include/mesh.hpp
	$(CC) $(CXXFLAGS) $(SRC) -o $(OBJ)

# Nettoyage
clean:
	rm -f $(OBJ)