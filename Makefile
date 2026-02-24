CXX = g++
CXXFLAGS = -Iinclude -Wall -Wextra -std=c++17 -O3
LDFLAGS = 

# Liste des fichiers sources (ajoute src/linear_sampling.cpp si tu en as un)
SRCS = src/main_temp.cpp src/fem.cpp src/mesh.cpp src/linear_sampling.cpp
# Génération automatique de la liste des objets (.o)
OBJS = $(SRCS:.cpp=.o)

OBJ = us_imaging.x

# Cible par défaut
all: $(OBJ) supp

# Lien final
$(OBJ): $(OBJS)
	$(CXX) $(CXXFLAGS) $(OBJS) -o $(OBJ) $(LDFLAGS)

# Règle pour compiler chaque .cpp en .o
%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

supp:
	rm -f $(OBJS)

# Nettoyage
clean:
	rm -f $(OBJS) $(OBJ)

.PHONY: all clean