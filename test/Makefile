CC = g++
FLAGS = -I ../include

# Cible par défaut
all: test_run test_run_mesh

# Compilation
test_run: test_math.cpp ../include/math.hpp
	$(CC) $(FLAGS) test_math.cpp -o test_math.x

test_run_mesh: test_mesh.cpp ../include/mesh.hpp
	$(CC) $(FLAGS) test_mesh.cpp -o test_mesh.x

test_run_lin_sampl: test_linear_sampling.cpp ../include/linear_sampling.hpp
	$(CC) $(FLAGS) test_linear_sampling.cpp -o test_linear_sampling.x

# Nettoyage
clean:
	rm -f *.x