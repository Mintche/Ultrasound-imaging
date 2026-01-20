CC = g++
FLAGS = -I ../include

# Cible par défaut
all: test_run

# Compilation
test_run: test_math.cpp ../include/math.hpp
	$(CC) $(FLAGS) test_math.cpp -o test_math

# Nettoyage
clean:
	rm -f test_run