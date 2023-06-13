# files are IISPH_solver.hpp, Vector.hpp, main.cpp

# Compiler
CC = g++

# Compiler flags
CFLAGS = -Wall -Wextra -Wpedantic -std=c++11 -O3 -fopenmp

# Linker flags
LDFLAGS =

# Executable name
EXEC = 3DIISPH

# Source files
SRC = main.cpp

# Object files
OBJ = $(SRC:.cpp=.o)

# Header files
HDR = IISPH_solver.hpp Vector.hpp particleInitialization.hpp

# Compile
$(EXEC): $(OBJ)
	$(CC) $(CFLAGS) $(OBJ) -o $(EXEC) $(LDFLAGS)

# Compile source files
%.o: %.cpp $(HDR)
	$(CC) $(CFLAGS) -c $< -o $@

# Clean
clean:
	rm -f $(OBJ) $(EXEC)

# Run
run:
	./$(EXEC)
