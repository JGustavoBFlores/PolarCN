.POSIX:

# Parameters:
#
#       FC          -   Fortran compiler.
#       SRC         -   Source file(s).
#       TARGET      -   Build target name (executable).

FC     = gfortran
SRC    = program.f
TARGET = simulator

.PHONY: all clean run

all: $(TARGET)

$(TARGET):
	$(FC) -o $(TARGET) *.f -llapack -lblas
# $(FC) -o $(TARGET) $(SRC) -llapack -lblas
clean:
	rm $(TARGET)

run:
	./$(TARGET)
