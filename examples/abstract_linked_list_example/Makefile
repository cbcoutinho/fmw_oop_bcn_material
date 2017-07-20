FC=gfortran
FFLAGS=-g -O0 
INC=-I.

EXECS=list_program

all: ${EXECS}

%.o: %.f90 
	$(FC) -c -o $@ $< $(FFLAGS) $(INC)

list_program.o: int_list.o 
list.o: node_list.o
int_list.o: list.o

OBJ=list.o node_list.o int_list.o list_program.o

list_program: $(OBJ)
	$(FC) -o $@ $^ $(FFLAGS) $(INC) 

clean:
	rm -f *.o *.mod *~ core $(EXECS)
