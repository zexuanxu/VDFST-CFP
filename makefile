PROGRAM = vdfst_cfp

# Define the Fortran compile rand its flags
F90= gfortran
F90FLAGS=

# Define the libraries
SYSLIBS=  -O2 -lc  -llapack -lblas -ltmglib 
USRLIB = 

# Define all object files which make up Modtools
OBJECTS= \
		global.o \
        start.o \
        utl.o \
		gwfmedia.o \
		gwfcon.o \
		transmedia.o \
		transcon.o \
        budget.o \
        vdfst_cfp.o \

# Define Task Function Program Modtools
all: $(PROGRAM)

# Define what Modtools is
$(PROGRAM): $(OBJECTS)
	-$(F90) $(F90FLAGS) -o $(PROGRAM) $(OBJECTS) $(USRLIB) $(SYSLIBS)

# Pth_Object codes of Modtools
%.o: %.f90
	$(F90) $(F90FLAGS) -c $<

# use make clean to remove .o and .mod files
clean:
	rm -f *.o *.mod


#  end
