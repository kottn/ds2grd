.SUFFIXES :
.SUFFIXES : .o .f90

PROG   = ../ds2grd.exe


#--- GNU
FC       = gfortran
#FFLAGS   =  -g -Wall -Wextra -fimplicit-none            \
            -static -Warray-temporaries                 \
            -ffpe-trap=invalid,zero,overflow,underflow  \
            -finit-real=nan -fbounds-check              \
            -ffree-line-length-0 -O0 -fbacktrace
#FFLAGS   = -O0
#FFLAGS   = -O2 -Wuninitialized
#FFLAGS   = -O3
#FFLAGS   = -Ofast -ffree-line-length-0
#FFLAGS  += -fopenmp

#--- Intel
#FC       = ifort
#FFLAGS   = -check all -warn all -CB -traceback
#FFLAGS   = -O3
#FFLAGS   = -fast
#FFLAGS  += -openmp

#--- Intel-MPI
#FC       = mpif90
#FFLAGS   = -check all -warn all -CB -traceback
#FFLAGS   = -O3 -axSSE4.1 -xSSE4.1
#FFLAGS   = -O3 -axAVX -xAVX
#FFLAGS   = -O3
#FFLAGS   = -fast


SRC     = \
ds2grd.f90






OBJS = $(SRC:%.f90=%.o)

.f90.o:
	$(FC) $(FFLAGS) $(FFLIB) -c $< -o $@

$(PROG) :$(OBJS)
	$(FC) $(OBJS) $(FFLAGS) $(FFLIB) -o $(PROG)

clean:
	@rm -f *.o *.mod *~ *__genmod.f90 $(PROG)

