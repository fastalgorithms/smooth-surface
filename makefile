
HOST = osx-gcc
#HOST = linux-gcc
#HOST = linux-gcc-openmp
HOST = osx-gcc-openmp
#HOST = osx-intel
#HOST = osx-intel-openmp
#HOST = linux-gfortran
#HOST = linux-gfortran-openmp
#HOST = amd-gfortran
#HOST = amd-gfortran-openmp
#HOST = linux-intel
#HOST = linux-intel-openmp


PROJECT = int2

GFORT = gfortran

ifeq ($(HOST),osx-gcc)
  FC = $(GFORT)
  FFLAGS = -O2 -w -fdefault-integer-8 -finteger-4-integer-8 \
              -fdefault-double-8 -fdefault-real-8 -freal-4-real-8 -std=legacy
  FLINK = $(GFORT) -w -fdefault-integer-8 -finteger-4-integer-8 \
             -fdefault-double-8 -fdefault-real-8 -freal-4-real-8 -std=legacy \
             -o $(PROJECT) -framework accelerate
endif

ifeq ($(HOST),osx-gcc-openmp)
  FC = $(GFORT)
  FFLAGS = -O2 -fopenmp -w -fdefault-integer-8 -finteger-4-integer-8 \
              -fdefault-double-8 -fdefault-real-8 -freal-4-real-8 -std=legacy
  FLINK = $(GFORT) -w  -fopenmp -fdefault-integer-8 -finteger-4-integer-8 \
             -fdefault-double-8 -fdefault-real-8 -freal-4-real-8 -std=legacy \
             -o $(PROJECT) -framework accelerate
endif


ifeq ($(HOST),osx-intel-openmp)
  FC = ifort
  FFLAGS = -i8 -r8 -O2 -w -qopenmp
  FLINK = ifort -i8 -r8 -w -mkl=parallel -qopenmp \
       -Wl,-stack_size,0x40000000 -o $(PROJECT)
  export OMP_STACKSIZE=2048M
endif

ifeq ($(HOST),osx-intel)
  FC = ifort
  FFLAGS = -i8 -r8 -O2 -w
  FLINK = ifort -i8 -r8 -w -mkl -o $(PROJECT)
endif


ifeq ($(HOST),linux-intel-openmp)
  FC = ifort
  FFLAGS = -i8 -r8 -O2 -w -qopenmp -xHost
  FLINK = ifort -i8 -r8 -w -mkl=parallel -qopenmp -o $(PROJECT)
endif


ifeq ($(HOST),linux-gcc)
   FC = gfortran
   FFLAGS = -O2 -g -w -fdefault-integer-8 -finteger-4-integer-8 -std=legacy
   FLINK = gfortran -g -w -fdefault-integer-8 -finteger-4-integer-8 -std=legacy \
       -o $(PROJECT) -lopenblas
endif

ifeq ($(HOST),linux-gcc-openmp)
   FC = gfortran 
   FFLAGS = -O2  -w --openmp -fdefault-integer-8 -finteger-4-integer-8
   FLINK = gfortran --openmp -fdefault-integer-8 -finteger-4-integer-8 \
       -w -o $(PROJECT) -llapack -lblas
endif

FFLAGS = -march=native -O2 -w -std=legacy -fopenmp
FLINK = gfortran -march=native -O2 -w -std=legacy -fopenmp -o $(PROJECT) -framework accelerate -lfmm3dbie -L/usr/local/lib



.PHONY: all clean list

SRC = src
EXM = examples

MOD_SOURCES = $(SRC)/Mod_TreeLRD.f90 \
 $(SRC)/ModType_Smooth_Surface.f90 \
 $(SRC)/Mod_Fast_Sigma.f90 \
 $(SRC)/Mod_Plot_Tools_sigma.f90 \
 $(SRC)/Mod_Feval.f90 \
 $(SRC)/Mod_Smooth_Surface.f90

SOURCES =  $(EXM)/test_surfsmooth.f90 \
 $(SRC)/koornexps.f90 \
 $(SRC)/cisurf_loadmsh.f90 \
 $(SRC)/cisurf_skeleton.f90 \
 $(SRC)/cisurf_plottools.f90 \
 $(SRC)/cisurf_tritools.f90 \
 $(SRC)/lapack_wrap.f90 \
 $(SRC)/pplot2.f \
 $(SRC)/tfmm_setsub.f \



MOD_OBJECTS = $(patsubst %.f,%.o,$(patsubst %.f90,%.o,$(MOD_SOURCES)))
OBJECTS = $(patsubst %.f,%.o,$(patsubst %.f90,%.o,$(SOURCES)))

#
# use only the file part of the filename, then manually specify
# the build location
#

%.o : %.f
	$(FC) $(FFLAGS) -c $< -o $@

%.o : %.f90
	$(FC) $(FFLAGS) -c $< -o $@


all: mods objects
	$(FLINK) $(OBJECTS) $(MOD_OBJECTS)
	./$(PROJECT)

mods: $(MOD_OBJECTS) 
objects: $(OBJECTS) 

clean:
	rm -f $(MOD_OBJECTS)
	rm -f *.mod
	rm -f $(OBJECTS)
	rm -f $(PROJECT)

list: $(SOURCES) $(MOD_SOURCES)
	$(warning Requires:  $^)



