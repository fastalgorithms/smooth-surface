
#HOST = osx-gcc
#HOST = osx-gcc-openmp
#HOST = osx-intel
#HOST = osx-intel-openmp
HOST = linux-gfortran
#HOST = amd-gfortran
#HOST = amd-gfortran-openmp
#HOST = linux-intel
#HOST = linux-intel-openmp


PROJECT = int2


ifeq ($(HOST),osx-gcc)
  FC = gfortran-8 -c -w
  FFLAGS = -O2
  FLINK = gfortran-8 -w -o $(PROJECT) -framework accelerate
endif

ifeq ($(HOST),osx-gcc-openmp)
  FC = gfortran-8 -c -fopenmp -w
  FFLAGS = -O2
  FLINK = gfortran-8 -fopenmp -w \
    -Wl,-stack_size,0x40000000 -o $(PROJECT) -framework accelerate
  export OMP_NUM_THREADS = 4
  export OMP_STACKSIZE = 2048M
endif



# ifeq ($(HOST),osx-intel)
#   FC = ifort -c -w 
#   FFLAGS = -O2
#   FLINK = ifort -mkl -o $(PROJECT)
# endif

# ifeq ($(HOST),osx-intel-openmp)
#   FC = ifort -c -w -qopenmp
#   FFLAGS = -O2
#   FLINK = ifort -w -mkl=parallel -qopenmp \
#     -Wl,-stack_size,0x40000000 -o $(PROJECT)
#   export OMP_NUM_THREADS=4
#   export OMP_STACKSIZE=2048M
# endif

ifeq ($(HOST),linux-gfortran)
   FC = gfortran -c -w
   FFLAGS = -O2 -mcmodel=medium 
   FLINK = gfortran -w -mcmodel=medium -o $(PROJECT) -llapack -lblas
endif

# ifeq ($(HOST),amd-gfortran)
#   FC = gfortran -c -w
#   FFLAGS = -O2 -mcmodel=medium 
#   #FLINK = gfortran -w -mcmodel=medium -o $(PROJECT) /usr/lib/liblapack.so.3 /usr/lib/libblas.so.3
#   FLINK = gfortran -w -mcmodel=medium -o $(PROJECT) -llapack -lblas
# endif

# ifeq ($(HOST),amd-gfortran-openmp)
#   FC = gfortran -c -w
#   FFLAGS = -O2 -mcmodel=medium -fopenmp
#   FLINK = gfortran -w -mcmodel=medium -fopenmp -o $(PROJECT) -llapack -lblas
#   export OMP_NUM_THREADS=8
#   export OMP_STACKSIZE=1024M
# endif

# ifeq ($(HOST),linux-intel)
#   FC = ifort -c -w
#   FFLAGS = -O2 -mcmodel=medium
#   FLINK = ifort -w -mcmodel=medium -o $(PROJECT) -mkl
# endif

# ifeq ($(HOST),linux-intel-openmp)
#   FC = ifort -c -w -qopenmp
#   FFLAGS = -O2 -mcmodel=medium
#   FLINK = ifort -w -mcmodel=medium -mkl=parallel -qopenmp -o $(PROJECT)
#   export OMP_NUM_THREADS=4
#   export OMP_STACKSIZE=1024M
# endif



.PHONY: all clean list

TFMM3D = ../lib/tfmm3d
SRC = ../src

MOD_SOURCES = $(SRC)/chebtarggridrouts_withoutfmm.f90

SOURCES =  testchebtargrouts.f90 \
 $(SRC)/lapack_wrap.f90 \
 $(TFMM3D)/prini.f \
 $(TFMM3D)/chebexps.f \
 $(TFMM3D)/d3hplratree.f \
 $(TFMM3D)/hkrand.f \
 $(TFMM3D)/dlaran.f \
 $(SRC)/pplot2.f \



TFMM3DLIB = ../lib/tfmm3d/tfmm3dlib.a

MOD_OBJECTS = $(patsubst %.f,%.o,$(patsubst %.f90,%.o,$(MOD_SOURCES)))
OBJECTS = $(patsubst %.f,%.o,$(patsubst %.f90,%.o,$(SOURCES)))

#
# use only the file part of the filename, then manually specify
# the build location
#

%.o : %.f
	$(FC) $(FFLAGS) $< -o $@

%.o : %.f90
	$(FC) $(FFLAGS) $< -o $@


all: mods objects
	$(FLINK) $(OBJECTS) $(MOD_OBJECTS)
	./$(PROJECT)

mods: $(MOD_OBJECTS) 
objects: $(OBJECTS) 

clean:
	rm -f $(MOD_OBJECTS)
	rm -f prefunrouts_withoutfmm.mod
	rm -f $(OBJECTS)
	rm -f $(PROJECT)

list: $(SOURCES) $(MOD_SOURCES)
	$(warning Requires:  $^)



