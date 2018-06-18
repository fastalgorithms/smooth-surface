
#HOST = osx-gcc-7
#HOST = osx-gcc-7-openmp
#HOST = osx-intel
#HOST = osx-intel-openmp
#HOST = linux-gfortran
#HOST = amd-gfortran
#HOST = amd-gfortran-openmp
HOST = linux-intel
#HOST = linux-intel-openmp
PROJECT = int2

ifeq ($(HOST),osx-gcc-7)
  FC = gfortran-7 -c -w
  FFLAGS = -O2
  FLINK = gfortran-7 -w -o $(PROJECT) -framework accelerate
endif

ifeq ($(HOST),osx-gcc-7-openmp)
  FC = gfortran-7 -c -fopenmp -w
  FFLAGS = -O2
  FLINK = gfortran-7 -fopenmp -w \
    -Wl,-stack_size,0x40000000 -o $(PROJECT) -framework accelerate
  export OMP_NUM_THREADS = 4
  export OMP_STACKSIZE = 2048M
endif



ifeq ($(HOST),osx-intel)
  FC = ifort -c -w 
  FFLAGS = -O2
  FLINK = ifort -mkl -o $(PROJECT)
endif

ifeq ($(HOST),osx-intel-openmp)
  FC = ifort -c -w -qopenmp
  FFLAGS = -O2
  FLINK = ifort -w -mkl=parallel -qopenmp \
    -Wl,-stack_size,0x40000000 -o $(PROJECT)
  export OMP_NUM_THREADS=4
  export OMP_STACKSIZE=2048M
endif

ifeq ($(HOST),linux-gfortran)
  FC = gfortran -c -w
  FFLAGS = -O2 -mcmodel=medium 
  FLINK = gfortran -w -mcmodel=medium -o $(PROJECT) -llapack -lblas
endif

ifeq ($(HOST),amd-gfortran)
  FC = gfortran -c -w
  FFLAGS = -O2 -mcmodel=medium 
  #FLINK = gfortran -w -mcmodel=medium -o $(PROJECT) /usr/lib/liblapack.so.3 /usr/lib/libblas.so.3
  FLINK = gfortran -w -mcmodel=medium -o $(PROJECT) -llapack -lblas
endif

ifeq ($(HOST),amd-gfortran-openmp)
  FC = gfortran -c -w
  FFLAGS = -O2 -mcmodel=medium -fopenmp
  FLINK = gfortran -w -mcmodel=medium -fopenmp -o $(PROJECT) -llapack -lblas
  export OMP_NUM_THREADS=8
  export OMP_STACKSIZE=1024M
endif

ifeq ($(HOST),linux-intel)
  FC = ifort -c -w
  FFLAGS = -O2 -mcmodel=medium
  FLINK = ifort -w -mcmodel=medium -o $(PROJECT) -mkl
endif

ifeq ($(HOST),linux-intel-openmp)
  FC = ifort -c -w -qopenmp
  FFLAGS = -O2 -mcmodel=medium
  FLINK = ifort -w -mcmodel=medium -mkl=parallel -qopenmp -o $(PROJECT)
  export OMP_NUM_THREADS=4
  export OMP_STACKSIZE=1024M
endif



.PHONY: lib driver clean list

SOURCES =  tfmm3dwrap.f tfmm3d.f l3dzero.f\
	treeplot.f prini.f l3dterms.f laprouts3d.f\
	Near_Interaction_code.f90 \
	l3dmpmpfinal4.f l3dloclocfinal4.f l3dmplocfinal4.f\
	prinm.f yrecursion.f ftophys.f phystof.f\
	legeexps.f d3hptree.f rotgen.f numthetafour.f numthetahalf2.f\
	lapweights.f pwrouts.f hkrand.f dlaran.f rotviarecur3.f rotgen2.f\
	l3dgqbxauxrouts.f lwtsexp_sep2.f 

SRC_DRIVER = tfmm3dwrap_dr.f


OBJECTS = $(patsubst %.f,%.o,$(patsubst %.f90,%.o,$(SOURCES)))

#
# use only the file part of the filename, then manually specify
# the build location
#

%.o : %.f
	$(FC) $(FFLAGS) $< -o $@

%.o : %.f90
	$(FC) $(FFLAGS) $< -o $@

lib: $(OBJECTS)
	ar rc tfmm3dlib.a $(OBJECTS)

driver: $(OBJECTS)
	$(FC) $(FFLAGS) tfmm3dwrap_dr.f
	$(FLINK) tfmm3dwrap_dr.o $(OBJECTS)
	./$(PROJECT)


clean:
	rm -f $(OBJECTS)
	rm -f $(PROJECT)
	rm -f tfmm3dlib.a

list: $(SOURCES)
	$(warning Requires:  $^)



