
#HOST = osx-gcc
#HOST = osx-gcc-openmp
#HOST = osx-intel
#HOST = osx-intel-openmp
#HOST = linux-gfortran
#HOST = linux-gfortran-openmp
#HOST = amd-gfortran
#HOST = amd-gfortran-openmp
#HOST = linux-intel
#HOST = linux-intel-openmp


PROJECT = int2


ifeq ($(HOST),osx-gcc)
  FC = gfortran-8
  FFLAGS = -O2 -w
  FLINK = gfortran-8 -w -o $(PROJECT) -framework accelerate
endif

ifeq ($(HOST),osx-gcc-openmp)
  FC = gfortran 
  FFLAGS = -O2 -fopenmp -w 
  FLINK = gfortran -fopenmp -w \
    -Wl,-stack_size,0x80000000 -o $(PROJECT) -framework accelerate
  export OMP_NUM_THREADS = 4
  export OMP_STACKSIZE = 2048M
endif


ifeq ($(HOST),osx-intel-openmp)
  FC = ifort
  FFLAGS = -i8 -r8 -O2 -w -qopenmp
  FLINK = ifort -i8 -r8 -w -mkl=parallel -qopenmp \
       -Wl,-stack_size,0x40000000 -o $(PROJECT)
  export OMP_NUM_THREADS=4
  export OMP_STACKSIZE=2048M
endif

ifeq ($(HOST),osx-intel)
  FC = ifort
  FFLAGS = -i8 -r8 -O2 -w
  FLINK = ifort -i8 -r8 -w -mkl -o $(PROJECT)
endif


ifeq ($(HOST),linux-intel-openmp)
  FC = ifort
  FFLAGS = -i8 -r8 -O2 -w -qopenmp
  FLINK = ifort -i8 -r8 -w -mkl=parallel -qopenmp -o $(PROJECT)
endif


ifeq ($(HOST),linux-gcc)
   FC = gfortran
   FFLAGS = -O2 -g -w -fdefault-integer-8 -finteger-4-integer-8
   FLINK = gfortran -g -w -fdefault-integer-8 -finteger-4-integer-8 \
       -o $(PROJECT) -llapack -lblas
endif

ifeq ($(HOST),linux-gcc-openmp)
   FC = gfortran 
   FFLAGS = -O2  -w --openmp -fdefault-integer-8 -finteger-4-integer-8
   FLINK = gfortran --openmp -fdefault-integer-8 -finteger-4-integer-8 \
       -w -o $(PROJECT) -llapack -lblas
endif




.PHONY: all clean list

TFMM3D = lib/tfmm3dlr
SRC = src
EXM = examples

MOD_SOURCES = $(SRC)/Mod_TreeLRD.f90 \
 $(SRC)/ModType_Smooth_Surface.f90 \
 $(SRC)/Mod_Fast_Sigma.f90 \
 $(SRC)/chebtarggridrouts.f90 \
 $(SRC)/Mod_Plot_Tools_sigma.f90 \
 $(SRC)/Mod_Feval.f90 \
 $(SRC)/Mod_Smooth_Surface.f90

SOURCES =  $(EXM)/test_surfsmooth.f90 \
 $(SRC)/koornexps.f90 \
 $(SRC)/ortho2eva.f \
 $(SRC)/ortho2eva_new.f90 \
 $(SRC)/ortho2exps.f \
 $(SRC)/orthom.f \
 $(SRC)/cisurf_loadmsh.f90 \
 $(SRC)/cisurf_skeleton.f90 \
 $(SRC)/cisurf_plottools.f90 \
 $(SRC)/cisurf_tritools.f90 \
 $(SRC)/lapack_wrap.f90 \
 $(SRC)/pplot2.f \
 $(TFMM3D)/tfmm3dlr_expout.f \
 $(TFMM3D)/tfmm3dlrwrap_expout.f \
 $(TFMM3D)/l3dzero.f \
 $(TFMM3D)/treeplot.f \
 $(TFMM3D)/prini.f \
 $(TFMM3D)/l3dterms.f \
 $(TFMM3D)/laprouts3d.f \
 $(TFMM3D)/nearfield.f90 \
 $(TFMM3D)/l3dmpmpfinal4.f \
 $(TFMM3D)/l3dloclocfinal4.f \
 $(TFMM3D)/l3dmplocfinal4.f \
 $(TFMM3D)/prinm.f \
 $(TFMM3D)/yrecursion.f \
 $(TFMM3D)/ftophys.f \
 $(TFMM3D)/phystof2.f \
 $(TFMM3D)/legeexps.f \
 $(TFMM3D)/d3hplratree.f \
 $(TFMM3D)/rotgen.f \
 $(TFMM3D)/numthetafour.f \
 $(TFMM3D)/numthetahalf2.f \
 $(TFMM3D)/lapweights.f \
 $(TFMM3D)/pwrouts2.f \
 $(TFMM3D)/hkrand.f \
 $(TFMM3D)/dlaran.f \
 $(TFMM3D)/rotviarecur3.f \
 $(TFMM3D)/rotgen2.f \
 $(TFMM3D)/l3dgqbxauxrouts2.f \
 $(TFMM3D)/lwtsexp_sep2.f \
 $(TFMM3D)/chebexps.f



#TFMM3DLIB = ../lib/tfmm3dlr/tfmm3dlib.a

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



