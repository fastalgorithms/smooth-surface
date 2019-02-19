
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
  FC = gfortran-8 -c -w
  FFLAGS = -O2
  FLINK = gfortran-8 -w -o $(PROJECT) -framework accelerate
endif

ifeq ($(HOST),osx-gcc-openmp)
  FC = gfortran -c -fopenmp -w
  FFLAGS = -O2
  FLINK = gfortran -fopenmp -w \
    -Wl,-stack_size,0x40000000 -o $(PROJECT) -framework accelerate
  export OMP_NUM_THREADS = 4
  export OMP_STACKSIZE = 2048M
endif



# ifeq ($(HOST),osx-intel)
#   FC = ifort -c -w 
#   FFLAGS = -O2
#   FLINK = ifort -mkl -o $(PROJECT)
# endif

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

ifeq ($(HOST),linux-gfortran-openmp)
   FC = gfortran -c -w --openmp
   FFLAGS = -O3  
   FLINK = gfortran --openmp -w -o $(PROJECT) -llapack -lblas
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

TFMM3D = lib/tfmm3dlr
SRC = src
EXM = examples

MOD_SOURCES = $(SRC)/Mod_Plot_Tools.f90 \
 $(SRC)/Mod_Tri_Tools.f90 \
 $(SRC)/Mod_GaussTri.f90 \
 $(SRC)/Mod_TreeLRD.f90 \
 $(SRC)/ModType_Smooth_Surface.f90 \
 $(SRC)/Mod_InOut.f90 \
 $(SRC)/Mod_Fast_Sigma.f90 \
 $(SRC)/Mod_Plot_Tools_sigma.f90 \
 $(SRC)/chebtarggridrouts.f90 \
 $(SRC)/Mod_Feval.f90 \
 $(SRC)/Mod_Smooth_Surface.f90 \

SOURCES =  $(EXM)/test_surfsmooth.f90 \
 src/koornexps.f90 \
 src/ortho2eva.f \
 src/ortho2eva_new.f90 \
 src/ortho2exps.f \
 src/orthom.f \
 $(SRC)/lapack_wrap.f90 \
 $(TFMM3D)/tfmm3dlr_expout.f \
 $(TFMM3D)/tfmm3dlrwrap_expout.f \
 $(TFMM3D)/l3dzero.f \
 $(TFMM3D)/treeplot.f \
 $(TFMM3D)/prini.f \
 $(TFMM3D)/l3dterms.f \
 $(TFMM3D)/laprouts3d.f \
 $(TFMM3D)/Near_Interaction_code.f90 \
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
 $(TFMM3D)/chebexps.f \
 $(SRC)/pplot2.f \



#TFMM3DLIB = ../lib/tfmm3dlr/tfmm3dlib.a

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



