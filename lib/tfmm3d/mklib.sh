gfortran -c -O3 tfmm3dwrap.f tfmm3d.f l3dzero.f\
	treeplot.f prini.f l3dterms.f laprouts3d.f\
	Near_Interaction_code.f90 \
	l3dmpmpfinal4.f l3dloclocfinal4.f l3dmplocfinal4.f\
	prinm.f yrecursion.f ftophys.f phystof.f\
	legeexps.f d3hptree.f rotgen.f numthetafour.f numthetahalf2.f\
	lapweights.f pwrouts.f hkrand.f dlaran.f rotviarecur3.f rotgen2.f\
	l3dgqbxauxrouts.f lwtsexp_sep2.f 
#
ar rc tfmm3dlib.a tfmm3dwrap.o tfmm3d.o l3dzero.o\
	treeplot.o prini.o l3dterms.o laprouts3d.o\
	Near_Interaction_code.o \
	l3dmpmpfinal4.o l3dloclocfinal4.o l3dmplocfinal4.o\
	prinm.o yrecursion.o ftophys.o phystof.o\
	legeexps.o d3hptree.o rotgen.o numthetafour.o numthetahalf2.o\
	lapweights.o pwrouts.o hkrand.o dlaran.o rotviarecur3.o rotgen2.o\
	l3dgqbxauxrouts.o lwtsexp_sep2.o 
#
