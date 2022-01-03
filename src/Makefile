

F90 = ifort -mkl



#modules
GW_modules = modmain.f90

# 3d-hubbard GW part
#GW_3dHubbard = Model_GW.f90
GW_3dHubbard = Model_GW.f90


# main subroutines and functions
GW_main = readinput.f90 LatticeConst.f90 parmeterR2K.f90 Grid.f90 NonInteracting.f90 Hartree.f90 Fock.f90 fourierT.f90 selfE_gw.f90 bosonic_gw.f90 chempot_wRenomalizedG_gw.f90 chempot_wRenomalizedG_hf.f90 chempot_wRenomalizedG_gwdmft.f90 PrintF.f90 CallocF.f90 ReadF.f90 convcheck.f90 GWDMFT.f90 

GW = $(GW_modules) $(GW_main) $(GW_3dHubbard)

OBJ = $(GW:.f90=.o)
#EXE = GW_3dHubbard
EXE = Model_GW

#-------------------------------------------------------------------------------
# Suffix rules
#-------------------------------------------------------------------------------
.SUFFIXES: .o .f90

3dHubbard_GW: $(OBJ)
	$(F90) $(OBJ) -o Model_GW
 
.f90.o:
	$(F90) -c $<


all: main

clean : 
	rm -rf *.o *.mod Model_GW 
