#C=mpif90
C=mpif90
#C=mpiifort
#L=-lblas -lscalapack -L /u/sa/br/pwu/myprogs/lib
#O=-O2
L=-L$(MKLROOT)/lib/intel64 -lmkl_scalapack_lp64 -lmkl_intel_lp64 -lmkl_sequential -lmkl_core -lmkl_blacs_openmpi_lp64 -lpthread -lm
I=-I$(MKLROOT)/include
O=
#D=-g -traceback -check bounds -gen-interfaces -warn interfaces
D=-g

ft_pdpotrf.o: ft_pdpotrf.f90
	$C -c ft_pdpotrf.f90

ftchol.o: ftchol.f90
	$C -c ftchol.f90
dsybc2_t: dsybc2_t.f90 ftchol.o
	$C dsybc2_t.f90 ftchol.o $I $L $O $D -o dsybc2_t

chk1_t: chk1_t.f90 ftchol.o
	$C chk1_t.f90 ftchol.o $I $L $O $D -o chk1_t

chk2_t: chk2_t.f90 ftchol.o
	$C chk2_t.f90 ftchol.o $I $L $O $D -o chk2_t

chk3_t: chk3_t.f90 ftchol.o
	$C chk3_t.f90 ftchol.o $I $L $O $D -o chk3_t

pdpotrf_t1: pdpotrf_t1.f90 ft_pdpotrf.o ftchol.o
	$C pdpotrf_t1.f90 ft_pdpotrf.o ftchol.o $I $L $O $D -o pdpotrf_t1
