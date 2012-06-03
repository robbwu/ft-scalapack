#C=mpif90
#C=mpif90
C=mpiifort
#L=-lblas -lscalapack -L /u/sa/br/pwu/myprogs/lib
#L=/opt/intel/Compiler/11.1/069/mkl/lib/em64t/libmkl_scalapack_lp64.a \
#/opt/intel/Compiler/11.1/069/mkl/lib/em64t/libmkl_intel_lp64.a \
#/opt/intel/Compiler/11.1/069/mkl/lib/em64t/libmkl_blacs_openmpi_lp64.a \
#/opt/intel/Compiler/11.1/069/mkl/lib/em64t/libmkl_sequential.a \
#/opt/intel/Compiler/11.1/069/mkl/lib/em64t/libmkl_core.a \
#/opt/intel/Compiler/11.1/069/mkl/lib/em64t/libmkl_sequential.a \
#/opt/intel/Compiler/11.1/069/mkl/lib/em64t/libmkl_core.a \
#-lpthread -i_dynamic \
#O=-O2
L=-L$(MKLROOT)/lib/intel64 -lmkl_scalapack_ilp64 -lmkl_intel_ilp64 -lmkl_sequential -lmkl_core -lmkl_blacs_intelmpi_ilp64 -lpthread -lm
I=-i8  -I$(MKLROOT)/include
O=
D=-g -traceback -check bounds -gen-interfaces -warn interfaces


dsybc2_t: dsybc2.f90 dsybc2_t.f90
	$C dsybc2_t.f90 dsybc2.f90 $I $L $O $D -o dsybc2_t
