.SUFFIXES: .f90 .f

#include ../build/make.sys
include ./make.sys

base = daxpy.o dcabs1.o dgedi.o dgefa.o dgemm.o dgemv.o dger.o dscal.o dswap.o F01QCF.o F01QEF.o F01RCF.o F01REF.o F06FBF.o  \
      F06FJF.o F06FRF.o F06HBF.o F06HRF.o F06KJF.o F06QHF.o F06THF.o idamax.o izamax.o lsame.o P01ABF.o P01ABW.o P01ABY.o   \
      P01ABZ.o X02AJF.o X04AAF.o X04BAF.o xerbla.o zaxpy.o zdscal.o zgedi.o zgefa.o zgemm.o zgemv.o zgerc.o zscal.o zswap.o
sub = HermDiag.o InvMatrx.o MatrCmpr.o MatrMult.o MyRandom.o ReSyDiag.o UDVdcmps.o UDVReOrtho.o s_matrix.o
objects = $(base) $(sub)

default: all

all: build-base build-sub build-lib

build-base: $(base)
build-sub: $(sub)
build-lib: $(objects)
	$(ARCHIVER) libMM.a $(objects)

.f.o:
	$(F90) $(FFLAGS) $*.f

.f90.o:
	$(F90) $(FFLAGS) $*.f90

clean:
	rm -f *.mod
	rm -f *.o
	rm -f libMM.a
