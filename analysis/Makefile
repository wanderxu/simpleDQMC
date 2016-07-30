# first delete all default suffixes
.SUFFIXES:
# then add your suffixes
.SUFFIXES: .f .f90 .o

FC= gfortran -O3
FLAGS= -c  -w
OBJS= m_variance.o cov_tau.o Files.o
LIBS= ../lib/libMM.a

default: all

all: cov_tau

cov_tau: $(OBJS)
	$(FC) $(OBJS) -o cov_tau.out $(LFLAGS) $(LIBS)
	
clean: 	
	rm *.mod *.o *.out

.f.o .f90.o:
	$(FC) $(SUFFIX) $(FLAGS) $<
