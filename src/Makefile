include make.flags
include modules.flags
all: ../bing/ring
../bing/ring:  basout.o  dlamch.o  greatr.o  hpdel.o  hpins.o  lsame.o twodq.o xerrwv.o grid.o integu.o ringrav.o param.o lib_array.o readriv.o
	$(FC) $(FFLAGS)  basout.o  dlamch.o  greatr.o  hpdel.o  hpins.o  lsame.o twodq.o xerrwv.o grid.o integu.o ringrav.o param.o lib_array.o readriv.o -o ../bin/ring.x
depend:
	mkdep90  -L /opt/intel/composerxe-2011.4.191/compiler/include/intel64/  -L/opt/intel/composerxe-2011.4.191/mkl/include/intel64/lp64/ -L/opt/intel/mkl/lib/intel64/  *.f90  >make.dep
	mkmoddep -f .compiler ./  >modules.flags

clean:
	rm -f  *.o *.mod *% *~ core 
# Dependencies
include make.dep

	
