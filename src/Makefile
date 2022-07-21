#####################
#  Scorpio Makefile #
#  Hsiang-Hsu Wang  #
#####################
EXEC=Scorpio
F90=gfortran
OBJECT=gridModule.o evolveAmbipolarDiffusion.o calcSG.o calcDT.o initSGWindows2D.o initSGWindows3D.o sgPlan.o DTPlan.o sgKernel.o inout.o initMPIWindows1D.o initMPIWindows2D.o initMPIWindows3D.o exchgBdryMPI.o testSuiteMPI.o limiterModule.o setCoordinates.o init1d.o init2d.o init3d.o setBdry1D.o setBdry2D.o setBdry3D.o dt1D.o dt2D.o dt3D.o source2D.o main.o riemannSolverModule.o rk2.o 

HOMEDIR=/home/hinny/Scorpio
ScorpioHOME=$(HOMEDIR)/Scorpio_v1.5/
HDF5HOME=/opt/hdf5-1.10.0/mpich/gcc-5.4.0
FFTWHOME=/opt/fftw-3.3.6/mpich/gcc-5.4.0
MPIHOME=/opt/mpich-3.2/gcc-5.4.0

LIB = -L$(FFTWHOME)/lib -lfftw3_mpi -lfftw3 -lm -L$(MPIHOME)/lib -lmpi -lz -L$(HDF5HOME)/lib -lhdf5_fortran -lhdf5hl_fortran -lhdf5_hl -lhdf5 -fPIC
INCL= -I$(HDF5HOME)/include -I$(MPIHOME)/include -I$(FFTWHOME)/include
VPATH==./:$(ScorpioHOME)
%.o:%.f03
#	$(F90) -O2 -Wall -fcheck=all -c $^ -o $@ $(INCL)
	$(F90) -O3 -c $^ -o $@ $(INCL)
Scorpio: $(OBJECT)
	$(F90) -o $(EXEC) $(OBJECT) $(LIB)
clean:
	rm *.o *.mod
cleaner:
	rm $(EXEC)
clean-data:
	rm *.h5
