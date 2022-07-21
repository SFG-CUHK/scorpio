program main
use mpi
use testSuiteMPI
implicit none
integer::np,rank,ierr


call MPI_INIT(ierr)
call MPI_COMM_SIZE(MPI_COMM_WORLD,np,ierr)

call setMPI(np)  
!!! zwg: set nprocs=np  where nprocs is the total number of processes
!!! zwg: this subroutine is in gridmodule.f30

call setTestOnOff(.true.) 
!!!! zwg : set logical variable testOnOff which is initialized to be .false. 
!!!! and then resetted by this sunroutine [subroutine setTestOnOff(.true.) ]

if(testOnOff) then
   
!!!! zwg: all these tests below are seted in the file of testSuitMPI.f90
  !call isoShockTube1D(gridID=1)
  !call SodShockTube1D(gridID=2)
  !call WCShockTube1D(gridID=3)
  !call isoShockTube2D(gridID=4)
  !call SodShockTube2D(gridID=5)
  !call WCShockTube2D(gridID=6)

  !call doubleMachReflection(gridID=7)

  !call strongRarefactionTest(gridID=8)
  !call ShuOsherShockTubeTest(gridID=9)
  !call BrioWuShockTube1D(gridID=10)
  !call RJ2aShockTube1D(gridID=11)
  !call BrioWuShockTube2D(gridID=12)

  !call OrsagTangVortex(gridID=13) 
  
  !call sphericalBlastWaveMHD2D(gridID=14)

  !call BFieldLoopTest2D(gridID=15)

  !call rotorMHD2D(gridID=16)

  !call currentSheetMHD2D(gridID=17)

  !call CPAlvenWave(gridID=18)
  !call linearSlowWaveTestMHD2D(gridID=19)
  !call linearFastWaveTestMHD2D(gridID=20)

  !call selfgravityTest2D(gridID=21)

  !call CShockTest1D(gridIDn=22,gridIDi=23)

  !call CShockTest2D(gridIDn=122,gridIDi=123)  !!! zwg: added by zwg

  !call WardleInstability(gridIDn=26,gridIDi=27)

  !call coreCollapseAD2D(gridIDn=28,gridIDi=29)

   !call SodShockTube3D(gridID=30)

  !call isoShockTube3D(gridID=31)
  !call BrioWuShockTube3D(gridID=32)
  !call FieldLoopAdvection(gridID=33)
  !call MHDBlastWave(gridID=34)

  call selfgravityTest3D(gridID=35)

  !call HDBlastWavePolar2D(gridID=36)

  !call HDBlastWaveCart2D(gridID=37)

  !call HDBlastWavePolarIso2D(gridID=38)
  !call HDBlastWaveCartIso2D(gridID=39)
  !call sphericalBlastWaveMHDPolar2D(gridID=40)
  !call forceBalanceMHDPolar2D(gridID=41)

  !call KelvinHelmholtzInstabilityHD2D(gridID=42)

  !call KelvinHelmholtzInstabilityMHD2D(gridID=43)

  !call RayleighTaylorInstabilityHD2D(gridID=44)

  !call RayleighTaylorInstabilityMHD2D(gridID=45)
  !call isoShockTubeMHD1D(gridID=46)
  !call OrsagTangVortexIso(gridID=47)

  !call isoShockTubeMHD3D(gridID=48)

  !call isoMHDsg3D(gridID=49)

  !call polyShockTubeHD1D(gridID=50)
  !call polyShockTubeHD2D(gridID=51)

   !call CShockTest3D(gridIDn=52,gridIDi=53)

   !call sphericalBlastWaveAD2D(gridIDn=54,gridIDi=55)

   !call TestDrivingTurbulence2DHD(gridID=399)

   !call TestDrivingTurbulence3DHD(gridID=400)

   !call TestDrivingTurbulence2DMHD(gridID=99)

  !call TestDrivingTurbulence3DMHD(gridID=100)

  !call AdiShockTubeMHD3D(gridID=199)

  !RICHTMYER¨CMESHKOV INSTABILITY
  !call HDRichtmyerMeshkovInstability2D(gridID=200)
  !call MHDRichtmyerMeshkovInstability2D(gridID=201)

  !call HDRichtmyerMeshkovInstability3D(gridID=202)

  !call MHDRichtmyerMeshkovInstability3D(gridID=203)
  !call HDRichtmyerMeshkovInstability2D_cylinder(gridID=300)
  !call MHDRichtmyerMeshkovInstability2D_cylinder(gridID=301)

  !call ADRichtmyerMeshkovInstability2D(gridIDn=210,gridIDi=210)
  !call ADRichtmyerMeshkovInstability3D(gridIDn=220,gridIDi=211)
endif



call MPI_FINALIZE(ierr)

end program main
