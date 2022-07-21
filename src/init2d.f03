subroutine init2d(this,q)
use gridModule
use testSuiteMPI
class(grid)::this
double precision,dimension(1-this%nbuf:this%nMesh(1)+this%nbuf,1-this%nbuf:this%nMesh(2)+this%nbuf,this%nvar)::q
double precision,dimension(1-this%nbuf:this%nMesh(1)+this%nbuf)::xc
double precision,dimension(1-this%nbuf:this%nMesh(2)+this%nbuf)::yc
integer::nx,nvar,nbuf,gridID,eosType
integer::i

gridID=this%gridID
xc=this%xc(1)%coords
yc=this%xc(2)%coords

if(testOnOff) then
select case (gridID)
   case (4)
     call initIsoShockTube2D(this,q)
   case (5)
     call initSodShockTube2D(this,q)
   case (6)
     call initWCShockTube2D(this,q)
   case (7)
     call initDoubleMachReflection(this,q)
   case (8)
     call initStrongRarefactionTest(this,q)
   case (9)
     call initShuOsherShocktubeTest(this,q)
   case (12)
     call initBrioWuShockTube2D(this,q)
   case (13)
     call initOrsagTangVortex(this,q)
   case (14)
     call initSphericalBlastWaveMHD2D(this,q)
   case (15)
     call initBFieldLoopTest2D(this,q)
   case (16)
     call initRotorMHD2D(this,q)
   case (17)
     call initCurrentSheetMHD2D(this,q)
   case (18)
     call initCPAlvenWave(this,q)
   case (19)
     call initLinearSlowWaveTestMHD2D(this,q)
   case (20)
     call initLinearFastWaveTestMHD2D(this,q)
   case (21)
     call initSelfgravityTest2D(this,q)
   case (26)
     call initWardleInstabilityn(this,q)
   case (27)
     call initWardleInstabilityi(this,q)
   case (122)
     call initCShockTest2Dn(this,q)
   case (123)
     call initCShockTest2Di(this,q)
   case (28)
     call initCoreCollapseAD2Dn(this,q)
   case (29)
     call initCoreCollapseAD2Di(this,q)
   case (36)
     call initHDBlastWavePolar2D(this,q)
   case (37)
     call initHDBlastWaveCart2D(this,q)
   case (38)
     call initHDBlastWavePolarIso2D(this,q)
   case (39)
     call initHDBlastWaveCartIso2D(this,q)
   case (40)
     call initSphericalBlastWaveMHDPolar2D(this,q)
   case (41)
     call initForceBalanceMHDPolar2D(this,q)
   case (42)
     call initKelvinHelmholtzInstabilityHD2D(this,q)
   case (43)
     call initKelvinHelmholtzInstabilityMHD2D(this,q)
   case (44)
     call initRayleighTaylorInstabilityHD2D(this,q)
   case (45)
     call initRayleighTaylorInstabilityMHD2D(this,q)
   case (47)
     call initOrsagTangVortexIso(this,q)
   case (51)
     call initPolyShockTubeHD2D(this,q)
   case (54)
     call initSphericalBlastWaveAD2Dn(this,q)
   case (55)
     call initSPhericalBlastWaveAD2Di(this,q)
   !case(99)
     !call initTestDrivingTurbulence2DMHD(this,q)

   case(200)
     call initHDRichtmyerMeshkovInstability2D(this,q) 
   case(201)
     call initMHDRichtmyerMeshkovInstability2D(this,q)
   case(399) !!! initialization  for 2DHD turbulence driving!!!!!
     call initTurbulenceDriving2DHD(this,q)


end select
endif

end subroutine init2d

subroutine init2d_for_FFTW(this,q,Den_ini)
use gridModule
use testSuiteMPI
class(grid)::this
integer:: Nx, Ny,Totalprocs

double precision,dimension(1-this%nbuf:this%nMesh(1)+this%nbuf,1-this%nbuf:this%nMesh(2)+this%nbuf,this%nvar)::q
double precision,dimension(1:this%nMesh(1),1:this%nMesh(2))::Den_ini

Nx=this%nMesh(1)
Ny=this%nMesh(2)
Totalprocs=nprocs

this%drivingWN_DT=8.0
this%Energy_DT=3.437068E6
this%zeta_DT=0.0


     !!!! zwg the following is added for check array(4,4)
     if(Totalprocs==1)then
       Den_ini(1,1)=1.0
       Den_ini(2,1)=2.0
       Den_ini(3,1)=5.0
       Den_ini(4,1)=0.0

       Den_ini(1,2)=3.0
       Den_ini(2,2)=2.0
       Den_ini(3,2)=5.0
       Den_ini(4,2)=0.0

       Den_ini(1,3)=1.0
       Den_ini(2,3)=2.0
       Den_ini(3,3)=5.0
       Den_ini(4,3)=0.0

       Den_ini(1,4)=3.0
       Den_ini(2,4)=2.0
       Den_ini(3,4)=5.0
       Den_ini(4,4)=0.0

     elseif(Totalprocs==2)then



       if(myid==0)then

       Den_ini(1,1)=1.0
       Den_ini(2,1)=2.0
       Den_ini(3,1)=5.0
       Den_ini(4,1)=0.0

       Den_ini(1,2)=3.0
       Den_ini(2,2)=2.0
       Den_ini(3,2)=5.0
       Den_ini(4,2)=0.0
       elseif(myid==1)then

       Den_ini(1,1)=1.0
       Den_ini(2,1)=2.0
       Den_ini(3,1)=5.0
       Den_ini(4,1)=0.0

       Den_ini(1,2)=3.0
       Den_ini(2,2)=2.0
       Den_ini(3,2)=5.0
       Den_ini(4,2)=0.0

       endif


     elseif(Totalprocs==4)then

       if(myid==0)then
       Den_ini(1,1)=1.0
       Den_ini(2,1)=2.0
       Den_ini(3,1)=5.0
       Den_ini(4,1)=0.0
       elseif(myid==1)then
       Den_ini(1,1)=3.0
       Den_ini(2,1)=2.0
       Den_ini(3,1)=5.0
       Den_ini(4,1)=0.0
       elseif(myid==2)then
       Den_ini(1,1)=1.0
       Den_ini(2,1)=2.0
       Den_ini(3,1)=5.0
       Den_ini(4,1)=0.0
       elseif(myid==3)then
       Den_ini(1,1)=3.0
       Den_ini(2,1)=2.0
       Den_ini(3,1)=5.0
       Den_ini(4,1)=0.0
       endif


     endif


end subroutine init2d_for_FFTW

