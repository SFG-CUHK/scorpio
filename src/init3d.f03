subroutine init3d(this,q)
use gridModule
use testSuiteMPI
class(grid)::this
double precision,dimension(1-this%nbuf:this%nMesh(1)+this%nbuf,1-this%nbuf:this%nMesh(2)+this%nbuf, & 
& 1-this%nbuf:this%nMesh(3)+this%nbuf, this%nvar)::q
double precision,dimension(1-this%nbuf:this%nMesh(1)+this%nbuf)::xc
double precision,dimension(1-this%nbuf:this%nMesh(2)+this%nbuf)::yc
double precision,dimension(1-this%nbuf:this%nMesh(3)+this%nbuf)::zc
integer::nx,ny,nz,nvar,nbuf,gridID,eosType
integer::i

gridID=this%gridID
xc=this%xc(1)%coords
yc=this%xc(2)%coords
zc=this%xc(3)%coords

if(testOnOff) then
select case (gridID)
   case (30)
     call initSodShockTube3D(this,q)
   case (31)
     call initIsoShockTube3D(this,q)
   case (32)
     call initBrioWuShockTube3D(this,q)
   case (33)
     call initFieldLoopAdvection(this,q)
   case (34)
     call initMHDBlastWave(this,q)
   case (35)
     call initselfgravityTest3D(this,q)
   case (48)
     call initIsoShockTubeMHD3D(this,q)
   case (49)
     call initIsoMHDsg3D(this,q)
   case (52)
     call initCShockTest3Dn(this,q)
   case (53)
     call initCShockTest3Di(this,q)
   case (199)
     call initAdiShockTubeMHD3D(this,q)
   case (202)
     call initHDRichtmyerMeshkovInstability3D(this,q)
   case (203)
     call initMHDRichtmyerMeshkovInstability3D(this,q)
   !case(100)
    !call initTestDrivingTurbulence3DMHD(this,q)
   case(400) !!! initialization for 3DHD turbulence driving!!!!!
     call initTurbulenceDriving3DHD(this,q)
end select
endif

end subroutine init3d



subroutine init3d_for_FFTW(this,q,Den_ini)
use gridModule
use testSuiteMPI
class(grid)::this
integer:: Nx, Ny,Nz,Totalprocs

double precision,dimension(1-this%nbuf:this%nMesh(1)+this%nbuf,1-this%nbuf:this%nMesh(2)+this%nbuf, &
                       &   1-this%nbuf:this%nMesh(3)+this%nbuf,this%nvar)::q
double precision,dimension(1:this%nMesh(1),1:this%nMesh(2),1:this%nMesh(3))::Den_ini

Nx=this%nMesh(1)
Ny=this%nMesh(2)
Ny=this%nMesh(3)

Totalprocs=nprocs

this%drivingWN_DT=8.0
this%Energy_DT=3.437068E6
this%zeta_DT=0.0


     !!!! zwg the following is added for check array(4,2,4)
     if(Totalprocs==1)then
       Den_ini(1,1,1)=1.0
       Den_ini(2,1,1)=2.0
       Den_ini(3,1,1)=5.0
       Den_ini(4,1,1)=0.0

       Den_ini(1,2,1)=3.0
       Den_ini(2,2,1)=2.0
       Den_ini(3,2,1)=5.0
       Den_ini(4,2,1)=0.0



       Den_ini(1,1,2)=1.0
       Den_ini(2,1,2)=2.0
       Den_ini(3,1,2)=5.0
       Den_ini(4,1,2)=0.0

       Den_ini(1,2,2)=3.0
       Den_ini(2,2,2)=2.0
       Den_ini(3,2,2)=5.0
       Den_ini(4,2,2)=0.0



       Den_ini(1,1,3)=1.0
       Den_ini(2,1,3)=2.0
       Den_ini(3,1,3)=5.0
       Den_ini(4,1,3)=0.0

       Den_ini(1,2,3)=3.0
       Den_ini(2,2,3)=2.0
       Den_ini(3,2,3)=5.0
       Den_ini(4,2,3)=0.0



       Den_ini(1,1,4)=1.0
       Den_ini(2,1,4)=2.0
       Den_ini(3,1,4)=5.0
       Den_ini(4,1,4)=0.0

       Den_ini(1,2,4)=3.0
       Den_ini(2,2,4)=2.0
       Den_ini(3,2,4)=5.0
       Den_ini(4,2,4)=0.0

     elseif(Totalprocs==2)then



       if(myid==0)then

       Den_ini(1,1,1)=1.0
       Den_ini(2,1,1)=2.0
       Den_ini(3,1,1)=5.0
       Den_ini(4,1,1)=0.0

       Den_ini(1,2,1)=3.0
       Den_ini(2,2,1)=2.0
       Den_ini(3,2,1)=5.0
       Den_ini(4,2,1)=0.0



       Den_ini(1,1,2)=1.0
       Den_ini(2,1,2)=2.0
       Den_ini(3,1,2)=5.0
       Den_ini(4,1,2)=0.0

       Den_ini(1,2,2)=3.0
       Den_ini(2,2,2)=2.0
       Den_ini(3,2,2)=5.0
       Den_ini(4,2,2)=0.0

       elseif(myid==1)then

       Den_ini(1,1,1)=1.0
       Den_ini(2,1,1)=2.0
       Den_ini(3,1,1)=5.0
       Den_ini(4,1,1)=0.0

       Den_ini(1,2,1)=3.0
       Den_ini(2,2,1)=2.0
       Den_ini(3,2,1)=5.0
       Den_ini(4,2,1)=0.0



       Den_ini(1,1,2)=1.0
       Den_ini(2,1,2)=2.0
       Den_ini(3,1,2)=5.0
       Den_ini(4,1,2)=0.0

       Den_ini(1,2,2)=3.0
       Den_ini(2,2,2)=2.0
       Den_ini(3,2,2)=5.0
       Den_ini(4,2,2)=0.0

       endif


     elseif(Totalprocs==4)then

       if(myid==0)then
       Den_ini(1,1,1)=1.0
       Den_ini(2,1,1)=2.0
       Den_ini(3,1,1)=5.0
       Den_ini(4,1,1)=0.0

       Den_ini(1,2,1)=3.0
       Den_ini(2,2,1)=2.0
       Den_ini(3,2,1)=5.0
       Den_ini(4,2,1)=0.0
       elseif(myid==1)then
       Den_ini(1,1,1)=1.0
       Den_ini(2,1,1)=2.0
       Den_ini(3,1,1)=5.0
       Den_ini(4,1,1)=0.0

       Den_ini(1,2,1)=3.0
       Den_ini(2,2,1)=2.0
       Den_ini(3,2,1)=5.0
       Den_ini(4,2,1)=0.0
       elseif(myid==2)then
       Den_ini(1,1,1)=1.0
       Den_ini(2,1,1)=2.0
       Den_ini(3,1,1)=5.0
       Den_ini(4,1,1)=0.0

       Den_ini(1,2,1)=3.0
       Den_ini(2,2,1)=2.0
       Den_ini(3,2,1)=5.0
       Den_ini(4,2,1)=0.0
       elseif(myid==3)then
       Den_ini(1,1,1)=1.0
       Den_ini(2,1,1)=2.0
       Den_ini(3,1,1)=5.0
       Den_ini(4,1,1)=0.0

       Den_ini(1,2,1)=3.0
       Den_ini(2,2,1)=2.0
       Den_ini(3,2,1)=5.0
       Den_ini(4,2,1)=0.0
       endif


     endif


end subroutine init3d_for_FFTW

