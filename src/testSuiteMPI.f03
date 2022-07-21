module testSuiteMPI
use gridModule
use mpi
implicit none
logical:: testOnOff=.false.

contains



subroutine setTestOnOff(switch)
logical::switch
testOnOff=switch
end subroutine



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!! Selfgravitating isothermal MHD 3D
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine isoMHDsg3D(gridID)
integer::gridID
type(grid)::g1
integer::ndim,nbuf,coordType,variable(8)
integer::nMesh(3),dims(3)
double precision::leftBdry(3),rightBdry(3)
logical::periods(3),reorder
double precision::t0,t1
integer::nstep,ierr

nstep=0
variable=0
ndim=3
nbuf=2
coordType=1
variable(1)=1  !! density
variable(2)=1  !! momx
variable(3)=1  !! momy
variable(4)=1  !! momz
variable(5)=1  !! bx
variable(6)=1  !! by
variable(7)=1  !! bz
variable(8)=1  !! ene, it is also required for isothermal MHD

nMesh(1)=64
nMesh(2)=64
nMesh(3)=64
leftBdry(1)=-0.5d0
leftBdry(2)=-0.5d0
leftBdry(3)=-0.5d0

rightBdry(1)=0.5d0
rightBdry(2)=0.5d0
rightBdry(3)=0.5d0

dims(1)=1
dims(2)=1
dims(3)=nprocs
periods(1)=.false.
periods(2)=.false.
periods(3)=.false.
reorder=.true.


call g1%enableSelfgravity()
call g1%setSgBdryType(sgBdryType=1) ! 0=isolated, 1=periodic
call g1%setGravConst(GravConst=1.d0)
call g1%setTopologyMPI(ndim,dims,periods,reorder)
call g1%setMesh(nMesh,leftBdry,rightBdry,nbuf,coordType,gridID)
call g1%setVariable(variable)
call g1%setSelfgravityKernel()
call g1%setSGMPIWindows()
call g1%setMPIWindows()
call g1%setEOS(eosType=1)
call g1%setSoundSpeed(snd=0.23d0)
call g1%setCFL(CFL=0.4d0)
call g1%setTime(fstart=0,tend=9.d0,dtout=0.3d0)
call g1%setSlopeLimiter(limiterType=3)
call g1%setSolverType(solverType=5)
call g1%setBoundaryType(boundaryType=1)
call g1%initVariable()
call g1%exchangeBdryMPI(g1%q,g1%winq)
call g1%setBoundary(g1%q)
call g1%calcSelfgravity(g1%q)
call g1%writeGrid()
call g1%writeGrid_MHD_vtk()
call g1%griddt()

g1%writeFlag=.false.

t0=MPI_WTIME()
do while (g1%t .lt. g1%tend)
   nstep=nstep+1
   call g1%griddt()
   call g1%evolveGridRK2()
   if(myid .eq. 0 ) then
     print *,"testSuite.f03: gridID=",g1%gridID,'t=',g1%t,"dt=",g1%dt
   endif
   if(g1%writeFlag .eqv. .true.) then
      t1=MPI_WTIME()
      if(myid .eq. 0) then
        print *,"testSuiteMPI.f03: elapse time=",t1-t0,"secs"
      endif
      call g1%writeGrid()
      call g1%writeGrid_MHD_vtk()
      g1%writeFlag=.false.
      t0=MPI_WTIME()
   endif
enddo !! end do while

end subroutine isoMHDsg3D




subroutine initIsoMHDsg3D(this,q)
class(grid)::this
double precision,dimension(1-this%nbuf:this%nMesh(1)+this%nbuf,1-this%nbuf:this%nMesh(2)+this%nbuf,&
                           1-this%nbuf:this%nMesh(3)+this%nbuf,this%nvar)::q
double precision,dimension(1-this%nbuf:this%nMesh(1)+this%nbuf)::xc
double precision,dimension(1-this%nbuf:this%nMesh(2)+this%nbuf)::yc
double precision,dimension(1-this%nbuf:this%nMesh(3)+this%nbuf)::zc
double precision::x,y,z,r,r2
double precision::pi,rho0,gam,ep,r0,Vrot,vx,vy,vz,Omega
integer::nx,ny,nz,nvar,nbuf,gridID,eosType
integer::i,j,k

pi=datan(1.d0)*4.d0
xc=this%xc(1)%coords
yc=this%xc(2)%coords
zc=this%xc(3)%coords
nx=this%nMesh(1)
ny=this%nMesh(2)
nz=this%nMesh(3)
gam=this%adiGamma
r0=0.05d0
rho0=1.d0
Omega=1.d0

do k=1,nz
  do j=1,ny
    do i=1,nx
         x=xc(i)
         y=yc(j)
         z=zc(k)
         r=dsqrt(x**2+y**2+z**2)
         r2=dsqrt(x**2+y**2)
         
         if(r .le. r0) then
           vx=-Omega*y  !! for a collapsing rotating sphere
           vy= Omega*x
           vz=0.d0
           q(i,j,k,1)=rho0 !rho0/r !rho0    !! for a collapsing rotating sphere
           q(i,j,k,2)=q(i,j,k,1)*vx
           q(i,j,k,3)=q(i,j,k,1)*vy
           q(i,j,k,4)=q(i,j,k,1)*vz
         else
           vx=0.d0
           vy=0.d0
           vz=0.d0
           q(i,j,k,1)=rho0/1000.d0 !/r0
           q(i,j,k,2)=0.d0
           q(i,j,k,3)=0.d0
           q(i,j,k,4)=0.d0
         endif
         q(i,j,k,5)=0.d0
         q(i,j,k,6)=0.d0
         q(i,j,k,7)=0.d0
         q(i,j,k,8)=0.d0
         q(i,j,k,9)=0.d0
         q(i,j,k,10)=0.d0
         q(i,j,k,11)=0.d0
    enddo
  enddo
enddo
end subroutine initIsoMHDsg3D

subroutine bdryIsoMHDsg3D(this,q)
class(grid)::this
double precision,dimension(1-this%nbuf:this%nMesh(1)+this%nbuf,1-this%nbuf:this%nMesh(2)+this%nbuf,&
& 1-this%nbuf:this%nMesh(3)+this%nbuf,this%nvar)::q
integer::nx,ny,nz,nbuf,i,j,k

nbuf=this%nbuf
nx=this%nMesh(1)
ny=this%nMesh(2)
nz=this%nMesh(3)

if(this%boundaryType .eq. 1) then  !! zero gradient boundary condition
!!!! left boundary
   if(this%left_mpi .lt. 0) then
     do i=1,nbuf
       q(1-i,:,:,:)=q(1,:,:,:)
     enddo
   endif
!!!! right boundary
   if(this%right_mpi .lt. 0) then
     do i=1,nbuf
       q(nx+i,:,:,:)=q(nx,:,:,:)
     enddo
   endif
!!!! up boundary
   if(this%up_mpi .lt. 0) then
     do j=1,nbuf
       q(:,1-j,:,:)=q(:,1,:,:)
     enddo
   endif
!!!! down boundary
   if(this%down_mpi .lt. 0) then
     do j=1,nbuf
       q(:,ny+j,:,:)=q(:,ny,:,:)
     enddo
   endif
!!!! top boundary
   if(this%top_mpi .lt. 0) then
     do k=1,nbuf
       q(:,:,1-k,:)=q(:,:,1,:)
     enddo
   endif
!!!! bottom boundary
   if(this%bottom_mpi .lt. 0) then
     do k=1,nbuf
       q(:,:,nz+k,:)=q(:,:,nz,:)
     enddo
   endif
!!!!
endif
end subroutine bdryIsoMHDsg3D




!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!! Isothermal MHD Shocktube (A. Mignone, 2007, JCP)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine isoShockTubeMHD3D(gridID)
integer::gridID
type(grid)::g1
integer::ndim,nbuf,coordType,variable(8)
integer::nMesh(3),dims(3)
double precision::leftBdry(3),rightBdry(3)
logical::periods(3),reorder
double precision::t0,t1
integer::nstep,ierr

nstep=0
variable=0
ndim=3
nbuf=2
coordType=1 !! Cartesian
variable(1)=1  !! rho
variable(2)=1  !! vx
variable(3)=1  !! vy
variable(4)=1  !! vz
variable(5)=1  !! bx
variable(6)=1  !! by
variable(7)=1  !! bz
variable(8)=1  !! energy

nMesh(1)=100
nMesh(2)=100
nMesh(3)=400
leftBdry(1)=0.d0
rightBdry(1)=0.25d0
leftBdry(2)=0.d0
rightBdry(2)=0.25d0
leftBdry(3)=0.d0
rightBdry(3)=1.d0
dims=(/2,1,4/)
call MPI_DIMS_CREATE(nprocs,ndim,dims,ierr)
!dims(1)=1
!dims(2)=1
!dims(3)=nprocs
periods(1)=.false.
periods(2)=.false.
periods(3)=.false.
reorder=.true.

call g1%setTopologyMPI(ndim,dims,periods,reorder)
call g1%setMesh(nMesh,leftBdry,rightBdry,nbuf,coordType,gridID)
call g1%setVariable(variable) !! den,vx,vy,vz,bx,by,bz,ene
call g1%setMPIWindows()
call g1%setEoS(eosType=1) !! isothermal
call g1%setSoundSpeed(snd=1.d0) !! sound speed
call g1%setCFL(CFL=0.8d0) !! courant number
call g1%setTime(fstart=0,tend=0.1d0,dtout=0.1d0)
call g1%setSlopeLimiter(limiterType=3)
call g1%setSolverType(solverType=5) !! 4=AdiHLLMHD, 5=AdiHLLDMHD
call g1%setBoundaryType(boundaryType=1)
call g1%initVariable()
call g1%setBoundary(g1%q)
call g1%exchangeBdryMPI(g1%q,g1%winq)
call g1%setBoundary(g1%q)
call g1%writeGrid()
call g1%writeGrid_MHD_vtk()

g1%writeFlag=.false.

t0=MPI_WTIME()
do while (g1%t .lt. g1%tend)
   nstep=nstep+1
   call g1%griddt()
   call g1%evolveGridRK2()
   if(myid .eq. 0 ) then
     print *,"testSuite.f03: gridID=",g1%gridID,'t=',g1%t,"dt=",g1%dt
   endif
   if(g1%writeFlag .eqv. .true.) then
      t1=MPI_WTIME()
      if(myid .eq. 0) then
        print *,"testSuiteMPI.f03: elapse time=",t1-t0,"secs"
      endif
      call g1%writeGrid()
      call g1%writeGrid_MHD_vtk()
      g1%writeFlag=.false.
      t0=MPI_WTIME()
   endif

enddo !! end do while


end subroutine isoShockTubeMHD3D

subroutine initIsoShockTubeMHD3D(this,q)
class(grid)::this
double precision,dimension(1-this%nbuf:this%nMesh(1)+this%nbuf,1-this%nbuf:this%nMesh(2)+this%nbuf,&
                           1-this%nbuf:this%nMesh(3)+this%nbuf,this%nvar)::q
double precision,dimension(1-this%nbuf:this%nMesh(1)+this%nbuf)::xc
double precision,dimension(1-this%nbuf:this%nMesh(2)+this%nbuf)::yc
double precision,dimension(1-this%nbuf:this%nMesh(3)+this%nbuf)::zc

integer::nx,ny,nz,nvar,nbuf,gridID,eosType
integer::i,j,k
double precision::gam,rho,vx,vy,vz,bxl,byl,bzl,bxr,byr,bzr,p
double precision::bxc,byc,bzc,pi,snd

xc=this%xc(1)%coords
yc=this%xc(2)%coords
zc=this%xc(3)%coords
nx=this%nMesh(1)
ny=this%nMesh(2)
nz=this%nMesh(3)
snd=this%snd
pi=4.d0*datan2(1.d0,1.d0)


do k=1,nz
  do j=1,ny
    do i=1,nx
      if(zc(k) .lt. 0.5d0) then
        rho=1.d0
         vx=0.d0
         vy=0.d0
         vz=0.d0
         bxl=0.d0
         bxr=bxl
         byl=5.d0/dsqrt(4.d0*pi)
         byr=byl
         bzl=3.d0/dsqrt(4.d0*pi)
         bzr=bzl
      else
         rho=0.1d0
          vx=0.d0
          vy=0.d0
          vz=0.d0
          bxl=0.d0
          bxr=bxl
          byl=2.d0/dsqrt(4.d0*pi)
          byr=byl
          bzl=3.d0/dsqrt(4.d0*pi)
          bzr=bzl
      endif
      bxc=0.5d0*(bxl+bxr)
      byc=0.5d0*(byl+byr)
      bzc=0.5d0*(bzl+bzr)

      q(i,j,k,1)=rho
      q(i,j,k,2)=rho*vx
      q(i,j,k,3)=rho*vy
      q(i,j,k,4)=rho*vz
      q(i,j,k,5)=bxl
      q(i,j,k,6)=byl
      q(i,j,k,7)=bzl
      q(i,j,k,8)=0.d0
      q(i,j,k,9)=bxr
      q(i,j,k,10)=byr
      q(i,j,k,11)=bzr
    enddo
  enddo
enddo

end subroutine initIsoShockTubeMHD3D

subroutine bdryIsoShockTubeMHD3D(this,q)
class(grid)::this
double precision,dimension(1-this%nbuf:this%nMesh(1)+this%nbuf,1-this%nbuf:this%nMesh(2)+this%nbuf,&
                           1-this%nbuf:this%nMesh(3)+this%nbuf,this%nvar)::q
integer::nx,ny,nz,nbuf,i,j,k

nbuf=this%nbuf
nx=this%nMesh(1)
ny=this%nMesh(2)
nz=this%nMesh(3)

if(this%boundaryType .eq. 1) then !! zero gradient
!!!!!! left boundary
  if(this%left_mpi .lt. 0) then
    do i=1,nbuf
       q(1-i,:,:,:)=q(1,:,:,:)
    enddo
  endif
!!!!!! right boundary
  if(this%right_mpi .lt. 0) then
    do i=1,nbuf
       q(nx+i,:,:,:)=q(nx,:,:,:)
    enddo
  endif
!!!!!! up boundary
  if(this%up_mpi .lt. 0) then
    do j=1,nbuf
       q(:,1-j,:,:)=q(:,1,:,:)
    enddo
  endif
!!!!!! down boundary
  if(this%down_mpi .lt. 0) then
    do j=1,nbuf
       q(:,ny+j,:,:)=q(:,ny,:,:)
    enddo
  endif
!!!! top boundary
   if(this%top_mpi .lt. 0) then
     do k=1,nbuf
       q(:,:,1-k,:)=q(:,:,1,:)
     enddo
   endif
!!!! bottom boundary
   if(this%bottom_mpi .lt. 0) then
     do k=1,nbuf
       q(:,:,nz+k,:)=q(:,:,nz,:)
     enddo
   endif
!!!!

endif

end subroutine bdryIsoShockTubeMHD3D



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!! Orszag & Tang isothermal MHD (A. Mignone, 2007, JCP)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine OrsagTangVortexIso(gridID)
integer::gridID
type(grid)::g1
integer::ndim,nbuf,coordType,variable(8)
integer::nMesh(2),dims(2)
double precision::leftBdry(2),rightBdry(2)
double precision::pi
logical::periods(2),reorder
double precision::t0,t1
integer::nstep,ierr


nstep=0
pi=4.d0*atan2(1.d0,1.d0)
variable=0
ndim=2
nbuf=2
coordType=1 !! Cartesian
variable(1)=1  !! rho
variable(2)=1  !! vx
variable(3)=1  !! vy
variable(4)=1  !! vz
variable(5)=1  !! bx
variable(6)=1  !! by
variable(7)=1  !! bz
variable(8)=1  !! energy, also required for isothermal simulation

nMesh(1)=400
nMesh(2)=400
leftBdry(1)=0.d0
rightBdry(1)=2.d0*pi
leftBdry(2)=0.d0
rightBdry(2)=2.d0*pi
dims=(/0,0/)
call MPI_DIMS_CREATE(nprocs,ndim,dims,ierr)
!dims(1)=1
!dims(2)=nprocs
periods(1)=.true.
periods(2)=.true.
reorder=.true.

call g1%setTopologyMPI(ndim,dims,periods,reorder)
call g1%setMesh(nMesh,leftBdry,rightBdry,nbuf,coordType,gridID)
call g1%setVariable(variable) !! den,vx,vy,vz,bx,by,bz,ene
call g1%setMPIWindows()
call g1%setEoS(eosType=1) !! isothermal
call g1%setSoundSpeed(snd=2.d0) !! isothermal sound speed
call g1%setCFL(CFL=0.4d0) !! courant number
call g1%setTime(fstart=0,tend=3.d0,dtout=0.5d0)
call g1%setSlopeLimiter(limiterType=3)
call g1%setSolverType(solverType=5) !! 4=AdiHLLMHD, 5=AdiHLLDMHD
call g1%setBoundaryType(boundaryType=3) !! 3=periodic
call g1%initVariable()
call g1%exchangeBdryMPI(g1%q,g1%winq)
call g1%setBoundary(g1%q)
call g1%writeGrid()
g1%writeFlag=.false.

t0=MPI_WTIME()
do while (g1%t .lt. g1%tend)
   call g1%griddt()
   call g1%evolveGridRK2()
   if(myid .eq. 0) then
     print *,"testSuite.f03: gridID=",g1%gridID,'t=',g1%t,"dt=",g1%dt
   endif
   if(g1%writeFlag .eqv. .true.) then
      t1=MPI_WTIME()
      if(myid .eq. 0) then
        print *,"testSuiteMPI.f03: elapse time=",t1-t0,"secs"
      endif
      call g1%writeGrid()
      g1%writeFlag=.false.
      t0=MPI_WTIME()
   endif
enddo !! end do while

end subroutine OrsagTangVortexIso

subroutine initOrsagTangVortexIso(this,q)
class(grid)::this
double precision,dimension(1-this%nbuf:this%nMesh(1)+this%nbuf,1-this%nbuf:this%nMesh(2)+this%nbuf,this%nvar)::q
double precision,dimension(1-this%nbuf:this%nMesh(1)+this%nbuf)::xc,xl,xr
double precision,dimension(1-this%nbuf:this%nMesh(2)+this%nbuf)::yc,yl,yr
integer::nx,ny,nvar,nbuf,gridID,eosType
integer::i,j
double precision::gam,rho,vx,vy,vz,bxl,byl,bzl,bxr,byr,bzr,p
double precision::bxc,byc,bzc
double precision::pi,b0,snd


pi=4.d0*atan2(1.d0,1.d0)
snd=this%snd
b0=snd*dsqrt(3.d0/5.d0)
xc=this%xc(1)%coords
yc=this%xc(2)%coords
xl=this%xl(1)%coords
yl=this%xl(2)%coords
xr=this%xr(1)%coords
yr=this%xr(2)%coords
nx=this%nMesh(1)
ny=this%nMesh(2)
gam=this%adiGamma

do j=1,ny
  do i=1,nx

    rho=1.d0
    vx=-snd*dsin(yc(j))
    vy= snd*dsin(xc(i))
    vz=0.d0
   bxl=-b0*dsin(yc(j))
   byl= b0*dsin(2.d0*xc(i))
   bzl=0.d0
   bxr=bxl
   byr=byl
   bzr=bzl

    bxc=0.5d0*(bxl+bxr)
    byc=0.5d0*(byl+byr)
    bzc=0.5d0*(bzl+bzr)

    q(i,j,1)=rho
    q(i,j,2)=rho*vx
    q(i,j,3)=rho*vy
    q(i,j,4)=rho*vz
    q(i,j,5)=bxl
    q(i,j,6)=byl
    q(i,j,7)=bzl
    q(i,j,8)=0.d0
    q(i,j,9)=bxr
    q(i,j,10)=byr
    q(i,j,11)=bzr
  enddo
enddo
end subroutine initOrsagTangVortexIso

subroutine bdryOrsagTangVortexIso(this,q)
class(grid)::this
double precision,dimension(1-this%nbuf:this%nMesh(1)+this%nbuf,1-this%nbuf:this%nMesh(2)+this%nbuf,this%nvar)::q
integer::nx,ny,nbuf,i,j

nbuf=this%nbuf
nx=this%nMesh(1)
ny=this%nMesh(2)

if(this%boundaryType .eq. 3) then !! periodic boundary condition
endif

end subroutine bdryOrsagTangVortexIso

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!! isothermal MHD Shocktube, (A. Mignone, 2007, JCP)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine isoShockTubeMHD1D(gridID)
integer::gridID
type(grid)::g1
integer::ndim,nbuf,coordType,variable(8)
integer::nMesh(1),dims(1)
double precision::leftBdry(1),rightBdry(1)
logical::periods(1),reorder
double precision::t0,t1
integer::nstep,ierr

nstep=0
variable=0
ndim=1
nbuf=2
coordType=1 !! Cartesian
variable(1)=1  !! rho
variable(2)=1  !! vx
variable(3)=1  !! vy
variable(4)=1  !! vz
variable(5)=1  !! bx
variable(6)=1  !! by
variable(7)=1  !! bz
variable(8)=1  !! energy, needs to be selected even for isothermal

nMesh(1)=800
leftBdry(1)=0.d0
rightBdry(1)=1.d0
dims(1)=nprocs
periods(1)=.false.
reorder=.true.

call g1%setTopologyMPI(ndim,dims,periods,reorder)
call g1%setMesh(nMesh,leftBdry,rightBdry,nbuf,coordType,gridID)
call g1%setVariable(variable) !! den,vx,vy,vz,bx,by,bz,ene
call g1%setMPIWindows()
call g1%setEoS(eosType=1) !! isothermal
call g1%setSoundSpeed(snd=1.d0) !! set sound speed
call g1%setCFL(CFL=0.8d0) !! courant number
call g1%setTime(fstart=0,tend=0.1d0,dtout=0.1d0)
call g1%setSlopeLimiter(limiterType=3)
call g1%setSolverType(solverType=4) !! 4=isoHLLMHD, 5=isoHLLDMHD
call g1%setBoundaryType(boundaryType=1)
call g1%initVariable()
call g1%exchangeBdryMPI(g1%q,g1%winq)

call g1%setBoundary(g1%q)
call g1%writeGrid()
g1%writeFlag=.false.

t0=MPI_WTIME()
do while (g1%t .lt. g1%tend)
   call g1%griddt()
   call g1%evolveGridRK2()
   nstep=nstep+1
   if(myid .eq. 0) then
     !print *,nstep," testSuite.f03: gridID=",g1%gridID,'t=',g1%t,"dt=",g1%dt
   endif
   if(g1%writeFlag .eqv. .true.) then
      t1=MPI_WTIME()
      if(myid .eq. 0) then
        print *,"testSuiteMPI.f03: elapse time=",t1-t0,"secs"
      endif
      call g1%writeGrid()
      g1%writeFlag=.false.
      t0=MPI_WTIME()
   endif
enddo !! end do while
end subroutine isoShockTubeMHD1D

subroutine initIsoShockTubeMHD1D(this,q)
class(grid)::this
double precision,dimension(1-this%nbuf:this%nMesh(1)+this%nbuf,this%nvar)::q
double precision,dimension(1-this%nbuf:this%nMesh(1)+this%nbuf)::xc
integer::nx,nvar,nbuf,gridID,eosType
integer::i
double precision::gam,rho,vx,vy,vz,bxl,byl,bzl,bxr,byr,bzr,p
double precision::bxc,byc,bzc,pi,b0

xc=this%xc(1)%coords
nx=this%nMesh(1)
gam=this%adiGamma
pi=4.d0*datan2(1.d0,1.d0)

do i=1,nx
    if(xc(i) .lt. 0.5d0) then
      rho=1.d0
       vx=0.d0
       vy=0.d0
       vz=0.d0
       bxl=3.d0/dsqrt(4.d0*pi)
       bxr=bxl
       byl=5.d0/dsqrt(4.d0*pi)
       byr=byl
       bzl=0.d0
       bzr=0.d0
    else
       rho=0.1d0
        vx=0.d0
        vy=0.d0
        vz=0.d0
        bxl=3.d0/dsqrt(4.d0*pi)
        bxr=bxl
        byl=2.d0/dsqrt(4.d0*pi)
        byr=byl
        bzl=0.d0
        bzr=0.d0
    endif
    bxc=0.5d0*(bxl+bxr)
    byc=0.5d0*(byl+byr)
    bzc=0.5d0*(bzl+bzr)

    q(i,1)=rho
    q(i,2)=rho*vx
    q(i,3)=rho*vy
    q(i,4)=rho*vz
    q(i,5)=bxl
    q(i,6)=byl
    q(i,7)=bzl
    q(i,8)=0.d0
    q(i,9)=bxr
    q(i,10)=byr
    q(i,11)=bzr
enddo

end subroutine initIsoShockTubeMHD1D

subroutine bdryIsoShockTubeMHD1D(this,q)
class(grid)::this
double precision,dimension(1-this%nbuf:this%nMesh(1)+this%nbuf,this%nvar)::q
integer::nbuf,nx,i

nbuf=this%nbuf
nx=this%nMesh(1)

if(this%boundaryType .eq. 1) then !! zero gradient
  if(this%left_mpi .lt. 0) then
    do i=1,nbuf
       q(i-nbuf,:)=q(nbuf-i+1,:)
    enddo
  endif
  if(this%right_mpi .lt. 0) then
    do i=1,nbuf
       q(nx+i,:)=q(nx-i+1,:)
    enddo
  endif
endif
end subroutine bdryIsoShockTubeMHD1D

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!! Rayleigh-Taylor instability MHD 2D
!!!!! http://www.astro.princeton.edu/~jstone/Athena/tests/kh/kh.html
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine RayleighTaylorInstabilityMHD2D(gridID)
integer::gridID
type(grid)::g1
integer::ndim,nbuf,coordType,variable(8)
integer::nMesh(2),dims(2)
double precision::leftBdry(2),rightBdry(2)
double precision::pi
logical::periods(2),reorder
double precision::t0,t1
integer::nstep,ierr

nstep=0
pi=4.d0*atan2(1.d0,1.d0)
variable=0
ndim=2
nbuf=2
coordType=1 !! Cartesian
variable(1)=1  !! rho
variable(2)=1  !! vx
variable(3)=1  !! vy
variable(4)=1  !! vz
variable(5)=1  !! bx
variable(6)=1  !! by
variable(7)=1  !! bz
variable(8)=1  !! ene

nMesh(1)=400
nMesh(2)=600
leftBdry(1)=-0.25d0
rightBdry(1)=0.25d0
leftBdry(2)=-0.375d0
rightBdry(2)=0.375d0
dims=(/0,0/)
call MPI_DIMS_CREATE(nprocs,ndim,dims,ierr)
dims(1)=1
dims(2)=nprocs
periods(1)=.true.
periods(2)=.false.
reorder=.true.

call g1%setTopologyMPI(ndim,dims,periods,reorder)
call g1%setMesh(nMesh,leftBdry,rightBdry,nbuf,coordType,gridID)
call g1%setVariable(variable) !! den,vx,vy,vz,bx,by,bz,ene
call g1%setMPIWindows()
call g1%setEoS(eosType=2) !! adiabatic
call g1%setadiGamma(gam=1.4d0) !! ratio of heat capacity
call g1%setCFL(CFL=0.4d0) !! courant number
call g1%setTime(fstart=0,tend=12.75d0,dtout=1.d0)
call g1%setSlopeLimiter(limiterType=3)
call g1%setSolverType(solverType=5) !! 4=AdiHLLMHD, 5=AdiHLLDMHD
call g1%setBoundaryType(boundaryType=0) !! 3=periodic
call g1%initVariable()
call g1%exchangeBdryMPI(g1%q,g1%winq)
call g1%setBoundary(g1%q)
call g1%writeGrid()
g1%writeFlag=.false.

t0=MPI_WTIME()
do while (g1%t .lt. g1%tend)
   call g1%griddt()
   call g1%evolveGridRK2()
   if(myid .eq. 0) then
     print *,"testSuite.f03: gridID=",g1%gridID,'t=',g1%t,"dt=",g1%dt
   endif
   if(g1%writeFlag .eqv. .true.) then
      t1=MPI_WTIME()
      if(myid .eq. 0) then
        print *,"testSuiteMPI.f03: elapse time=",t1-t0,"secs"
      endif
      call g1%writeGrid()
      g1%writeFlag=.false.
      t0=MPI_WTIME()
   endif
enddo !! end do while
end subroutine RayleighTaylorInstabilityMHD2D

subroutine initRayleighTaylorInstabilityMHD2D(this,q)
class(grid)::this
double precision,dimension(1-this%nbuf:this%nMesh(1)+this%nbuf,1-this%nbuf:this%nMesh(2)+this%nbuf,this%nvar)::q
double precision,dimension(1-this%nbuf:this%nMesh(1)+this%nbuf)::xc
double precision,dimension(1-this%nbuf:this%nMesh(2)+this%nbuf)::yc
integer::nx,ny,nvar,nbuf,gridID,eosType
integer::i,j
double precision::gam,rho,vx,vy,vz,p,p0,pi,rand1,rand2,A,bx0,by0,bz0

  nx=this%nMesh(1)
  ny=this%nMesh(2)
 gam=this%adiGamma
  xc=this%xc(1)%coords
  yc=this%xc(2)%coords
  p0=2.5d0
  pi=datan2(1.d0,1.d0)*4.d0
  call init_random_seed()
  rand1=0.01d0*(rand()-0.5d0)
  rand2=0.01d0*(rand()-0.5d0)
  bx0=0.0125d0
  by0=0.d0
  bz0=0.d0

  do j=1-nbuf,ny+nbuf
    do i=1-nbuf,nx+nbuf
       if(yc(j) .gt. 0.d0) then
         A=rand1
       else
         A=rand2
       endif
       vx=0.d0
       vz=0.d0
       vy=A*(rand()-0.5d0)*(1+dcos(8.d0*pi*yc(j)/3.d0))
       !vy=0.01d0*(1.d0+dcos(4.d0*pi*xc(i)))*(1.d0+dcos(3.d0*pi*yc(j)))/4.d0
      if(yc(j) .gt. 0.d0) then
        rho=2.d0
      else
        rho=1.d0
      endif
          p=p0-0.1d0*rho*yc(j)

       q(i,j,1)=rho
       q(i,j,2)=rho*vx
       q(i,j,3)=rho*vy
       q(i,j,4)=rho*vz
       q(i,j,5)=bx0
       q(i,j,6)=by0
       q(i,j,7)=bz0
       q(i,j,8)=p/(gam-1.d0)+0.5d0*rho*(vx**2+vy**2+vz**2)+0.5d0*(bx0**2+by0**2+bz0**2)
       q(i,j,9)=bx0
       q(i,j,10)=by0
       q(i,j,11)=bz0
    enddo
  enddo
end subroutine initRayleighTaylorInstabilityMHD2D

subroutine bdryRayleighTaylorInstabilityMHD2D(this,q)
class(grid)::this
double precision,dimension(1-this%nbuf:this%nMesh(1)+this%nbuf,1-this%nbuf:this%nMesh(2)+this%nbuf,this%nvar)::q
integer::nx,ny,nbuf,i,j

nbuf=this%nbuf
nx=this%nMesh(1)
ny=this%nMesh(2)

if(this%up_mpi .lt. 0)then
  do j=1-nbuf,0
    do i=1-nbuf,nx+nbuf
         q(i,j,1)=q(i,-j+1,1)
         q(i,j,2)=q(i,-j+1,2)
         q(i,j,3)=-q(i,-j+1,3)
         q(i,j,4)=q(i,-j+1,4)
         q(i,j,5)=q(i,-j+1,5)
         q(i,j,6)=q(i,-j+1,6)
         q(i,j,7)=q(i,-j+1,7)
         q(i,j,8)=q(i,-j+1,8)
         q(i,j,9)=q(i,-j+1,9)
         q(i,j,10)=q(i,-j+1,10)
         q(i,j,11)=q(i,-j+1,11)
    enddo
  enddo
endif

if(this%down_mpi .lt. 0)then
  do j=1,nbuf
    do i=1-nbuf,nx+nbuf
         q(i,ny+j,1)=q(i,ny-j,1)
         q(i,ny+j,2)=q(i,ny-j,2)
         q(i,ny+j,3)=-q(i,ny-j,3)
         q(i,ny+j,4)=q(i,ny-j,4)
         q(i,ny+j,5)=q(i,ny-j,5)
         q(i,ny+j,6)=q(i,ny-j,6)
         q(i,ny+j,7)=q(i,ny-j,7)
         q(i,ny+j,8)=q(i,ny-j,8)
         q(i,ny+j,9)=q(i,ny-j,9)
         q(i,ny+j,10)=q(i,ny-j,10)
         q(i,ny+j,11)=q(i,ny-j,11)
    enddo
  enddo
endif
end subroutine bdryRayleighTaylorInstabilityMHD2D

subroutine sourceRayleighTaylorInstabilityMHD2D(this,q,q1)
class(grid)::this
double precision,dimension(1-this%nbuf:this%nMesh(1)+this%nbuf,1-this%nbuf:this%nMesh(2)+this%nbuf,this%nvar)::q,q1
integer::nx,ny,nbuf,i,j

nbuf=this%nbuf
nx=this%nMesh(1)
ny=this%nMesh(2)

do j=1,ny
  do i=1,nx
    q1(i,j,3)=q1(i,j,3)-0.1d0*q(i,j,1)*this%dt
    q1(i,j,8)=q1(i,j,8)-q(i,j,3)*0.1d0*this%dt
  enddo
enddo

end subroutine sourceRayleighTaylorInstabilityMHD2D

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!! Rayleigh-Taylor instability HD 2D
!!!!! http://www.astro.princeton.edu/~jstone/Athena/tests/kh/kh.html
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine RayleighTaylorInstabilityHD2D(gridID)
integer::gridID
type(grid)::g1
integer::ndim,nbuf,coordType,variable(8)
integer::nMesh(2),dims(2)
double precision::leftBdry(2),rightBdry(2)
double precision::pi
logical::periods(2),reorder
double precision::t0,t1
integer::nstep,ierr

nstep=0
pi=4.d0*atan2(1.d0,1.d0)
variable=0
ndim=2
nbuf=2
coordType=1 !! Cartesian
variable(1)=1  !! rho
variable(2)=1  !! vx
variable(3)=1  !! vy
variable(4)=1  !! energy

nMesh(1)=400
nMesh(2)=600
leftBdry(1)=-0.25d0
rightBdry(1)=0.25d0
leftBdry(2)=-0.375d0
rightBdry(2)=0.375d0
dims=(/0,0/)
call MPI_DIMS_CREATE(nprocs,ndim,dims,ierr)
periods(1)=.true.
periods(2)=.false.
reorder=.true.

call g1%setTopologyMPI(ndim,dims,periods,reorder)
call g1%setMesh(nMesh,leftBdry,rightBdry,nbuf,coordType,gridID)
call g1%setVariable(variable) !! den,vx,vy,vz,bx,by,bz,ene
call g1%setMPIWindows()
call g1%setEoS(eosType=2) !! adiabatic
call g1%setadiGamma(gam=1.4d0) !! ratio of heat capacity
call g1%setCFL(CFL=0.4d0) !! courant number
call g1%setTime(fstart=0,tend=12.75d0,dtout=0.05d0)
call g1%setSlopeLimiter(limiterType=3)
call g1%setSolverType(solverType=3) !! 4=AdiHLLMHD, 5=AdiHLLDMHD
call g1%setBoundaryType(boundaryType=0) !! 3=periodic
call g1%initVariable()
call g1%setBoundary(g1%q)
call g1%exchangeBdryMPI(g1%q,g1%winq)
call g1%writeGrid()
g1%writeFlag=.false.

t0=MPI_WTIME()
do while (g1%t .lt. g1%tend)
   call g1%griddt()
   call g1%evolveGridRK2()
   if(myid .eq. 0) then
     print *,"testSuite.f03: gridID=",g1%gridID,'t=',g1%t,"dt=",g1%dt
   endif
   if(g1%writeFlag .eqv. .true.) then
      t1=MPI_WTIME()
      if(myid .eq. 0) then
        print *,"testSuiteMPI.f03: elapse time=",t1-t0,"secs"
      endif
      call g1%writeGrid()
      g1%writeFlag=.false.
      t0=MPI_WTIME()
   endif
enddo !! end do while
end subroutine RayleighTaylorInstabilityHD2D

subroutine initRayleighTaylorInstabilityHD2D(this,q)
class(grid)::this
double precision,dimension(1-this%nbuf:this%nMesh(1)+this%nbuf,1-this%nbuf:this%nMesh(2)+this%nbuf,this%nvar)::q
double precision,dimension(1-this%nbuf:this%nMesh(1)+this%nbuf)::xc
double precision,dimension(1-this%nbuf:this%nMesh(2)+this%nbuf)::yc
integer::nx,ny,nvar,nbuf,gridID,eosType
integer::i,j
double precision::gam,rho,vx,vy,p,p0,pi,rand1,rand2,A

  nx=this%nMesh(1)
  ny=this%nMesh(2)
 gam=this%adiGamma
  xc=this%xc(1)%coords
  yc=this%xc(2)%coords
  p0=2.5d0
  pi=datan2(1.d0,1.d0)*4.d0
  call init_random_seed()
  rand1=0.01d0*(rand()-0.5d0)
  rand2=0.01d0*(rand()-0.5d0)

  do j=1-nbuf,ny+nbuf
    do i=1-nbuf,nx+nbuf
       if(yc(j) .gt. 0.d0) then
         A=rand1
       else
         A=rand2
       endif
       vx=0.d0
       vy=A*(1+dcos(8.d0*pi*yc(j)/3.d0))/2.d0!0.01d0*(1.d0+dcos(4.d0*pi*xc(i)))*(1.d0+dcos(3.d0*pi*yc(j)))/4.d0
      if(yc(j) .gt. 0.d0) then
        rho=2.d0
      else
        rho=1.d0
      endif
          p=p0-0.1d0*rho*yc(j)

       q(i,j,1)=rho
       q(i,j,2)=rho*vx
       q(i,j,3)=rho*vy
       q(i,j,4)=p/(gam-1.d0)+0.5d0*rho*(vx**2+vy**2)
    enddo
  enddo
end subroutine initRayleighTaylorInstabilityHD2D


subroutine init_random_seed()
      INTEGER :: i, n, clock
      INTEGER, DIMENSION(:), ALLOCATABLE :: seed

      CALL RANDOM_SEED(size = n)
      ALLOCATE(seed(n))

      CALL SYSTEM_CLOCK(COUNT=clock)

      seed = clock + 37 * (/ (i - 1, i = 1, n) /)*(myid+1)
      CALL RANDOM_SEED(PUT = seed)

      DEALLOCATE(seed)
end subroutine init_random_seed

subroutine bdryRayleighTaylorInstabilityHD2D(this,q)
class(grid)::this
double precision,dimension(1-this%nbuf:this%nMesh(1)+this%nbuf,1-this%nbuf:this%nMesh(2)+this%nbuf,this%nvar)::q
integer::nx,ny,nbuf,i,j

nbuf=this%nbuf
nx=this%nMesh(1)
ny=this%nMesh(2)

if(this%up_mpi .lt. 0)then
  do j=1-nbuf,0
    do i=1-nbuf,nx+nbuf
         q(i,j,1)=q(i,-j+1,1)
         q(i,j,2)=q(i,-j+1,2)
         q(i,j,3)=-q(i,-j+1,3)
         q(i,j,4)=q(i,-j+1,4)
    enddo
  enddo
endif

if(this%down_mpi .lt. 0)then
  do j=1,nbuf
    do i=1-nbuf,nx+nbuf
         q(i,ny+j,1)=q(i,ny-j,1)
         q(i,ny+j,2)=q(i,ny-j,2)
         q(i,ny+j,3)=-q(i,ny-j,3)
         q(i,ny+j,4)=q(i,ny-j,4)
    enddo
  enddo
endif
end subroutine bdryRayleighTaylorInstabilityHD2D

subroutine sourceRayleighTaylorInstabilityHD2D(this,q,q1)
class(grid)::this
double precision,dimension(1-this%nbuf:this%nMesh(1)+this%nbuf,1-this%nbuf:this%nMesh(2)+this%nbuf,this%nvar)::q,q1
integer::nx,ny,nbuf,i,j

nbuf=this%nbuf
nx=this%nMesh(1)
ny=this%nMesh(2)

do j=1,ny
  do i=1,nx
    q1(i,j,3)=q1(i,j,3)-0.1d0*q(i,j,1)*this%dt
    q1(i,j,4)=q1(i,j,4)-q(i,j,3)*0.1d0*this%dt
  enddo
enddo

end subroutine sourceRayleighTaylorInstabilityHD2D

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!! Kelvin-Helmholtz instability MHD 2D
!!!!! http://www.astro.princeton.edu/~jstone/Athena/tests/kh/kh.html
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine KelvinHelmholtzInstabilityMHD2D(gridID)
integer::gridID
type(grid)::g1
integer::ndim,nbuf,coordType,variable(8)
integer::nMesh(2),dims(2)
double precision::leftBdry(2),rightBdry(2)
double precision::pi
logical::periods(2),reorder
double precision::t0,t1
integer::nstep,ierr

nstep=0
pi=4.d0*atan2(1.d0,1.d0)
variable=0
ndim=2
nbuf=2
coordType=1    !! Cartesian
variable(1)=1  !! rho
variable(2)=1  !! vx
variable(3)=1  !! vy
variable(4)=1  !! vz
variable(5)=1  !! bx
variable(6)=1  !! by
variable(7)=1  !! bz
variable(8)=1  !! energy

nMesh(1)=512
nMesh(2)=512
leftBdry(1)=-0.5d0
rightBdry(1)=0.5d0
leftBdry(2)=-0.5d0
rightBdry(2)=0.5d0
dims=(/0,0/)
call MPI_DIMS_CREATE(nprocs,ndim,dims,ierr)
periods(1)=.true.
periods(2)=.true.
reorder=.true.

call g1%setTopologyMPI(ndim,dims,periods,reorder)
call g1%setMesh(nMesh,leftBdry,rightBdry,nbuf,coordType,gridID)
call g1%setVariable(variable) !! den,vx,vy,vz,bx,by,bz,ene
call g1%setMPIWindows()
call g1%setEoS(eosType=2) !! adiabatic
call g1%setadiGamma(gam=1.4d0) !! ratio of heat capacity
call g1%setCFL(CFL=0.4d0) !! courant number
call g1%setTime(fstart=0,tend=3.d0,dtout=0.2d0)
call g1%setSlopeLimiter(limiterType=3)
call g1%setSolverType(solverType=5) !! 4=AdiHLLMHD, 5=AdiHLLDMHD
call g1%setBoundaryType(boundaryType=3) !! 3=periodic
call g1%initVariable()
call g1%exchangeBdryMPI(g1%q,g1%winq)
call g1%setBoundary(g1%q)
call g1%writeGrid()
g1%writeFlag=.false.

t0=MPI_WTIME()
do while (g1%t .lt. g1%tend)
   call g1%griddt()
   call g1%evolveGridRK2()
   if(myid .eq. 0) then
     print *,"testSuite.f03: gridID=",g1%gridID,'t=',g1%t,"dt=",g1%dt
   endif
   if(g1%writeFlag .eqv. .true.) then
      t1=MPI_WTIME()
      if(myid .eq. 0) then
        print *,"testSuiteMPI.f03: elapse time=",t1-t0,"secs"
      endif
      call g1%writeGrid()
      g1%writeFlag=.false.
      t0=MPI_WTIME()
   endif
enddo !! end do while
end subroutine KelvinHelmholtzInstabilityMHD2D

subroutine initKelvinHelmholtzInstabilityMHD2D(this,q)
class(grid)::this
double precision,dimension(1-this%nbuf:this%nMesh(1)+this%nbuf,1-this%nbuf:this%nMesh(2)+this%nbuf,this%nvar)::q
double precision,dimension(1-this%nbuf:this%nMesh(1)+this%nbuf)::xc
double precision,dimension(1-this%nbuf:this%nMesh(2)+this%nbuf)::yc
integer::nx,ny,nvar,nbuf,gridID,eosType
integer::i,j
double precision::gam,rho,vx,vy,p,b0,vz


  nx=this%nMesh(1)
  ny=this%nMesh(2)
 gam=this%adiGamma
  xc=this%xc(1)%coords
  yc=this%xc(2)%coords
  b0=0.25d0

  do j=1,ny
    do i=1,nx

       vy=0.d0
      if(dabs(yc(j)) .gt. 0.25d0) then
         vx=-0.5d0
        rho=1.d0
      else
         vx=0.5d0
        rho=2.d0
      endif
         p=2.5d0
        vz=0.d0

       q(i,j,1)=rho
       q(i,j,2)=rho*(vx+(rand()-0.5d0)*0.01d0)
       q(i,j,3)=rho*(vy+(rand()-0.5d0)*0.01d0)
       q(i,j,4)=rho*vz
       q(i,j,5)=b0
       q(i,j,6)=0.d0
       q(i,j,7)=0.d0
       q(i,j,8)=p/(gam-1.d0)+0.5d0*rho*(vx**2+vy**2+vz**2)+0.5d0*(b0**2)
       q(i,j,9)=b0
       q(i,j,10)=0.d0
       q(i,j,11)=0.d0
    enddo
  enddo
end subroutine initKelvinHelmholtzInstabilityMHD2D

subroutine bdryKelvinHelmholtzInstabilityMHD2D(this,q)
class(grid)::this
double precision,dimension(1-this%nbuf:this%nMesh(1)+this%nbuf,1-this%nbuf:this%nMesh(2)+this%nbuf,this%nvar)::q
integer::nx,ny,nbuf,i,j

nbuf=this%nbuf
nx=this%nMesh(1)
ny=this%nMesh(2)

if(this%boundaryType .eq. 3) then !! periodic boundary condition
endif
end subroutine bdryKelvinHelmholtzInstabilityMHD2D


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!! Kelvin-Helmholtz instability HD 2D
!!!!! http://www.astro.princeton.edu/~jstone/Athena/tests/kh/kh.html
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine KelvinHelmholtzInstabilityHD2D(gridID)
integer::gridID
type(grid)::g1
integer::ndim,nbuf,coordType,variable(8)
integer::nMesh(2),dims(2)
double precision::leftBdry(2), rightBdry(2)
logical::periods(2),reorder
double precision::t0,t1
integer::nstep,ierr,Nx,Ny

nstep=0
variable=0
ndim=2
nbuf=2
coordType=1    !! cartesian coordinates
variable(1)=1  !! density
variable(2)=1  !! momentum in x
variable(3)=1  !! momentum in y
variable(8)=1  !! energy

Nx=512
Ny=512
nMesh(1)=Nx
nMesh(2)=Ny
leftBdry(1)=-0.5d0
leftBdry(2)=-0.5d0
rightBdry(1)=0.5d0
rightBdry(2)=0.5d0
dims=(/0,0/)
call MPI_DIMS_CREATE(nprocs,ndim,dims,ierr)
periods(1)=.true.
periods(2)=.true.
reorder=.true.

call g1%setTopologyMPI(ndim,dims,periods,reorder)
call g1%setGridID(gridID=gridID)
call g1%setTime(fstart=0,tend=3.d0,dtout=0.01d0)
if (g1%fstart .eq. 0) then
  call g1%setMesh(nMesh,leftBdry,rightBdry,nbuf,coordType,gridID)
  call g1%setVariable(variable) !! den,vx,vy,vz,bx,by,bz,ene
  call g1%setMPIWindows()
  call g1%setEoS(eosType=2) !! adiabatic
  call g1%setadiGamma(gam=1.4d0) !! ratio of heat capacity
  call g1%setCFL(CFL=0.4d0) !! courant number
  call g1%setSlopeLimiter(limiterType=3)  !! 3=minmod
  call g1%setSolverType(solverType=3) !! 3=HLLC, 4=AdiHLLMHD, 5=AdiHLLDMHD
  call g1%setBoundaryType(boundaryType=3) !! 1=zero gradient, 3=periodic
  call g1%initVariable()
  call g1%exchangeBdryMPI(g1%q,g1%winq)
  call g1%setBoundary(g1%q)
  call g1%writeGrid()
  call g1%writeGrid_HD_vtk()
endif
g1%writeFlag=.false.

t0=MPI_WTIME()
do while (g1%t .lt. g1%tend)
   nstep=nstep+1
   call g1%griddt()
   call g1%evolveGridRK2()
   if(myid .eq. 0 ) then
     print *,"testSuite.f03: gridID=",g1%gridID,'t=',g1%t,"dt=",g1%dt
   endif
   if(g1%writeFlag .eqv. .true.) then
      t1=MPI_WTIME()
      if(myid .eq. 0) then
        print *,"testSuiteMPI.f03: elapse time=",t1-t0,"secs"
      endif
      call g1%writeGrid()
      call g1%writeGrid_HD_vtk()
      g1%writeFlag=.false.
      t0=MPI_WTIME()
   endif
enddo !! end do while
end subroutine KelvinHelmholtzInstabilityHD2D

subroutine initKelvinHelmholtzInstabilityHD2D(this,q)
class(grid)::this
double precision,dimension(1-this%nbuf:this%nMesh(1)+this%nbuf,1-this%nbuf:this%nMesh(2)+this%nbuf,this%nvar)::q
double precision,dimension(1-this%nbuf:this%nMesh(1)+this%nbuf)::xc
double precision,dimension(1-this%nbuf:this%nMesh(2)+this%nbuf)::yc
integer::nx,ny,nvar,nbuf,gridID,eosType
integer::i,j
double precision::gam,rho,vx,vy,p


  nx=this%nMesh(1)
  ny=this%nMesh(2)
 gam=this%adiGamma
  xc=this%xc(1)%coords
  yc=this%xc(2)%coords

  do j=1,ny
    do i=1,nx
       vy=0.d0
      if(dabs(yc(j)) .gt. 0.25d0) then
         vx=-0.5d0
        rho=1.d0
      else
         vx=0.5d0
        rho=2.d0 
      endif
         p=2.5d0

       q(i,j,1)=rho
       q(i,j,2)=rho*(vx+(rand()-0.5d0)*0.01d0)
       q(i,j,3)=rho*(vy+(rand()-0.5d0)*0.01d0)
       q(i,j,4)=p/(gam-1.d0)+0.5d0*rho*(vx**2+vy**2)
    enddo
  enddo
end subroutine initKelvinHelmholtzInstabilityHD2D

subroutine bdryKelvinHelmholtzInstabilityHD2D(this,q)
class(grid)::this
double precision,dimension(1-this%nbuf:this%nMesh(1)+this%nbuf,1-this%nbuf:this%nMesh(2)+this%nbuf,this%nvar)::q
integer::nx,ny,nbuf,i,j

nbuf=this%nbuf
nx=this%nMesh(1)
ny=this%nMesh(2)

if(this%boundaryType .eq. 3) then !! periodic boundary condition
endif

end subroutine bdryKelvinHelmholtzInstabilityHD2D
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!! force balance test MHD for polar 2D
!!!!! Skinner & Ostriker, 2010, ApJS, 188, 290
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine forceBalanceMHDPolar2D(gridID)
integer::gridID
type(grid)::g1
integer::ndim,nbuf,coordType,variable(8)
integer::nMesh(2),dims(2)
double precision::leftBdry(2),rightBdry(2)
double precision::pi
logical::periods(2),reorder
double precision::t0,t1
integer::nstep,ierr

nstep=0
pi=4.d0*atan2(1.d0,1.d0)
variable=0
ndim=2
nbuf=2
coordType=3 !! polar coordinates, uniform in r
variable(1)=1  !! rho
variable(2)=1  !! vx
variable(3)=1  !! vy
variable(4)=1  !! vz
variable(5)=1  !! bx
variable(6)=1  !! by
variable(7)=1  !! bz
variable(8)=1  !! energy

nMesh(1)=200*8
nMesh(2)=150*8
leftBdry(1)=1.d0
rightBdry(1)=2.d0
leftBdry(2)=-0.25d0
rightBdry(2)=0.25d0
dims=(/0,0/)
call MPI_DIMS_CREATE(nprocs,ndim,dims,ierr)
periods(1)=.false.
periods(2)=.true.
reorder=.true.

call g1%setTopologyMPI(ndim,dims,periods,reorder)
call g1%setMesh(nMesh,leftBdry,rightBdry,nbuf,coordType,gridID)
call g1%setVariable(variable) !! den,vx,vy,vz,bx,by,bz,ene
call g1%setMPIWindows()
call g1%setEoS(eosType=2) !! adiabatic
call g1%setadiGamma(gam=5.d0/3.d0) !! ratio of heat capacity
call g1%setCFL(CFL=0.4d0) !! courant number
call g1%setTime(fstart=0,tend=10.d0,dtout=1.d0)
call g1%setSlopeLimiter(limiterType=3)
call g1%setSolverType(solverType=5) !! 4=AdiHLLMHD, 5=AdiHLLDMHD
call g1%setBoundaryType(boundaryType=0) !! 3=periodic
call g1%initVariable()
call g1%exchangeBdryMPI(g1%q,g1%winq)
call g1%setBoundary(g1%q)
call g1%writeGrid()
g1%writeFlag=.false.

t0=MPI_WTIME()
do while (g1%t .lt. g1%tend)
   call g1%griddt()
   call g1%evolveGridRK2()
   if(myid .eq. 0) then
     print *,"testSuite.f03: gridID=",g1%gridID,'t=',g1%t,"dt=",g1%dt
   endif
   if(g1%writeFlag .eqv. .true.) then
      t1=MPI_WTIME()
      if(myid .eq. 0) then
        print *,"testSuiteMPI.f03: elapse time=",t1-t0,"secs"
      endif
      call g1%writeGrid()
      g1%writeFlag=.false.
      t0=MPI_WTIME()
   endif
enddo !! end do while

end subroutine forceBalanceMHDPolar2D

subroutine initForceBalanceMHDPolar2D(this,q)
class(grid)::this
double precision,dimension(1-this%nbuf:this%nMesh(1)+this%nbuf,1-this%nbuf:this%nMesh(2)+this%nbuf,this%nvar)::q
double precision,dimension(1-this%nbuf:this%nMesh(1)+this%nbuf)::xc,xl,xr
double precision,dimension(1-this%nbuf:this%nMesh(2)+this%nbuf)::yc,yl,yr
integer::nx,ny,nvar,nbuf,gridID,eosType
integer::i,j
double precision::gam,rho,vx,vy,vz,bxl,byl,bzl,bxr,byr,bzr,p
double precision::bxc,byc,bzc
double precision::pi,b0
double precision::rc,x,y,z,dr,dphi


pi=4.d0*atan2(1.d0,1.d0)
b0=1.d0!dsqrt(4.d0*pi/2.d0)
xc=this%xc(1)%coords
yc=this%xc(2)%coords
xl=this%xl(1)%coords
yl=this%xl(2)%coords
xr=this%xr(1)%coords
yr=this%xr(2)%coords
nx=this%nMesh(1)
ny=this%nMesh(2)
gam=this%adiGamma

do j=1,ny
  do i=1,nx
    rho=1.d0
    vx=0.d0
    vy=0.d0
    vz=0.d0
   bxl=0.d0
   byl=b0/xc(i)
   bzl=0.d0
   bxr=0.d0
   byr=b0/xc(i)
   bzr=0.d0
     p=1.d0

    bxc=0.d0
    byc=0.5d0*(byl+byr)
    bzc=0.d0

    q(i,j,1)=rho
    q(i,j,2)=rho*vx
    q(i,j,3)=rho*vy
    q(i,j,4)=rho*vz
    q(i,j,5)=bxl
    q(i,j,6)=byl
    q(i,j,7)=bzl
    q(i,j,8)=p/(gam-1.d0)+0.5d0*rho*(vx**2.d0+vy**2.d0+vz**2.d0)+0.5d0*(bxc**2.d0+byc**2.d0+bzc**2.d0)
    q(i,j,9)=bxr
    q(i,j,10)=byr
    q(i,j,11)=bzr
  enddo
enddo

end subroutine initForceBalanceMHDPolar2D

subroutine bdryForceBalanceMHDPolar2D(this,q)
class(grid)::this
double precision,dimension(1-this%nbuf:this%nMesh(1)+this%nbuf,1-this%nbuf:this%nMesh(2)+this%nbuf,this%nvar)::q
double precision,dimension(1-this%nbuf:this%nMesh(1)+this%nbuf)::xc,xl,xr
double precision,dimension(1-this%nbuf:this%nMesh(2)+this%nbuf)::yc,yl,yr
integer::nx,ny,nbuf,i,j
double precision::rho,vx,vy,vz,bxl,byl,bzl,bxr,byr,bzr,bxc,byc,bzc,p,gam,b0

nbuf=this%nbuf
b0=1.d0!dsqrt(4.d0*pi/2.d0)
xc=this%xc(1)%coords
yc=this%xc(2)%coords
xl=this%xl(1)%coords
yl=this%xl(2)%coords
xr=this%xr(1)%coords
yr=this%xr(2)%coords
nx=this%nMesh(1)
ny=this%nMesh(2)
gam=this%adiGamma

!!!!!!!! specify boundary condition in r-boundary
if(this%left_mpi .lt. 0) then 
do j=1-nbuf,ny+nbuf
  do i=1,nbuf
    rho=1.d0
    vx=0.d0
    vy=0.d0
    vz=0.d0
   bxl=0.d0
   byl=b0/xc(1-i)
   bzl=0.d0
   bxr=0.d0
   byr=b0/xc(1-i)
   bzr=0.d0
     p=1.d0

    bxc=0.d0
    byc=0.5d0*(byl+byr)
    bzc=0.d0

    q(1-i,j,1)=rho
    q(1-i,j,2)=rho*vx
    q(1-i,j,3)=rho*vy
    q(1-i,j,4)=rho*vz
    q(1-i,j,5)=bxl
    q(1-i,j,6)=byl
    q(1-i,j,7)=bzl
    q(1-i,j,8)=p/(gam-1.d0)+0.5d0*rho*(vx**2.d0+vy**2.d0+vz**2.d0)+0.5d0*(bxc**2.d0+byc**2.d0+bzc**2.d0)
    q(1-i,j,9)=bxr
    q(1-i,j,10)=byr
    q(1-i,j,11)=bzr    
  enddo
enddo
endif

if(this%right_mpi .lt. 0) then
do j=1-nbuf,ny+nbuf
  do i=1,nbuf
    rho=1.d0
    vx=0.d0
    vy=0.d0
    vz=0.d0
   bxl=0.d0
   byl=b0/xc(nx+i)
   bzl=0.d0
   bxr=0.d0
   byr=b0/xc(nx+i)
   bzr=0.d0
     p=1.d0

    bxc=0.d0
    byc=0.5d0*(byl+byr)
    bzc=0.d0

    q(nx+i,j,1)=rho
    q(nx+i,j,2)=rho*vx
    q(nx+i,j,3)=rho*vy
    q(nx+i,j,4)=rho*vz
    q(nx+i,j,5)=bxl
    q(nx+i,j,6)=byl
    q(nx+i,j,7)=bzl
    q(nx+i,j,8)=p/(gam-1.d0)+0.5d0*rho*(vx**2.d0+vy**2.d0+vz**2.d0)+0.5d0*(bxc**2.d0+byc**2.d0+bzc**2.d0)
    q(nx+i,j,9)=bxr
    q(nx+i,j,10)=byr
    q(nx+i,j,11)=bzr
  enddo
enddo
endif

end subroutine bdryForceBalanceMHDPolar2D
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!! Spherical Blast Wave Test (Zachary, Malagoli & Colella,
!!!!! (SIAM J. Sci. Comp, 1994, 15, 263)
!!!!! Balsara, D., & Spicer, D., JCP 149, 270 (1999);
!!!!! Londrillo, P. & Del Zanna, L., ApJ 530, 508 (2000).
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine sphericalBlastWaveMHDPolar2D(gridID)
integer::gridID
type(grid)::g1
integer::ndim,nbuf,coordType,variable(8)
integer::nMesh(2),dims(2)
double precision::leftBdry(2),rightBdry(2)
double precision::pi
logical::periods(2),reorder
double precision::t0,t1
integer::nstep,ierr

nstep=0
pi=4.d0*atan2(1.d0,1.d0)
variable=0
ndim=2
nbuf=2
coordType=3 !! polar coordinates, uniform in r
variable(1)=1  !! rho
variable(2)=1  !! vx
variable(3)=1  !! vy
variable(4)=1  !! vz
variable(5)=1  !! bx
variable(6)=1  !! by
variable(7)=1  !! bz
variable(8)=1  !! energy

nMesh(1)=400
nMesh(2)=600
leftBdry(1)=1.d0
rightBdry(1)=2.d0
leftBdry(2)=-0.5d0
rightBdry(2)=0.5d0
dims=(/0,0/)
call MPI_DIMS_CREATE(nprocs,ndim,dims,ierr)
periods(1)=.true.
periods(2)=.true.
reorder=.true.

call g1%setTopologyMPI(ndim,dims,periods,reorder)
call g1%setMesh(nMesh,leftBdry,rightBdry,nbuf,coordType,gridID)
call g1%setVariable(variable) !! den,vx,vy,vz,bx,by,bz,ene
call g1%setMPIWindows()
call g1%setEoS(eosType=2) !! adiabatic
call g1%setadiGamma(gam=5.d0/3.d0) !! ratio of heat capacity
call g1%setCFL(CFL=0.4d0) !! courant number
call g1%setTime(fstart=0,tend=0.1d0,dtout=0.1d0)
call g1%setSlopeLimiter(limiterType=3)
call g1%setSolverType(solverType=5) !! 4=AdiHLLMHD, 5=AdiHLLDMHD
call g1%setBoundaryType(boundaryType=3) !! 3=periodic
call g1%initVariable()
call g1%exchangeBdryMPI(g1%q,g1%winq)
call g1%setBoundary(g1%q)
call g1%writeGrid()
g1%writeFlag=.false.

t0=MPI_WTIME()
do while (g1%t .lt. g1%tend)
   call g1%griddt()
   call g1%evolveGridRK2()
   if(myid .eq. 0) then
     print *,"testSuite.f03: gridID=",g1%gridID,'t=',g1%t,"dt=",g1%dt
   endif
   if(g1%writeFlag .eqv. .true.) then
      t1=MPI_WTIME()
      if(myid .eq. 0) then
        print *,"testSuiteMPI.f03: elapse time=",t1-t0,"secs"
      endif
      call g1%writeGrid()
      g1%writeFlag=.false.
      t0=MPI_WTIME()
   endif
enddo !! end do while

end subroutine sphericalBlastWaveMHDPolar2D

subroutine initSphericalBlastWaveMHDPolar2D(this,q)
class(grid)::this
double precision,dimension(1-this%nbuf:this%nMesh(1)+this%nbuf,1-this%nbuf:this%nMesh(2)+this%nbuf,this%nvar)::q
double precision,dimension(1-this%nbuf:this%nMesh(1)+this%nbuf)::xc,xl,xr
double precision,dimension(1-this%nbuf:this%nMesh(2)+this%nbuf)::yc,yl,yr
integer::nx,ny,nvar,nbuf,gridID,eosType
integer::i,j
double precision::gam,rho,vx,vy,vz,bxl,byl,bzl,bxr,byr,bzr,p
double precision::bxc,byc,bzc
double precision::pi,b0
double precision::rc,x,y,z,dr,dphi


pi=4.d0*atan2(1.d0,1.d0)
b0=1.d0!dsqrt(4.d0*pi/2.d0)
xc=this%xc(1)%coords
yc=this%xc(2)%coords
xl=this%xl(1)%coords
yl=this%xl(2)%coords
xr=this%xr(1)%coords
yr=this%xr(2)%coords
nx=this%nMesh(1)
ny=this%nMesh(2)
gam=this%adiGamma

do j=1,ny
  do i=1,nx
    x=xc(i)*dcos(yc(j))
    y=xc(i)*dsin(yc(j))

   rc=dsqrt((x-1.5d0)**2+y**2)

    rho=1.d0
    vx=0.d0
    vy=0.d0
    vz=0.d0
   bxl=b0/xl(i)*xc(i)
   byl=b0
   bzl=0.d0
   bxr=b0/xr(i)*xc(i)
   byr=byl
   bzr=bzl
     p=0.1d0
    if(rc .lt. 0.1d0) then
     p=10.d0
    endif

    bxc=0.5d0*(xl(i)*bxl+xr(i)*bxr)/xc(i)
    byc=0.5d0*(byl+byr)
    bzc=0.5d0*(bzl+bzr)

    q(i,j,1)=rho
    q(i,j,2)=rho*vx
    q(i,j,3)=rho*vy
    q(i,j,4)=rho*vz
    q(i,j,5)=bxl
    q(i,j,6)=byl
    q(i,j,7)=bzl
    q(i,j,8)=p/(gam-1.d0)+0.5d0*rho*(vx**2.d0+vy**2.d0+vz**2.d0)+0.5d0*(bxc**2.d0+byc**2.d0+bzc**2.d0)
    q(i,j,9)=bxr
    q(i,j,10)=byr
    q(i,j,11)=bzr
  enddo
enddo

end subroutine initSphericalBlastWaveMHDPolar2D

subroutine bdrySphericalBlastWaveMHDPolar2D(this,q)
class(grid)::this
double precision,dimension(1-this%nbuf:this%nMesh(1)+this%nbuf,1-this%nbuf:this%nMesh(2)+this%nbuf,this%nvar)::q
integer::nx,ny,nbuf,i,j

nbuf=this%nbuf
nx=this%nMesh(1)
ny=this%nMesh(2)

if(this%boundaryType .eq. 3) then !! periodic boundary condition
endif

end subroutine bdrySphericalBlastWaveMHDPolar2D
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!! Blast Wave in 2D Cartesian coordinates for isothermal gas
!!!!! Skinner & Ostriker, 2010, ApJS, 188, 290
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine HDBlastWaveCartIso2D(gridID)
integer::gridID
type(grid)::g1
integer::ndim,nbuf,coordType,variable(8)
integer::nMesh(2),dims(2)
double precision::leftBdry(2), rightBdry(2)
logical::periods(2),reorder
double precision::t0,t1
integer::nstep,ierr,Nx,Ny

nstep=0
variable=0
ndim=2
nbuf=2
coordType=1    !! cartesian coordinates
variable(1)=1  !! density
variable(2)=1  !! momentum in x
variable(3)=1  !! momentum in y

Nx=800
Ny=1200
nMesh(1)=Nx
nMesh(2)=Ny
leftBdry(1)=-0.5d0
leftBdry(2)=-0.75d0
rightBdry(1)=0.5d0
rightBdry(2)=0.75d0
dims=(/0,0/)
call MPI_DIMS_CREATE(nprocs,ndim,dims,ierr)
periods(1)=.false.
periods(2)=.false.
reorder=.true.

call g1%setTopologyMPI(ndim,dims,periods,reorder)
call g1%setGridID(gridID=gridID)
call g1%setTime(fstart=0,tend=0.15d0,dtout=0.05d0)
if (g1%fstart .eq. 0) then
  call g1%setMesh(nMesh,leftBdry,rightBdry,nbuf,coordType,gridID)
  call g1%setVariable(variable) !! den,vx,vy,vz,bx,by,bz,ene
  call g1%setMPIWindows()
  call g1%setEoS(eosType=1) !! isothermal
  call g1%setSoundSpeed(snd=1.d0)
  call g1%setCFL(CFL=0.3d0) !! courant number
  call g1%setSlopeLimiter(limiterType=3)  !! 3=minmod
  call g1%setSolverType(solverType=1) !! 1=exact Riemann solver 
  call g1%setBoundaryType(boundaryType=1) !! 1=zero gradient
  call g1%initVariable()
  call g1%exchangeBdryMPI(g1%q,g1%winq)
  call g1%setBoundary(g1%q)
  call g1%writeGrid()
endif
g1%writeFlag=.false.

t0=MPI_WTIME()
do while (g1%t .lt. g1%tend)
   nstep=nstep+1
   call g1%griddt()
   call g1%evolveGridRK2()
   if(myid .eq. 0 ) then
     print *,"testSuite.f03: gridID=",g1%gridID,'t=',g1%t,"dt=",g1%dt
   endif
   if(g1%writeFlag .eqv. .true.) then
      t1=MPI_WTIME()
      if(myid .eq. 0) then
        print *,"testSuiteMPI.f03: elapse time=",t1-t0,"secs"
      endif
      call g1%writeGrid()
      g1%writeFlag=.false.
      t0=MPI_WTIME()
   endif
enddo !! end do while

end subroutine HDBlastWaveCartIso2D

subroutine initHDBlastWaveCartIso2D(this,q)
class(grid)::this
double precision,dimension(1-this%nbuf:this%nMesh(1)+this%nbuf,1-this%nbuf:this%nMesh(2)+this%nbuf,this%nvar)::q
double precision,dimension(1-this%nbuf:this%nMesh(1)+this%nbuf)::xc
double precision,dimension(1-this%nbuf:this%nMesh(2)+this%nbuf)::yc
integer::nx,ny,nvar,nbuf,gridID,eosType
integer::i,j
double precision::gam,rho,vx,vy,p,r,R0


  nx=this%nMesh(1)
  ny=this%nMesh(2)
 gam=this%adiGamma
  xc=this%xc(1)%coords
  yc=this%xc(2)%coords
  R0=0.1d0

  do j=1,ny
    do i=1,nx
       vx=0.d0
       vy=0.d0
       r=dsqrt(xc(i)**2+yc(j)**2)

       if(r .le. R0) then
         rho=10.d0
       else 
         rho=0.1d0
       endif

       q(i,j,1)=rho
       q(i,j,2)=rho*vx
       q(i,j,3)=rho*vy
    enddo
  enddo
end subroutine initHDBlastWaveCartIso2D


subroutine bdryHDBlastWaveCartIso2D(this,q)
class(grid)::this
double precision,dimension(1-this%nbuf:this%nMesh(1)+this%nbuf,1-this%nbuf:this%nMesh(2)+this%nbuf,this%nvar)::q
integer::nx,ny,nbuf,i,j

nbuf=this%nbuf
nx=this%nMesh(1)
ny=this%nMesh(2)

if(this%boundaryType .eq. 1) then !! zero gradient
!!!!!! left boundary
  if(this%left_mpi .lt. 0) then
    do i=1,nbuf
       q(1-i,:,:)=q(1,:,:)
    enddo
  endif
!!!!!! right boundary
  if(this%right_mpi .lt. 0) then
    do i=1,nbuf
       q(nx+i,:,:)=q(nx,:,:)
    enddo
  endif
!!!!!! up boundary
  if(this%up_mpi .lt. 0) then
    do j=1,nbuf
       q(:,1-j,:)=q(:,1,:)
    enddo
  endif
!!!!!! down boundary
  if(this%down_mpi .lt. 0) then
    do j=1,nbuf
       q(:,ny+j,:)=q(:,ny,:)
    enddo
  endif
endif
end subroutine bdryHDBlastWaveCartIso2D
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!! Blast Wave in 2D polar coordinates for isothermal gas
!!!!! Skinner & Ostriker, 2010, ApJS, 188, 290
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine HDBlastWavePolarIso2D(gridID)
integer::gridID
type(grid)::g1
integer::ndim,nbuf,coordType,variable(8)
integer::nMesh(2),dims(2)
double precision::leftBdry(2), rightBdry(2)
logical::periods(2),reorder
double precision::t0,t1
integer::nstep,ierr,Nr,Nphi

nstep=0
variable=0
ndim=2
nbuf=2
coordType=3   !! 3=uniform mesh in r, 2= log in r
variable(1)=1  !! density
variable(2)=1  !! momentum in r
variable(3)=1  !! momentum in phi

Nr=800
Nphi=1200
nMesh(1)=Nr
nMesh(2)=Nphi
leftBdry(1)=1.d0
leftBdry(2)=-0.5d0
rightBdry(1)=2.d0
rightBdry(2)=0.5d0
dims=(/0,0/)
call MPI_DIMS_CREATE(nprocs,ndim,dims,ierr)
periods(1)=.false.
periods(2)=.false.
reorder=.true.

call g1%setTopologyMPI(ndim,dims,periods,reorder)
call g1%setGridID(gridID=gridID)
call g1%setTime(fstart=0,tend=0.15d0,dtout=0.05d0)
if (g1%fstart .eq. 0) then
  call g1%setMesh(nMesh,leftBdry,rightBdry,nbuf,coordType,gridID)
  call g1%setVariable(variable) !! den,vx,vy,vz,bx,by,bz,ene
  call g1%setMPIWindows()
  call g1%setEoS(eosType=1) !! isothermal
  call g1%setSoundSpeed(snd=1.d0) !! isothermal sound speed
  call g1%setCFL(CFL=0.4d0) !! courant number
  call g1%setSlopeLimiter(limiterType=3)  !! 3=minmod
  call g1%setSolverType(solverType=1) !! 1=exact Riemann solver
  call g1%setBoundaryType(boundaryType=1) !! 1=zero gradient
  call g1%initVariable()
  call g1%exchangeBdryMPI(g1%q,g1%winq)
  call g1%setBoundary(g1%q)
  call g1%writeGrid()
endif
g1%writeFlag=.false.

t0=MPI_WTIME()
do while (g1%t .lt. g1%tend)
   nstep=nstep+1
   call g1%griddt()
   call g1%evolveGridRK2()
   if(myid .eq. 0 ) then
     print *,"testSuite.f03: gridID=",g1%gridID,'t=',g1%t,"dt=",g1%dt
   endif
   if(g1%writeFlag .eqv. .true.) then
      t1=MPI_WTIME()
      if(myid .eq. 0) then
        print *,"testSuiteMPI.f03: elapse time=",t1-t0,"secs"
      endif
      call g1%writeGrid()
      g1%writeFlag=.false.
      t0=MPI_WTIME()
   endif
enddo !! end do while

end subroutine HDBlastWavePolarIso2D

subroutine initHDBlastWavePolarIso2D(this,q)
class(grid)::this
double precision,dimension(1-this%nbuf:this%nMesh(1)+this%nbuf,1-this%nbuf:this%nMesh(2)+this%nbuf,this%nvar)::q
double precision,dimension(1-this%nbuf:this%nMesh(1)+this%nbuf)::rc
double precision,dimension(1-this%nbuf:this%nMesh(2)+this%nbuf)::phic
integer::nr,nphi,nvar,nbuf,gridID,eosType
integer::i,j
double precision::gam,rho,vr,vphi,p
double precision::R0,r,xc,yc,rmid

  rc=this%xc(1)%coords
phic=this%xc(2)%coords


  nr=this%nMesh(1)
nphi=this%nMesh(2)
 gam=this%adiGamma
  R0=0.1d0
  rmid=1.5d0

  do j=1,nphi
    do i=1,nr
       vr=0.d0
       vphi=0.d0
       xc=rc(i)*dcos(phic(j))
       yc=rc(i)*dsin(phic(j))
       r=dsqrt((xc-rmid)**2+yc**2)

       if(r .le. R0) then
         rho=10.d0
       else 
         rho=0.1d0
       endif

       q(i,j,1)=rho
       q(i,j,2)=rho*vr
       q(i,j,3)=rho*vphi
    enddo
  enddo
end subroutine initHDBlastWavePolarIso2D


subroutine bdryHDBlastWavePolarIso2D(this,q)
class(grid)::this
double precision,dimension(1-this%nbuf:this%nMesh(1)+this%nbuf,1-this%nbuf:this%nMesh(2)+this%nbuf,this%nvar)::q
integer::nx,ny,nbuf,i,j

nbuf=this%nbuf
nx=this%nMesh(1)
ny=this%nMesh(2)

if(this%boundaryType .eq. 1) then !! zero gradient
!!!!!! left boundary
  if(this%left_mpi .lt. 0) then
    do i=1,nbuf
       q(1-i,:,:)=q(1,:,:)
    enddo
  endif
!!!!!! right boundary
  if(this%right_mpi .lt. 0) then
    do i=1,nbuf
       q(nx+i,:,:)=q(nx,:,:)
    enddo
  endif
!!!!!! up boundary
  if(this%up_mpi .lt. 0) then
    do j=1,nbuf
       q(:,1-j,:)=q(:,1,:)
    enddo
  endif
!!!!!! down boundary
  if(this%down_mpi .lt. 0) then
    do j=1,nbuf
       q(:,ny+j,:)=q(:,ny,:)
    enddo
  endif
endif
end subroutine bdryHDBlastWavePolarIso2D
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!! Blast Wave in 2D Cartesian coordinates
!!!!! Skinner & Ostriker, 2010, ApJS, 188, 290
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine HDBlastWaveCart2D(gridID)
integer::gridID
type(grid)::g1
integer::ndim,nbuf,coordType,variable(8)
integer::nMesh(2),dims(2)
double precision::leftBdry(2), rightBdry(2)
logical::periods(2),reorder
double precision::t0,t1
integer::nstep,ierr,Nx,Ny

nstep=0
variable=0
ndim=2
nbuf=2
coordType=1    !! cartesian coordinates
variable(1)=1  !! density
variable(2)=1  !! momentum in x
variable(3)=1  !! momentum in y
variable(8)=1  !! energy

Nx=800
Ny=1200
nMesh(1)=Nx
nMesh(2)=Ny
leftBdry(1)=-0.5d0
leftBdry(2)=-0.75d0
rightBdry(1)=0.5d0
rightBdry(2)=0.75d0
dims=(/0,0/)
call MPI_DIMS_CREATE(nprocs,ndim,dims,ierr)
periods(1)=.false.
periods(2)=.false.
reorder=.true.

call g1%setTopologyMPI(ndim,dims,periods,reorder)
call g1%setGridID(gridID=gridID)
call g1%setTime(fstart=0,tend=0.2d0,dtout=0.05d0)
if (g1%fstart .eq. 0) then
  call g1%setMesh(nMesh,leftBdry,rightBdry,nbuf,coordType,gridID)
  call g1%setVariable(variable) !! den,vx,vy,vz,bx,by,bz,ene
  call g1%setMPIWindows()
  call g1%setEoS(eosType=2) !! adiabatic
  call g1%setadiGamma(gam=5.0/3.d0) !! ratio of heat capacity
  call g1%setCFL(CFL=0.4d0) !! courant number
  call g1%setSlopeLimiter(limiterType=3)  !! 3=minmod
  call g1%setSolverType(solverType=3) !! 3=HLLC, 4=AdiHLLMHD, 5=AdiHLLDMHD
  call g1%setBoundaryType(boundaryType=1) !! 1=zero gradient
  call g1%initVariable()
  call g1%exchangeBdryMPI(g1%q,g1%winq)
  call g1%setBoundary(g1%q)
  call g1%writeGrid()
  call g1%writeGrid_HD_vtk()
endif
g1%writeFlag=.false.

t0=MPI_WTIME()
do while (g1%t .lt. g1%tend)
   nstep=nstep+1
   call g1%griddt()
   call g1%evolveGridRK2()
   if(myid .eq. 0 ) then
     print *,"testSuite.f03: gridID=",g1%gridID,'t=',g1%t,"dt=",g1%dt
   endif
   if(g1%writeFlag .eqv. .true.) then
      t1=MPI_WTIME()
      if(myid .eq. 0) then
        print *,"testSuiteMPI.f03: elapse time=",t1-t0,"secs"
      endif
      call g1%writeGrid()
      call g1%writeGrid_HD_vtk()
      g1%writeFlag=.false.
      t0=MPI_WTIME()
   endif
enddo !! end do while

end subroutine HDBlastWaveCart2D

subroutine initHDBlastWaveCart2D(this,q)
class(grid)::this
double precision,dimension(1-this%nbuf:this%nMesh(1)+this%nbuf,1-this%nbuf:this%nMesh(2)+this%nbuf,this%nvar)::q
double precision,dimension(1-this%nbuf:this%nMesh(1)+this%nbuf)::xc
double precision,dimension(1-this%nbuf:this%nMesh(2)+this%nbuf)::yc
integer::nx,ny,nvar,nbuf,gridID,eosType
integer::i,j
double precision::gam,rho,vx,vy,p,r,R0


  nx=this%nMesh(1)
  ny=this%nMesh(2)
 gam=this%adiGamma
  xc=this%xc(1)%coords
  yc=this%xc(2)%coords
  R0=0.1d0

  do j=1,ny
    do i=1,nx
      rho=1.d0
       vx=0.d0
       vy=0.d0
       r=dsqrt(xc(i)**2+yc(j)**2)

       if(r .le. R0) then
         p=10.d0
       else 
         p=0.1d0
       endif

       q(i,j,1)=rho
       q(i,j,2)=rho*vx
       q(i,j,3)=rho*vy
       q(i,j,4)=p/(gam-1.d0)+0.5d0*rho*(vx**2+vy**2)
    enddo
  enddo
end subroutine initHDBlastWaveCart2D


subroutine bdryHDBlastWaveCart2D(this,q)
class(grid)::this
double precision,dimension(1-this%nbuf:this%nMesh(1)+this%nbuf,1-this%nbuf:this%nMesh(2)+this%nbuf,this%nvar)::q
integer::nx,ny,nbuf,i,j

nbuf=this%nbuf
nx=this%nMesh(1)
ny=this%nMesh(2)

if(this%boundaryType .eq. 1) then !! zero gradient
!!!!!! left boundary
  if(this%left_mpi .lt. 0) then
    do i=1,nbuf
       q(1-i,:,:)=q(1,:,:)
    enddo
  endif
!!!!!! right boundary
  if(this%right_mpi .lt. 0) then
    do i=1,nbuf
       q(nx+i,:,:)=q(nx,:,:)
    enddo
  endif
!!!!!! up boundary
  if(this%up_mpi .lt. 0) then
    do j=1,nbuf
       q(:,1-j,:)=q(:,1,:)
    enddo
  endif
!!!!!! down boundary
  if(this%down_mpi .lt. 0) then
    do j=1,nbuf
       q(:,ny+j,:)=q(:,ny,:)
    enddo
  endif
endif
end subroutine bdryHDBlastWaveCart2D
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!! Blast Wave in 2D polar coordinates
!!!!! Skinner & Ostriker, 2010, ApJS, 188, 290
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine HDBlastWavePolar2D(gridID)
integer::gridID
type(grid)::g1
integer::ndim,nbuf,coordType,variable(8)
integer::nMesh(2),dims(2)
double precision::leftBdry(2), rightBdry(2)
logical::periods(2),reorder
double precision::t0,t1
integer::nstep,ierr,Nr,Nphi

nstep=0
variable=0
ndim=2
nbuf=2
coordType=2   !! 3=uniform mesh in r, 2= log in r
variable(1)=1  !! density
variable(2)=1  !! momentum in r
variable(3)=1  !! momentum in phi
variable(8)=1  !! energy

Nr=800
Nphi=1200
nMesh(1)=Nr
nMesh(2)=Nphi
leftBdry(1)=1.d0
leftBdry(2)=-0.5d0
rightBdry(1)=2.d0
rightBdry(2)=0.5d0
dims=(/0,0/)
call MPI_DIMS_CREATE(nprocs,ndim,dims,ierr)
periods(1)=.false.
periods(2)=.false.
reorder=.true.

call g1%setTopologyMPI(ndim,dims,periods,reorder)
call g1%setGridID(gridID=gridID)
call g1%setTime(fstart=0,tend=0.2d0,dtout=0.1d0)
if (g1%fstart .eq. 0) then
  call g1%setMesh(nMesh,leftBdry,rightBdry,nbuf,coordType,gridID)
  call g1%setVariable(variable) !! den,vx,vy,vz,bx,by,bz,ene
  call g1%setMPIWindows()
  call g1%setEoS(eosType=2) !! adiabatic
  call g1%setadiGamma(gam=5.0/3.d0) !! ratio of heat capacity
  call g1%setCFL(CFL=0.4d0) !! courant number
  call g1%setSlopeLimiter(limiterType=3)  !! 3=minmod
  call g1%setSolverType(solverType=3) !! 3=HLLC, 4=AdiHLLMHD, 5=AdiHLLDMHD
  call g1%setBoundaryType(boundaryType=1) !! 1=zero gradient
  call g1%initVariable()
  call g1%exchangeBdryMPI(g1%q,g1%winq)
  call g1%setBoundary(g1%q)
  call g1%writeGrid()
endif
g1%writeFlag=.false.

t0=MPI_WTIME()
do while (g1%t .lt. g1%tend)
   nstep=nstep+1
   call g1%griddt()
   call g1%evolveGridRK2()
   if(myid .eq. 0 ) then
     print *,"testSuite.f03: gridID=",g1%gridID,'t=',g1%t,"dt=",g1%dt
   endif
   if(g1%writeFlag .eqv. .true.) then
      t1=MPI_WTIME()
      if(myid .eq. 0) then
        print *,"testSuiteMPI.f03: elapse time=",t1-t0,"secs"
      endif
      call g1%writeGrid()
      g1%writeFlag=.false.
      t0=MPI_WTIME()
   endif
enddo !! end do while

end subroutine HDBlastWavePolar2D

subroutine initHDBlastWavePolar2D(this,q)
class(grid)::this
double precision,dimension(1-this%nbuf:this%nMesh(1)+this%nbuf,1-this%nbuf:this%nMesh(2)+this%nbuf,this%nvar)::q
double precision,dimension(1-this%nbuf:this%nMesh(1)+this%nbuf)::rc
double precision,dimension(1-this%nbuf:this%nMesh(2)+this%nbuf)::phic
integer::nr,nphi,nvar,nbuf,gridID,eosType
integer::i,j
double precision::gam,rho,vr,vphi,p
double precision::R0,r,xc,yc,rmid

  rc=this%xc(1)%coords
phic=this%xc(2)%coords


  nr=this%nMesh(1)
nphi=this%nMesh(2)
 gam=this%adiGamma
  R0=0.1d0
  rmid=1.5d0

  do j=1,nphi
    do i=1,nr
      rho=1.d0
       vr=0.d0
       vphi=0.d0
       xc=rc(i)*dcos(phic(j))
       yc=rc(i)*dsin(phic(j))
       r=dsqrt((xc-rmid)**2+yc**2)

       if(r .le. R0) then
         p=10.d0
       else 
         p=0.1d0
       endif

       q(i,j,1)=rho
       q(i,j,2)=rho*vr
       q(i,j,3)=rho*vphi
       q(i,j,4)=p/(gam-1.d0)+0.5d0*rho*(vr**2+vphi**2)
    enddo
  enddo
end subroutine initHDBlastWavePolar2D


subroutine bdryHDBlastWavePolar2D(this,q)
class(grid)::this
double precision,dimension(1-this%nbuf:this%nMesh(1)+this%nbuf,1-this%nbuf:this%nMesh(2)+this%nbuf,this%nvar)::q
integer::nx,ny,nbuf,i,j

nbuf=this%nbuf
nx=this%nMesh(1)
ny=this%nMesh(2)

if(this%boundaryType .eq. 1) then !! zero gradient
!!!!!! left boundary
  if(this%left_mpi .lt. 0) then
    do i=1,nbuf
       q(1-i,:,:)=q(1,:,:)
    enddo
  endif
!!!!!! right boundary
  if(this%right_mpi .lt. 0) then
    do i=1,nbuf
       q(nx+i,:,:)=q(nx,:,:)
    enddo
  endif
!!!!!! up boundary
  if(this%up_mpi .lt. 0) then
    do j=1,nbuf
       q(:,1-j,:)=q(:,1,:)
    enddo
  endif
!!!!!! down boundary
  if(this%down_mpi .lt. 0) then
    do j=1,nbuf
       q(:,ny+j,:)=q(:,ny,:)
    enddo
  endif
endif
end subroutine bdryHDBlastWavePolar2D
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!! Selfgravity Test 3D
!!!!! PASJ: Publ. Astron. Soc. Japan 62, 301-314, 2010 April 25
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine selfgravityTest3D(gridID)
integer::gridID
type(grid)::g1
integer::ndim,nbuf,coordType,variable(8)
integer::nMesh(3),dims(3)
double precision::leftBdry(3),rightBdry(3)
logical::periods(3),reorder
double precision::t0,t1
integer::nstep,ierr

nstep=0
variable=0
ndim=3
nbuf=2
coordType=1

variable(1)=1  !! density
variable(2)=1  !! momx
variable(3)=1  !! momy
variable(4)=1  !! momz
variable(8)=1  !! ene

nMesh(1)=128
nMesh(2)=128
nMesh(3)=128
leftBdry(1)=-2.d0
leftBdry(2)=-2.d0
leftBdry(3)=-2.d0

rightBdry(1)=2.d0
rightBdry(2)=2.d0
rightBdry(3)=2.d0

dims(1)=1
dims(2)=1
dims(3)=nprocs ! FFTW sg

periods(1)=.false.
periods(2)=.false.
periods(3)=.false.
reorder=.true.


call g1%enableSelfgravity()
call g1%setGravConst(GravConst=1.d0)
call g1%setTopologyMPI(ndim,dims,periods,reorder)
call g1%setMesh(nMesh,leftBdry,rightBdry,nbuf,coordType,gridID)
call g1%setVariable(variable)
call g1%setSelfgravityKernel()
call g1%setSGMPIWindows()
call g1%setMPIWindows()
call g1%setEOS(eosType=2)
call g1%setadiGamma(gam=5.d0/3.d0) !! 1.1d0 for a collapsing rotating sphere
call g1%setCFL(CFL=0.4d0)
call g1%setTime(fstart=0,tend=3.d0,dtout=0.3d0)
call g1%setSlopeLimiter(limiterType=3)
call g1%setSolverType(solverType=3)
call g1%setBoundaryType(boundaryType=1)
call g1%initVariable()
call g1%exchangeBdryMPI(g1%q,g1%winq)
call g1%setBoundary(g1%q)
call g1%calcSelfgravity(g1%q)
call g1%writeGrid()
! call g1%writeGrid_HD_vtk()
call g1%griddt()
g1%writeFlag=.false.

t0=MPI_WTIME()
do while (g1%t .lt. g1%tend)
   nstep=nstep+1
   call g1%griddt()
   call g1%evolveGridRK2()
   if(myid .eq. 0 ) then
     print *,"testSuite.f03: gridID=",g1%gridID,'t=',g1%t,"dt=",g1%dt
   endif
   if(g1%writeFlag .eqv. .true.) then
      t1=MPI_WTIME()
      if(myid .eq. 0) then
        print *,"testSuiteMPI.f03: elapse time=",t1-t0,"secs"
      endif
      call g1%writeGrid()
      call g1%writeGrid_HD_vtk()
      g1%writeFlag=.false.
      t0=MPI_WTIME()
      if(g1%dt .lt. 1d-8) then
        stop
      endif
   endif
   if(g1%dt .lt. 1d-8) then
     g1%writeFlag=.true.
     g1%fnum=g1%fnum+1
   endif

enddo !! end do while

end subroutine selfgravityTest3D

subroutine initSelfgravityTest3D(this,q)
class(grid)::this
double precision,dimension(1-this%nbuf:this%nMesh(1)+this%nbuf,1-this%nbuf:this%nMesh(2)+this%nbuf,&
                           1-this%nbuf:this%nMesh(3)+this%nbuf,this%nvar)::q
double precision,dimension(1-this%nbuf:this%nMesh(1)+this%nbuf)::xc
double precision,dimension(1-this%nbuf:this%nMesh(2)+this%nbuf)::yc
double precision,dimension(1-this%nbuf:this%nMesh(3)+this%nbuf)::zc
double precision::x,y,z,r,r2
double precision::pi,rho0,gam,ep,r0,Vrot,vx,vy,vz,Omega
integer::nx,ny,nz,nvar,nbuf,gridID,eosType
integer::i,j,k

pi=datan(1.d0)*4.d0
xc=this%xc(1)%coords
yc=this%xc(2)%coords
zc=this%xc(3)%coords
nx=this%nMesh(1)
ny=this%nMesh(2)
nz=this%nMesh(3)
gam=this%adiGamma
r0=1.d0
rho0=1.d0/(2.d0*pi*r0)
Omega=0.3d0

do k=1,nz
  do j=1,ny
    do i=1,nx
         x=xc(i)
         y=yc(j)
         z=zc(k)
         r=dsqrt(x**2+y**2+z**2)
         r2=dsqrt(x**2+y**2)
         
         if(r .le. r0) then
           vx=0.d0 !-Omega*y  !! for a collapsing rotating sphere
           vy=0.d0 ! Omega*x
           vz=0.d0
           q(i,j,k,1)=rho0/r !rho0    !! for a collapsing rotating sphere
           q(i,j,k,2)=0.d0 !q(i,j,k,1)*vx
           q(i,j,k,3)=0.d0 !q(i,j,k,1)*vy
           q(i,j,k,4)=0.d0 !q(i,j,k,1)*vz
         else
           vx=0.d0
           vy=0.d0
           vz=0.d0
           q(i,j,k,1)=rho0/1000.d0!/r0
           q(i,j,k,2)=0.d0
           q(i,j,k,3)=0.d0
           q(i,j,k,4)=0.d0
         endif
         q(i,j,k,5)=0.05d0*q(i,j,k,1)+0.5d0*q(i,j,k,1)*(vx**2+vy**2+vz**2)
    enddo
  enddo
enddo
end subroutine initSelfgravityTest3D

subroutine bdrySelfgravityTest3D(this,q)
class(grid)::this
double precision,dimension(1-this%nbuf:this%nMesh(1)+this%nbuf,1-this%nbuf:this%nMesh(2)+this%nbuf,&
& 1-this%nbuf:this%nMesh(3)+this%nbuf,this%nvar)::q
integer::nx,ny,nz,nbuf,i,j,k

nbuf=this%nbuf
nx=this%nMesh(1)
ny=this%nMesh(2)
nz=this%nMesh(3)

if(this%boundaryType .eq. 1) then  !! zero gradient boundary condition
!!!! left boundary
   if(this%left_mpi .lt. 0) then
     do i=1,nbuf
       q(1-i,:,:,:)=q(1,:,:,:)
     enddo
   endif
!!!! right boundary
   if(this%right_mpi .lt. 0) then
     do i=1,nbuf
       q(nx+i,:,:,:)=q(nx,:,:,:)
     enddo
   endif
!!!! up boundary
   if(this%up_mpi .lt. 0) then
     do j=1,nbuf
       q(:,1-j,:,:)=q(:,1,:,:)
     enddo
   endif
!!!! down boundary
   if(this%down_mpi .lt. 0) then
     do j=1,nbuf
       q(:,ny+j,:,:)=q(:,ny,:,:)
     enddo
   endif
!!!! top boundary
   if(this%top_mpi .lt. 0) then
     do k=1,nbuf
       q(:,:,1-k,:)=q(:,:,1,:)
     enddo
   endif
!!!! bottom boundary
   if(this%bottom_mpi .lt. 0) then
     do k=1,nbuf
       q(:,:,nz+k,:)=q(:,:,nz,:)
     enddo
   endif
!!!!
endif
end subroutine bdrySelfgravityTest3D
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!! MHD Blast Wave (Gardiner & Stone, 2008, JCP, 227, 4123)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine MHDBlastWave(gridID)
integer::gridID
type(grid)::g1
integer::ndim,nbuf,coordType,variable(8)
integer::nMesh(3),dims(3)
double precision::leftBdry(3), rightBdry(3)
logical::periods(3),reorder
double precision::t0,t1
integer::nstep,ierr,N

nstep=0
variable=0
ndim=3
nbuf=2
coordType=1
variable(1)=1
variable(2)=1
variable(3)=1
variable(4)=1
variable(5)=1
variable(6)=1
variable(7)=1
variable(8)=1

N=200
nMesh(1)=N
nMesh(2)=N
nMesh(3)=N
leftBdry(1)=-0.5d0
leftBdry(2)=-0.5d0
leftBdry(3)=-0.5d0
rightBdry(1)=0.5d0
rightBdry(2)=0.5d0
rightBdry(3)=0.5d0
dims=(/0,0,0/)
call MPI_DIMS_CREATE(nprocs,ndim,dims,ierr)
periods(1)=.true.
periods(2)=.true.
periods(3)=.true.
reorder=.true.

call g1%setTopologyMPI(ndim,dims,periods,reorder)
call g1%setGridID(gridID=gridID)
call g1%setTime(fstart=0,tend=0.02d0,dtout=0.01d0)
if (g1%fstart .eq. 0) then
!call g1%setTopologyMPI(ndim,dims,periods,reorder)
  call g1%setMesh(nMesh,leftBdry,rightBdry,nbuf,coordType,gridID)
  call g1%setVariable(variable) !! den,vx,vy,vz,bx,by,bz,ene
  call g1%setMPIWindows()
  call g1%setEoS(eosType=2) !! adiabatic
  call g1%setadiGamma(gam=5.d0/3.d0) !! ratio of heat capacity
  call g1%setCFL(CFL=0.4d0) !! courant number
!call g1%setTime(fstart=1,tend=0.02d0,dtout=0.01d0)
  call g1%setSlopeLimiter(limiterType=3)  !! minmod
  call g1%setSolverType(solverType=5) !! 4=AdiHLLMHD, 5=AdiHLLDMHD
  call g1%setBoundaryType(boundaryType=3) !! periodic
  call g1%initVariable()
  call g1%exchangeBdryMPI(g1%q,g1%winq)
  call g1%setBoundary(g1%q)
  call g1%writeGrid()
endif
g1%writeFlag=.false.

t0=MPI_WTIME()
do while (g1%t .lt. g1%tend)
   nstep=nstep+1
   call g1%griddt()
   call g1%evolveGridRK2()
   if(myid .eq. 0 ) then
     print *,"testSuite.f03: gridID=",g1%gridID,'t=',g1%t,"dt=",g1%dt
   endif
   if(g1%writeFlag .eqv. .true.) then
      t1=MPI_WTIME()
      if(myid .eq. 0) then
        print *,"testSuiteMPI.f03: elapse time=",t1-t0,"secs"
      endif
      call g1%writeGrid()
      g1%writeFlag=.false.
      t0=MPI_WTIME()
   endif

enddo !! end do while

end subroutine MHDBlastWave


subroutine initMHDBlastWave(this,q)
class(grid)::this
double precision,dimension(1-this%nbuf:this%nMesh(1)+this%nbuf,1-this%nbuf:this%nMesh(2)+this%nbuf,&
                           1-this%nbuf:this%nMesh(3)+this%nbuf,this%nvar)::q
double precision,dimension(1-this%nbuf:this%nMesh(1)+this%nbuf)::xc
double precision,dimension(1-this%nbuf:this%nMesh(2)+this%nbuf)::yc
double precision,dimension(1-this%nbuf:this%nMesh(3)+this%nbuf)::zc
integer::nx,ny,nz,nvar,nbuf,gridID,eosType
integer::i,j,k
double precision::gam,rho,vx,vy,vz,p
double precision::R0,bx,by,bz,r

xc=this%xc(1)%coords
yc=this%xc(2)%coords
zc=this%xc(3)%coords


nx=this%nMesh(1)
ny=this%nMesh(2)
nz=this%nMesh(3)
gam=this%adiGamma
R0=0.125d0


do k=1,nz
  do j=1,ny
    do i=1,nx
      rho=1.d0
       vx=0.d0
       vy=0.d0
       vz=0.d0
       bx=10.d0/dsqrt(2.d0)
       by=0.d0
       bz=10.d0/dsqrt(2.d0)
        r=dsqrt(xc(i)**2+yc(j)**2+zc(k)**2)
       if(r .le. R0) then
         p=100.d0
       else 
         p=1.d0
       endif

       q(i,j,k,1)=rho
       q(i,j,k,2)=rho*vx
       q(i,j,k,3)=rho*vy
       q(i,j,k,4)=rho*vz
       q(i,j,k,5)=bx
       q(i,j,k,6)=by
       q(i,j,k,7)=bz
       q(i,j,k,8)=p/(gam-1.d0)+0.5d0*rho*(vx**2+vy**2+vz**2)+0.5d0*(bx**2+by**2+bz**2)
       q(i,j,k,9)=bx
       q(i,j,k,10)=by
       q(i,j,k,11)=bz
    enddo
  enddo
enddo
end subroutine initMHDBlastWave


subroutine bdryMHDBlastWave(this,q)
class(grid)::this
double precision,dimension(1-this%nbuf:this%nMesh(1)+this%nbuf,1-this%nbuf:this%nMesh(2)+this%nbuf,&
                           1-this%nbuf:this%nMesh(3)+this%nbuf,this%nvar)::q
integer::i,j,k,nx,ny,nz,nbuf

nbuf=this%nbuf
nx=this%nMesh(1)
ny=this%nMesh(2)
nz=this%nMesh(3)

if(this%boundaryType .eq. 3) then  !! periodic boundary condition
  !! MPI takes care the data exchange...
endif

end subroutine bdryMHDBlastWave
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!! Field Loop Advection (Gardiner & Stone, 2008, JCP, 227, 4123)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine FieldLoopAdvection(gridID)
integer::gridID
type(grid)::g1
integer::ndim,nbuf,coordType,variable(8)
integer::nMesh(3),dims(3)
double precision::leftBdry(3), rightBdry(3)
logical::periods(3),reorder
double precision::t0,t1
integer::nstep,ierr,N

nstep=0
variable=0
ndim=3
nbuf=2
coordType=1
variable(1)=1
variable(2)=1
variable(3)=1
variable(4)=1
variable(5)=1
variable(6)=1
variable(7)=1
variable(8)=1

N=128
nMesh(1)=N
nMesh(2)=N
nMesh(3)=2*N
leftBdry(1)=-0.5d0
leftBdry(2)=-0.5d0
leftBdry(3)= -1.d0
rightBdry(1)=0.5d0
rightBdry(2)=0.5d0
rightBdry(3)= 1.d0
dims=(/0,0,0/)
call MPI_DIMS_CREATE(nprocs,ndim,dims,ierr)
periods(1)=.true.
periods(2)=.true.
periods(3)=.true.
reorder=.true.

call g1%setTopologyMPI(ndim,dims,periods,reorder)
call g1%setMesh(nMesh,leftBdry,rightBdry,nbuf,coordType,gridID)
call g1%setVariable(variable) !! den,vx,vy,vz,bx,by,bz,ene
call g1%setMPIWindows()
call g1%setEoS(eosType=2) !! adiabatic
call g1%setadiGamma(gam=5.d0/3.d0) !! ratio of heat capacity
call g1%setCFL(CFL=0.4d0) !! courant number
call g1%setTime(fstart=0,tend=1.d0,dtout=0.01d0)
call g1%setSlopeLimiter(limiterType=1)  !! minmod
call g1%setSolverType(solverType=5) !! 4=AdiHLLMHD, 5=AdiHLLDMHD
call g1%setBoundaryType(boundaryType=3) !! periodic
call g1%initVariable()
call g1%exchangeBdryMPI(g1%q,g1%winq)
call g1%setBoundary(g1%q)
call g1%writeGrid()
g1%writeFlag=.false.

t0=MPI_WTIME()
do while (g1%t .lt. g1%tend)
   nstep=nstep+1
   call g1%griddt()
   call g1%evolveGridRK2()
   if(myid .eq. 0 ) then
     print *,"testSuite.f03: gridID=",g1%gridID,'t=',g1%t,"dt=",g1%dt
   endif
   if(g1%writeFlag .eqv. .true.) then
      t1=MPI_WTIME()
      if(myid .eq. 0) then
        print *,"testSuiteMPI.f03: elapse time=",t1-t0,"secs"
      endif
      call g1%writeGrid()
      g1%writeFlag=.false.
      t0=MPI_WTIME()
   endif

enddo !! end do while

end subroutine FieldLoopAdvection

subroutine initFieldLoopAdvection(this,q)
class(grid)::this
double precision,dimension(1-this%nbuf:this%nMesh(1)+this%nbuf,1-this%nbuf:this%nMesh(2)+this%nbuf,&
                           1-this%nbuf:this%nMesh(3)+this%nbuf,this%nvar)::q
double precision,dimension(1-this%nbuf:this%nMesh(1)+this%nbuf,1-this%nbuf:this%nMesh(2)+this%nbuf,&
                           1-this%nbuf:this%nMesh(3)+this%nbuf)::Ax,Ay,Az
double precision,dimension(1-this%nbuf:this%nMesh(1)+this%nbuf)::xc,xl,xr
double precision,dimension(1-this%nbuf:this%nMesh(2)+this%nbuf)::yc,yl,yr
double precision,dimension(1-this%nbuf:this%nMesh(3)+this%nbuf)::zc,zl,zr
integer::nx,ny,nz,nvar,nbuf,gridID,eosType
integer::i,j,k
double precision::gam,rho,vx,vy,vz,bxl,byl,bzl,bxr,byr,bzr,p
double precision::bxc,byc,bzc
double precision::x1,x2,x3,r,B0,R0,x,y,z,dx,dy,dz

xc=this%xc(1)%coords
yc=this%xc(2)%coords
zc=this%xc(3)%coords

xl=this%xl(1)%coords
yl=this%xl(2)%coords
zl=this%xl(3)%coords

xr=this%xr(1)%coords
yr=this%xr(2)%coords
zr=this%xr(3)%coords

nx=this%nMesh(1)
ny=this%nMesh(2)
nz=this%nMesh(3)
gam=this%adiGamma
R0=0.3d0
B0=1d-3

do k=1,nz+1
  do j=1,ny+1
    do i=1,nx+1
       x=xl(i)
       y=yl(j)
       z=zl(k) 
      x1=(2.d0*x+z)/dsqrt(5.d0)
      x2=y
      x3=(2.d0*z-x)/dsqrt(5.d0)
      ! r=dsqrt(x1**2+x2**2)
      !if(r .le. R0) then
      if(dsqrt(x1**2+x2**2) .le. R0) then
        r=dsqrt(x1**2+x2**2)
        Ax(i,j,k)=B0*(R0-r)*(-1.d0/dsqrt(5.d0)) 
        Ay(i,j,k)=0.d0
        Az(i,j,k)=B0*(R0-r)*( 2.d0/dsqrt(5.d0)) 
      elseif(dsqrt((x1-2.d0/dsqrt(5.d0))**2+x2**2) .lt. R0) then
        r=dsqrt((x1-2.d0/dsqrt(5.d0))**2+x2**2)
        Ax(i,j,k)=B0*(R0-r)*(-1.d0/dsqrt(5.d0))
        Ay(i,j,k)=0.d0
        Az(i,j,k)=B0*(R0-r)*( 2.d0/dsqrt(5.d0))
      elseif(dsqrt((x1+2.d0/dsqrt(5.d0))**2+x2**2) .lt. R0) then
        r=dsqrt((x1+2.d0/dsqrt(5.d0))**2+x2**2)
        Ax(i,j,k)=B0*(R0-r)*(-1.d0/dsqrt(5.d0))
        Ay(i,j,k)=0.d0
        Az(i,j,k)=B0*(R0-r)*( 2.d0/dsqrt(5.d0))        
      else
        Ax(i,j,k)=0.d0
        Ay(i,j,k)=0.d0
        Az(i,j,k)=0.d0
      endif
    enddo
  enddo
enddo

do k=1,nz
  do j=1,ny
    do i=1,nx
      dx=this%dx(1)%coords(i)
      dy=this%dx(2)%coords(j)
      dz=this%dx(3)%coords(k)
      rho=1.d0
       vx=1.d0
       vy=1.d0
       vz=2.d0
        p=1.d0

       bxl= (    Az(i,j+1,k)-Az(i,j,k))/dy
       byl= (    Ax(i,j,k+1)-Ax(i,j,k))/dz-(    Az(i+1,j,k)-Az(i,j,k))/dx
       bzl=-(    Ax(i,j+1,k)-Ax(i,j,k))/dy

       bxr= (Az(i+1,j+1,k)-Az(i+1,j,k))/dy
       byr= (Ax(i,j+1,k+1)-Ax(i,j+1,k))/dz-(Az(i+1,j+1,k)-Az(i,j+1,k))/dx
       bzr=-(Ax(i,j+1,k+1)-Ax(i,j,k+1))/dy

       bxc=0.5d0*(bxl+bxr)
       byc=0.5d0*(byl+byr)
       bzc=0.5d0*(bzl+bzr)

       q(i,j,k,1)=rho
       q(i,j,k,2)=rho*vx
       q(i,j,k,3)=rho*vy
       q(i,j,k,4)=rho*vz
       q(i,j,k,5)=bxl
       q(i,j,k,6)=byl
       q(i,j,k,7)=bzl
       q(i,j,k,8)=p/(gam-1.d0)+0.5d0*rho*(vx**2+vy**2+vz**2)+0.5d0*(bxc**2+byc**2+bzc**2)
       q(i,j,k,9)=bxr
       q(i,j,k,10)=byr
       q(i,j,k,11)=bzr
    enddo
  enddo
enddo
end subroutine initFieldLoopAdvection

subroutine bdryFieldLoopAdvection(this,q)
class(grid)::this
double precision,dimension(1-this%nbuf:this%nMesh(1)+this%nbuf,1-this%nbuf:this%nMesh(2)+this%nbuf,&
                           1-this%nbuf:this%nMesh(3)+this%nbuf,this%nvar)::q
integer::i,j,k,nx,ny,nz,nbuf

nbuf=this%nbuf
nx=this%nMesh(1)
ny=this%nMesh(2)
nz=this%nMesh(3)

if(this%boundaryType .eq. 3) then  !! periodic boundary condition
  !! MPI takes care the data exchange...
endif

end subroutine bdryFieldLoopAdvection

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!! Brio-Wu Shocktube 3D (1988, JCP, 75, 400)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine BrioWuShockTube3D(gridID)
integer::gridID
type(grid)::g1
integer::ndim,nbuf,coordType,variable(8)
integer::nMesh(3),dims(3)
double precision::leftBdry(3),rightBdry(3)
logical::periods(3),reorder
double precision::t0,t1
integer::nstep,ierr

nstep=0
variable=0
ndim=3
nbuf=2
coordType=1 !! Cartesian
variable(1)=1  !! rho
variable(2)=1  !! vx
variable(3)=1  !! vy
variable(4)=1  !! vz
variable(5)=1  !! bx
variable(6)=1  !! by
variable(7)=1  !! bz
variable(8)=1  !! energy

nMesh(1)=4
nMesh(2)=400
nMesh(3)=4
leftBdry(1)=0.d0
rightBdry(1)=0.01d0
leftBdry(2)=0.d0
rightBdry(2)=1.d0
leftBdry(3)=0.d0
rightBdry(3)=0.01d0
dims=(/0,0,0/)
call MPI_DIMS_CREATE(nprocs,ndim,dims,ierr)
!dims(1)=1
!dims(2)=1
!dims(3)=nprocs
periods(1)=.false.
periods(2)=.false.
periods(3)=.false.
reorder=.true.

call g1%setTopologyMPI(ndim,dims,periods,reorder)
call g1%setMesh(nMesh,leftBdry,rightBdry,nbuf,coordType,gridID)
call g1%setVariable(variable) !! den,vx,vy,vz,bx,by,bz,ene
call g1%setMPIWindows()
call g1%setEoS(eosType=2) !! adiabatic
call g1%setadiGamma(gam=2.d0) !! ratio of heat capacity
call g1%setCFL(CFL=0.8d0) !! courant number
call g1%setTime(fstart=0,tend=0.1d0,dtout=0.1d0)
call g1%setSlopeLimiter(limiterType=3)
call g1%setSolverType(solverType=5) !! 4=AdiHLLMHD, 5=AdiHLLDMHD
call g1%setBoundaryType(boundaryType=1)
call g1%initVariable()
call g1%setBoundary(g1%q)
call g1%exchangeBdryMPI(g1%q,g1%winq)
!call g1%setBoundary(g1%q)
call g1%writeGrid()
g1%writeFlag=.false.

t0=MPI_WTIME()
do while (g1%t .lt. g1%tend)
   nstep=nstep+1
   call g1%griddt()
   call g1%evolveGridRK2()
   if(myid .eq. 0 ) then
     print *,"testSuite.f03: gridID=",g1%gridID,'t=',g1%t,"dt=",g1%dt
   endif
   if(g1%writeFlag .eqv. .true.) then
      t1=MPI_WTIME()
      if(myid .eq. 0) then
        print *,"testSuiteMPI.f03: elapse time=",t1-t0,"secs"
      endif
      call g1%writeGrid()
      g1%writeFlag=.false.
      t0=MPI_WTIME()
   endif

enddo !! end do while


end subroutine BrioWuShockTube3D

subroutine initBrioWuShockTube3D(this,q)
class(grid)::this
double precision,dimension(1-this%nbuf:this%nMesh(1)+this%nbuf,1-this%nbuf:this%nMesh(2)+this%nbuf,&
                           1-this%nbuf:this%nMesh(3)+this%nbuf,this%nvar)::q
double precision,dimension(1-this%nbuf:this%nMesh(1)+this%nbuf)::xc
double precision,dimension(1-this%nbuf:this%nMesh(2)+this%nbuf)::yc
double precision,dimension(1-this%nbuf:this%nMesh(3)+this%nbuf)::zc

integer::nx,ny,nz,nvar,nbuf,gridID,eosType
integer::i,j,k
double precision::gam,rho,vx,vy,vz,bxl,byl,bzl,bxr,byr,bzr,p
double precision::bxc,byc,bzc

xc=this%xc(1)%coords
yc=this%xc(2)%coords
zc=this%xc(3)%coords
nx=this%nMesh(1)
ny=this%nMesh(2)
nz=this%nMesh(3)
gam=this%adiGamma

do k=1,nz
  do j=1,ny
    do i=1,nx
      if(yc(j) .lt. 0.5d0) then
        rho=1.d0
         vx=0.d0
         vy=0.d0
         vz=0.d0
         bxl=1.d0
         bxr=1.d0
         byl=0.75d0
         byr=0.75d0
         bzl=0.d0
         bzr=0.d0
           p=1.d0
      else
         rho=0.125d0
          vx=0.d0
          vy=0.d0
          vz=0.d0
          bxl=-1.d0
          bxr=-1.d0
          byl=0.75d0
          byr=0.75d0
          bzl=0.d0
          bzr=0.d0
            p=0.1d0
      endif
      bxc=0.5d0*(bxl+bxr)
      byc=0.5d0*(byl+byr)
      bzc=0.5d0*(bzl+bzr)

      q(i,j,k,1)=rho
      q(i,j,k,2)=rho*vx
      q(i,j,k,3)=rho*vy
      q(i,j,k,4)=rho*vz
      q(i,j,k,5)=bxl
      q(i,j,k,6)=byl
      q(i,j,k,7)=bzl
      q(i,j,k,8)=p/(gam-1.d0)+0.5d0*rho*(vx**2.d0+vy**2.d0+vz**2.d0)+0.5d0*(bxc**2.d0+byc**2.d0+bzc**2.d0)
      q(i,j,k,9)=bxr
      q(i,j,k,10)=byr
      q(i,j,k,11)=bzr
    enddo
  enddo
enddo

end subroutine initBrioWuShockTube3D

subroutine bdryBrioWuShockTube3D(this,q)
class(grid)::this
double precision,dimension(1-this%nbuf:this%nMesh(1)+this%nbuf,1-this%nbuf:this%nMesh(2)+this%nbuf,&
                           1-this%nbuf:this%nMesh(3)+this%nbuf,this%nvar)::q
integer::nx,ny,nz,nbuf,i,j,k

nbuf=this%nbuf
nx=this%nMesh(1)
ny=this%nMesh(2)
nz=this%nMesh(3)

if(this%boundaryType .eq. 1) then !! zero gradient
!!!!!! left boundary
  if(this%left_mpi .lt. 0) then
    do i=1,nbuf
       q(1-i,:,:,:)=q(1,:,:,:)
    enddo
  endif
!!!!!! right boundary
  if(this%right_mpi .lt. 0) then
    do i=1,nbuf
       q(nx+i,:,:,:)=q(nx,:,:,:)
    enddo
  endif
!!!!!! up boundary
  if(this%up_mpi .lt. 0) then
    do j=1,nbuf
       q(:,1-j,:,:)=q(:,1,:,:)
    enddo
  endif
!!!!!! down boundary
  if(this%down_mpi .lt. 0) then
    do j=1,nbuf
       q(:,ny+j,:,:)=q(:,ny,:,:)
    enddo
  endif
!!!! top boundary
   if(this%top_mpi .lt. 0) then
     do k=1,nbuf
       q(:,:,1-k,:)=q(:,:,1,:)
     enddo
   endif
!!!! bottom boundary
   if(this%bottom_mpi .lt. 0) then
     do k=1,nbuf
       q(:,:,nz+k,:)=q(:,:,nz,:)
     enddo
   endif
!!!!

endif

end subroutine bdryBrioWuShockTube3D
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!! isothermal Shocktube problem 3D
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine isoShockTube3D(gridID)
integer::gridID
type(grid)::g1
integer::ndim,nbuf,coordType,variable(8)
integer::nMesh(3),dims(3)
double precision::leftBdry(3),rightBdry(3)
logical::periods(3),reorder
double precision::t0,t1
integer::nstep,ierr

nstep=0
variable=0
ndim=3
nbuf=2
coordType=1
variable(1)=1
variable(2)=1
variable(3)=1
variable(4)=1

nMesh(1)=20
nMesh(2)=400
nMesh(3)=20
leftBdry(1)=-0.025d0
leftBdry(2)=-0.5d0
leftBdry(3)=-0.025d0
rightBdry(1)=0.025d0
rightBdry(2)=0.5d0
rightBdry(3)=0.025d0
dims=(/0,0,0/)
call MPI_DIMS_CREATE(nprocs,ndim,dims,ierr)
periods(1)=.false.
periods(2)=.false.
periods(3)=.false.
reorder=.true.

call g1%setTopologyMPI(ndim,dims,periods,reorder)
call g1%setMesh(nMesh,leftBdry,rightBdry,nbuf,coordType,gridID)
call g1%setVariable(variable) !! den,vx,vy,vz
call g1%setMPIWindows()
call g1%setEoS(eosType=1) !! isothermal
call g1%setSoundSpeed(snd=1.d0) !! isothermal sound speed
call g1%setCFL(CFL=0.6d0) !! Courant number
call g1%setTime(fstart=0,tend=0.2d0,dtout=0.2d0)
call g1%setSlopeLimiter(limiterType=1)
call g1%setSolverType(solverType=2)
call g1%setBoundaryType(boundaryType=1)
call g1%initVariable()
call g1%exchangeBdryMPI(g1%q,g1%winq)
call g1%setBoundary(g1%q)
call g1%writeGrid()
g1%writeFlag=.false.

t0=MPI_WTIME()
do while (g1%t .lt. g1%tend)
   call g1%griddt()
   call g1%evolveGridRK2()
   if(myid .eq. 0) then
     print *,"testSuite.f03: gridID=",g1%gridID,'t=',g1%t,"dt=",g1%dt
   endif
   if(g1%writeFlag .eqv. .true.) then
      t1=MPI_WTIME()
      if(myid .eq. 0) then
        print *,"testSuiteMPI.f03: elapse time=",t1-t0,"secs"
      endif
      call g1%writeGrid()
      g1%writeFlag=.false.
      t0=MPI_WTIME()
   endif
enddo !! end do while

end subroutine isoShockTube3D

subroutine initIsoShockTube3D(this,q)
class(grid)::this
double precision,dimension(1-this%nbuf:this%nMesh(1)+this%nbuf,1-this%nbuf:this%nMesh(2)+this%nbuf,&
& 1-this%nbuf:this%nMesh(3)+this%nbuf,this%nvar)::q
double precision,dimension(1-this%nbuf:this%nMesh(1)+this%nbuf)::xc
double precision,dimension(1-this%nbuf:this%nMesh(2)+this%nbuf)::yc
double precision,dimension(1-this%nbuf:this%nMesh(3)+this%nbuf)::zc
integer::nx,ny,nz,nvar,nbuf,gridID,eosType
integer::i,j,k

xc=this%xc(1)%coords
yc=this%xc(2)%coords
zc=this%xc(3)%coords
nx=this%nMesh(1)
ny=this%nMesh(2)
nz=this%nMesh(3)

 do k=1,nz
   do j=1,ny
     do i=1,nx
       if(yc(j) .lt. 0.d0) then
          q(i,j,k,1)=1.d0
       else
          q(i,j,k,1)=0.25d0
       endif
          q(i,j,k,2)=0.d0
          q(i,j,k,3)=0.d0
          q(i,j,k,4)=0.d0
     enddo
   enddo
 enddo
end subroutine initIsoShockTube3D

subroutine bdryIsoShockTube3D(this,q)
class(grid)::this
double precision,dimension(1-this%nbuf:this%nMesh(1)+this%nbuf,1-this%nbuf:this%nMesh(2)+this%nbuf,&
& 1-this%nbuf:this%nMesh(3)+this%nbuf,this%nvar)::q
integer::nx,ny,nz,nbuf,i,j,k

nbuf=this%nbuf
nx=this%nMesh(1)
ny=this%nMesh(2)
nz=this%nMesh(3)

if(this%boundaryType .eq. 1) then
!!!! left boundary
   if(this%left_mpi .lt. 0) then
     do i=1,nbuf
       q(1-i,:,:,:)=q(1,:,:,:)
     enddo
   endif
!!!! right boundary
   if(this%right_mpi .lt. 0) then
     do i=1,nbuf
       q(nx+i,:,:,:)=q(nx,:,:,:)
     enddo
   endif
!!!! up boundary
   if(this%up_mpi .lt. 0) then
     do j=1,nbuf
       q(:,1-j,:,:)=q(:,1,:,:)
     enddo
   endif
!!!! down boundary
   if(this%down_mpi .lt. 0) then
     do j=1,nbuf
       q(:,ny+j,:,:)=q(:,ny,:,:)
     enddo
   endif
!!!! top boundary
   if(this%top_mpi .lt. 0) then
     do k=1,nbuf
       q(:,:,1-k,:)=q(:,:,1,:)
     enddo
   endif
!!!! bottom boundary
   if(this%bottom_mpi .lt. 0) then
     do k=1,nbuf
       q(:,:,nz+k,:)=q(:,:,nz,:)
     enddo
   endif
!!!!
endif
end subroutine bdryIsoShockTube3D

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!! Sod shocktube problem 3D (Sod, 1978, JCP, 27, 1)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine SodShockTube3D(gridID)
integer::gridID
type(grid)::g1
integer::ndim,nbuf,coordType,variable(8)
integer::nMesh(3),dims(3)
double precision::leftBdry(3),rightBdry(3)
logical::periods(3),reorder
double precision::t0,t1
integer::nstep,ierr

nstep=0
variable=0
ndim=3
nbuf=2
coordType=1
variable(1)=1 !! density
variable(2)=1 !! momx
variable(3)=1 !! momy
variable(4)=1 !! momz
variable(8)=1 !! ene

nMesh(1)=20
nMesh(2)=20
nMesh(3)=400
leftBdry(1)=-0.025d0
leftBdry(2)=-0.025d0
leftBdry(3)=-0.5d0
rightBdry(1)=0.025d0
rightBdry(2)=0.025d0
rightBdry(3)=0.5d0
dims=(/0,0,0/)
call MPI_DIMS_CREATE(nprocs,ndim,dims,ierr)
dims(1)=1
dims(2)=1
dims(3)=nprocs
periods(1)=.false.
periods(2)=.false.
periods(3)=.false.
reorder=.true.

call g1%setTopologyMPI(ndim,dims,periods,reorder)
call g1%setMesh(nMesh,leftBdry,rightBdry,nbuf,coordType,gridID)
call g1%setVariable(variable)
call g1%setMPIWindows()
call g1%setEoS(eosType=2) !! adiabatic
call g1%setadiGamma(gam=1.4d0)
call g1%setCFL(CFL=0.6d0)
call g1%setTime(fstart=0,tend=0.2d0,dtout=0.2d0)
call g1%setSlopeLimiter(limiterType=1)
call g1%setSolverType(solverType=3)
call g1%setBoundaryType(boundaryType=1)
call g1%initVariable()
call g1%setBoundary(g1%q)
call g1%exchangeBdryMPI(g1%q,g1%winq)
call g1%setBoundary(g1%q)
call g1%writeGrid()
call g1%writeGrid_HD_vtk()
g1%writeFlag=.false.

t0=MPI_WTIME()
do while (g1%t .lt. g1%tend)
   call g1%griddt()
   call g1%evolveGridRK2()
   if(myid .eq. 0) then
     print *,"testSuite.f03: gridID=",g1%gridID,'t=',g1%t,"dt=",g1%dt
   endif 
   if(g1%writeFlag .eqv. .true.) then
      t1=MPI_WTIME()
      if(myid .eq. 0) then
        print *,"testSuiteMPI.f03: elapse time=",t1-t0,"secs"
      endif
      call g1%writeGrid()
      call g1%writeGrid_HD_vtk()
      g1%writeFlag=.false.
      t0=MPI_WTIME()
   endif

enddo !! end do while


end subroutine SodShockTube3D

subroutine initSodShockTube3D(this,q)
class(grid)::this
double precision,dimension(1-this%nbuf:this%nMesh(1)+this%nbuf,1-this%nbuf:this%nMesh(2)+this%nbuf,&
& 1-this%nbuf:this%nMesh(3)+this%nbuf,this%nvar)::q
double precision,dimension(1-this%nbuf:this%nMesh(1)+this%nbuf)::xc
double precision,dimension(1-this%nbuf:this%nMesh(2)+this%nbuf)::yc
double precision,dimension(1-this%nbuf:this%nMesh(3)+this%nbuf)::zc
integer::nx,ny,nz,nvar,nbuf,gridID,eosType
integer::i,j,k

xc=this%xc(1)%coords
yc=this%xc(2)%coords
zc=this%xc(3)%coords
nx=this%nMesh(1)
ny=this%nMesh(2)
nz=this%nMesh(3)

 do k=1,nz
   do j=1,ny
     do i=1,nx
       if(zc(k) .lt. 0.d0) then
          q(i,j,k,1)=1.d0
          q(i,j,k,5)=2.5d0
       else
          q(i,j,k,1)=0.125d0
          q(i,j,k,5)=0.25d0
       endif
          q(i,j,k,2)=0.d0
          q(i,j,k,3)=0.d0
          q(i,j,k,4)=0.d0
     enddo
   enddo
 enddo
end subroutine initSodShockTube3D

subroutine bdrySodShockTube3D(this,q)
class(grid)::this
double precision,dimension(1-this%nbuf:this%nMesh(1)+this%nbuf,1-this%nbuf:this%nMesh(2)+this%nbuf,&
& 1-this%nbuf:this%nMesh(3)+this%nbuf,this%nvar)::q
integer::nx,ny,nz,nbuf,i,j,k

nbuf=this%nbuf
nx=this%nMesh(1)
ny=this%nMesh(2)
nz=this%nMesh(3)

if(this%boundaryType .eq. 1) then
!!!! left boundary
   if(this%left_mpi .lt. 0) then
     do i=1,nbuf
       q(1-i,:,:,:)=q(1,:,:,:)
     enddo
   endif
!!!! right boundary
   if(this%right_mpi .lt. 0) then
     do i=1,nbuf
       q(nx+i,:,:,:)=q(nx,:,:,:)
     enddo
   endif
!!!! up boundary
   if(this%up_mpi .lt. 0) then
     do j=1,nbuf
       q(:,1-j,:,:)=q(:,1,:,:)
     enddo
   endif 
!!!! down boundary
   if(this%down_mpi .lt. 0) then
     do j=1,nbuf
       q(:,ny+j,:,:)=q(:,ny,:,:)
     enddo
   endif
!!!! top boundary
   if(this%top_mpi .lt. 0) then
     do k=1,nbuf
       q(:,:,1-k,:)=q(:,:,1,:)
     enddo
   endif
!!!! bottom boundary
   if(this%bottom_mpi .lt. 0) then
     do k=1,nbuf
       q(:,:,nz+k,:)=q(:,:,nz,:)
     enddo
   endif
!!!!
endif
end subroutine bdrySodShockTube3D

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!! Core collapse AD 2D (for the demonstration in the proposal)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine coreCollapseAD2D(gridIDn,gridIdi)
implicit none
integer::gridIDn,gridIDi
type(grid)::gn,gi
integer::ndim,nbuf,coordType,nvariable(8),ivariable(8)
integer::nMesh(2),dims(2)
double precision::leftBdry(2),rightBdry(2)
logical::periods(2),reorder
double precision::t0,t1,dt,tend
integer::nstep,ierr
!!! free fall time of a thin disk: c_s/(sqrt(pi)*Gravconst*\Sigma)
!!! Jeans length: c_s^2/(GravConst*\Sigma)
!!! Jeans mass: c_s^4/(GravConst^2*\Sigma)
!!! GravConst=4.3011e-3 [km^2*pc/(s^2*Msun)]
!!! [M]=Msun
!!! [L]=pc
!!! [T]=pc/km*s = 3.086e13 sec = 9.78e5 yrs
!!! critical mass-to-flux ratio: M/(\Phi_B)=1/(2*pi*sqrt(GravConst))=2.4271
!!! sound speed of neutrals: 0.344 km/s 
!!! sound spedd of ions: 0.096 km/s
!!! free fall time: 1.73e6 yr = 1.772 code unit
!!! B field = 15 uG = 5.13 code unit
nstep=0
nvariable=0
ivariable=0
ndim=2
nbuf=2
coordType=1
nvariable(1)=1
nvariable(2)=1
nvariable(3)=1
nvariable(8)=1

ivariable=1

nMesh(1)=1024
nMesh(2)=1024
leftBdry(1)=-1.d0
leftBdry(2)=-1.d0
rightBdry(1)=1.d0
rightBdry(2)=1.d0
dims(1)=1
dims(2)=nprocs
periods(1)=.true.
periods(2)=.true.
reorder=.true.
tend=1.772d0*2.d0/2.0

call gn%enableSelfgravity()
call gi%enableSelfgravity()
call gn%setGravConst(GravConst=4.3011d-3)
call gi%setGravConst(GravConst=4.3011d-3)
call gn%enableAD(enable_ad=.true.)
call gi%enableAD(enable_ad=.true.)
call gn%setADparams(mu_ad=7.d0/3.d0,alpha_ad=7.32d4)
call gi%setADparams(mu_ad=30.d0,alpha_ad=7.32d4)
call gn%setTopologyMPI(ndim,dims,periods,reorder)
call gi%setTopologyMPI(ndim,dims,periods,reorder)
call gn%setMesh(nMesh,leftBdry,rightBdry,nbuf,coordType,gridIDn)
call gi%setMesh(nMesh,leftBdry,rightBdry,nbuf,coordType,gridIDi)
call gn%setVariable(nvariable)
call gi%setVariable(ivariable)
call gn%setSelfgravityKernel()
call gi%setSelfgravityKernel()
call gn%setSGMPIWindows()
call gi%setSGMPIWindows()
call gn%setMPIWindows()
call gi%setMPIWIndows()
call gn%setEOS(eosType=2)
call gi%setEOS(eosType=2)
call gn%setadiGamma(gam=1.001d0)
call gi%setadiGamma(gam=5.d0/3.d0)
call gn%setCFL(CFL=0.4d0)
call gi%setCFL(CFL=0.4d0)
call gn%setTime(fstart=0,tend=tend,dtout=tend/50.d0)
call gi%setTime(fstart=0,tend=tend,dtout=tend/50.d0)
call gn%setSlopeLimiter(limiterType=3)
call gi%setSlopeLimiter(limiterType=3)
call gn%setSolverType(solverType=3)
call gi%setSolverType(solverType=5)
call gn%setBoundaryType(boundaryType=0)
call gi%setBoundaryType(boundaryType=0)
call gn%initVariable()
call gi%initVariable()
call gn%exchangeBdryMPI(gn%q,gn%winq)
call gi%exchangeBdryMPI(gi%q,gi%winq)
call gn%setBoundary(gn%q)
call gi%setBoundary(gi%q)
call gn%calcSelfgravity(gn%q)
call gi%calcSelfgravity(gi%q)
call gn%writeGrid()
call gi%writeGrid()
call gn%writeGrid_HD_vtk()
call gi%writeGrid_MHD_vtk()

t0=MPI_WTIME()
do while (gn%t .lt. gn%tend)
   call MPI_BARRIER(MPI_COMM_WORLD,ierr)
   call gn%griddt()
   call gi%griddt()
   call MPI_BARRIER(MPI_COMM_WORLD,ierr)
   dt=dmin1(gn%dt,gi%dt)
   if(dt .gt. gn%toutput-gn%t) then
     dt=gn%toutput-gn%t
     gn%writeFlag=.true.
     gi%writeFlag=.true.
     gn%toutput=gn%toutput+gn%dtout
     gi%toutput=gi%toutput+gi%dtout
     gn%fnum=gn%fnum+1
     gi%fnum=gi%fnum+1
   elseif (dt .gt. gn%tend-gn%t) then
     dt=gn%tend-gn%t
     gn%writeFlag=.true.
     gi%writeFlag=.true.
     gn%fnum=gn%fnum+1
     gi%fnum=gi%fnum+1
   endif
   gn%dt=dt
   gi%dt=dt
   call MPI_BARRIER(MPI_COMM_WORLD,ierr)
   !call gn%evolveGridRK2()
   !call gi%evolveGridRK2()
   call rk2ADsg_2D(gn,gn%q,gn%q1,gn%q2,gi,gi%q,gi%q1,gi%q2)
   call MPI_BARRIER(MPI_COMM_WORLD,ierr)
   !call evolveAD2D(gn,gi,gn%q,gi%q)
   call gn%exchangeBdryMPI(gn%q,gn%winq)
   call gi%exchangeBdryMPI(gi%q,gi%winq)
   call gn%setBoundary(gn%q)
   call gi%setBoundary(gi%q)
   if(myid .eq. 0) then
     print *,"testSuite.f03: gridID=",gn%gridID,'t=',gn%t,"dt=",gn%dt
   endif
   if(gn%writeFlag .eqv. .true.) then
     if(gi%writeFlag .eqv. .true.) then
       t1=MPI_WTIME()
       if(myid .eq. 0) then
          print *,"testSuiteMPI.f03: elapse time=",t1-t0,"secs"
       endif
       call gn%writeGrid()
       call gi%writeGrid()
       call gn%writeGrid_HD_vtk()
       call gi%writeGrid_MHD_vtk()
       gn%writeFlag=.false.
       gi%writeFlag=.false.
       t0=MPI_WTIME()
     endif
   endif
enddo !! end do while

end subroutine coreCollapseAD2D

subroutine initCoreCollapseAD2Dn(this,q)
implicit none
class(grid)::this
double precision,dimension(1-this%nbuf:this%nMesh(1)+this%nbuf,1-this%nbuf:this%nMesh(2)+this%nbuf,this%nvar)::q
integer::nx,ny,i,j
double precision::x,y,r,denc,rflat,snd,p

nx=this%nMesh(1)
ny=this%nMesh(2)
denc=25.465d0*3.d0 !! Msun/pc^2
rflat=0.7d0  !! pc
snd=0.344d0  !! km/s

do j=1,ny
  do i=1,nx
    x=this%xc(1)%coords(i)
    y=this%xc(2)%coords(j)
    r=dsqrt(x**2+y**2)
    q(i,j,1)=denc/(1.d0+(r/rflat)**2)
    q(i,j,2)=0.d0
    q(i,j,3)=0.d0
    p=snd**2*q(i,j,1)/this%adiGamma
    q(i,j,4)=p/(this%adiGamma-1.d0)
  enddo
enddo
end subroutine initCoreCollapseAD2Dn

subroutine initCoreCollapseAD2Di(this,q)
implicit none
class(grid)::this
double precision,dimension(1-this%nbuf:this%nMesh(1)+this%nbuf,1-this%nbuf:this%nMesh(2)+this%nbuf,this%nvar)::q
integer::nx,ny,i,j
double precision::x,y,r,denc,rflat,snd,p,b0

nx=this%nMesh(1)
ny=this%nMesh(2)
denc=25.465d0*3.d0/2.33d0*30.d0*1d-6 !! Msun/pc^2
rflat=0.7d0  !! pc
snd=0.096d0  !! km/s
b0=5.13d0/5.d0  !! code unit

do j=1,ny
  do i=1,nx
    x=this%xc(1)%coords(i)
    y=this%xc(2)%coords(j)
    r=dsqrt(x**2+y**2)
    q(i,j,1)=denc/(1.d0+(r/rflat)**2)
    q(i,j,2)=0.d0
    q(i,j,3)=0.d0
    q(i,j,4)=0.d0
    q(i,j,5)=0.d0
    q(i,j,6)=0.d0
    q(i,j,7)=b0*dsqrt(q(i,j,1)/denc)
    p=snd**2*q(i,j,1)/this%adiGamma
    q(i,j,8)=p/(this%adiGamma-1.d0)+0.5d0*(q(i,j,7)**2)
    q(i,j,9)=0.d0
    q(i,j,10)=0.d0
    q(i,j,11)=q(i,j,7)
  enddo
enddo
end subroutine initCoreCollapseAD2Di

subroutine bdryCoreCollapseAD2Dn(this,q)
implicit none
class(grid)::this
double precision,dimension(1-this%nbuf:this%nMesh(1)+this%nbuf,1-this%nbuf:this%nMesh(2)+this%nbuf,this%nvar)::q
integer::nx,ny,i,j,nbuf
double precision::x,y,r,denc,rflat,snd,p

nbuf=this%nbuf
nx=this%nMesh(1)
ny=this%nMesh(2)

denc=25.465d0*3.d0 !! Msun/pc^2
rflat=0.7d0  !! pc
snd=0.344d0  !! km/s

  !!!!!! left boundary
  if(this%left_mpi .lt. 0) then
    do j=1-nbuf,ny+nbuf
      do i=1,nbuf
        x=this%xc(1)%coords(1-i)
        y=this%xc(2)%coords(j)
        r=dsqrt(x**2+y**2)
        q(1-i,j,1)=denc/(1.d0+(r/rflat)**2)
        q(1-i,j,2)=0.d0
        q(1-i,j,3)=0.d0
        p=snd**2*q(1-i,j,1)/this%adiGamma
        q(1-i,j,4)=p/(this%adiGamma-1.d0)
      enddo
    enddo
  endif
!!!!!! right boundary
  if(this%right_mpi .lt. 0) then
    do j=1-nbuf,ny+nbuf
      do i=1,nbuf
        x=this%xc(1)%coords(nx+i)
        y=this%xc(2)%coords(j)
        r=dsqrt(x**2+y**2)
        q(nx+i,j,1)=denc/(1.d0+(r/rflat)**2)
        q(nx+i,j,2)=0.d0
        q(nx+i,j,3)=0.d0
        p=snd**2*q(nx+i,j,1)/this%adiGamma
        q(nx+i,j,4)=p/(this%adiGamma-1.d0)
      enddo
    enddo
  endif
!!!!!! up boundary
  if(this%up_mpi .lt. 0) then
    do j=1,nbuf
      do i=1-nbuf,nx+nbuf
        x=this%xc(1)%coords(i)
        y=this%xc(2)%coords(1-j)
        r=dsqrt(x**2+y**2)
        q(i,1-j,1)=denc/(1.d0+(r/rflat)**2)
        q(i,1-j,2)=0.d0
        q(i,1-j,3)=0.d0
        p=snd**2*q(i,1-j,1)/this%adiGamma
        q(i,1-j,4)=p/(this%adiGamma-1.d0)
      enddo
    enddo
  endif
!!!!!! down boundary
  if(this%down_mpi .lt. 0) then
    do j=1,nbuf
      do i=1-nbuf,nx+nbuf
        x=this%xc(1)%coords(i)
        y=this%xc(2)%coords(ny+j)
        r=dsqrt(x**2+y**2)
        q(i,ny+j,1)=denc/(1.d0+(r/rflat)**2)
        q(i,ny+j,2)=0.d0
        q(i,ny+j,3)=0.d0
        p=snd**2*q(i,ny+j,1)/this%adiGamma
        q(i,ny+j,4)=p/(this%adiGamma-1.d0)
      enddo
    enddo
  endif

end subroutine bdryCoreCollapseAD2Dn

subroutine bdryCoreCollapseAD2Di(this,q)
implicit none
class(grid)::this
double precision,dimension(1-this%nbuf:this%nMesh(1)+this%nbuf,1-this%nbuf:this%nMesh(2)+this%nbuf,this%nvar)::q
integer::nx,ny,i,j,nbuf
double precision::x,y,r,denc,rflat,snd,p,b0

nbuf=this%nbuf
nx=this%nMesh(1)
ny=this%nMesh(2)
denc=25.465d0*3.d0/2.33d0*30.d0*1d-6 !! Msun/pc^2
rflat=0.7d0  !! pc
snd=0.096d0  !! km/s
b0=5.13d0/5.d0  !! code unit

  !!!!!! left boundary
  if(this%left_mpi .lt. 0) then
    do j=1-nbuf,ny+nbuf
      do i=1,nbuf
        x=this%xc(1)%coords(1-i)
        y=this%xc(2)%coords(j)
        r=dsqrt(x**2+y**2)
        q(1-i,j,1)=denc/(1.d0+(r/rflat)**2)
        q(1-i,j,2)=0.d0
        q(1-i,j,3)=0.d0
        q(1-i,j,4)=0.d0
        q(1-i,j,5)=0.d0
        q(1-i,j,6)=0.d0
        q(1-i,j,7)=b0*dsqrt(q(1-i,j,1)/denc)
        p=snd**2*q(1-i,j,1)/this%adiGamma
        q(1-i,j,8)=p/(this%adiGamma-1.d0)+0.5d0*(q(1-i,j,7)**2)
        q(1-i,j,9)=0.d0
        q(1-i,j,10)=0.d0
        q(1-i,j,11)=q(1-i,j,7)
      enddo
    enddo
  endif
!!!!!! right boundary
  if(this%right_mpi .lt. 0) then
    do j=1-nbuf,ny+nbuf
      do i=1,nbuf
        x=this%xc(1)%coords(nx+i)
        y=this%xc(2)%coords(j)
        r=dsqrt(x**2+y**2)
        q(nx+i,j,1)=denc/(1.d0+(r/rflat)**2)
        q(nx+i,j,2)=0.d0
        q(nx+i,j,3)=0.d0
        q(nx+i,j,4)=0.d0
        q(nx+i,j,5)=0.d0
        q(nx+i,j,6)=0.d0
        q(nx+i,j,7)=b0*dsqrt(q(nx+i,j,1)/denc)
        p=snd**2*q(nx+i,j,1)/this%adiGamma
        q(nx+i,j,8)=p/(this%adiGamma-1.d0)+0.5d0*(q(nx+i,j,7)**2)
        q(nx+i,j,9)=0.d0
        q(nx+i,j,10)=0.d0
        q(nx+i,j,11)=q(nx+i,j,7)
      enddo
    enddo
  endif
!!!!!! up boundary
  if(this%up_mpi .lt. 0) then
    do j=1,nbuf
      do i=1-nbuf,nx+nbuf
        x=this%xc(1)%coords(i)
        y=this%xc(2)%coords(1-j)
        r=dsqrt(x**2+y**2)
        q(i,1-j,1)=denc/(1.d0+(r/rflat)**2)
        q(i,1-j,2)=0.d0
        q(i,1-j,3)=0.d0
        q(i,1-j,4)=0.d0
        q(i,1-j,5)=0.d0
        q(i,1-j,6)=0.d0
        q(i,1-j,7)=b0*dsqrt(q(i,1-j,1)/denc)
        p=snd**2*q(i,1-j,1)/this%adiGamma
        q(i,1-j,8)=p/(this%adiGamma-1.d0)+0.5d0*(q(i,1-j,7)**2)
        q(i,1-j,9)=0.d0
        q(i,1-j,10)=0.d0
        q(i,1-j,11)=q(i,1-j,7)
      enddo
    enddo
  endif
!!!!!! down boundary
  if(this%down_mpi .lt. 0) then
    do j=1,nbuf
      do i=1-nbuf,nx+nbuf
        x=this%xc(1)%coords(i)
        y=this%xc(2)%coords(ny+j)
        r=dsqrt(x**2+y**2)
        q(i,ny+j,1)=denc/(1.d0+(r/rflat)**2)
        q(i,ny+j,2)=0.d0
        q(i,ny+j,3)=0.d0
        q(i,ny+j,4)=0.d0
        q(i,ny+j,5)=0.d0
        q(i,ny+j,6)=0.d0
        q(i,ny+j,7)=b0*dsqrt(q(i,ny+j,1)/denc)
        p=snd**2*q(i,ny+j,1)/this%adiGamma
        q(i,ny+j,8)=p/(this%adiGamma-1.d0)+0.5d0*(q(i,ny+j,7)**2)
        q(i,ny+j,9)=0.d0
        q(i,ny+j,10)=0.d0
        q(i,ny+j,11)=q(i,ny+j,7)
      enddo
    enddo
  endif


end subroutine bdryCoreCollapseAD2Di

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!! Wardle instability
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine WardleInstability(gridIDn,gridIDi)
implicit none
integer::gridIDn,gridIDi
type(grid)::gn2,gi2,gn1,gi1
integer::ndim,nbuf,coordType,nvariable(8),ivariable(8)
integer::nMesh(2),dims(2)
double precision::leftBdry(2),rightBdry(2)
logical::periods(2),reorder
double precision::t0,t1,dt,tend
integer::nstep,ierr
integer::j,k,istart,iend,offset,icopy

call WardleInstabilityInit(gn1,gi1,gridIDn=24,gridIdi=25)

nstep=0
nvariable=0
ivariable=0
ndim=2
nbuf=2
coordType=1

nvariable(1)=1
nvariable(2)=1
nvariable(3)=1
nvariable(8)=1

ivariable=1

nMesh(1)=768
nMesh(2)=64
leftBdry(1)=0.d0
leftBdry(2)=0.d0
rightBdry(1)=0.00742d0
rightBdry(2)=dble(nMesh(2))/dble(nMesh(1))*rightBdry(1)
call MPI_DIMS_CREATE(nprocs,ndim,dims,ierr)
dims(1)=nprocs
dims(2)=1
periods(1)=.false.
periods(2)=.true.
reorder=.true.
tend=1000.d0/9.78d5

call gn2%enableAD(enable_ad=.true.)
call gi2%enableAD(enable_ad=.true.)
call gn2%setADparams(mu_ad=7.d0/3.d0,alpha_ad=7.31d4)
call gi2%setADparams(mu_ad=10.d0,alpha_ad=7.31d4)
call gn2%setTopologyMPI(ndim,dims,periods,reorder)
call gi2%setTopologyMPI(ndim,dims,periods,reorder)
call gn2%setMesh(nMesh,leftBdry,rightBdry,nbuf,coordType,gridIDn)
call gi2%setMesh(nMesh,leftBdry,rightBdry,nbuf,coordType,gridIDi)
call gn2%setVariable(nvariable)
call gi2%setVariable(ivariable)
call gn2%setMPIWindows()
call gi2%setMPIWindows()
call gn2%setEOS(eosType=2)
call gi2%setEOS(eosType=2)
call gn2%setadiGamma(gam=1.01d0)
call gi2%setadiGamma(gam=5.d0/3.d0)
call gn2%setCFL(CFL=0.4d0)
call gi2%setCFL(CFL=0.4d0)
call gn2%setTime(fstart=0,tend=tend,dtout=tend/100.d0)
call gi2%setTime(fstart=0,tend=tend,dtout=tend/100.d0)
call gn2%setSlopeLimiter(limiterType=3)
call gi2%setSlopeLimiter(limiterType=3)
call gn2%setSolverType(solverType=3)
call gi2%setSolverType(solverType=5)
call gn2%setBoundaryType(boundaryType=0)
call gi2%setBoundaryType(boundaryType=0)

!!!!!!!!!!   Neutrals !!!!!!!!!!!!!!!!
!!!! copy density
offset=0
icopy=300
do j=1,gn2%nMesh(2)+2*nbuf
  istart=1+(gn2%nMesh(1)+2*nbuf)*(j-1)+offset
  iend=(gn2%nMesh(1)+2*nbuf)*j+offset
  gn2%q(istart:iend)=gn1%q(1:gn1%nMesh(1)+2*nbuf)
  !gn2%q(istart:istart+icopy)=gn2%q(istart+icopy+1)
enddo
!!! copy momx
offset=(gn2%nMesh(1)+2*nbuf)*(gn2%nMesh(2)+2*nbuf)
do j=1,gn2%nMesh(2)+2*nbuf
  istart=1+(gn2%nMesh(1)+2*nbuf)*(j-1)+offset
  iend=(gn2%nMesh(1)+2*nbuf)*j+offset
  gn2%q(istart:iend)=gn1%q((gn1%nMesh(1)+2*nbuf)+1:2*(gn1%nMesh(1)+2*nbuf))
  !gn2%q(istart:istart+icopy)=gn2%q(istart+icopy+1)
enddo
!!! copy momy
offset=2*(gn2%nMesh(1)+2*nbuf)*(gn2%nMesh(2)+2*nbuf)
do j=1,gn2%nMesh(2)+2*nbuf
  istart=1+(gn2%nMesh(1)+2*nbuf)*(j-1)+offset
  iend=(gn2%nMesh(1)+2*nbuf)*j+offset
  gn2%q(istart:iend)=0.d0
enddo
!!! copy ene
offset=3*(gn2%nMesh(1)+2*nbuf)*(gn2%nMesh(2)+2*nbuf)
do j=1,gn2%nMesh(2)+2*nbuf
  istart=1+(gn2%nMesh(1)+2*nbuf)*(j-1)+offset
  iend=(gn2%nMesh(1)+2*nbuf)*j+offset
  gn2%q(istart:iend)=gn1%q(2*(gn1%nMesh(1)+2*nbuf)+1:3*(gn1%nMesh(1)+2*nbuf))
  !gn2%q(istart:istart+icopy)=gn2%q(istart+icopy+1)
enddo
!!!!!!!!!!!!   ions   !!!!!!!!!!!!!!!!!
do k=1,gi2%nvar
  offset=(k-1)*(gn2%nMesh(1)+2*nbuf)*(gn2%nMesh(2)+2*nbuf)
  do j=1,gn2%nMesh(2)+2*nbuf
    istart=1+(gn2%nMesh(1)+2*nbuf)*(j-1)+offset
    iend=(gn2%nMesh(1)+2*nbuf)*j+offset
    gi2%q(istart:iend)=gi1%q((k-1)*(gn1%nMesh(1)+2*nbuf)+1:k*(gn1%nMesh(1)+2*nbuf))
    !gi2%q(istart:istart+icopy)=gi2%q(istart+icopy+1)
  enddo
enddo


call gn2%initVariable()
call gi2%initVariable()
call gn2%exchangeBdryMPI(gn2%q,gn2%winq)
call gi2%exchangeBdryMPI(gi2%q,gi2%winq)
call gn2%setBoundary(gn2%q)
call gi2%setBoundary(gi2%q)
call gn2%writeGrid()
call gi2%writeGrid()

t0=MPI_WTIME()
do while (gn2%t .lt. gn2%tend)
   call MPI_BARRIER(MPI_COMM_WORLD,ierr)
   call gn2%griddt()
   call gi2%griddt()
   call MPI_BARRIER(MPI_COMM_WORLD,ierr)
   dt=dmin1(gn2%dt,gi2%dt)
   if(dt .gt. gn2%toutput-gn2%t) then
     dt=gn2%toutput-gn2%t
     gn2%writeFlag=.true.
     gi2%writeFlag=.true.
     gn2%toutput=gn2%toutput+gn2%dtout
     gi2%toutput=gi2%toutput+gi2%dtout
     gn2%fnum=gn2%fnum+1
     gi2%fnum=gi2%fnum+1
   elseif (dt .gt. gn2%tend-gn2%t) then
     dt=gn2%tend-gn2%t
     gn2%writeFlag=.true.
     gi2%writeFlag=.true.
     gn2%fnum=gn2%fnum+1
     gi2%fnum=gi2%fnum+1
   endif
   gn2%dt=dt
   gi2%dt=dt
   call MPI_BARRIER(MPI_COMM_WORLD,ierr)
   !call gn2%evolveGridRK2()
   !call gi2%evolveGridRK2()
   call rk2AD_2D(gn2,gn2%q,gn2%q1,gn2%q2,gi2,gi2%q,gi2%q1,gi2%q2)
   call MPI_BARRIER(MPI_COMM_WORLD,ierr)
   !call evolveAD2D(gn2,gi2,gn2%q,gi2%q)
   call gn2%exchangeBdryMPI(gn2%q,gn2%winq)
   call gi2%exchangeBdryMPI(gi2%q,gi2%winq)
   call gn2%setBoundary(gn2%q)
   call gi2%setBoundary(gi2%q)
   if(myid .eq. 0) then
     print *,"testSuite.f03: gridID=",gn2%gridID,'t=',gn2%t,"dt=",gn2%dt
   endif
   if(gn2%writeFlag .eqv. .true.) then
     if(gi2%writeFlag .eqv. .true.) then
       t1=MPI_WTIME()
       if(myid .eq. 0) then
          print *,"testSuiteMPI.f03: elapse time=",t1-t0,"secs"
       endif
       call gn2%writeGrid()
       call gi2%writeGrid()
       gn2%writeFlag=.false.
       gi2%writeFlag=.false.
       t0=MPI_WTIME()
     endif
   endif
enddo !! end do while
end subroutine WardleInstability

subroutine initWardleInstabilityn(this,q)
implicit none
class(grid)::this
double precision,dimension(1-this%nbuf:this%nMesh(1)+this%nbuf,1-this%nbuf:this%nMesh(2)+this%nbuf,this%nvar)::q
double precision::A,ky,x,y,x0,den,vy
integer::nx,ny,i,j

ky=10.167d3
nx=this%nMesh(1)
ny=this%nMesh(2)
A=1d-3
x0=6.53d-3

do j=1,ny
  do i=1,nx
     x=this%xc(1)%coords(i)
     y=this%xc(2)%coords(j)
     if(ky*(x-x0) .gt. 0 .and. ky*(x-x0) .lt. 3.1416d0) then
       den=q(i,j,1)
        vy=q(i,j,3)/den
        vy=vy+A*dcos(ky*y)*dsin(ky*(x-x0))
        q(i,j,3)=den*vy
     endif
  enddo
enddo
end subroutine initWardleInstabilityn

subroutine initWardleInstabilityi(this,q)
implicit none
class(grid)::this
double precision,dimension(1-this%nbuf:this%nMesh(1)+this%nbuf,1-this%nbuf:this%nMesh(2)+this%nbuf,this%nvar)::q
double precision::A,ky,x,y,x0,den,vy
integer::nx,ny,i,j

ky=10.167d3
nx=this%nMesh(1)
ny=this%nMesh(2)
A=1d-3
x0=6.53d-3

do j=1,ny
  do i=1,nx
     x=this%xc(1)%coords(i)
     y=this%xc(2)%coords(j)
     if(ky*(x-x0) .gt. 0 .and. ky*(x-x0) .lt. 3.1416d0) then
       den=q(i,j,1)
        vy=q(i,j,3)/den
        vy=vy+A*dcos(ky*y)*dsin(ky*(x-x0))
        q(i,j,3)=den*vy
     endif
  enddo
enddo
end subroutine initWardleInstabilityi

subroutine bdryWardleInstabilityn(this,q)
implicit none
class(grid)::this
double precision,dimension(1-this%nbuf:this%nMesh(1)+this%nbuf,1-this%nbuf:this%nMesh(2)+this%nbuf,this%nvar)::q
integer::i,j,nx,ny,nbuf
double precision::vx,vy,p,snd,den

nx=this%nMesh(1)
ny=this%nMesh(2)
nbuf=this%nbuf
den=3.45d0
vx=-7.79266d0
vy=0.d0
snd=0.344d0
p=snd**2*den/this%adiGamma

if(this%left_mpi .lt. 0) then
  do i=1,nbuf
    q(i-nbuf,:,:)=q(nbuf-i+1,:,:)
  enddo

  do i=1,nbuf
    !q(i-nbuf,:,2)=-q(nbuf-i+1,:,2)
  enddo
endif

if(this%right_mpi .lt. 0) then
  do i=1,nbuf
    q(nx+i,:,1)=den
    q(nx+i,:,2)=den*vx
    q(nx+i,:,3)=0.d0
    q(nx+i,:,4)=p/(this%adiGamma-1.d0)+0.5d0*den*(vx**2+vy**2)
  enddo
endif

if(this%up_mpi .lt. 0) then
  do j=1,nbuf
    q(:,j-nbuf,:)=q(:,1,:)
  enddo
endif

if(this%down_mpi .lt. 0) then
  do j=1,nbuf
    q(:,ny+j,:)=q(:,ny,:)
  enddo
endif
end subroutine bdryWardleInstabilityn

subroutine bdryWardleInstabilityi(this,q)
implicit none
class(grid)::this
double precision,dimension(1-this%nbuf:this%nMesh(1)+this%nbuf,1-this%nbuf:this%nMesh(2)+this%nbuf,this%nvar)::q
integer::i,j,nx,ny,nbuf
double precision::vx,vy,vz,den,bxl,byl,bzl,bxr,byr,bzr,bxc,byc,bzc,snd,p,pi

pi=datan2(1.d0,1.d0)*4.d0
nx=this%nMesh(1)
ny=this%nMesh(2)
nbuf=this%nbuf
den=7.4d-3
vx=-7.79266d0
vy=0.d0
vz=0.d0
bxl=0.d0
bzl=0.d0
byl=0.4949d0
bxr=0.d0
bzr=0.d0
byr=0.4949d0
snd=9.6d-2
p=snd**2*den/this%adiGamma
bxc=0.5d0*(bxl+bxr)
byc=0.5d0*(byl+byr)
bzc=0.5d0*(bzl+bzr)

if(this%left_mpi .lt. 0) then
  do i=1,nbuf
    q(i-nbuf,:,:)=q(nbuf-i+1,:,:)
  enddo
  do i=1,nbuf
    !q(i-nbuf,:,2)=-q(nbuf-i+1,:,2)
    !q(i-nbuf,:,5)=-q(nbuf-i+1,:,5)
    !q(i-nbuf,:,9)=-q(nbuf-i+1,:,9)
  enddo
endif

if(this%right_mpi .lt. 0) then
  do i=1,nbuf
    q(nx+i,:,1)=den
    q(nx+i,:,2)=den*vx
    q(nx+i,:,3)=den*vy
    q(nx+i,:,4)=den*vz
    q(nx+i,:,5)=bxl
    q(nx+i,:,6)=byl
    q(nx+i,:,7)=bzl
    q(nx+i,:,8)=p/(this%adiGamma-1.d0)+0.5d0*den*(vx**2+vy**2+vz**2)+0.5d0*(bxc**2+byc**2+bzc**2)
    q(nx+i,:,9)=bxr
    q(nx+i,:,10)=byr
    q(nx+i,:,11)=bzr
  enddo
endif


if(this%up_mpi .lt. 0) then
  do j=1,nbuf
    q(:,j-nbuf,:)=q(:,1,:)
  enddo
endif

if(this%down_mpi .lt. 0) then
  do j=1,nbuf
    q(:,ny+j,:)=q(:,ny,:)
  enddo
endif

end subroutine bdryWardleInstabilityi


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!! prepare 1D initial condition for Wardle instability
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine WardleInstabilityInit(gn,gi,gridIDn,gridIDi)
integer::gridIDn,gridIDi
type(grid)::gn,gi
integer::ndim,nbuf,coordType,nvariable(8),ivariable(8)
integer::nMesh(1),dims(1)
double precision::leftBdry(1),rightBdry(1)
logical::periods(1),reorder
double precision::t0,t1,dt,tend,tend0
integer::nstep,ierr

nstep=0
nvariable=0
ivariable=0
ndim=1
nbuf=2
coordType=1
nvariable(1)=1
nvariable(2)=1
nvariable(8)=1

ivariable=1

nMesh(1)=768
leftBdry(1)=0.d0  !! pc
rightBdry(1)=0.00742d0 !! pc

dims(1)=nprocs
periods(1)=.false.
reorder=.true.
tend0=600.d0/9.78d5
tend=0.01636d0*0.0037d0/0.00329d0

call gn%enableAD(enable_ad=.true.)
call gi%enableAD(enable_ad=.true.)
call gn%setADparams(mu_ad=7.d0/3.d0,alpha_ad=7.31d4)
call gi%setADparams(mu_ad=10.d0 ,alpha_ad=7.31d4)
call gn%setTopologyMPI(ndim,dims,periods,reorder)
call gi%setTopologyMPI(ndim,dims,periods,reorder)
call gn%setMesh(nMesh,leftBdry,rightBdry,nbuf,coordType,gridIDn)
call gi%setMesh(nMesh,leftBdry,rightBdry,nbuf,coordType,gridIDi)
call gn%setVariable(nvariable)
call gi%setVariable(ivariable)
call gn%setMPIWindows()
call gi%setMPIWIndows()
call gn%setEOS(eosType=2)
call gi%setEOS(eosType=2)
call gn%setadiGamma(gam=1.01d0)
call gi%setadiGamma(gam=5.d0/3.d0)
call gn%setCFL(CFL=0.4d0)
call gi%setCFL(CFL=0.4d0)
call gn%setTime(fstart=0,tend=tend,dtout=tend)
call gi%setTime(fstart=0,tend=tend,dtout=tend)
call gn%setSlopeLimiter(limiterType=3)
call gi%setSlopeLimiter(limiterType=3)
call gn%setSolverType(solverType=3)
call gi%setSolverType(solverType=5)
call gn%setBoundaryType(boundaryType=0)
call gi%setBoundaryType(boundaryType=0)
call gn%initVariable()
call gi%initVariable()
call gn%exchangeBdryMPI(gn%q,gn%winq)
call gi%exchangeBdryMPI(gi%q,gi%winq)
call gn%setBoundary(gn%q)
call gi%setBoundary(gi%q)
call gn%writeGrid()
call gi%writeGrid()

gn%writeFlag=.false.
gi%writeFlag=.false.

t0=MPI_WTIME()
do while (gn%t .lt. gn%tend)
   call MPI_BARRIER(MPI_COMM_WORLD,ierr)
   call gn%griddt()
   call gi%griddt()
   call MPI_BARRIER(MPI_COMM_WORLD,ierr)
   dt=dmin1(gn%dt,gi%dt)
   if(dt .gt. gn%toutput-gn%t) then
     dt=gn%toutput-gn%t
     gn%writeFlag=.true.
     gi%writeFlag=.true.
     gn%toutput=gn%toutput+gn%dtout
     gi%toutput=gi%toutput+gi%dtout
     gn%fnum=gn%fnum+1 
     gi%fnum=gi%fnum+1    
   elseif (dt .gt. gn%tend-gn%t) then
     dt=gn%tend-gn%t
     gn%writeFlag=.true.
     gi%writeFlag=.true.
     gn%fnum=gn%fnum+1
     gi%fnum=gi%fnum+1
   endif
   gn%dt=dt
   gi%dt=dt
   call MPI_BARRIER(MPI_COMM_WORLD,ierr)
   !call gn%evolveGridRK2()
   !call gi%evolveGridRK2()
   call rk2AD_1D(gn,gn%q,gn%q1,gn%q2,gi,gi%q,gi%q1,gi%q2)
   call MPI_BARRIER(MPI_COMM_WORLD,ierr)
   !call evolveAD1D(gn,gi,gn%q,gi%q)
   call gn%exchangeBdryMPI(gn%q,gn%winq)
   call gi%exchangeBdryMPI(gi%q,gi%winq)
   call gn%setBoundary(gn%q)
   call gi%setBoundary(gi%q)
   if(myid .eq. 0) then
     print *,"testSuite.f03: gridID=",gn%gridID,'t=',gn%t,"dt=",gn%dt
   endif
   if(gn%writeFlag .eqv. .true.) then
     if(gi%writeFlag .eqv. .true.) then
       t1=MPI_WTIME()
       if(myid .eq. 0) then
          print *,"testSuiteMPI.f03: elapse time=",t1-t0,"secs"
       endif
       call gn%writeGrid()
       call gi%writeGrid()
       gn%writeFlag=.false.
       gi%writeFlag=.false.
       t0=MPI_WTIME()
     endif
   endif
enddo !! end do while

end subroutine WardleInstabilityInit

subroutine initWardleInstabilityInitn(this,q)
implicit none
class(grid)::this
double precision,dimension(1-this%nbuf:this%nMesh(1)+this%nbuf,this%nvar)::q
integer::i,nx
double precision::vx,p,snd,den

nx=this%nMesh(1)
den=3.45d0 !! Msun/pc^3
vx=-7.79266d0 !! km/s
snd=0.344d0 !! km/s
p=snd**2*den/this%adiGamma
do i=1,nx
  q(i,1)=3.45d0
  q(i,2)=q(i,1)*vx
  q(i,3)=p/(this%adiGamma-1.d0)+0.5d0*den*(vx**2)
enddo
end subroutine initWardleInstabilityInitn

subroutine initWardleInstabilityIniti(this,q)
implicit none
class(grid)::this
double precision,dimension(1-this%nbuf:this%nMesh(1)+this%nbuf,this%nvar)::q
integer::i,nx
double precision::vx,vy,vz,p,snd,den,bxl,byl,bzl,bxr,byr,bzr
double precision::pi
double precision::bxc,byc,bzc
pi=4.d0*datan2(1.d0,1.d0)
nx=this%nMesh(1)
den=7.4d-3 !! Msun/pc^3
vx=-7.79266d0
vy=0.d0
vz=0.d0
bxl=0.d0
bzl=0.d0
byl=0.4949d0
bxr=0.d0
bzr=0.d0
byr=0.4949d0
snd=9.6d-2 !! km/s
p=snd**2*den/this%adiGamma
bxc=0.5d0*(bxl+bxr)
byc=0.5d0*(byl+byr)
bzc=0.5d0*(bzl+bzr)
do i=1,nx
  q(i,1)=den
  q(i,2)=den*vx
  q(i,3)=den*vy
  q(i,4)=den*vz
  q(i,5)=bxl
  q(i,6)=byl
  q(i,7)=bzl
  q(i,8)=p/(this%adiGamma-1.d0)+0.5d0*den*(vx**2+vy**2+vz**2)+0.5d0*(bxc**2+byc**2+bzc**2)
  q(i,9)=bxr
  q(i,10)=byr
  q(i,11)=bzr
enddo
end subroutine initWardleInstabilityIniti

subroutine bdryWardleInstabilityInitn(this,q)
implicit none
class(grid)::this
double precision,dimension(1-this%nbuf:this%nMesh(1)+this%nbuf,this%nvar)::q
integer::i,nx,nbuf
double precision::vx,p,snd,den

nx=this%nMesh(1)
nbuf=this%nbuf
den=3.45d0 
vx=-7.79266d0
snd=0.344d0
p=snd**2*den/this%adiGamma

if(this%left_mpi .lt. 0) then
  do i=1,nbuf
    q(i-nbuf,:)=q(nbuf-i+1,:)
  enddo
  
  do i=1,nbuf
    q(i-nbuf,2)=-q(nbuf-i+1,2)
  enddo
endif

if(this%right_mpi .lt. 0) then
  do i=1,nbuf
    q(nx+i,1)=den
    q(nx+i,2)=den*vx
    q(nx+i,3)=p/(this%adiGamma-1.d0)+0.5d0*den*(vx**2)
  enddo
endif
end subroutine bdryWardleInstabilityInitn

subroutine bdryWardleInstabilityIniti(this,q)
implicit none
class(grid)::this
double precision,dimension(1-this%nbuf:this%nMesh(1)+this%nbuf,this%nvar)::q
integer::i,nx,nbuf
double precision::vx,vy,vz,den,bxl,byl,bzl,bxr,byr,bzr,bxc,byc,bzc,snd,p,pi

pi=datan2(1.d0,1.d0)*4.d0
nx=this%nMesh(1)
nbuf=this%nbuf
den=7.4d-3
vx=-7.79266d0
vy=0.d0
vz=0.d0
bxl=0.d0
bzl=0.d0
byl=0.4949d0
bxr=0.d0
bzr=0.d0
byr=0.4949d0
snd=9.6d-2
p=snd**2*den/this%adiGamma
bxc=0.5d0*(bxl+bxr)
byc=0.5d0*(byl+byr)
bzc=0.5d0*(bzl+bzr)

if(this%left_mpi .lt. 0) then
  do i=1,nbuf
    q(i-nbuf,:)=q(nbuf-i+1,:)
  enddo
  do i=1,nbuf
    q(i-nbuf,2)=-q(nbuf-i+1,2)
    q(i-nbuf,5)=-q(nbuf-i+1,5)
    q(i-nbuf,9)=-q(nbuf-i+1,9)
  enddo
endif

if(this%right_mpi .lt. 0) then
  do i=1,nbuf
    q(nx+i,1)=den
    q(nx+i,2)=den*vx
    q(nx+i,3)=den*vy
    q(nx+i,4)=den*vz
    q(nx+i,5)=bxl
    q(nx+i,6)=byl
    q(nx+i,7)=bzl
    q(nx+i,8)=p/(this%adiGamma-1.d0)+0.5d0*den*(vx**2+vy**2+vz**2)+0.5d0*(bxc**2+byc**2+bzc**2)
    q(nx+i,9)=bxr
    q(nx+i,10)=byr
    q(nx+i,11)=bzr
  enddo
endif

end subroutine bdryWardleInstabilityIniti




!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!! C-Shock test 1D
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine CShockTest1D(gridIDn,gridIDi)
integer::gridIDn,gridIDi
type(grid)::gn,gi
integer::ndim,nbuf,coordType,nvariable(8),ivariable(8)
integer::nMesh(1),dims(1)
double precision::leftBdry(1),rightBdry(1)
logical::periods(1),reorder
double precision::t0,t1,dt,tend
integer::nstep,ierr

nstep=0
nvariable=0
ivariable=0
ndim=1
nbuf=2
coordType=1
nvariable(1)=1
nvariable(2)=1
nvariable(8)=1

ivariable=1

nMesh(1)=512
leftBdry(1)=0.d0  !! pc
rightBdry(1)=0.02d0 !! pc

dims(1)=nprocs
periods(1)=.false.
reorder=.true.
tend=0.01636d0

call gn%enableAD(enable_ad=.true.)
call gi%enableAD(enable_ad=.true.)
call gn%setADparams(mu_ad=2.33d0,alpha_ad=7.31d4)
call gi%setADparams(mu_ad=30.d0 ,alpha_ad=7.31d4)
call gn%setTopologyMPI(ndim,dims,periods,reorder)
call gi%setTopologyMPI(ndim,dims,periods,reorder)
call gn%setMesh(nMesh,leftBdry,rightBdry,nbuf,coordType,gridIDn)
call gi%setMesh(nMesh,leftBdry,rightBdry,nbuf,coordType,gridIDi)
call gn%setVariable(nvariable)
call gi%setVariable(ivariable)
call gn%setMPIWindows()
call gi%setMPIWIndows()
call gn%setEOS(eosType=2)
call gi%setEOS(eosType=2)
call gn%setadiGamma(gam=1.0001d0)
call gi%setadiGamma(gam=5.d0/3.d0)
call gn%setCFL(CFL=0.2d0)
call gi%setCFL(CFL=0.2d0)
call gn%setTime(fstart=0,tend=tend,dtout=tend/10.d0)
call gi%setTime(fstart=0,tend=tend,dtout=tend/10.d0)
call gn%setSlopeLimiter(limiterType=3)
call gi%setSlopeLimiter(limiterType=3)
call gn%setSolverType(solverType=3)
call gi%setSolverType(solverType=5)
call gn%setBoundaryType(boundaryType=0)
call gi%setBoundaryType(boundaryType=0)
call gn%initVariable()
call gi%initVariable()
call gn%exchangeBdryMPI(gn%q,gn%winq)
call gi%exchangeBdryMPI(gi%q,gi%winq)
call gn%setBoundary(gn%q)
call gi%setBoundary(gi%q)
call gn%writeGrid()
call gi%writeGrid()

gn%writeFlag=.false.
gi%writeFlag=.false.

t0=MPI_WTIME()
do while (gn%t .lt. gn%tend)
   call MPI_BARRIER(MPI_COMM_WORLD,ierr)
   call gn%griddt()
   call gi%griddt()
   call MPI_BARRIER(MPI_COMM_WORLD,ierr)
   dt=dmin1(gn%dt,gi%dt)
   if(dt .gt. gn%toutput-gn%t) then
     dt=gn%toutput-gn%t
     gn%writeFlag=.true.
     gi%writeFlag=.true.
     gn%toutput=gn%toutput+gn%dtout
     gi%toutput=gi%toutput+gi%dtout
     gn%fnum=gn%fnum+1 
     gi%fnum=gi%fnum+1    
   elseif (dt .gt. gn%tend-gn%t) then
     dt=gn%tend-gn%t
     gn%writeFlag=.true.
     gi%writeFlag=.true.
     gn%fnum=gn%fnum+1
     gi%fnum=gi%fnum+1
   endif
   gn%dt=dt
   gi%dt=dt
   call MPI_BARRIER(MPI_COMM_WORLD,ierr)
   !call gn%evolveGridRK2()
   !call gi%evolveGridRK2()
   call rk2AD_1D(gn,gn%q,gn%q1,gn%q2,gi,gi%q,gi%q1,gi%q2)
   call MPI_BARRIER(MPI_COMM_WORLD,ierr)
   !call evolveAD1D(gn,gi,gn%q,gi%q)
   call gn%exchangeBdryMPI(gn%q,gn%winq)
   call gi%exchangeBdryMPI(gi%q,gi%winq)
   call gn%setBoundary(gn%q)
   call gi%setBoundary(gi%q)
   if(myid .eq. 0) then
     print *,"testSuite.f03: gridID=",gn%gridID,'t=',gn%t,"dt=",gn%dt
   endif
   if(gn%writeFlag .eqv. .true.) then
     if(gi%writeFlag .eqv. .true.) then
       t1=MPI_WTIME()
       if(myid .eq. 0) then
          print *,"testSuiteMPI.f03: elapse time=",t1-t0,"secs"
       endif
       call gn%writeGrid()
       call gi%writeGrid()
       gn%writeFlag=.false.
       gi%writeFlag=.false.
       t0=MPI_WTIME()
     endif
   endif
enddo !! end do while

end subroutine CShockTest1D

subroutine initCShockTest1Dn(this,q)
implicit none
class(grid)::this
double precision,dimension(1-this%nbuf:this%nMesh(1)+this%nbuf,this%nvar)::q
integer::i,nx
double precision::vx,p,snd,den

nx=this%nMesh(1)
den=3.45d0 !! Msun/pc^3
vx=-2.2605d0 !! km/s
snd=0.344d0 !! km/s
p=snd**2*den/this%adiGamma
do i=1,nx
  q(i,1)=3.45d0
  q(i,2)=q(i,1)*vx
  q(i,3)=p/(this%adiGamma-1.d0)+0.5d0*den*(vx**2)
enddo
end subroutine initCShockTest1Dn

subroutine initCShockTest1Di(this,q)
implicit none
class(grid)::this
double precision,dimension(1-this%nbuf:this%nMesh(1)+this%nbuf,this%nvar)::q
integer::i,nx
double precision::vx,vy,vz,p,snd,den,bxl,byl,bzl,bxr,byr,bzr
double precision::pi
double precision::bxc,byc,bzc
pi=4.d0*datan2(1.d0,1.d0)
nx=this%nMesh(1)
den=7.4d-3 !! Msun/pc^3
vx=-2.2605d0
vy=0.d0
vz=0.d0
bxl=0.d0
bzl=0.d0
byl=1.71d0
bxr=0.d0
bzr=0.d0
byr=1.71d0
snd=9.6d-2 !! km/s
p=snd**2*den/this%adiGamma
bxc=0.5d0*(bxl+bxr)
byc=0.5d0*(byl+byr)
bzc=0.5d0*(bzl+bzr)
do i=1,nx
  q(i,1)=den
  q(i,2)=den*vx
  q(i,3)=den*vy
  q(i,4)=den*vz
  q(i,5)=bxl
  q(i,6)=byl
  q(i,7)=bzl
  q(i,8)=p/(this%adiGamma-1.d0)+0.5d0*den*(vx**2+vy**2+vz**2)+0.5d0*(bxc**2+byc**2+bzc**2)
  q(i,9)=bxr
  q(i,10)=byr
  q(i,11)=bzr
enddo
end subroutine initCShockTest1Di

subroutine bdryCShockTest1Dn(this,q)
implicit none
class(grid)::this
double precision,dimension(1-this%nbuf:this%nMesh(1)+this%nbuf,this%nvar)::q
integer::i,nx,nbuf
double precision::vx,p,snd,den

nx=this%nMesh(1)
nbuf=this%nbuf
den=3.45d0 
vx=-2.2605d0
snd=0.344d0
p=snd**2*den/this%adiGamma

if(this%left_mpi .lt. 0) then
  do i=1,nbuf
    q(i-nbuf,:)=q(nbuf-i+1,:)
  enddo
  
  do i=1,nbuf
    q(i-nbuf,2)=-q(nbuf-i+1,2)
  enddo
endif

if(this%right_mpi .lt. 0) then
  do i=1,nbuf
    q(nx+i,1)=den
    q(nx+i,2)=den*vx
    q(nx+i,3)=p/(this%adiGamma-1.d0)+0.5d0*den*(vx**2)
  enddo
endif
end subroutine bdryCShockTest1Dn

subroutine bdryCShockTest1Di(this,q)
implicit none
class(grid)::this
double precision,dimension(1-this%nbuf:this%nMesh(1)+this%nbuf,this%nvar)::q
integer::i,nx,nbuf
double precision::vx,vy,vz,den,bxl,byl,bzl,bxr,byr,bzr,bxc,byc,bzc,snd,p,pi

pi=datan2(1.d0,1.d0)*4.d0
nx=this%nMesh(1)
nbuf=this%nbuf
den=7.4d-3
vx=-2.2605d0
vy=0.d0
vz=0.d0
bxl=0.d0
bzl=0.d0
byl=1.71d0
bxr=0.d0
bzr=0.d0
byr=1.71d0
snd=9.6d-2
p=snd**2*den/this%adiGamma
bxc=0.5d0*(bxl+bxr)
byc=0.5d0*(byl+byr)
bzc=0.5d0*(bzl+bzr)

if(this%left_mpi .lt. 0) then
  do i=1,nbuf
    q(i-nbuf,:)=q(nbuf-i+1,:)
  enddo
  do i=1,nbuf
    q(i-nbuf,2)=-q(nbuf-i+1,2)
    !q(i-nbuf,5)=-q(nbuf-i+1,5)
    !q(i-nbuf,9)=-q(nbuf-i+1,9)
  enddo
endif

if(this%right_mpi .lt. 0) then
  do i=1,nbuf
    q(nx+i,1)=den
    q(nx+i,2)=den*vx
    q(nx+i,3)=den*vy
    q(nx+i,4)=den*vz
    q(nx+i,5)=bxl
    q(nx+i,6)=byl
    q(nx+i,7)=bzl
    q(nx+i,8)=p/(this%adiGamma-1.d0)+0.5d0*den*(vx**2+vy**2+vz**2)+0.5d0*(bxc**2+byc**2+bzc**2)
    q(nx+i,9)=bxr
    q(nx+i,10)=byr
    q(nx+i,11)=bzr
  enddo
endif

end subroutine bdryCShockTest1Di




!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!! C-Shock test 2D
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine CShockTest2D(gridIDn,gridIDi)
implicit none
integer::gridIDn,gridIDi
type(grid)::gn,gi
integer::ndim,nbuf,coordType,nvariable(8),ivariable(8)
integer::nMesh(2),dims(2)
double precision::leftBdry(2),rightBdry(2)
logical::periods(2),reorder
double precision::t0,t1,dt,tend
integer::nstep,ierr

nstep=0
nvariable=0
ivariable=0
ndim=2
nbuf=2
coordType=1
nvariable(1)=1
nvariable(2)=1
nvariable(3)=1
nvariable(8)=1

ivariable=1

nMesh(1)=512
nMesh(2)=4
leftBdry(1)=0.d0  !! pc
rightBdry(1)=0.02d0 !! pc
leftBdry(2)=0.d0  !! pc
rightBdry(2)=0.02d0/128.0 !! pc

dims(1)=nprocs
dims(2)=1

periods(1)=.false.
periods(2)=.true.

reorder=.true.
tend=0.01636d0

call gn%enableAD(enable_ad=.true.)
call gi%enableAD(enable_ad=.true.)

call gn%setADparams(mu_ad=2.33d0,alpha_ad=7.31d4)
call gi%setADparams(mu_ad=30.d0 ,alpha_ad=7.31d4)

call gn%setTopologyMPI(ndim,dims,periods,reorder)
call gi%setTopologyMPI(ndim,dims,periods,reorder)

call gn%setMesh(nMesh,leftBdry,rightBdry,nbuf,coordType,gridIDn)
call gi%setMesh(nMesh,leftBdry,rightBdry,nbuf,coordType,gridIDi)

call gn%setVariable(nvariable)
call gi%setVariable(ivariable)

call gn%setMPIWindows()
call gi%setMPIWIndows()

call gn%setEOS(eosType=2)
call gi%setEOS(eosType=2)

call gn%setadiGamma(gam=1.001d0)
call gi%setadiGamma(gam=5.d0/3.d0)

call gn%setCFL(CFL=0.2d0)
call gi%setCFL(CFL=0.2d0)

call gn%setTime(fstart=0,tend=tend,dtout=tend/1.d0)
call gi%setTime(fstart=0,tend=tend,dtout=tend/1.d0)

call gn%setSlopeLimiter(limiterType=1)
call gi%setSlopeLimiter(limiterType=1)

call gn%setSolverType(solverType=3)
call gi%setSolverType(solverType=5)

call gn%setBoundaryType(boundaryType=0)
call gi%setBoundaryType(boundaryType=0)

call gn%initVariable()
call gi%initVariable()

call gn%exchangeBdryMPI(gn%q,gn%winq)
call gi%exchangeBdryMPI(gi%q,gi%winq)

call gn%setBoundary(gn%q)
call gi%setBoundary(gi%q)

call gn%writeGrid()
call gi%writeGrid()

call gn%writeGrid_HD_vtk()
call gi%writeGrid_MHD_vtk()

gn%writeFlag=.false.
gi%writeFlag=.false.

t0=MPI_WTIME()
do while (gn%t .lt. gn%tend)
   call MPI_BARRIER(MPI_COMM_WORLD,ierr)
   call gn%griddt()
   call gi%griddt()
   call MPI_BARRIER(MPI_COMM_WORLD,ierr)
   dt=dmin1(gn%dt,gi%dt)
   if(dt .gt. gn%toutput-gn%t) then
     dt=gn%toutput-gn%t
     gn%writeFlag=.true.
     gi%writeFlag=.true.
     gn%toutput=gn%toutput+gn%dtout
     gi%toutput=gi%toutput+gi%dtout
     gn%fnum=gn%fnum+1 
     gi%fnum=gi%fnum+1    
   elseif (dt .gt. gn%tend-gn%t) then
     dt=gn%tend-gn%t
     gn%writeFlag=.true.
     gi%writeFlag=.true.
     gn%fnum=gn%fnum+1
     gi%fnum=gi%fnum+1
   endif
   gn%dt=dt
   gi%dt=dt
   call MPI_BARRIER(MPI_COMM_WORLD,ierr)
   !call gn%evolveGridRK2()
   !call gi%evolveGridRK2()
   call rk2AD_2D_HSHSMD(gn,gn%q,gn%q1,gn%q2,gi,gi%q,gi%q1,gi%q2)
   call MPI_BARRIER(MPI_COMM_WORLD,ierr)
   !call evolveAD1D(gn,gi,gn%q,gi%q)
   call gn%exchangeBdryMPI(gn%q,gn%winq)
   call gi%exchangeBdryMPI(gi%q,gi%winq)
   call gn%setBoundary(gn%q)
   call gi%setBoundary(gi%q)
   if(myid .eq. 0) then
     print *,"testSuite.f03: gridID=",gn%gridID,'t=',gn%t,"dt=",gn%dt
   endif
   if(gn%writeFlag .eqv. .true.) then
     if(gi%writeFlag .eqv. .true.) then
       t1=MPI_WTIME()
       if(myid .eq. 0) then
          print *,"testSuiteMPI.f03: elapse time=",t1-t0,"secs"
       endif
       call gn%writeGrid()
       call gi%writeGrid()

       call gn%writeGrid_HD_vtk()
       call gi%writeGrid_MHD_vtk()

       gn%writeFlag=.false.
       gi%writeFlag=.false.
       t0=MPI_WTIME()
     endif
   endif
enddo !! end do while

end subroutine CShockTest2D

subroutine initCShockTest2Dn(this,q)
implicit none
class(grid)::this
double precision,dimension(1-this%nbuf:this%nMesh(1)+this%nbuf,1-this%nbuf:this%nMesh(2)+this%nbuf,this%nvar)::q
integer::i,j,nx,ny
double precision::vx,p,snd,den

nx=this%nMesh(1)
ny=this%nMesh(2)

den=3.45d0 !! Msun/pc^3
vx=-2.2605d0 !! km/s

snd=0.344d0 !! km/s
 !snd=0.2667  !!!km/s

p=snd**2*den/this%adiGamma

do j=1,ny
do i=1,nx
  q(i,j,1)=3.45d0
  q(i,j,2)=q(i,j,1)*vx
  q(i,j,3)=0.0
  q(i,j,4)=p/(this%adiGamma-1.d0)+0.5d0*den*(vx**2)
enddo
enddo
end subroutine initCShockTest2Dn

subroutine initCShockTest2Di(this,q)
implicit none
class(grid)::this
double precision,dimension(1-this%nbuf:this%nMesh(1)+this%nbuf,1-this%nbuf:this%nMesh(2)+this%nbuf,this%nvar)::q
integer::nx,ny,i,j

double precision::vx,vy,vz,p,snd,den,bxl,byl,bzl,bxr,byr,bzr
double precision::pi
double precision::bxc,byc,bzc
pi=4.d0*datan2(1.d0,1.d0)

nx=this%nMesh(1)
ny=this%nMesh(2)
den=7.405d-3 !! Msun/pc^3
vx=-2.2605d0
vy=0.d0
vz=0.d0
bxl=0.d0
bzl=0.d0
byl=1.71d0
bxr=0.d0
bzr=0.d0
byr=1.71d0

snd=9.6d-2 !! km/s
 !snd=7.432d-2 !! km/s

p=snd**2*den/this%adiGamma
bxc=0.5d0*(bxl+bxr)
byc=0.5d0*(byl+byr)
bzc=0.5d0*(bzl+bzr)

do j=1,ny
do i=1,nx
  q(i,j,1)=den
  q(i,j,2)=den*vx
  q(i,j,3)=den*vy
  q(i,j,4)=den*vz
  q(i,j,5)=bxl
  q(i,j,6)=byl
  q(i,j,7)=bzl
  q(i,j,8)=p/(this%adiGamma-1.d0)+0.5d0*den*(vx**2+vy**2+vz**2)+0.5d0*(bxc**2+byc**2+bzc**2)
  q(i,j,9)=bxr
  q(i,j,10)=byr
  q(i,j,11)=bzr
enddo
enddo
end subroutine initCShockTest2Di

subroutine bdryCShockTest2Dn(this,q)
implicit none
class(grid)::this
double precision,dimension(1-this%nbuf:this%nMesh(1)+this%nbuf,1-this%nbuf:this%nMesh(2)+this%nbuf,this%nvar)::q
integer::i,j,nx,ny,nbuf
double precision::vx,p,snd,den

nx=this%nMesh(1)
ny=this%nMesh(2)

!!!!!!!!!!!!!!!!!!!!!!
nbuf=this%nbuf
den=3.45d0 
vx=-2.2605d0
snd=0.344d0
p=snd**2*den/this%adiGamma

if(this%left_mpi .lt. 0) then
  do i=1,nbuf
    q(i-nbuf,:,:)=q(nbuf-i+1,:,:)
  enddo
  
  do i=1,nbuf
    q(i-nbuf,:,2)=-q(nbuf-i+1,:,2)
  enddo
endif

if(this%right_mpi .lt. 0) then
  do i=1,nbuf
    q(nx+i,:,1)=den
    q(nx+i,:,2)=den*vx
    q(nx+i,:,3)=0.0
    q(nx+i,:,4)=p/(this%adiGamma-1.d0)+0.5d0*den*(vx**2)
  enddo
endif
end subroutine bdryCShockTest2Dn

subroutine bdryCShockTest2Di(this,q)
implicit none
class(grid)::this
double precision,dimension(1-this%nbuf:this%nMesh(1)+this%nbuf,1-this%nbuf:this%nMesh(2)+this%nbuf,this%nvar)::q
integer::i,j,nx,ny,nbuf
double precision::vx,vy,vz,den,bxl,byl,bzl,bxr,byr,bzr,bxc,byc,bzc,snd,p,pi

pi=datan2(1.d0,1.d0)*4.d0
nx=this%nMesh(1)
ny=this%nMesh(2)
nbuf=this%nbuf
den=7.4d-3
vx=-2.2605d0
vy=0.d0
vz=0.d0
bxl=0.d0
bzl=0.d0
byl=1.71d0
bxr=0.d0
bzr=0.d0
byr=1.71d0

!snd=9.6d-2
  snd=7.432d-2 !! km/s
p=snd**2*den/this%adiGamma
bxc=0.5d0*(bxl+bxr)
byc=0.5d0*(byl+byr)
bzc=0.5d0*(bzl+bzr)

if(this%left_mpi .lt. 0) then
  do i=1,nbuf
    q(i-nbuf,:,:)=q(nbuf-i+1,:,:)
  enddo
  do i=1,nbuf
    q(i-nbuf,:,2)=-q(nbuf-i+1,:,2)
  enddo
endif

if(this%right_mpi .lt. 0) then
  do i=1,nbuf
    q(nx+i,:,1)=den
    q(nx+i,:,2)=den*vx
    q(nx+i,:,3)=den*vy
    q(nx+i,:,4)=den*vz
    q(nx+i,:,5)=bxl
    q(nx+i,:,6)=byl
    q(nx+i,:,7)=bzl
    q(nx+i,:,8)=p/(this%adiGamma-1.d0)+0.5d0*den*(vx**2+vy**2+vz**2)+0.5d0*(bxc**2+byc**2+bzc**2)
    q(nx+i,:,9)=bxr
    q(nx+i,:,10)=byr
    q(nx+i,:,11)=bzr
  enddo
endif

end subroutine bdryCShockTest2Di


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!! Selfgravity Test 2D
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine selfgravityTest2D(gridID)
integer::gridID
type(grid)::g1
integer::ndim,nbuf,coordType,variable(8)
integer::nMesh(2),dims(2)
double precision::leftBdry(2),rightBdry(2)
logical::periods(2),reorder
double precision::t0,t1
integer::nstep,ierr

nstep=0
variable=0
ndim=2
nbuf=2
coordType=1
variable(1)=1
variable(2)=1
variable(3)=1
variable(8)=1

nMesh(1)=128
nMesh(2)=128
leftBdry(1)=0.d0
leftBdry(2)=0.d0
rightBdry(1)=1.d0
rightBdry(2)=1.d0

dims(1)=1
dims(2)=nprocs
periods(1)=.true.
periods(2)=.true.
reorder=.true.


call g1%enableSelfgravity()
call g1%setGravConst(GravConst=1.d0)
call g1%setTopologyMPI(ndim,dims,periods,reorder)
call g1%setMesh(nMesh,leftBdry,rightBdry,nbuf,coordType,gridID)
call g1%setVariable(variable)
call g1%setSelfgravityKernel()
call g1%setSGMPIWindows()
call g1%setMPIWindows()
call g1%setEOS(eosType=2)
call g1%setadiGamma(gam=1.4d0)
call g1%setCFL(CFL=0.9d0)
call g1%setTime(fstart=0,tend=0.1d0,dtout=0.1d0)
call g1%setSlopeLimiter(limiterType=3)
call g1%setSolverType(solverType=3)
call g1%setBoundaryType(boundaryType=3)
call g1%initVariable()
call g1%exchangeBdryMPI(g1%q,g1%winq)
call g1%setBoundary(g1%q)
call g1%calcSelfgravity(g1%q)
call g1%writeGrid()

end subroutine selfgravityTest2D

subroutine initSelfgravityTest2D(this,q)
class(grid)::this
double precision,dimension(1-this%nbuf:this%nMesh(1)+this%nbuf,1-this%nbuf:this%nMesh(2)+this%nbuf,this%nvar)::q
double precision,dimension(1-this%nbuf:this%nMesh(1)+this%nbuf)::xc
double precision,dimension(1-this%nbuf:this%nMesh(2)+this%nbuf)::yc
integer::nx,ny,nvar,nbuf,gridID,eosType
integer::i,j

xc=this%xc(1)%coords
yc=this%xc(2)%coords
nx=this%nMesh(1)
ny=this%nMesh(2)

 do j=1,ny
   do i=1,nx
         q(i,j,1)=dble(myid)
         q(i,j,2)=0.d0
         q(i,j,3)=0.d0
         q(i,j,4)=1.d0
   enddo
 enddo
end subroutine initSelfgravityTest2D

subroutine bdrySelfgravityTest2D(this,q)
class(grid)::this
double precision,dimension(1-this%nbuf:this%nMesh(1)+this%nbuf,1-this%nbuf:this%nMesh(2)+this%nbuf,this%nvar)::q

if(this%boundaryType .eq. 3) then
endif
end subroutine bdrySelfgravityTest2D

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!! Woodward-Colella blast wave problem (1984, JCP, 54, 115)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine WCShockTube2D(gridID)
integer::gridID
type(grid)::g1
integer::ndim,nbuf,coordType,variable(8)
integer::nMesh(2),dims(2)
double precision::leftBdry(2),rightBdry(2)
logical::periods(2),reorder
double precision::t0,t1
integer::nstep,ierr

nstep=0
variable=0
ndim=2
nbuf=2
coordType=1
variable(1)=1
variable(2)=1
variable(3)=1
variable(8)=1

nMesh(1)=7
nMesh(2)=2400
leftBdry(1)=0.d0
leftBdry(2)=0.d0
rightBdry(1)=dble(nMesh(1))/dble(nMesh(2))
rightBdry(2)=1.d0
dims=(/0,0/)
call MPI_DIMS_CREATE(nprocs,ndim,dims,ierr)
periods(1)=.false.
periods(2)=.false.
reorder=.true.

call g1%setTopologyMPI(ndim,dims,periods,reorder)
call g1%setMesh(nMesh,leftBdry,rightBdry,nbuf,coordType,gridID)
call g1%setVariable(variable) !! den,vx,vy,ene
call g1%setMPIWindows()
call g1%setEoS(eosType=2) !! adiabatic
call g1%setadiGamma(gam=1.4d0) !! ratio of heat capacity Cp/Cv
call g1%setCFL(CFL=0.9d0) !! Courant number
call g1%setTime(fstart=0,tend=0.038d0,dtout=0.038d0)
call g1%setSlopeLimiter(limiterType=3)
call g1%setSolverType(solverType=3)
call g1%setBoundaryType(boundaryType=2)
call g1%initVariable()
call g1%exchangeBdryMPI(g1%q,g1%winq)
call g1%setBoundary(g1%q)
call g1%writeGrid()
g1%writeFlag=.false.

t0=MPI_WTIME()
do while (g1%t .lt. g1%tend)
   call g1%griddt()
   call g1%evolveGridRK2()
   if(myid .eq. 0) then
     print *,"testSuite.f03: gridID=",g1%gridID,'t=',g1%t,"dt=",g1%dt
   endif
   if(g1%writeFlag .eqv. .true.) then
      t1=MPI_WTIME()
      if(myid .eq. 0) then
         print *,"testSuiteMPI.f03: elapse time=",t1-t0,"secs"
      endif
      call g1%writeGrid()
      g1%writeFlag=.false.
      t0=MPI_WTIME()
   endif
enddo !! end do while

end subroutine WCShockTube2D

subroutine initWCShockTube2D(this,q)
class(grid)::this
double precision,dimension(1-this%nbuf:this%nMesh(1)+this%nbuf,1-this%nbuf:this%nMesh(2)+this%nbuf,this%nvar)::q
double precision,dimension(1-this%nbuf:this%nMesh(1)+this%nbuf)::xc
double precision,dimension(1-this%nbuf:this%nMesh(2)+this%nbuf)::yc
integer::nx,ny,nvar,nbuf,gridID,eosType
integer::i,j

xc=this%xc(1)%coords
yc=this%xc(2)%coords
nx=this%nMesh(1)
ny=this%nMesh(2)

 do j=1,ny
   do i=1,nx
      if(yc(j) .lt. 0.1d0) then
         q(i,j,4)=2500.d0
      elseif(yc(j) .gt. 0.1d0 .and. yc(j) .lt. 0.9d0) then
         q(i,j,4)=0.025d0
      else
         q(i,j,4)=250.d0
      endif
         q(i,j,1)=1.d0
         q(i,j,2)=0.d0
         q(i,j,3)=0.d0
   enddo
 enddo

end subroutine initWCShockTube2D

subroutine bdryWCShockTube2D(this,q)
class(grid)::this
double precision,dimension(1-this%nbuf:this%nMesh(1)+this%nbuf,1-this%nbuf:this%nMesh(2)+this%nbuf,this%nvar)::q
integer::nx,ny,nbuf,i,j

nbuf=this%nbuf
nx=this%nMesh(1)
ny=this%nMesh(2)

if(this%boundaryType .eq. 2) then !! reflective boundary condition

if(this%left_mpi .lt. 0) then
  do j=1,ny
    do i=1,nbuf
       q(i-nbuf,j,:)=q(nbuf-i+1,j,:)
    enddo
  enddo

  do j=1,ny
    do i=1,nbuf
       q(i-nbuf,j,2)=-q(nbuf-i+1,j,2)
    enddo
  enddo
endif

if(this%right_mpi .lt. 0) then
  do j=1,ny
    do i=1,nbuf
       q(nx+i,j,:)=q(nx-i+1,j,:)
    enddo
  enddo

  do j=1,ny
    do i=1,nbuf
       q(nx+i,j,2)=-q(nx-i+1,j,2)
    enddo
  enddo
endif

if(this%up_mpi .lt. 0) then
  do j=1,nbuf
    do i=1,nx
       q(i,j-nbuf,:)=q(i,nbuf-j+1,:)
    enddo
  enddo

  do j=1,nbuf
    do i=1,nx
       q(i,j-nbuf,3)=-q(i,nbuf-j+1,3)
    enddo
  enddo
endif

if(this%down_mpi .lt. 0) then
  do j=1,nbuf
    do i=1,nx
       q(i,ny+j,:)=q(i,ny-j+1,:)
    enddo
  enddo

  do j=1,nbuf
    do i=1,nx
       q(i,ny+j,3)=-q(i,ny-j+1,3)
    enddo
  enddo
endif

endif

end subroutine bdryWCShockTube2D



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!! Sod shocktube problem 2D (Sod, 1978, JCP, 27, 1)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine SodShockTube2D(gridID)
integer::gridID
type(grid)::g1
integer::ndim,nbuf,coordType,variable(8)
integer::nMesh(2),dims(2)
double precision::leftBdry(2),rightBdry(2)
logical::periods(2),reorder
double precision::t0,t1
integer::nstep,ierr

nstep=0
variable=0
ndim=2
nbuf=2
coordType=1
variable(1)=1
variable(2)=1
variable(3)=1
variable(8)=1

nMesh(1)=20
nMesh(2)=200
leftBdry(1)=-0.05d0
leftBdry(2)=-0.5d0
rightBdry(1)=0.05d0
rightBdry(2)=0.5d0
dims=(/0,0/)
call MPI_DIMS_CREATE(nprocs,ndim,dims,ierr)
dims(1)=nprocs
dims(2)=1
periods(1)=.false.
periods(2)=.false.
reorder=.true.

call g1%setTopologyMPI(ndim,dims,periods,reorder)
call g1%setMesh(nMesh,leftBdry,rightBdry,nbuf,coordType,gridID)
call g1%setVariable(variable) !! den,vx
call g1%setMPIWindows()
call g1%setEoS(eosType=2) !! isothermal
call g1%setadiGamma(gam=1.4d0) !! isothermal sound speed
call g1%setCFL(CFL=0.6d0) !! Courant number
call g1%setTime(fstart=0,tend=0.2d0,dtout=0.2d0)
call g1%setSlopeLimiter(limiterType=1)
call g1%setSolverType(solverType=2)
call g1%setBoundaryType(boundaryType=1)
call g1%initVariable()
call g1%exchangeBdryMPI(g1%q,g1%winq)
call g1%setBoundary(g1%q)
call g1%writeGrid()
g1%writeFlag=.false.

t0=MPI_WTIME()
do while (g1%t .lt. g1%tend)
   call g1%griddt()
   call g1%evolveGridRK2()
   if(myid .eq. 0 ) then
     print *,"testSuite.f03: gridID=",g1%gridID,'t=',g1%t,"dt=",g1%dt
   endif
   if(g1%writeFlag .eqv. .true.) then
      t1=MPI_WTIME()
      if(myid .eq. 0) then
        print *,"testSuiteMPI.f03: elapse time=",t1-t0,"secs"
      endif
      call g1%writeGrid()
      g1%writeFlag=.false.
      t0=MPI_WTIME()
   endif
enddo !! end do while

end subroutine SodShockTube2D

subroutine initSodShockTube2D(this,q)
class(grid)::this
double precision,dimension(1-this%nbuf:this%nMesh(1)+this%nbuf,1-this%nbuf:this%nMesh(2)+this%nbuf,this%nvar)::q
double precision,dimension(1-this%nbuf:this%nMesh(1)+this%nbuf)::xc
double precision,dimension(1-this%nbuf:this%nMesh(2)+this%nbuf)::yc
integer::nx,ny,nvar,nbuf,gridID,eosType
integer::i,j

xc=this%xc(1)%coords
yc=this%xc(2)%coords
nx=this%nMesh(1)
ny=this%nMesh(2)

 do j=1,ny
   do i=1,nx
      if(yc(j) .lt. 0.d0) then
         q(i,j,1)=1.d0
         q(i,j,4)=2.5d0
      else
         q(i,j,1)=0.125d0
         q(i,j,4)=0.25d0
      endif
         q(i,j,2)=0.d0
         q(i,j,3)=0.d0
   enddo
 enddo

end subroutine initSodShockTube2D

subroutine bdrySodShockTube2D(this,q)
class(grid)::this
double precision,dimension(1-this%nbuf:this%nMesh(1)+this%nbuf,1-this%nbuf:this%nMesh(2)+this%nbuf,this%nvar)::q
integer::nx,ny,nbuf,i,j

nbuf=this%nbuf
nx=this%nMesh(1)
ny=this%nMesh(2)

if(this%boundaryType .eq. 1) then !! zero gradient
!!!!!! left boundary
  if(this%left_mpi .lt. 0) then
    do i=1,nbuf
       q(1-i,:,:)=q(1,:,:)
    enddo
  endif
!!!!!! right boundary
  if(this%right_mpi .lt. 0) then
    do i=1,nbuf
       q(nx+i,:,:)=q(nx,:,:)
    enddo
  endif
!!!!!! up boundary
  if(this%up_mpi .lt. 0) then
    do j=1,nbuf
       q(:,1-j,:)=q(:,1,:)
    enddo
  endif
!!!!!! down boundary
  if(this%down_mpi .lt. 0) then
    do j=1,nbuf
       q(:,ny+j,:)=q(:,ny,:)
    enddo
  endif
endif
end subroutine bdrySodShockTube2D


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!! isothermal Shocktube problem 2D
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine isoShockTube2D(gridID)
integer::gridID
type(grid)::g1
integer::ndim,nbuf,coordType,variable(8)
integer::nMesh(2),dims(2)
double precision::leftBdry(2),rightBdry(2)
logical::periods(2),reorder
double precision::t0,t1
integer::nstep,ierr

nstep=0
variable=0
ndim=2
nbuf=2
coordType=1
variable(1)=1
variable(2)=1
variable(3)=1

nMesh(1)=5
nMesh(2)=500
leftBdry(1)=-0.005d0
leftBdry(2)=-0.5d0
rightBdry(1)=0.005d0
rightBdry(2)=0.5d0
dims=(/0,0/)
call MPI_DIMS_CREATE(nprocs,ndim,dims,ierr)
periods(1)=.false.
periods(2)=.false.
reorder=.true.

call g1%setTopologyMPI(ndim,dims,periods,reorder)
call g1%setMesh(nMesh,leftBdry,rightBdry,nbuf,coordType,gridID)
call g1%setVariable(variable) !! den,vx
call g1%setMPIWindows()
call g1%setEoS(eosType=1) !! isothermal
call g1%setSoundSpeed(snd=1.d0) !! isothermal sound speed
call g1%setCFL(CFL=0.6d0) !! Courant number
call g1%setTime(fstart=0,tend=0.2d0,dtout=0.2d0)
call g1%setSlopeLimiter(limiterType=1)
call g1%setSolverType(solverType=2)
call g1%setBoundaryType(boundaryType=1)
call g1%initVariable()
call g1%exchangeBdryMPI(g1%q,g1%winq)
call g1%setBoundary(g1%q)
call g1%writeGrid()
g1%writeFlag=.false.

t0=MPI_WTIME()
do while (g1%t .lt. g1%tend)
   call g1%griddt()
   call g1%evolveGridRK2()
   if(myid .eq. 0) then
     print *,"testSuite.f03: gridID=",g1%gridID,'t=',g1%t,"dt=",g1%dt
   endif
   if(g1%writeFlag .eqv. .true.) then
      t1=MPI_WTIME()
      if(myid .eq. 0) then
        print *,"testSuiteMPI.f03: elapse time=",t1-t0,"secs"
      endif
      call g1%writeGrid()
      g1%writeFlag=.false.
      t0=MPI_WTIME()
   endif
enddo !! end do while

end subroutine isoShockTube2D

subroutine initIsoShockTube2D(this,q)
class(grid)::this
double precision,dimension(1-this%nbuf:this%nMesh(1)+this%nbuf,1-this%nbuf:this%nMesh(2)+this%nbuf,this%nvar)::q
double precision,dimension(1-this%nbuf:this%nMesh(1)+this%nbuf)::xc
double precision,dimension(1-this%nbuf:this%nMesh(2)+this%nbuf)::yc
integer::nx,ny,nvar,nbuf,gridID,eosType
integer::i,j

xc=this%xc(1)%coords
yc=this%xc(2)%coords
nx=this%nMesh(1)
ny=this%nMesh(2)

 do j=1,ny
   do i=1,nx
      if(yc(j) .lt. 0.d0) then
         q(i,j,1)=1.d0
      else
         q(i,j,1)=0.25d0
      endif
         q(i,j,2)=0.d0
         q(i,j,3)=0.d0
   enddo
 enddo
end subroutine initIsoShockTube2D

subroutine bdryIsoShockTube2D(this,q)
class(grid)::this
double precision,dimension(1-this%nbuf:this%nMesh(1)+this%nbuf,1-this%nbuf:this%nMesh(2)+this%nbuf,this%nvar)::q
integer::nx,ny,nbuf,i,j

nbuf=this%nbuf
nx=this%nMesh(1)
ny=this%nMesh(2)

if(this%boundaryType .eq. 1) then !! zero gradient
!!!!!! left boundary
  if(this%left_mpi .lt. 0) then
    do i=1,nbuf
       q(1-i,:,:)=q(1,:,:)
    enddo
  endif
!!!!!! right boundary
  if(this%right_mpi .lt. 0) then
    do i=1,nbuf
       q(nx+i,:,:)=q(nx,:,:)
    enddo
  endif
!!!!!! up boundary
  if(this%up_mpi .lt. 0) then
    do j=1,nbuf
       q(:,1-j,:)=q(:,1,:)
    enddo
  endif
!!!!!! down boundary
  if(this%down_mpi .lt. 0) then
    do j=1,nbuf
       q(:,ny+j,:)=q(:,ny,:)
    enddo
  endif
endif

end subroutine bdryIsoShockTube2D

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!! polytropic Shocktube problem 2D
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine polyShockTubeHD2D(gridID)
integer::gridID
type(grid)::g1
integer::ndim,nbuf,coordType,variable(8)
integer::nMesh(2),dims(2)
double precision::leftBdry(2),rightBdry(2)
logical::periods(2),reorder
double precision::t0,t1
integer::nstep,ierr

nstep=0
variable=0
ndim=2
nbuf=2
coordType=1
variable(1)=1
variable(2)=1
variable(3)=1

nMesh(1)=500
nMesh(2)=5
leftBdry(1)=-2.d0
leftBdry(2)=-0.02d0
rightBdry(1)=2.d0
rightBdry(2)=0.02d0
dims=(/0,0/)
call MPI_DIMS_CREATE(nprocs,ndim,dims,ierr)
periods(1)=.false.
periods(2)=.false.
reorder=.true.
g1%polyK=0.5d0
g1%polyGamma=2.d0

call g1%setTopologyMPI(ndim,dims,periods,reorder)
call g1%setMesh(nMesh,leftBdry,rightBdry,nbuf,coordType,gridID)
call g1%setVariable(variable) !! den,vx
call g1%setMPIWindows()
call g1%setEoS(eosType=3) !! polytropic gas
!call g1%setSoundSpeed(snd=1.d0) !! isothermal sound speed
call g1%setCFL(CFL=0.6d0) !! Courant number
call g1%setTime(fstart=0,tend=0.5d0,dtout=0.5d0)
call g1%setSlopeLimiter(limiterType=1)
call g1%setSolverType(solverType=6)
call g1%setBoundaryType(boundaryType=1)
call g1%initVariable()
call g1%exchangeBdryMPI(g1%q,g1%winq)
call g1%setBoundary(g1%q)
call g1%writeGrid()
g1%writeFlag=.false.

t0=MPI_WTIME()
do while (g1%t .lt. g1%tend)
   call g1%griddt()
   call g1%evolveGridRK2()
   if(myid .eq. 0) then
     print *,"testSuite.f03: gridID=",g1%gridID,'t=',g1%t,"dt=",g1%dt
   endif
   if(g1%writeFlag .eqv. .true.) then
      t1=MPI_WTIME()
      if(myid .eq. 0) then
        print *,"testSuiteMPI.f03: elapse time=",t1-t0,"secs"
      endif
      call g1%writeGrid()
      g1%writeFlag=.false.
      t0=MPI_WTIME()
   endif
enddo !! end do while

end subroutine polyShockTubeHD2D

subroutine initPolyShockTubeHD2D(this,q)
class(grid)::this
double precision,dimension(1-this%nbuf:this%nMesh(1)+this%nbuf,1-this%nbuf:this%nMesh(2)+this%nbuf,this%nvar)::q
double precision,dimension(1-this%nbuf:this%nMesh(1)+this%nbuf)::xc
double precision,dimension(1-this%nbuf:this%nMesh(2)+this%nbuf)::yc
integer::nx,ny,nvar,nbuf,gridID,eosType
integer::i,j

xc=this%xc(1)%coords
yc=this%xc(2)%coords
nx=this%nMesh(1)
ny=this%nMesh(2)

 do j=1,ny
   do i=1,nx
      if(xc(i) .lt. 0.d0) then
         q(i,j,1)=3.d0
      else
         q(i,j,1)=1.d0
      endif
         q(i,j,2)=0.d0
         q(i,j,3)=0.d0
   enddo
 enddo
end subroutine initPolyShockTubeHD2D

subroutine bdryPolyShockTubeHD2D(this,q)
class(grid)::this
double precision,dimension(1-this%nbuf:this%nMesh(1)+this%nbuf,1-this%nbuf:this%nMesh(2)+this%nbuf,this%nvar)::q
integer::nx,ny,nbuf,i,j

nbuf=this%nbuf
nx=this%nMesh(1)
ny=this%nMesh(2)

if(this%boundaryType .eq. 1) then !! zero gradient
!!!!!! left boundary
  if(this%left_mpi .lt. 0) then
    do i=1,nbuf
       q(1-i,:,:)=q(1,:,:)
    enddo
  endif
!!!!!! right boundary
  if(this%right_mpi .lt. 0) then
    do i=1,nbuf
       q(nx+i,:,:)=q(nx,:,:)
    enddo
  endif
!!!!!! up boundary
  if(this%up_mpi .lt. 0) then
    do j=1,nbuf
       q(:,1-j,:)=q(:,1,:)
    enddo
  endif
!!!!!! down boundary
  if(this%down_mpi .lt. 0) then
    do j=1,nbuf
       q(:,ny+j,:)=q(:,ny,:)
    enddo
  endif
endif

end subroutine bdryPolyShockTubeHD2D

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!! Brio-Wu Shocktube 2D (1988, JCP, 75, 400)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine BrioWuShockTube2D(gridID)
integer::gridID
type(grid)::g1
integer::ndim,nbuf,coordType,variable(8)
integer::nMesh(2),dims(2)
double precision::leftBdry(2),rightBdry(2)
logical::periods(2),reorder
double precision::t0,t1
integer::nstep,ierr

nstep=0
variable=0
ndim=2
nbuf=2
coordType=1 !! Cartesian
variable(1)=1  !! rho
variable(2)=1  !! vx
variable(3)=1  !! vy
variable(4)=1  !! vz
variable(5)=1  !! bx
variable(6)=1  !! by
variable(7)=1  !! bz
variable(8)=1  !! energy

nMesh(1)=4
nMesh(2)=400
leftBdry(1)=0.d0
rightBdry(1)=0.01d0
leftBdry(2)=0.d0
rightBdry(2)=1.d0
dims=(/0,0/)
call MPI_DIMS_CREATE(nprocs,ndim,dims,ierr)
dims(1)=1
dims(2)=nprocs
periods(1)=.false.
periods(2)=.false.
reorder=.true.

call g1%setTopologyMPI(ndim,dims,periods,reorder)
call g1%setMesh(nMesh,leftBdry,rightBdry,nbuf,coordType,gridID)
call g1%setVariable(variable) !! den,vx,vy,vz,bx,by,bz,ene
call g1%setMPIWindows()
call g1%setEoS(eosType=2) !! adiabatic
call g1%setadiGamma(gam=2.d0) !! ratio of heat capacity
call g1%setCFL(CFL=0.8d0) !! courant number
call g1%setTime(fstart=0,tend=0.1d0,dtout=0.1d0)
call g1%setSlopeLimiter(limiterType=3)
call g1%setSolverType(solverType=4) !! 4=AdiHLLMHD, 5=AdiHLLDMHD
call g1%setBoundaryType(boundaryType=1)
call g1%initVariable()
call g1%exchangeBdryMPI(g1%q,g1%winq)
call g1%setBoundary(g1%q)
call g1%writeGrid()
g1%writeFlag=.false.

t0=MPI_WTIME()
do while (g1%t .lt. g1%tend)
   call g1%griddt()
   call g1%evolveGridRK2()
   if(myid .eq. 0 ) then
     print *,"testSuite.f03: gridID=",g1%gridID,'t=',g1%t,"dt=",g1%dt
   endif
   if(g1%writeFlag .eqv. .true.) then
      t1=MPI_WTIME()
      if(myid .eq. 0) then
        print *,"testSuiteMPI.f03: elapse time=",t1-t0,"secs"
      endif
      call g1%writeGrid()
      g1%writeFlag=.false.
      t0=MPI_WTIME()
   endif
enddo !! end do while
end subroutine BrioWuShockTube2D

subroutine initBrioWuShockTube2D(this,q)
class(grid)::this
double precision,dimension(1-this%nbuf:this%nMesh(1)+this%nbuf,1-this%nbuf:this%nMesh(2)+this%nbuf,this%nvar)::q
double precision,dimension(1-this%nbuf:this%nMesh(1)+this%nbuf)::xc
double precision,dimension(1-this%nbuf:this%nMesh(2)+this%nbuf)::yc
integer::nx,ny,nvar,nbuf,gridID,eosType
integer::i,j
double precision::gam,rho,vx,vy,vz,bxl,byl,bzl,bxr,byr,bzr,p
double precision::bxc,byc,bzc

xc=this%xc(1)%coords
yc=this%xc(2)%coords
nx=this%nMesh(1)
ny=this%nMesh(2)
gam=this%adiGamma

do j=1,ny
  do i=1,nx
    if(yc(j) .lt. 0.5d0) then
      rho=1.d0
       vx=0.d0
       vy=0.d0
       vz=0.d0
       bxl=1.d0!0.75d0
       bxr=1.d0!0.75d0
       byl=0.75d0!1.d0
       byr=0.75d0!1.d0
       bzl=0.d0
       bzr=0.d0
         p=1.d0
    else
       rho=0.125d0
        vx=0.d0
        vy=0.d0
        vz=0.d0
        bxl=-1.d0!0.75d0
        bxr=-1.d0!0.75d0
        byl=0.75d0!-1.d0
        byr=0.75d0!-1.d0
        bzl=0.d0
        bzr=0.d0
          p=0.1d0
    endif
    bxc=0.5d0*(bxl+bxr)
    byc=0.5d0*(byl+byr)
    bzc=0.5d0*(bzl+bzr)

    q(i,j,1)=rho
    q(i,j,2)=rho*vx
    q(i,j,3)=rho*vy
    q(i,j,4)=rho*vz
    q(i,j,5)=bxl
    q(i,j,6)=byl
    q(i,j,7)=bzl
    q(i,j,8)=p/(gam-1.d0)+0.5d0*rho*(vx**2.d0+vy**2.d0+vz**2.d0)+0.5d0*(bxc**2.d0+byc**2.d0+bzc**2.d0)
    q(i,j,9)=bxr
    q(i,j,10)=byr
    q(i,j,11)=bzr
  enddo
enddo

end subroutine initBrioWuShockTube2D

subroutine bdryBrioWuShockTube2D(this,q)
class(grid)::this
double precision,dimension(1-this%nbuf:this%nMesh(1)+this%nbuf,1-this%nbuf:this%nMesh(2)+this%nbuf,this%nvar)::q
integer::nx,ny,nbuf,i,j

nbuf=this%nbuf
nx=this%nMesh(1)
ny=this%nMesh(2)

if(this%boundaryType .eq. 1) then !! zero gradient
!!!!!! left boundary
  if(this%left_mpi .lt. 0) then
    do i=1,nbuf
       q(1-i,:,:)=q(1,:,:)
    enddo
  endif
!!!!!! right boundary
  if(this%right_mpi .lt. 0) then
    do i=1,nbuf
       q(nx+i,:,:)=q(nx,:,:)
    enddo
  endif
!!!!!! up boundary
  if(this%up_mpi .lt. 0) then
    do j=1,nbuf
       q(:,1-j,:)=q(:,1,:)
    enddo
  endif
!!!!!! down boundary
  if(this%down_mpi .lt. 0) then
    do j=1,nbuf
       q(:,ny+j,:)=q(:,ny,:)
    enddo
  endif
endif

end subroutine bdryBrioWuShockTube2D

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!! Shu Osher Shocktube test (1989, JCP, 83,32)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!! nMesh(1)= 200 cannot capture small k oscillation
!!!!! well. A third-order limiter is requited for
!!!!! resolving short wavelength at low resolution.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine ShuOsherShocktubeTest(gridID)
integer::gridID
type(grid)::g1
integer::ndim,nbuf,coordType,variable(8)
integer::nMesh(2),dims(2)
double precision::leftBdry(2),rightBdry(2)
logical::periods(2),reorder
double precision::t0,t1
integer::nstep,ierr

nstep=0
variable=0
ndim=2
nbuf=2
coordType=1 !! Cartesian
variable(1)=1
variable(2)=1
variable(3)=1
variable(8)=1

nMesh(1)=400
nMesh(2)=4
leftBdry(1)=-1.d0
leftBdry(2)=-4.d0/dble(nMesh(1))
rightBdry(1)=1.d0
rightBdry(2)=4.d0/dble(nMesh(1))
dims=(/0,0/)
call MPI_DIMS_CREATE(nprocs,ndim,dims,ierr)
periods(1)=.false.
periods(2)=.false.
reorder=.true.

call g1%setTopologyMPI(ndim,dims,periods,reorder)
call g1%setMesh(nMesh,leftBdry,rightBdry,nbuf,coordType,gridID)
call g1%setVariable(variable) !! den,vx,vy,ene
call g1%setMPIWindows()
call g1%setEoS(eosType=2) !! adiabatic
call g1%setadiGamma(gam=1.4d0) !! ratio of heat capacity Cp/Cv
call g1%setCFL(CFL=0.8d0) !! Courant number
call g1%setTime(fstart=0,tend=0.47d0,dtout=0.47d0)
call g1%setSlopeLimiter(limiterType=2)
call g1%setSolverType(solverType=3)
call g1%setBoundaryType(boundaryType=1)
call g1%initVariable()
call g1%exchangeBdryMPI(g1%q,g1%winq)
call g1%setBoundary(g1%q)
call g1%writeGrid()
g1%writeFlag=.false.

t0=MPI_WTIME()
do while (g1%t .lt. g1%tend)
   call g1%griddt()
   call g1%evolveGridRK2()
   if(myid .eq. 0) then
     print *,"testSuite.f03: gridID=",g1%gridID,'t=',g1%t,"dt=",g1%dt
   endif
   if(g1%writeFlag .eqv. .true.) then
      t1=MPI_WTIME()
      if(myid .eq. 0) then
        print *,"testSuiteMPI.f03: elapse time=",t1-t0,"secs"
      endif
      call g1%writeGrid()
      g1%writeFlag=.false.
      t0=MPI_WTIME()
   endif
enddo !! end do while
end subroutine ShuOsherShocktubeTest

subroutine initShuOsherShocktubeTest(this,q)
class(grid)::this
double precision,dimension(1-this%nbuf:this%nMesh(1)+this%nbuf,1-this%nbuf:this%nMesh(2)+this%nbuf,this%nvar)::q
double precision,dimension(1-this%nbuf:this%nMesh(1)+this%nbuf)::xc
double precision,dimension(1-this%nbuf:this%nMesh(2)+this%nbuf)::yc
integer::nx,ny,nvar,nbuf,gridID,eosType
integer::i,j
double precision::gam,vx,rho,pi,p

xc=this%xc(1)%coords
yc=this%xc(2)%coords
nx=this%nMesh(1)
ny=this%nMesh(2)
gam=this%adiGamma
pi = datan(1.d0)*4.d0

 do j=1,ny
   do i=1,nx
      if(xc(i) .lt. -0.8d0) then
        rho = 3.857143d0
         vx = 2.629369d0
          p = 10.33333d0
        q(i,j,1)=rho
        q(i,j,2)=rho*vx
        q(i,j,3)=0.d0
        q(i,j,4)=p/(gam-1.d0)+0.5*rho*vx**2
      else
        rho = 1.d0+0.2d0*dsin(5.d0*pi*xc(i))
         vx = 0.d0
          p = 1.d0
        q(i,j,1)=rho
        q(i,j,2)=rho*vx
        q(i,j,3)=0.d0
        q(i,j,4)=p/(gam-1.d0)+0.5*rho*vx**2
      endif
   enddo
 enddo

end subroutine initShuOsherShocktubeTest

subroutine bdryShuOsherShocktubeTest(this,q)
class(grid)::this
double precision,dimension(1-this%nbuf:this%nMesh(1)+this%nbuf,1-this%nbuf:this%nMesh(2)+this%nbuf,this%nvar)::q
integer::nx,ny,nbuf,i,j

nbuf=this%nbuf
nx=this%nMesh(1)
ny=this%nMesh(2)

if(this%boundaryType .eq. 1) then !! zero gradient
!!!!!! left boundary
  if(this%left_mpi .lt. 0) then
    do i=1,nbuf
       q(1-i,:,:)=q(1,:,:)
    enddo
  endif
!!!!!! right boundary
  if(this%right_mpi .lt. 0) then
    do i=1,nbuf
       q(nx+i,:,:)=q(nx,:,:)
    enddo
  endif
!!!!!! up boundary
  if(this%up_mpi .lt. 0) then
    do j=1,nbuf
       q(:,1-j,:)=q(:,1,:)
    enddo
  endif
!!!!!! down boundary
  if(this%down_mpi .lt. 0) then
    do j=1,nbuf
       q(:,ny+j,:)=q(:,ny,:)
    enddo
  endif
endif

end subroutine bdryShuOsherShocktubeTest


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!! Strong Rarefaction tests (Einfeldt et al. 1991, JCP, 92, 273)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine strongRarefactionTest(gridID)
integer::gridID
type(grid)::g1
integer::ndim,nbuf,coordType,variable(8)
integer::nMesh(2),dims(2)
double precision::leftBdry(2),rightBdry(2)
logical::periods(2),reorder
double precision::t0,t1
integer::nstep,ierr

nstep=0
variable=0
ndim=2
nbuf=2
coordType=1 !! Cartesian
variable(1)=1
variable(2)=1
variable(3)=1
variable(8)=1

nMesh(1)=200
nMesh(2)=8
leftBdry(1)=-0.45d0
leftBdry(2)=-8.d0/200.d0*0.45d0
rightBdry(1)=0.45d0
rightBdry(2)=8.d0/200.d0*0.45d0
dims=(/0,0/)
call MPI_DIMS_CREATE(nprocs,ndim,dims,ierr)
periods(1)=.false.
periods(2)=.false.
reorder=.true.

call g1%setTopologyMPI(ndim,dims,periods,reorder)
call g1%setMesh(nMesh,leftBdry,rightBdry,nbuf,coordType,gridID)
call g1%setVariable(variable) !! den,vx,vy,ene
call g1%setMPIWindows()
call g1%setEoS(eosType=2) !! adiabatic
call g1%setadiGamma(gam=1.4d0) !! ratio of heat capacity Cp/Cv
call g1%setCFL(CFL=0.8d0) !! Courant number
call g1%setTime(fstart=0,tend=0.1d0,dtout=0.1d0)
call g1%setSlopeLimiter(limiterType=3)
call g1%setSolverType(solverType=3)
call g1%setBoundaryType(boundaryType=1)
call g1%initVariable()
call g1%exchangeBdryMPI(g1%q,g1%winq)
call g1%setBoundary(g1%q)
call g1%writeGrid()

g1%writeFlag=.false.

t0=MPI_WTIME()
do while (g1%t .lt. g1%tend)
   call g1%griddt()
   call g1%evolveGridRK2()
   if(myid .eq. 0) then
     print *,"testSuite.f03: gridID=",g1%gridID,'t=',g1%t,"dt=",g1%dt
   endif
   if(g1%writeFlag .eqv. .true.) then
      t1=MPI_WTIME()
      if(myid .eq. 0) then
        print *,"testSuiteMPI.f03: elapse time=",t1-t0,"secs"
      endif
      call g1%writeGrid()
      g1%writeFlag=.false.
      t0=MPI_WTIME()
   endif
enddo !! end do while

end subroutine strongRareFactionTest

subroutine initStrongRareFactionTest(this,q)
class(grid)::this
double precision,dimension(1-this%nbuf:this%nMesh(1)+this%nbuf,1-this%nbuf:this%nMesh(2)+this%nbuf,this%nvar)::q
double precision,dimension(1-this%nbuf:this%nMesh(1)+this%nbuf)::xc
double precision,dimension(1-this%nbuf:this%nMesh(2)+this%nbuf)::yc
integer::nx,ny,nvar,nbuf,gridID,eosType
integer::i,j

xc=this%xc(1)%coords
yc=this%xc(2)%coords
nx=this%nMesh(1)
ny=this%nMesh(2)

 do j=1,ny
   do i=1,nx
      if(xc(i) .lt. 0.d0) then
        q(i,j,2)=-2.d0
      else
        q(i,j,2)=2.d0
      endif
        q(i,j,1)=1.d0
        q(i,j,3)=0.d0
        q(i,j,4)=3.d0
   enddo
 enddo

end subroutine initStrongRarefactionTest

subroutine bdryStrongRarefactionTest(this,q)
class(grid)::this
double precision,dimension(1-this%nbuf:this%nMesh(1)+this%nbuf,1-this%nbuf:this%nMesh(2)+this%nbuf,this%nvar)::q
integer::nx,ny,nbuf,i,j

nbuf=this%nbuf
nx=this%nMesh(1)
ny=this%nMesh(2)

if(this%boundaryType .eq. 1) then !! zero gradient
!!!!!! left boundary
  if(this%left_mpi .lt. 0) then
    do i=1,nbuf
       q(1-i,:,:)=q(1,:,:)
    enddo
  endif
!!!!!! right boundary
  if(this%right_mpi .lt. 0) then
    do i=1,nbuf
       q(nx+i,:,:)=q(nx,:,:)
    enddo
  endif
!!!!!! up boundary
  if(this%up_mpi .lt. 0) then
    do j=1,nbuf
       q(:,1-j,:)=q(:,1,:)
    enddo
  endif
!!!!!! down boundary
  if(this%down_mpi .lt. 0) then
    do j=1,nbuf
       q(:,ny+j,:)=q(:,ny,:)
    enddo
  endif
endif

end subroutine bdryStrongRarefactionTest
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!! Linear Fast Wave Test MHD 2D
!!!!! Gardiner & Stone, JCP, 2005, 205, 509
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine linearFastWaveTestMHD2D(gridID)
integer::gridID
integer::N
type(grid)::g1
integer::ndim,nbuf,coordType,variable(8)
integer::nMesh(2),dims(2)
double precision::leftBdry(2),rightBdry(2)
double precision::pi,snd,wavelen
logical::periods(2),reorder
double precision::t0,t1
integer::nstep,ierr

nstep=0
pi=4.d0*atan2(1.d0,1.d0)
variable=0
ndim=2
nbuf=2
wavelen=1.d0
coordType=1 !! Cartesian
variable(1)=1  !! rho
variable(2)=1  !! vx
variable(3)=1  !! vy
variable(4)=1  !! vz
variable(5)=1  !! bx
variable(6)=1  !! by
variable(7)=1  !! bz
variable(8)=1  !! energy

snd=2.d0 ! 0.5=slow mode, 1=alfven wave, 2=fast mode
N=256
nMesh(1)=2*N
nMesh(2)=N
leftBdry(1)=0.d0
rightBdry(1)=dsqrt(5.d0)
leftBdry(2)=0.d0
rightBdry(2)=dsqrt(5.d0)/2.d0
dims=(/0,0/)
call MPI_DIMS_CREATE(nprocs,ndim,dims,ierr)
periods(1)=.true.
periods(2)=.true.
reorder=.true.

call g1%setTopologyMPI(ndim,dims,periods,reorder)
call g1%setGridID(gridID=gridID)
call g1%setTime(fstart=0,tend=wavelen/snd,dtout=wavelen/snd/10.d0)
if (g1%fstart .eq. 0) then
   call g1%setMesh(nMesh,leftBdry,rightBdry,nbuf,coordType,gridID)
   call g1%setVariable(variable) !! den,vx,vy,vz,bx,by,bz,ene
   call g1%setMPIWindows()
   call g1%setEoS(eosType=2) !! adiabatic
   call g1%setadiGamma(gam=5.d0/3.d0) !! ratio of heat capacity
   call g1%setCFL(CFL=0.4d0) !! courant number
   call g1%setSlopeLimiter(limiterType=2)
   call g1%setSolverType(solverType=5) !! 4=AdiHLLMHD, 5=AdiHLLDMHD
   call g1%setBoundaryType(boundaryType=3) !! 1=zero gradient, 3=periodic
   call g1%initVariable()
   call g1%exchangeBdryMPI(g1%q,g1%winq)
   call g1%setBoundary(g1%q)
   call g1%writeGrid()
endif
g1%writeFlag=.false.

t0=MPI_WTIME()
do while (g1%t .lt. g1%tend)
   call g1%griddt()
   call g1%evolveGridRK2()
   if(myid .eq. 0) then
     print *,"testSuite.f03: gridID=",g1%gridID,'t=',g1%t,"dt=",g1%dt
   endif
   if(g1%writeFlag .eqv. .true.) then
      t1=MPI_WTIME()
      if(myid .eq. 0) then
        print *,"testSuiteMPI.f03: elapse time=",t1-t0,"secs"
      endif
      call g1%writeGrid()
      g1%writeFlag=.false.
      t0=MPI_WTIME()
   endif
enddo !! end do while
end subroutine linearFastWaveTestMHD2D

subroutine initLinearFastWaveTestMHD2D(this,q)
class(grid)::this
double precision,dimension(1-this%nbuf:this%nMesh(1)+this%nbuf,1-this%nbuf:this%nMesh(2)+this%nbuf,this%nvar)::q
double precision,dimension(1-this%nbuf:this%nMesh(1)+this%nbuf,1-this%nbuf:this%nMesh(2)+this%nbuf)::Az
double precision,dimension(1-this%nbuf:this%nMesh(1)+this%nbuf)::xc,xl,xr
double precision,dimension(1-this%nbuf:this%nMesh(2)+this%nbuf)::yc,yl,yr
integer::nx,ny,nvar,nbuf,gridID,eosType
integer::i,j
double precision::gam,rho,vx,vy,vz,bxl,byl,bzl,bxr,byr,bzr,p
double precision::bxc,byc,bzc
double precision::pi,b0
double precision::rc
double precision::ang,x1,x2
double precision::b1,b2,b3,v1,v2,v3
double precision::x,y,dx,dy,bx,by,bz
double precision::eps

pi=4.d0*atan2(1.d0,1.d0)
b0=1.d0
xc=this%xc(1)%coords
yc=this%xc(2)%coords
xl=this%xl(1)%coords
yl=this%xl(2)%coords
xr=this%xr(1)%coords
yr=this%xr(2)%coords
nx=this%nMesh(1)
ny=this%nMesh(2)
gam=this%adiGamma
ang=datan2(2.d0,1.d0)
eps=1.d-6

do j=1,ny+1
  do i=1,nx+1
    x=xl(i)
    y=yl(j)
   x1= x*dcos(ang)+y*dsin(ang)
   x2=-x*dsin(ang)+y*dcos(ang)
   Az(i,j)=-eps*2.d0/(3.d0*pi)*dsqrt(2.d0/5.d0)*dsin(2*pi*x1)
  enddo
enddo

do j=1,ny
  do i=1,nx
    dx=this%dx(1)%coords(i)
    dy=this%dx(2)%coords(j)
    x1=xc(i)*dcos(ang)+yc(j)*dsin(ang)
    x2=-xc(i)*dsin(ang)+yc(j)*dcos(ang)

    rho=1.d0+eps*6.d0/(6.d0*dsqrt(5.d0))*dcos(2*pi*x1)
    v1=(0.d0+eps*12.d0/(6.d0*dsqrt(5.d0))*dcos(2*pi*x1))/rho !! 1.d0 for entropy mode test
    v2=(0.d0-eps*4.d0*dsqrt(2.d0)/(6.d0*dsqrt(5.d0))*dcos(2*pi*x1))/rho
    v3=(0.d0-eps*2.d0/(6.d0*dsqrt(5.d0))*dcos(2*pi*x1))/rho

    b1=1.d0
    b2=dsqrt(2.d0)
    b3=0.5d0

    vx=-v2*dsin(ang)+v1*dcos(ang)
    vy= v2*dcos(ang)+v1*dsin(ang)
    vz= v3

    bx=-b2*dsin(ang)+b1*dcos(ang)
    by= b2*dcos(ang)+b1*dsin(ang)
    bz= b3

   bxl=bx+(Az(i,j+1)-Az(i,j))/dy
   byl=by-(Az(i+1,j)-Az(i,j))/dx
   bzl=bz+4.d0/(6.d0*dsqrt(5.d0))*eps*dcos(2*pi*x1)
   bxr=bx+(Az(i+1,j+1)-Az(i+1,j))/dy
   byr=by-(Az(i+1,j+1)-Az(i,j+1))/dx
   bzr=bzl
     p=1.d0/gam

    bxc=0.5d0*(bxl+bxr)
    byc=0.5d0*(byl+byr)
    bzc=0.5d0*(bzl+bzr)

    q(i,j,1)=rho
    q(i,j,2)=rho*vx
    q(i,j,3)=rho*vy
    q(i,j,4)=rho*vz
    q(i,j,5)=bxl
    q(i,j,6)=byl
    q(i,j,7)=bzl
    q(i,j,8)=p/(gam-1.d0)+0.5d0*(b1**2.d0+b2**2.d0+b3**2.d0)+eps*27.d0/(6.d0*dsqrt(5.d0))*dcos(2*pi*x1)
    q(i,j,9)=bxr
    q(i,j,10)=byr
    q(i,j,11)=bzr
  enddo
enddo
end subroutine initLinearFastWaveTestMHD2D

subroutine bdryLinearFastWaveTestMHD2D(this,q)
class(grid)::this
double precision,dimension(1-this%nbuf:this%nMesh(1)+this%nbuf,1-this%nbuf:this%nMesh(2)+this%nbuf,this%nvar)::q
integer::nx,ny,nbuf,i,j

nbuf=this%nbuf
nx=this%nMesh(1)
ny=this%nMesh(2)

if(this%boundaryType .eq. 3) then !! periodic boundary condition
endif
end subroutine bdryLinearFastWaveTestMHD2D
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!! Linear Slow Wave Test MHD 2D
!!!!! Gardiner & Stone, JCP, 2005, 205, 509
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine linearSlowWaveTestMHD2D(gridID)
integer::gridID
integer::N
type(grid)::g1
integer::ndim,nbuf,coordType,variable(8)
integer::nMesh(2),dims(2)
double precision::leftBdry(2),rightBdry(2)
double precision::pi,snd,wavelen
logical::periods(2),reorder
double precision::t0,t1
integer::nstep,ierr

nstep=0
pi=4.d0*atan2(1.d0,1.d0)
variable=0
ndim=2
nbuf=2
wavelen=1.d0
coordType=1 !! Cartesian
variable(1)=1  !! rho
variable(2)=1  !! vx
variable(3)=1  !! vy
variable(4)=1  !! vz
variable(5)=1  !! bx
variable(6)=1  !! by
variable(7)=1  !! bz
variable(8)=1  !! energy

snd=0.5d0 ! 0.5=slow mode, 1=alfven wave, 2=fast mode
N=256
nMesh(1)=2*N
nMesh(2)=N
leftBdry(1)=0.d0
rightBdry(1)=dsqrt(5.d0)
leftBdry(2)=0.d0
rightBdry(2)=dsqrt(5.d0)/2.d0
dims=(/0,0/)
call MPI_DIMS_CREATE(nprocs,ndim,dims,ierr)
periods(1)=.true.
periods(2)=.true.
reorder=.true.

call g1%setTopologyMPI(ndim,dims,periods,reorder)
call g1%setGridID(gridID=gridID)
call g1%setTime(fstart=0,tend=wavelen/snd,dtout=wavelen/snd/10.d0)
if (g1%fstart .eq. 0) then
   call g1%setMesh(nMesh,leftBdry,rightBdry,nbuf,coordType,gridID)
   call g1%setVariable(variable) !! den,vx,vy,vz,bx,by,bz,ene
   call g1%setMPIWindows()
   call g1%setEoS(eosType=2) !! adiabatic
   call g1%setadiGamma(gam=5.d0/3.d0) !! ratio of heat capacity
   call g1%setCFL(CFL=0.4d0) !! courant number
   call g1%setSlopeLimiter(limiterType=2)
   call g1%setSolverType(solverType=5) !! 4=AdiHLLMHD, 5=AdiHLLDMHD
   call g1%setBoundaryType(boundaryType=3) !! 1=zero gradient, 3=periodic
   call g1%initVariable()
   call g1%exchangeBdryMPI(g1%q,g1%winq)
   call g1%setBoundary(g1%q)
   call g1%writeGrid()
endif
g1%writeFlag=.false.

t0=MPI_WTIME()
do while (g1%t .lt. g1%tend)
   call g1%griddt()
   call g1%evolveGridRK2()
   if(myid .eq. 0) then
     print *,"testSuite.f03: gridID=",g1%gridID,'t=',g1%t,"dt=",g1%dt
   endif
   if(g1%writeFlag .eqv. .true.) then
      t1=MPI_WTIME()
      if(myid .eq. 0) then
        print *,"testSuiteMPI.f03: elapse time=",t1-t0,"secs"
      endif
      call g1%writeGrid()
      g1%writeFlag=.false.
      t0=MPI_WTIME()
   endif
enddo !! end do while
end subroutine linearSlowWaveTestMHD2D

subroutine initLinearSlowWaveTestMHD2D(this,q)
class(grid)::this
double precision,dimension(1-this%nbuf:this%nMesh(1)+this%nbuf,1-this%nbuf:this%nMesh(2)+this%nbuf,this%nvar)::q
double precision,dimension(1-this%nbuf:this%nMesh(1)+this%nbuf,1-this%nbuf:this%nMesh(2)+this%nbuf)::Az
double precision,dimension(1-this%nbuf:this%nMesh(1)+this%nbuf)::xc,xl,xr
double precision,dimension(1-this%nbuf:this%nMesh(2)+this%nbuf)::yc,yl,yr
integer::nx,ny,nvar,nbuf,gridID,eosType
integer::i,j
double precision::gam,rho,vx,vy,vz,bxl,byl,bzl,bxr,byr,bzr,p
double precision::bxc,byc,bzc
double precision::pi,b0
double precision::rc
double precision::ang,x1,x2
double precision::b1,b2,b3,v1,v2,v3
double precision::x,y,dx,dy,bx,by,bz
double precision::eps

pi=4.d0*atan2(1.d0,1.d0)
b0=1.d0
xc=this%xc(1)%coords
yc=this%xc(2)%coords
xl=this%xl(1)%coords
yl=this%xl(2)%coords
xr=this%xr(1)%coords
yr=this%xr(2)%coords
nx=this%nMesh(1)
ny=this%nMesh(2)
gam=this%adiGamma
ang=datan2(2.d0,1.d0)
eps=1.d-6

do j=1,ny+1
  do i=1,nx+1
    x=xl(i)
    y=yl(j)
   x1= x*dcos(ang)+y*dsin(ang)
   x2=-x*dsin(ang)+y*dcos(ang)
   Az(i,j)=eps*1.d0/(3.d0*pi)*dsqrt(2.d0/5.d0)*dsin(2*pi*x1)
  enddo
enddo

do j=1,ny
  do i=1,nx
    dx=this%dx(1)%coords(i)
    dy=this%dx(2)%coords(j)
    x1=xc(i)*dcos(ang)+yc(j)*dsin(ang)
    x2=-xc(i)*dsin(ang)+yc(j)*dcos(ang)

    rho=1.d0+eps*12.d0/(6.d0*dsqrt(5.d0))*dcos(2*pi*x1)
    v1=(0.d0+eps*6.d0/(6.d0*dsqrt(5.d0))*dcos(2*pi*x1))/rho !! 1.d0 for entropy mode test
    v2=(0.d0+eps*8.d0*dsqrt(2.d0)/(6.d0*dsqrt(5.d0))*dcos(2*pi*x1))/rho
    v3=(0.d0+eps*4.d0/(6.d0*dsqrt(5.d0))*dcos(2*pi*x1))/rho

    b1=1.d0
    b2=dsqrt(2.d0)
    b3=0.5d0

    vx=-v2*dsin(ang)+v1*dcos(ang)
    vy= v2*dcos(ang)+v1*dsin(ang)
    vz= v3

    bx=-b2*dsin(ang)+b1*dcos(ang)
    by= b2*dcos(ang)+b1*dsin(ang)
    bz= b3

   bxl=bx+(Az(i,j+1)-Az(i,j))/dy
   byl=by-(Az(i+1,j)-Az(i,j))/dx
   bzl=bz-2.d0/(6.d0*dsqrt(5.d0))*eps*dcos(2*pi*x1)
   bxr=bx+(Az(i+1,j+1)-Az(i+1,j))/dy
   byr=by-(Az(i+1,j+1)-Az(i,j+1))/dx
   bzr=bzl
     p=1.d0/gam

    bxc=0.5d0*(bxl+bxr)
    byc=0.5d0*(byl+byr)
    bzc=0.5d0*(bzl+bzr)

    q(i,j,1)=rho
    q(i,j,2)=rho*vx
    q(i,j,3)=rho*vy
    q(i,j,4)=rho*vz
    q(i,j,5)=bxl
    q(i,j,6)=byl
    q(i,j,7)=bzl
    q(i,j,8)=p/(gam-1.d0)+0.5d0*(b1**2.d0+b2**2.d0+b3**2.d0)+eps*9.d0/(6.d0*dsqrt(5.d0))*dcos(2*pi*x1)
    q(i,j,9)=bxr
    q(i,j,10)=byr
    q(i,j,11)=bzr
  enddo
enddo
end subroutine initLinearSlowWaveTestMHD2D

subroutine bdryLinearSlowWaveTestMHD2D(this,q)
class(grid)::this
double precision,dimension(1-this%nbuf:this%nMesh(1)+this%nbuf,1-this%nbuf:this%nMesh(2)+this%nbuf,this%nvar)::q
integer::nx,ny,nbuf,i,j

nbuf=this%nbuf
nx=this%nMesh(1)
ny=this%nMesh(2)

if(this%boundaryType .eq. 3) then !! periodic boundary condition
endif

end subroutine bdryLinearSlowWaveTestMHD2D
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!! Circularly Polarized Alven Wave
!!!!! Toth, 2000, JCP, 161, 605
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine CPAlvenWave(gridID)
integer::gridID
integer::N
type(grid)::g1
integer::ndim,nbuf,coordType,variable(8)
integer::nMesh(2),dims(2)
double precision::leftBdry(2),rightBdry(2)
double precision::pi
logical::periods(2),reorder
double precision::t0,t1
integer::nstep,ierr

nstep=0
pi=4.d0*atan2(1.d0,1.d0)
variable=0
ndim=2
nbuf=2
coordType=1 !! Cartesian
variable(1)=1  !! rho
variable(2)=1  !! vx
variable(3)=1  !! vy
variable(4)=1  !! vz
variable(5)=1  !! bx
variable(6)=1  !! by
variable(7)=1  !! bz
variable(8)=1  !! energy

N=32
nMesh(1)=2*N
nMesh(2)=N
leftBdry(1)=0.d0
rightBdry(1)=dsqrt(5.d0)
leftBdry(2)=0.d0
rightBdry(2)=dsqrt(5.d0)/2.d0
dims=(/0,0/)
call MPI_DIMS_CREATE(nprocs,ndim,dims,ierr)
periods(1)=.true.
periods(2)=.true.
reorder=.true.

call g1%setTopologyMPI(ndim,dims,periods,reorder)
call g1%setGridID(gridID=gridID)
call g1%setTime(fstart=0,tend=1.d0,dtout=1.d0)
if (g1%fstart .eq. 0) then
  call g1%setMesh(nMesh,leftBdry,rightBdry,nbuf,coordType,gridID)
  call g1%setVariable(variable) !! den,vx,vy,vz,bx,by,bz,ene
  call g1%setMPIWindows()
  call g1%setEoS(eosType=2) !! adiabatic
  call g1%setadiGamma(gam=5.d0/3.d0) !! ratio of heat capacity
  call g1%setCFL(CFL=0.4d0) !! courant number
  call g1%setSlopeLimiter(limiterType=2)
  call g1%setSolverType(solverType=5) !! 4=AdiHLLMHD, 5=AdiHLLDMHD
  call g1%setBoundaryType(boundaryType=3) !! 1=zero gradient, 3=periodic
  call g1%initVariable()
  call g1%exchangeBdryMPI(g1%q,g1%winq)
  call g1%setBoundary(g1%q)
  call g1%writeGrid()
endif
g1%writeFlag=.false.

t0=MPI_WTIME()
do while (g1%t .lt. g1%tend)
   call g1%griddt()
   call g1%evolveGridRK2()
   if(myid .eq. 0) then
     print *,"testSuite.f03: gridID=",g1%gridID,'t=',g1%t,"dt=",g1%dt
   endif
   if(g1%writeFlag .eqv. .true.) then
      t1=MPI_WTIME()
      if(myid .eq. 0) then
        print *,"testSuiteMPI.f03: elapse time=",t1-t0,"secs"
      endif
      call g1%writeGrid()
      g1%writeFlag=.false.
      t0=MPI_WTIME()
   endif
enddo !! end do while
end subroutine CPAlvenWave

subroutine initCPAlvenWave(this,q)
class(grid)::this
double precision,dimension(1-this%nbuf:this%nMesh(1)+this%nbuf,1-this%nbuf:this%nMesh(2)+this%nbuf,this%nvar)::q
double precision,dimension(1-this%nbuf:this%nMesh(1)+this%nbuf,1-this%nbuf:this%nMesh(2)+this%nbuf)::Az
double precision,dimension(1-this%nbuf:this%nMesh(1)+this%nbuf)::xc,xl,xr
double precision,dimension(1-this%nbuf:this%nMesh(2)+this%nbuf)::yc,yl,yr
integer::nx,ny,nvar,nbuf,gridID,eosType
integer::i,j
double precision::gam,rho,vx,vy,vz,bxl,byl,bzl,bxr,byr,bzr,p
double precision::bxc,byc,bzc
double precision::pi,b0
double precision::rc
double precision::ang,x1,x2
double precision::b1,b2,b3,v1,v2,v3
double precision::x,y,dx,dy

pi=4.d0*atan2(1.d0,1.d0)
b0=1.d0
xc=this%xc(1)%coords
yc=this%xc(2)%coords
xl=this%xl(1)%coords
yl=this%xl(2)%coords
xr=this%xr(1)%coords
yr=this%xr(2)%coords
nx=this%nMesh(1)
ny=this%nMesh(2)
gam=this%adiGamma
ang=datan2(2.d0,1.d0)

do j=1,ny+1
  do i=1,nx+1
    x=xl(i)
    y=yl(j)
    x1= x*dcos(ang)+y*dsin(ang)
    x2=-x*dsin(ang)+y*dcos(ang)
    Az(i,j)=0.1d0/(2.d0*pi)*dcos(2.d0*pi*x1)+x2
  enddo
enddo

do j=1,ny
  do i=1,nx
    dx=this%dx(1)%coords(i)
    dy=this%dx(2)%coords(j)
    x1=xc(i)*dcos(ang)+yc(j)*dsin(ang)
    x2=-xc(i)*dsin(ang)+yc(j)*dcos(ang)

    rho=1.d0
    v1=0.d0
    v2=0.1d0*dsin(2.d0*pi*x1)
    v3=0.1d0*dcos(2.d0*pi*x1)

    vx=-v2*dsin(ang)+v1*dcos(ang)
    vy= v2*dcos(ang)+v1*dsin(ang)
    vz= v3

   bxl=(Az(i,j+1)-Az(i,j))/dy
   byl=-(Az(i+1,j)-Az(i,j))/dx
   bzl=0.1d0*dcos(2.d0*pi*x1)
   bxr=(Az(i+1,j+1)-Az(i+1,j))/dy
   byr=-(Az(i+1,j+1)-Az(i,j+1))/dx
   bzr=bzl
     p=0.1d0

    bxc=0.5d0*(bxl+bxr)
    byc=0.5d0*(byl+byr)
    bzc=0.5d0*(bzl+bzr)

    q(i,j,1)=rho
    q(i,j,2)=rho*vx
    q(i,j,3)=rho*vy
    q(i,j,4)=rho*vz
    q(i,j,5)=bxl
    q(i,j,6)=byl
    q(i,j,7)=bzl
    q(i,j,8)=p/(gam-1.d0)+0.5d0*rho*(vx**2.d0+vy**2.d0+vz**2.d0)+0.5d0*(bxc**2.d0+byc**2.d0+bzc**2.d0)
    q(i,j,9)=bxr
    q(i,j,10)=byr
    q(i,j,11)=bzr
  enddo
enddo
end subroutine initCPAlvenWave

subroutine bdryCPAlvenWave(this,q)
class(grid)::this
double precision,dimension(1-this%nbuf:this%nMesh(1)+this%nbuf,1-this%nbuf:this%nMesh(2)+this%nbuf,this%nvar)::q
integer::nx,ny,nbuf,i,j

nbuf=this%nbuf
nx=this%nMesh(1)
ny=this%nMesh(2)

if(this%boundaryType .eq. 3) then !! periodic boundary condition
endif

end subroutine bdryCPAlvenWave


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!! Current Sheet Test
!!!!! (Hawley, J.F., & Stone, J.M., Comp. Phys. Comm. 89, 127 (1995))
!!!!! Scorpio will crash when \beta<0.002, i.e., p=0.002
!!!!! Scorpio is more robust than the Athena code for this test
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine currentSheetMHD2D(gridID)
integer::gridID
type(grid)::g1
integer::ndim,nbuf,coordType,variable(8)
integer::nMesh(2),dims(2)
double precision::leftBdry(2),rightBdry(2)
double precision::pi
logical::periods(2),reorder
double precision::t0,t1
integer::nstep,ierr

nstep=0
pi=4.d0*atan2(1.d0,1.d0)
variable=0
ndim=2
nbuf=2
coordType=1 !! Cartesian
variable(1)=1  !! rho
variable(2)=1  !! vx
variable(3)=1  !! vy
variable(4)=1  !! vz
variable(5)=1  !! bx
variable(6)=1  !! by
variable(7)=1  !! bz
variable(8)=1  !! energy

nMesh(1)=256
nMesh(2)=256
leftBdry(1)=0.d0
rightBdry(1)=2.d0
leftBdry(2)=0.d0
rightBdry(2)=2.d0
dims=(/0,0/)
call MPI_DIMS_CREATE(nprocs,ndim,dims,ierr)
periods(1)=.true.
periods(2)=.true.
reorder=.true.

call g1%setTopologyMPI(ndim,dims,periods,reorder)
call g1%setGridID(gridID=gridID)
call g1%setTime(fstart=0,tend=10.d0,dtout=1.d0)

if (g1%fstart .eq. 0) then
  call g1%setMesh(nMesh,leftBdry,rightBdry,nbuf,coordType,gridID)
  call g1%setVariable(variable) !! den,vx,vy,vz,bx,by,bz,ene
  call g1%setMPIWindows()
  call g1%setEoS(eosType=2) !! adiabatic
  call g1%setadiGamma(gam=5.d0/3.d0) !! ratio of heat capacity
  call g1%setCFL(CFL=0.4d0) !! courant number
  call g1%setSlopeLimiter(limiterType=3)
  call g1%setSolverType(solverType=5) !! 4=AdiHLLMHD, 5=AdiHLLDMHD
  call g1%setBoundaryType(boundaryType=3) !! 1=zero gradient, 3=periodic
  call g1%initVariable()
  call g1%exchangeBdryMPI(g1%q,g1%winq)
  call g1%setBoundary(g1%q)
  call g1%writeGrid()
  call g1%writeGrid_MHD_vtk()
endif
g1%writeFlag=.false.

t0=MPI_WTIME()
do while (g1%t .lt. g1%tend)
   call g1%griddt()
   call g1%evolveGridRK2()

   if(myid .eq. 0) then
     print *,"testSuite.f03: gridID=",g1%gridID,'t=',g1%t,"dt=",g1%dt
   endif

   if(g1%writeFlag .eqv. .true.) then
      t1=MPI_WTIME()
      if(myid .eq. 0) then
        print *,"testSuiteMPI.f03: elapse time=",t1-t0,"secs"
      endif
      call g1%writeGrid()
      call g1%writeGrid_MHD_vtk()
      g1%writeFlag=.false.
      t0=MPI_WTIME()
   endif
enddo !! end do while
end subroutine currentSheetMHD2D

subroutine initCurrentSheetMHD2D(this,q)
class(grid)::this
double precision,dimension(1-this%nbuf:this%nMesh(1)+this%nbuf,1-this%nbuf:this%nMesh(2)+this%nbuf,this%nvar)::q
double precision,dimension(1-this%nbuf:this%nMesh(1)+this%nbuf)::xc,xl,xr
double precision,dimension(1-this%nbuf:this%nMesh(2)+this%nbuf)::yc,yl,yr
integer::nx,ny,nvar,nbuf,gridID,eosType
integer::i,j
double precision::gam,rho,vx,vy,vz,bxl,byl,bzl,bxr,byr,bzr,p
double precision::bxc,byc,bzc
double precision::pi,b0
double precision::rc


pi=4.d0*atan2(1.d0,1.d0)
b0=1.d0
xc=this%xc(1)%coords
yc=this%xc(2)%coords
xl=this%xl(1)%coords
yl=this%xl(2)%coords
xr=this%xr(1)%coords
yr=this%xr(2)%coords
nx=this%nMesh(1)
ny=this%nMesh(2)
gam=this%adiGamma

do j=1,ny
  do i=1,nx

    rho=1.d0

    !vx=0.1d0*dsin(2*pi*yc(j))
    !vy=0.d0
    !vz=0.d0

    vx=0.0
    vy=0.1d0*dsin(2*pi*xc(i))
    vz=0.d0

   bxl=0.d0
   byl=b0
   if(xc(i) .ge. 0.5d0 .and. xc(i) .lt. 1.5d0) then
     byl=-b0
   endif



        if (yc(i)<1.0)then
         bxl =  b0*tanh((yc(j)-0.5)/(1.0/30.0))
        else
         bxl = -b0*tanh((yc(j)-1.5)/(1.0/30.0))
        endif

        byl=0.d0

   bzl=0.d0
   bxr=bxl
   byr=byl
   bzr=bzl

     
    p=0.002d0

    bxc=0.5d0*(bxl+bxr)
    byc=0.5d0*(byl+byr)
    bzc=0.5d0*(bzl+bzr)

    q(i,j,1)=rho
    q(i,j,2)=rho*vx
    q(i,j,3)=rho*vy
    q(i,j,4)=rho*vz
    q(i,j,5)=bxl
    q(i,j,6)=byl
    q(i,j,7)=bzl
    q(i,j,8)=p/(gam-1.d0)+0.5d0*rho*(vx**2.d0+vy**2.d0+vz**2.d0)+0.5d0*(bxc**2.d0+byc**2.d0+bzc**2.d0)
    q(i,j,9)=bxr
    q(i,j,10)=byr
    q(i,j,11)=bzr
  enddo
enddo
end subroutine initCurrentSheetMHD2D

subroutine bdryCurrentSheetMHD2D(this,q)
class(grid)::this
double precision,dimension(1-this%nbuf:this%nMesh(1)+this%nbuf,1-this%nbuf:this%nMesh(2)+this%nbuf,this%nvar)::q
integer::nx,ny,nbuf,i,j

nbuf=this%nbuf
nx=this%nMesh(1)
ny=this%nMesh(2)

if(this%boundaryType .eq. 3) then !! periodic boundary condition
endif

end subroutine bdryCurrentSheetMHD2D

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!! MHD rotor 2D (Balsara and Spicer, 1999)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine rotorMHD2D(gridID)
integer::gridID
type(grid)::g1
integer::ndim,nbuf,coordType,variable(8)
integer::nMesh(2),dims(2)
double precision::leftBdry(2),rightBdry(2)
double precision::pi
logical::periods(2),reorder
double precision::t0,t1
integer::nstep,ierr

nstep=0
pi=4.d0*atan2(1.d0,1.d0)
variable=0
ndim=2
nbuf=2
coordType=1 !! Cartesian
variable(1)=1  !! rho
variable(2)=1  !! vx
variable(3)=1  !! vy
variable(4)=1  !! vz
variable(5)=1  !! bx
variable(6)=1  !! by
variable(7)=1  !! bz
variable(8)=1  !! energy

nMesh(1)=201
nMesh(2)=201
leftBdry(1)=0.d0
rightBdry(1)=1.d0
leftBdry(2)=0.d0
rightBdry(2)=1.d0
dims=(/0,0/)
call MPI_DIMS_CREATE(nprocs,ndim,dims,ierr)
periods(1)=.false.
periods(2)=.false.
reorder=.true.

call g1%setTopologyMPI(ndim,dims,periods,reorder)
call g1%setMesh(nMesh,leftBdry,rightBdry,nbuf,coordType,gridID)
call g1%setVariable(variable) !! den,vx,vy,vz,bx,by,bz,ene
call g1%setMPIWindows()
call g1%setEoS(eosType=2) !! adiabatic
call g1%setadiGamma(gam=1.4d0) !! ratio of heat capacity
call g1%setCFL(CFL=0.4d0) !! courant number
call g1%setTime(fstart=0,tend=0.15d0,dtout=0.15d0)
call g1%setSlopeLimiter(limiterType=1)
call g1%setSolverType(solverType=5) !! 4=AdiHLLMHD, 5=AdiHLLDMHD
call g1%setBoundaryType(boundaryType=1) !! 1=zero gradient, 3=periodic
call g1%initVariable()
call g1%exchangeBdryMPI(g1%q,g1%winq)
call g1%setBoundary(g1%q)
call g1%writeGrid()
call g1%writeGrid_MHD_vtk()
g1%writeFlag=.false.

t0=MPI_WTIME()
do while (g1%t .lt. g1%tend)
   call g1%griddt()
   call g1%evolveGridRK2()
   if(myid .eq. 0) then
     print *,"testSuite.f03: gridID=",g1%gridID,'t=',g1%t,"dt=",g1%dt
   endif

   if(g1%writeFlag .eqv. .true.) then
      t1=MPI_WTIME()
      if(myid .eq. 0) then
        print *,"testSuiteMPI.f03: elapse time=",t1-t0,"secs"
      endif
      call g1%writeGrid()
      call g1%writeGrid_MHD_vtk()
      g1%writeFlag=.false.
      t0=MPI_WTIME()
   endif
enddo !! end do while

end subroutine rotorMHD2D

subroutine initRotorMHD2D(this,q)
class(grid)::this
double precision,dimension(1-this%nbuf:this%nMesh(1)+this%nbuf,1-this%nbuf:this%nMesh(2)+this%nbuf,this%nvar)::q
double precision,dimension(1-this%nbuf:this%nMesh(1)+this%nbuf)::xc,xl,xr
double precision,dimension(1-this%nbuf:this%nMesh(2)+this%nbuf)::yc,yl,yr
integer::nx,ny,nvar,nbuf,gridID,eosType
integer::i,j
double precision::gam,rho,vx,vy,vz,bxl,byl,bzl,bxr,byr,bzr,p
double precision::bxc,byc,bzc
double precision::pi,b0
double precision::r,r0,r1,u0,fr


pi=4.d0*atan2(1.d0,1.d0)
xc=this%xc(1)%coords
yc=this%xc(2)%coords
xl=this%xl(1)%coords
yl=this%xl(2)%coords
xr=this%xr(1)%coords
yr=this%xr(2)%coords
nx=this%nMesh(1)
ny=this%nMesh(2)
gam=this%adiGamma
r0=0.1d0
r1=0.115d0
u0=2.d0
b0=5.d0/dsqrt(4.d0*pi)

do j=1,ny
  do i=1,nx
   r=dsqrt((xc(i)-0.5d0)**2+(yc(j)-0.5d0)**2)
   fr=(r1-r)/(r1-r0)
   if (r .le. r0) then
     rho=10.d0
      vx=-u0*(yc(j)-0.5d0)/r0
      vy= u0*(xc(i)-0.5d0)/r0
   elseif(r .gt. r0 .and. r .lt. r1) then
     rho=1.d0+9.d0*fr
      vx=-fr*u0*(yc(j)-0.5d0)/r
      vy= fr*u0*(xc(i)-0.5d0)/r
   else
     rho=1.d0
      vx=0.d0
      vy=0.d0
   endif


    vz=0.d0
   bxl=b0
   byl=0.d0
   bzl=0.d0
   bxr=bxl
   byr=byl
   bzr=bzl
     p=1.d0

    bxc=0.5d0*(bxl+bxr)
    byc=0.5d0*(byl+byr)
    bzc=0.5d0*(bzl+bzr)

    q(i,j,1)=rho
    q(i,j,2)=rho*vx
    q(i,j,3)=rho*vy
    q(i,j,4)=rho*vz
    q(i,j,5)=bxl
    q(i,j,6)=byl
    q(i,j,7)=bzl
    q(i,j,8)=p/(gam-1.d0)+0.5d0*rho*(vx**2.d0+vy**2.d0+vz**2.d0)+0.5d0*(bxc**2.d0+byc**2.d0+bzc**2.d0)
    q(i,j,9)=bxr
    q(i,j,10)=byr
    q(i,j,11)=bzr
  enddo
enddo
end subroutine initRotorMHD2D

subroutine bdryRotorMHD2D(this,q)
class(grid)::this
double precision,dimension(1-this%nbuf:this%nMesh(1)+this%nbuf,1-this%nbuf:this%nMesh(2)+this%nbuf,this%nvar)::q
integer::nx,ny,nbuf,i,j

nbuf=this%nbuf
nx=this%nMesh(1)
ny=this%nMesh(2)

if(this%boundaryType .eq. 1) then !! zero gradient
!!!!!! left boundary
  if(this%left_mpi .lt. 0) then
    do i=1,nbuf
       q(1-i,:,:)=q(1,:,:)
    enddo
  endif
!!!!!! right boundary
  if(this%right_mpi .lt. 0) then
    do i=1,nbuf
       q(nx+i,:,:)=q(nx,:,:)
    enddo
  endif
!!!!!! up boundary
  if(this%up_mpi .lt. 0) then
    do j=1,nbuf
       q(:,1-j,:)=q(:,1,:)
    enddo
  endif
!!!!!! down boundary
  if(this%down_mpi .lt. 0) then
    do j=1,nbuf
       q(:,ny+j,:)=q(:,ny,:)
    enddo
  endif
endif

end subroutine bdryRotorMHD2D
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!! Magnetic Field Loop Test
!!!!! Toth, G. & Odstrcil, D (1996, JCP, 128,82)
!!!!! DeVore (1991, JCP, 92, 142)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine BFieldLoopTest2D(gridID)
integer::gridID
type(grid)::g1
integer::ndim,nbuf,coordType,variable(8)
integer::nMesh(2),dims(2)
double precision::leftBdry(2),rightBdry(2)
double precision::pi
logical::periods(2),reorder
double precision::t0,t1
integer::nstep,ierr

nstep=0
pi=4.d0*atan2(1.d0,1.d0)
variable=0
ndim=2
nbuf=2
coordType=1 !! Cartesian
variable(1)=1  !! rho
variable(2)=1  !! vx
variable(3)=1  !! vy
variable(4)=1  !! vz
variable(5)=1  !! bx
variable(6)=1  !! by
variable(7)=1  !! bz
variable(8)=1  !! energy

nMesh(1)=256
nMesh(2)=148
leftBdry(1)=-1.d0
rightBdry(1)=1.d0
leftBdry(2)= -1.d0/256.d0*148.d0 !-1.d0/(2.d0*dcos(pi/6.d0))
rightBdry(2)=1/256.d0*148.d0!1.d0/(2.d0*dcos(pi/6.d0))
dims=(/0,0/)
call MPI_DIMS_CREATE(nprocs,ndim,dims,ierr)
periods(1)=.true.
periods(2)=.true.
reorder=.true.

call g1%setTopologyMPI(ndim,dims,periods,reorder)
call g1%setMesh(nMesh,leftBdry,rightBdry,nbuf,coordType,gridID)
call g1%setVariable(variable) !! den,vx,vy,vz,bx,by,bz,ene
call g1%setMPIWindows()
call g1%setEoS(eosType=2) !! adiabatic
call g1%setadiGamma(gam=5.d0/3.d0) !! ratio of heat capacity
call g1%setCFL(CFL=0.4d0) !! courant number
call g1%setTime(fstart=0,tend=16.d0/dsqrt(3.d0),dtout=16.d0/dsqrt(3.d0)/1.d0)
call g1%setSlopeLimiter(limiterType=3)
call g1%setSolverType(solverType=5) !! 4=AdiHLLMHD, 5=AdiHLLDMHD
call g1%setBoundaryType(boundaryType=3) !! 3=periodic
call g1%initVariable()
call g1%exchangeBdryMPI(g1%q,g1%winq)
call g1%setBoundary(g1%q)
call g1%writeGrid()
call g1%writeGrid_MHD_vtk()
g1%writeFlag=.false.

t0=MPI_WTIME()
do while (g1%t .lt. g1%tend)
   call g1%griddt()
   call g1%evolveGridRK2()
   if(myid .eq. 0) then
     !print *,"testSuite.f03: gridID=",g1%gridID,'t=',g1%t,"dt=",g1%dt
   endif
   if(g1%writeFlag .eqv. .true.) then
      t1=MPI_WTIME()
      if(myid .eq. 0)then
         print *,"testSuiteMPI.f03: elapse time=",t1-t0,"secs"
      endif
      call g1%writeGrid()
      call g1%writeGrid_MHD_vtk()
      g1%writeFlag=.false.
      t0=MPI_WTIME()
   endif
enddo !! end do while

end subroutine BFieldLoopTest2D

subroutine initBFieldLoopTest2D(this,q)
class(grid)::this
double precision,dimension(1-this%nbuf:this%nMesh(1)+this%nbuf,1-this%nbuf:this%nMesh(2)+this%nbuf,this%nvar)::q
double precision,dimension(1-this%nbuf:this%nMesh(1)+this%nbuf,1-this%nbuf:this%nMesh(2)+this%nbuf)::Az
double precision,dimension(1-this%nbuf:this%nMesh(1)+this%nbuf)::xc,xl,xr
double precision,dimension(1-this%nbuf:this%nMesh(2)+this%nbuf)::yc,yl,yr
integer::nx,ny,nvar,nbuf,gridID,eosType
integer::i,j
double precision::gam,rho,vx,vy,vz,bxl,byl,bzl,bxr,byr,bzr,p
double precision::bxc,byc,bzc
double precision::pi,A
double precision::r,x,y,dx,dy


pi=4.d0*atan2(1.d0,1.d0)
A=1.d-3
xc=this%xc(1)%coords
yc=this%xc(2)%coords
xl=this%xl(1)%coords
yl=this%xl(2)%coords
xr=this%xr(1)%coords
yr=this%xr(2)%coords
nx=this%nMesh(1)
ny=this%nMesh(2)
gam=this%adiGamma

do j=1,ny+1
  do i=1,nx+1
    x=xl(i)
    y=yl(j)
    Az(i,j)=dmax1(A*(0.3d0-dsqrt(x**2+y**2)),0.d0)
  enddo
enddo

do j=1,ny
  do i=1,nx
    dx=this%dx(1)%coords(i)
    dy=this%dx(2)%coords(j)
    rho=1.d0
    vx=dsin(pi/3.d0)
    vy=dcos(pi/3.d0)
    vz=0.d0
    bzl=0.d0
    bzr=0.d0

     bxl=(Az(i,j+1)-Az(i,j))/dy
     byl=-(Az(i+1,j)-Az(i,j))/dx

     bxr=(Az(i+1,j+1)-Az(i+1,j))/dy
     byr=-(Az(i+1,j+1)-Az(i,j+1))/dx

     p=1.d0

    bxc=0.5d0*(bxl+bxr)
    byc=0.5d0*(byl+byr)
    bzc=0.5d0*(bzl+bzr)

    q(i,j,1)=rho
    q(i,j,2)=rho*vx
    q(i,j,3)=rho*vy
    q(i,j,4)=rho*vz
    q(i,j,5)=bxl
    q(i,j,6)=byl
    q(i,j,7)=bzl
    q(i,j,8)=p/(gam-1.d0)+0.5d0*rho*(vx**2.d0+vy**2.d0+vz**2.d0)+0.5d0*(bxc**2.d0+byc**2.d0+bzc**2.d0)
    q(i,j,9)=bxr
    q(i,j,10)=byr
    q(i,j,11)=bzr
  enddo
enddo

end subroutine initBFieldLoopTest2D

subroutine bdryBFieldLoopTest2D(this,q)
class(grid)::this
double precision,dimension(1-this%nbuf:this%nMesh(1)+this%nbuf,1-this%nbuf:this%nMesh(2)+this%nbuf,this%nvar)::q
integer::nx,ny,nbuf,i,j

nbuf=this%nbuf
nx=this%nMesh(1)
ny=this%nMesh(2)

if(this%boundaryType .eq. 3) then !! periodic boundary condition
endif

end subroutine bdryBFieldLoopTest2D
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!! Spherical Blast Wave Test (Zachary, Malagoli & Colella,
!!!!! (SIAM J. Sci. Comp, 1994, 15, 263)
!!!!! Balsara, D., & Spicer, D., JCP 149, 270 (1999);
!!!!! Londrillo, P. & Del Zanna, L., ApJ 530, 508 (2000).
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine sphericalBlastWaveMHD2D(gridID)
integer::gridID
type(grid)::g1
integer::ndim,nbuf,coordType,variable(8)
integer::nMesh(2),dims(2)
double precision::leftBdry(2),rightBdry(2)
double precision::pi
logical::periods(2),reorder
double precision::t0,t1
integer::nstep,ierr

nstep=0
pi=4.d0*atan2(1.d0,1.d0)
variable=0
ndim=2
nbuf=2
coordType=1 !! Cartesian
variable(1)=1  !! rho
variable(2)=1  !! vx
variable(3)=1  !! vy
variable(4)=1  !! vz
variable(5)=1  !! bx
variable(6)=1  !! by
variable(7)=1  !! bz
variable(8)=1  !! energy

nMesh(1)=400
nMesh(2)=600
leftBdry(1)=-0.5d0
rightBdry(1)=0.5d0
leftBdry(2)=-0.75d0
rightBdry(2)=0.75d0
dims=(/0,0/)
call MPI_DIMS_CREATE(nprocs,ndim,dims,ierr)
periods(1)=.true.
periods(2)=.true.
reorder=.true.

call g1%setTopologyMPI(ndim,dims,periods,reorder)
call g1%setMesh(nMesh,leftBdry,rightBdry,nbuf,coordType,gridID)
call g1%setVariable(variable) !! den,vx,vy,vz,bx,by,bz,ene
call g1%setMPIWindows()
call g1%setEoS(eosType=2) !! adiabatic
call g1%setadiGamma(gam=5.d0/3.d0) !! ratio of heat capacity
call g1%setCFL(CFL=0.4d0) !! courant number
call g1%setTime(fstart=0,tend=0.1d0,dtout=0.1d0)
call g1%setSlopeLimiter(limiterType=3)
call g1%setSolverType(solverType=5) !! 4=AdiHLLMHD, 5=AdiHLLDMHD
call g1%setBoundaryType(boundaryType=3) !! 3=periodic
call g1%initVariable()
call g1%exchangeBdryMPI(g1%q,g1%winq)
call g1%setBoundary(g1%q)
call g1%writeGrid()
call g1%writeGrid_MHD_vtk()
g1%writeFlag=.false.

t0=MPI_WTIME()
do while (g1%t .lt. g1%tend)
   call g1%griddt()
   call g1%evolveGridRK2()
   if(myid .eq. 0) then
     print *,"testSuite.f03: gridID=",g1%gridID,'t=',g1%t,"dt=",g1%dt
   endif
   if(g1%writeFlag .eqv. .true.) then
      t1=MPI_WTIME()
      if(myid .eq. 0) then
        print *,"testSuiteMPI.f03: elapse time=",t1-t0,"secs"
      endif
      call g1%writeGrid()
      call g1%writeGrid_MHD_vtk()
      g1%writeFlag=.false.
      t0=MPI_WTIME()
   endif
enddo !! end do while

end subroutine sphericalBlastWaveMHD2D

subroutine initSphericalBlastWaveMHD2D(this,q)
class(grid)::this
double precision,dimension(1-this%nbuf:this%nMesh(1)+this%nbuf,1-this%nbuf:this%nMesh(2)+this%nbuf,this%nvar)::q
double precision,dimension(1-this%nbuf:this%nMesh(1)+this%nbuf)::xc,xl,xr
double precision,dimension(1-this%nbuf:this%nMesh(2)+this%nbuf)::yc,yl,yr
integer::nx,ny,nvar,nbuf,gridID,eosType
integer::i,j
double precision::gam,rho,vx,vy,vz,bxl,byl,bzl,bxr,byr,bzr,p
double precision::bxc,byc,bzc
double precision::pi,b0
double precision::rc


pi=4.d0*atan2(1.d0,1.d0)
b0=1.d0!dsqrt(4.d0*pi/2.d0)
xc=this%xc(1)%coords
yc=this%xc(2)%coords
xl=this%xl(1)%coords
yl=this%xl(2)%coords
xr=this%xr(1)%coords
yr=this%xr(2)%coords
nx=this%nMesh(1)
ny=this%nMesh(2)
gam=this%adiGamma

do j=1,ny
  do i=1,nx
   rc=dsqrt(xc(i)**2+yc(j)**2)

    rho=1.d0
    vx=0.d0
    vy=0.d0
    vz=0.d0
   bxl=b0
   byl=b0
   bzl=0.d0
   bxr=bxl
   byr=byl
   bzr=bzl
     p=0.1d0
    if(rc .lt. 0.1d0) then
     p=10.d0
    endif

    bxc=0.5d0*(bxl+bxr)
    byc=0.5d0*(byl+byr)
    bzc=0.5d0*(bzl+bzr)

    q(i,j,1)=rho
    q(i,j,2)=rho*vx
    q(i,j,3)=rho*vy
    q(i,j,4)=rho*vz
    q(i,j,5)=bxl
    q(i,j,6)=byl
    q(i,j,7)=bzl
    q(i,j,8)=p/(gam-1.d0)+0.5d0*rho*(vx**2.d0+vy**2.d0+vz**2.d0)+0.5d0*(bxc**2.d0+byc**2.d0+bzc**2.d0)
    q(i,j,9)=bxr
    q(i,j,10)=byr
    q(i,j,11)=bzr
  enddo
enddo

end subroutine initSphericalBlastWaveMHD2D

subroutine bdrySphericalBlastWaveMHD2D(this,q)
class(grid)::this
double precision,dimension(1-this%nbuf:this%nMesh(1)+this%nbuf,1-this%nbuf:this%nMesh(2)+this%nbuf,this%nvar)::q
integer::nx,ny,nbuf,i,j

nbuf=this%nbuf
nx=this%nMesh(1)
ny=this%nMesh(2)

if(this%boundaryType .eq. 3) then !! periodic boundary condition
endif

end subroutine bdrySphericalBlastWaveMHD2D

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!! Double Mach Reflection (1984, JCP, 54, 115)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine doubleMachReflection(gridID)
integer::gridID
type(grid)::g1
integer::ndim,nbuf,coordType,variable(8)
integer::nMesh(2),dims(2)
double precision::leftBdry(2),rightBdry(2)
logical::periods(2),reorder
double precision::t0,t1
integer::nstep,ierr

variable=0
ndim=2
nbuf=2
coordType=1 !! Cartesian
variable(1)=1
variable(2)=1
variable(3)=1
variable(8)=1

nMesh(1)=520
nMesh(2)=160
leftBdry(1)=0.d0
leftBdry(2)=0.d0
rightBdry(1)=3.25d0
rightBdry(2)=1.d0
dims=(/0,0/)
call MPI_DIMS_CREATE(nprocs,ndim,dims,ierr)
periods(1)=.false.
periods(2)=.false.
reorder=.true.

call g1%setTopologyMPI(ndim,dims,periods,reorder)
call g1%setMesh(nMesh,leftBdry,rightBdry,nbuf,coordType,gridID)
call g1%setVariable(variable) !! den,vx,vy,ene
call g1%setMPIWindows()
call g1%setEoS(eosType=2) !! adiabatic
call g1%setadiGamma(gam=1.4d0) !! ratio of heat capacity Cp/Cv
call g1%setCFL(CFL=0.7d0) !! Courant number
call g1%setTime(fstart=0,tend=0.2d0,dtout=0.2d0)
call g1%setSlopeLimiter(limiterType=3)
call g1%setSolverType(solverType=3)
call g1%setBoundaryType(boundaryType=0)
call g1%initVariable()
call g1%exchangeBdryMPI(g1%q,g1%winq)
call g1%setBoundary(g1%q)
call g1%writeGrid()
call g1%writeGrid_HD_vtk()
g1%writeFlag=.false.

t0=MPI_WTIME()
do while (g1%t .lt. g1%tend)
   call g1%griddt()
   call g1%evolveGridRK2()
   if(myid .eq. 0) then
     !print *,"testSuite.f03: gridID=",g1%gridID,'t=',g1%t,"dt=",g1%dt
   endif
   if(g1%writeFlag .eqv. .true.) then
      t1=MPI_WTIME()
      if(myid .eq. 0) then
        print *,"testSuiteMPI.f03: elapse time=",t1-t0,"secs"
      endif
      call g1%writeGrid()
      call g1%writeGrid_HD_vtk
      g1%writeFlag=.false.
      t0=MPI_WTIME()
   endif
enddo !! end do while

end subroutine doubleMachReflection

subroutine initDoubleMachReflection(this,q)
class(grid)::this
double precision,dimension(1-this%nbuf:this%nMesh(1)+this%nbuf,1-this%nbuf:this%nMesh(2)+this%nbuf,this%nvar)::q
double precision,dimension(1-this%nbuf:this%nMesh(1)+this%nbuf)::xc
double precision,dimension(1-this%nbuf:this%nMesh(2)+this%nbuf)::yc
integer::nx,ny,nvar,nbuf,gridID,eosType
integer::i,j
double precision::d0,e0,u0,v0,xshock

xc=this%xc(1)%coords
yc=this%xc(2)%coords
nx=this%nMesh(1)
ny=this%nMesh(2)

d0=8.d0
e0=291.25d0
u0= 8.25d0*dsqrt(3.d0)/2.d0
v0=-8.25d0*0.5d0

do j=1,ny
  do i=1,nx
     xshock=0.1666666666+yc(j)/dsqrt(3.d0)
     q(i,j,1)=1.4d0
     q(i,j,2)=0.d0
     q(i,j,3)=0.d0
     q(i,j,4)=2.5d0
     if(xc(i) .le. xshock) then
        q(i,j,1)=d0
        q(i,j,2)=d0*u0
        q(i,j,3)=d0*v0
        q(i,j,4)=e0+0.5d0*d0*(u0**2+v0**2)
     endif
  enddo
enddo
end subroutine initDoubleMachReflection

subroutine bdryDoubleMachReflection(this,q)
class(grid)::this
double precision,dimension(1-this%nbuf:this%nMesh(1)+this%nbuf,1-this%nbuf:this%nMesh(2)+this%nbuf,this%nvar)::q
double precision,dimension(1-this%nbuf:this%nMesh(1)+this%nbuf)::xc
double precision,dimension(1-this%nbuf:this%nMesh(2)+this%nbuf)::yc
integer::nx,ny,nbuf,i,j
double precision::d0,e0,u0,v0,xshock

d0=8.d0
e0=291.25d0
u0= 8.25d0*dsqrt(3.d0)/2.d0
v0=-8.25d0*0.5d0
nbuf=this%nbuf
nx=this%nMesh(1)
ny=this%nMesh(2)
xc=this%xc(1)%coords
yc=this%xc(2)%coords

!!!!!!!!! left boundary !!!!!!!!!!
if(this%left_mpi .lt. 0) then
  do j=1,ny
    do i=1-nbuf,0
       q(i,j,1)=d0
       q(i,j,2)=d0*u0
       q(i,j,3)=d0*v0
       q(i,j,4)=e0+0.5d0*d0*(u0**2+v0**2)
    enddo
  enddo
endif
!!!!!!!! up boundary !!!!!!!!!
if(this%up_mpi .lt. 0)then
  do j=1-nbuf,0
    do i=1,nx
       if(xc(i) .le. 0.1666666666) then
         q(i,j,1)=d0
         q(i,j,2)=d0*u0
         q(i,j,3)=d0*v0
         q(i,j,4)=e0+0.5d0*d0*(u0**2+v0**2)
       else
         q(i,j,1)=q(i,-j+1,1)
         q(i,j,2)=q(i,-j+1,2)
         q(i,j,3)=-q(i,-j+1,3)
         q(i,j,4)=q(i,-j+1,4)
       endif
    enddo
  enddo
endif
!!!!!!!! bottom boundary !!!!!!!!!
if(this%down_mpi .lt. 0) then
  xshock=0.1666666666+(1.d0+20.d0*this%t)/dsqrt(3.d0)

  do j=ny+1,ny+nbuf
    do i=1,nx
       if(xc(i) .le. xshock) then
         q(i,j,1)=d0
         q(i,j,2)=d0*u0
         q(i,j,3)=d0*v0
         q(i,j,4)=e0+0.5d0*d0*(u0**2+v0**2)
       else
         q(i,j,1)=1.4d0
         q(i,j,2)=0.d0
         q(i,j,3)=0.d0
         q(i,j,4)=2.5d0
       endif
    enddo
  enddo
endif
!!!!!!! right boundary !!!!!!!
if(this%right_mpi .lt. 0) then
  do j=1,ny
    do i=1,nbuf
       q(nx+i,j,:)=q(nx-i+1,j,:)
    enddo
  enddo
endif
end subroutine bdryDoubleMachReflection



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!! Orszag & Tang (J. Fluid Mech., 90, 129, 1998)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine OrsagTangVortex(gridID)  
!!!!! zwg : gridID=13 for this ADI OrsagTangVortex problem


integer::gridID  
!!!!zwg: gridID=13 for OrsagTangVortex
 
type(grid)::g1
  
!!!! zwg : the definition of type grid is in the file of gridModule.f03
!!!! which include a lot informatuons

integer::ndim,nbuf,coordType,variable(8)
!!! zwg ndim=2, nbuf=2,
!!! zwg: coordType=1 means using Cartesian coordinate
 
integer::nMesh(2),dims(2)

!!!! zwg: nmesh(2) numbers of grid points in x and y directions for 2D problems
!!!! zwg: for 1D plblems nmesh should be nmesh(1), and for 3D problems nmesh should be nmesh(3)
!!!! zwg: but what is dims(2)???? to be defined

double precision::leftBdry(2),rightBdry(2)
!!!!! zwg: leftmost positions and rightmost positions
!!!!! zwg moreover, left means the minmum positions and right means the maxum positions

double precision::pi

logical::periods(2),reorder

double precision::t0,t1

integer::nstep,ierr,idims


nstep=0
pi=4.d0*atan2(1.d0,1.d0)

variable=0
ndim=2
nbuf=2

coordType=1 !! Cartesian
variable(1)=1  !! rho
variable(2)=1  !! vx
variable(3)=1  !! vy
variable(4)=1  !! vz
variable(5)=1  !! bx
variable(6)=1  !! by
variable(7)=1  !! bz
variable(8)=1  !! energy

nMesh(1)=401
nMesh(2)=403

leftBdry(1)=0.d0
rightBdry(1)=1.d0
leftBdry(2)=0.d0
rightBdry(2)=1.d0

dims=(/2,4/) 
!!! zwg: initialize dims,which represents the cores in each direction
!!! zwg: 

call MPI_DIMS_CREATE(nprocs,ndim,dims,ierr)
!!!! zwg : to get processes in each diractions
!!!! zwg : for 1D problem, dims must be dims(1)
!!!! zwg : for 2D problem, dims must be dims(2)
!!!! zwg : for 3D problem, dims must be dims(3) 

!!!! zwg : also, you can reset dims here!
!dims(1)=1  
!dims(2)=nprocs






periods(1)=.true.
periods(2)=.true.


reorder=.true.

!!! zwg: what is reorder used for????
!!! zwg: usded for MPI_CART_CREATE which is included (or called ) in subroutine setTopologyMPI
!!! setTopologyMPI


call g1%setTopologyMPI(ndim,dims,periods,reorder)
!!!print*, "set Topology MPI done!!!"

!!!!!! zwg added for *.pvts output formate
!if(myid==0) then
!do idims=1, ndim

  !if(  mod(nMesh(idims), dims(idims)).ne. 0) then 
   
     !print *, "testSuiteMPI.f90: not suitable partition of whole domain for totalprocess "

     !stop
  !endif

!enddo
!endif  !!! end of if(myid==0) then

call MPI_BARRIER(MPI_COMM_WORLD,ierr)


call g1%setMesh(nMesh,leftBdry,rightBdry,nbuf,coordType,gridID)

call g1%setVariable(variable) !! den,vx,vy,vz,bx,by,bz,ene

!!!! zwg: remains to be further investigation!!!!!
call g1%setMPIWindows()

call g1%setEoS(eosType=2) !! adiabatic
call g1%setadiGamma(gam=5.d0/3.d0) !! ratio of heat capacity
call g1%setCFL(CFL=0.4d0) !! courant number
call g1%setTime(fstart=0,tend=10.0d0,dtout=0.5d0)
call g1%setSlopeLimiter(limiterType=3)
call g1%setSolverType(solverType=5) !! 4=AdiHLLMHD, 5=AdiHLLDMHD
call g1%setBoundaryType(boundaryType=3) !! 3=periodic
call g1%initVariable()
call g1%exchangeBdryMPI(g1%q,g1%winq)
call g1%setBoundary(g1%q)
call g1%writeGrid()
call g1%writeGrid_MHD_vtk()
g1%writeFlag=.false.

t0=MPI_WTIME()
do while (g1%t .lt. g1%tend)
   call g1%griddt()
   call g1%evolveGridRK2()
   if(myid .eq. 0) then
     print *,"testSuite.f03: gridID=",g1%gridID,'t=',g1%t,"dt=",g1%dt
   endif
   if(g1%writeFlag .eqv. .true.) then
      t1=MPI_WTIME()
      if(myid .eq. 0) then
        print *,"testSuiteMPI.f03: elapse time=",t1-t0,"secs"
      endif
      call g1%writeGrid()
      call g1%writeGrid_MHD_vtk()
      g1%writeFlag=.false.
      t0=MPI_WTIME()
   endif
enddo !! end do while

end subroutine OrsagTangVortex

subroutine initOrsagTangVortex(this,q)
class(grid)::this
double precision,dimension(1-this%nbuf:this%nMesh(1)+this%nbuf,1-this%nbuf:this%nMesh(2)+this%nbuf,this%nvar)::q
double precision,dimension(1-this%nbuf:this%nMesh(1)+this%nbuf)::xc,xl,xr
double precision,dimension(1-this%nbuf:this%nMesh(2)+this%nbuf)::yc,yl,yr
integer::nx,ny,nvar,nbuf,gridID,eosType
integer::i,j
double precision::gam,rho,vx,vy,vz,bxl,byl,bzl,bxr,byr,bzr,p
double precision::bxc,byc,bzc
double precision::pi,b0


pi=4.d0*atan2(1.d0,1.d0)
b0=1.d0/dsqrt(4.d0*pi)
xc=this%xc(1)%coords
yc=this%xc(2)%coords
xl=this%xl(1)%coords
yl=this%xl(2)%coords
xr=this%xr(1)%coords
yr=this%xr(2)%coords
nx=this%nMesh(1)
ny=this%nMesh(2)
gam=this%adiGamma

do j=1,ny
  do i=1,nx

    rho=25.d0/(36.d0*pi)
    vx=-dsin(2*pi*yc(j))
    vy= dsin(2*pi*xc(i))
    vz=0.d0
   bxl=-b0*dsin(2.d0*pi*yc(j))
   byl= b0*dsin(4.d0*pi*xc(i))
   bzl=0.d0
   bxr=bxl
   byr=byl
   bzr=bzl
     p=5.d0/(12.d0*pi)

    bxc=0.5d0*(bxl+bxr)
    byc=0.5d0*(byl+byr)
    bzc=0.5d0*(bzl+bzr)

    q(i,j,1)=rho
    q(i,j,2)=rho*vx
    q(i,j,3)=rho*vy
    q(i,j,4)=rho*vz
    q(i,j,5)=bxl
    q(i,j,6)=byl
    q(i,j,7)=bzl
    q(i,j,8)=p/(gam-1.d0)+0.5d0*rho*(vx**2.d0+vy**2.d0+vz**2.d0)+0.5d0*(bxc**2.d0+byc**2.d0+bzc**2.d0)
    q(i,j,9)=bxr
    q(i,j,10)=byr
    q(i,j,11)=bzr
  enddo
enddo
end subroutine initOrsagTangVortex

subroutine bdryOrsagTangVortex(this,q)
class(grid)::this
double precision,dimension(1-this%nbuf:this%nMesh(1)+this%nbuf,1-this%nbuf:this%nMesh(2)+this%nbuf,this%nvar)::q
integer::nx,ny,nbuf,i,j

nbuf=this%nbuf
nx=this%nMesh(1)
ny=this%nMesh(2)

if(this%boundaryType .eq. 3) then !! periodic boundary condition
endif

end subroutine bdryOrsagTangVortex




!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!! Ryu & Jones, (1995, AJ, 442, 228), Figure 2a
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine RJ2aShockTube1D(gridID)
integer::gridID
type(grid)::g1
integer::ndim,nbuf,coordType,variable(8)
integer::nMesh(1),dims(1)
double precision::leftBdry(1),rightBdry(1)
logical::periods(1),reorder
double precision::t0,t1
integer::nstep,ierr

nstep=0
variable=0
ndim=1
nbuf=2
coordType=1 !! Cartesian
variable(1)=1  !! rho
variable(2)=1  !! vx
variable(3)=1  !! vy
variable(4)=1  !! vz
variable(5)=1  !! bx
variable(6)=1  !! by
variable(7)=1  !! bz
variable(8)=1  !! energy

nMesh(1)=1600
leftBdry(1)=0.d0
rightBdry(1)=1.d0
dims(1)=nprocs
periods(1)=.false.
reorder=.false.

call g1%setTopologyMPI(ndim,dims,periods,reorder)
call g1%setMesh(nMesh,leftBdry,rightBdry,nbuf,coordType,gridID)
call g1%setVariable(variable) !! den,vx,vy,vz,bx,by,bz,ene
call g1%setMPIWindows()
call g1%setEoS(eosType=2) !! adiabatic
call g1%setadiGamma(gam=5.d0/3.d0) !! ratio of heat capacity
call g1%setCFL(CFL=0.8d0) !! courant number
call g1%setTime(fstart=0,tend=0.2d0,dtout=0.2d0)
call g1%setSlopeLimiter(limiterType=3)
call g1%setSolverType(solverType=5) !! 4=AdiHLLMHD, 5=AdiHLLDMHD
call g1%setBoundaryType(boundaryType=1)
call g1%initVariable()
call g1%exchangeBdryMPI(g1%q,g1%winq)
call g1%setBoundary(g1%q)
call g1%writeGrid()
g1%writeFlag=.false.

t0=MPI_WTIME()
do while (g1%t .lt. g1%tend)
   call g1%griddt()
   call g1%evolveGridRK2()
   nstep=nstep+1
   !print *,nstep," testSuite.f03: gridID=",g1%gridID,'t=',g1%t,"dt=",g1%dt

   if(g1%writeFlag .eqv. .true.) then
      t1=MPI_WTIME()
      if(myid .eq. 0) then
        print *,"testSuiteMPI.f03: elapse time=",t1-t0,"secs"
      endif
      call g1%writeGrid()
      g1%writeFlag=.false.
      t0=MPI_WTIME()
   endif
enddo !! end do while
end subroutine RJ2aShockTube1D

subroutine initRJ2aShockTube1D(this,q)
class(grid)::this
double precision,dimension(1-this%nbuf:this%nMesh(1)+this%nbuf,this%nvar)::q
double precision,dimension(1-this%nbuf:this%nMesh(1)+this%nbuf)::xc
integer::nx,nvar,nbuf,gridID,eosType
integer::i
double precision::gam,rho,vx,vy,vz,bxl,byl,bzl,bxr,byr,bzr,p
double precision::bxc,byc,bzc
double precision::pi

pi=atan2(1.d0,1.d0)*4.d0
xc=this%xc(1)%coords
nx=this%nMesh(1)
gam=this%adiGamma

do i=1,nx
    if(xc(i) .lt. 0.5d0) then
      rho=1.08d0
       vx=1.2d0
       vy=0.01d0
       vz=0.5d0
       bxl=2.d0/dsqrt(4.d0*pi)
       bxr=2.d0/dsqrt(4.d0*pi)
       byl=3.6d0/dsqrt(4.d0*pi)
       byr=3.6d0/dsqrt(4.d0*pi)
       bzl=2.d0/dsqrt(4.d0*pi)
       bzr=2.d0/dsqrt(4.d0*pi)
         p=0.95d0
    else
       rho=1.d0
        vx=0.d0
        vy=0.d0
        vz=0.d0
        bxl=2.d0/dsqrt(4.d0*pi)
        bxr=2.d0/dsqrt(4.d0*pi)
        byl=4.d0/dsqrt(4.d0*pi)
        byr=4.d0/dsqrt(4.d0*pi)
        bzl=2.d0/dsqrt(4.d0*pi)
        bzr=2.d0/dsqrt(4.d0*pi)
          p=1.d0
    endif
    bxc=0.5d0*(bxl+bxr)
    byc=0.5d0*(byl+byr)
    bzc=0.5d0*(bzl+bzr)

    q(i,1)=rho
    q(i,2)=rho*vx
    q(i,3)=rho*vy
    q(i,4)=rho*vz
    q(i,5)=bxl
    q(i,6)=byl
    q(i,7)=bzl
    q(i,8)=p/(gam-1.d0)+0.5d0*rho*(vx**2.d0+vy**2.d0+vz**2.d0)+0.5d0*(bxc**2.d0+byc**2.d0+bzc**2.d0)
    q(i,9)=bxr
    q(i,10)=byr
    q(i,11)=bzr
enddo
end subroutine initRJ2aShockTube1D

subroutine bdryRJ2aShockTube1D(this,q)
class(grid)::this
double precision,dimension(1-this%nbuf:this%nMesh(1)+this%nbuf,this%nvar)::q
integer::nbuf,nx,i

nbuf=this%nbuf
nx=this%nMesh(1)

if(this%boundaryType .eq. 1) then !! zero gradient
  if(this%left_mpi .lt. 0) then
    do i=1,nbuf
       q(i-nbuf,:)=q(nbuf-i+1,:)
    enddo
  endif
  if(this%right_mpi .lt. 0) then
    do i=1,nbuf
       q(nx+i,:)=q(nx-i+1,:)
    enddo
  endif
endif
end subroutine bdryRJ2aShockTube1D
 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!! Brio-Wu Shocktube, (1988, JCP, 75, 400)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine BrioWuShockTube1D(gridID)
integer::gridID
type(grid)::g1
integer::ndim,nbuf,coordType,variable(8)
integer::nMesh(1),dims(1)
double precision::leftBdry(1),rightBdry(1)
logical::periods(1),reorder
double precision::t0,t1
integer::nstep,ierr

nstep=0
variable=0
ndim=1
nbuf=2
coordType=1 !! Cartesian
variable(1)=1  !! rho
variable(2)=1  !! vx
variable(3)=1  !! vy
variable(4)=1  !! vz
variable(5)=1  !! bx
variable(6)=1  !! by
variable(7)=1  !! bz
variable(8)=1  !! energy

nMesh(1)=800
leftBdry(1)=0.d0
rightBdry(1)=1.d0
dims(1)=nprocs
periods(1)=.false.
reorder=.true.

call g1%setTopologyMPI(ndim,dims,periods,reorder)
call g1%setMesh(nMesh,leftBdry,rightBdry,nbuf,coordType,gridID)
call g1%setVariable(variable) !! den,vx,vy,vz,bx,by,bz,ene
call g1%setMPIWindows()
call g1%setEoS(eosType=2) !! adiabatic
call g1%setadiGamma(gam=2.d0) !! ratio of heat capacity
call g1%setCFL(CFL=0.8d0) !! courant number
call g1%setTime(fstart=0,tend=0.1d0,dtout=0.1d0/200.d0)
call g1%setSlopeLimiter(limiterType=3)
call g1%setSolverType(solverType=5) !! 4=AdiHLLMHD, 5=AdiHLLDMHD
call g1%setBoundaryType(boundaryType=1)
call g1%initVariable()
call g1%exchangeBdryMPI(g1%q,g1%winq)

call g1%setBoundary(g1%q)
call g1%writeGrid()
g1%writeFlag=.false.

t0=MPI_WTIME()
do while (g1%t .lt. g1%tend)
   call g1%griddt()
   call g1%evolveGridRK2()
   nstep=nstep+1
   if(myid .eq. 0) then
     !print *,nstep," testSuite.f03: gridID=",g1%gridID,'t=',g1%t,"dt=",g1%dt
   endif
   if(g1%writeFlag .eqv. .true.) then
      t1=MPI_WTIME()
      if(myid .eq. 0) then
        print *,"testSuiteMPI.f03: elapse time=",t1-t0,"secs"
      endif
      call g1%writeGrid()
      g1%writeFlag=.false.
      t0=MPI_WTIME()
   endif
enddo !! end do while
end subroutine BrioWuShockTube1D

subroutine initBrioWuShockTube1D(this,q)
class(grid)::this
double precision,dimension(1-this%nbuf:this%nMesh(1)+this%nbuf,this%nvar)::q
double precision,dimension(1-this%nbuf:this%nMesh(1)+this%nbuf)::xc
integer::nx,nvar,nbuf,gridID,eosType
integer::i
double precision::gam,rho,vx,vy,vz,bxl,byl,bzl,bxr,byr,bzr,p
double precision::bxc,byc,bzc

xc=this%xc(1)%coords
nx=this%nMesh(1)
gam=this%adiGamma

do i=1,nx
    if(xc(i) .lt. 0.5d0) then
      rho=1.d0
       vx=0.d0
       vy=0.d0
       vz=0.d0
       bxl=0.75d0
       bxr=0.75d0
       byl=1.d0
       byr=1.d0
       bzl=0.d0
       bzr=0.d0
         p=1.d0
    else
       rho=0.125d0
        vx=0.d0
        vy=0.d0
        vz=0.d0
        bxl=0.75d0
        bxr=0.75d0
        byl=-1.d0
        byr=-1.d0
        bzl=0.d0
        bzr=0.d0
          p=0.1d0
    endif
    bxc=0.5d0*(bxl+bxr)
    byc=0.5d0*(byl+byr)
    bzc=0.5d0*(bzl+bzr)

    q(i,1)=rho
    q(i,2)=rho*vx
    q(i,3)=rho*vy
    q(i,4)=rho*vz
    q(i,5)=bxl
    q(i,6)=byl
    q(i,7)=bzl
    q(i,8)=p/(gam-1.d0)+0.5d0*rho*(vx**2.d0+vy**2.d0+vz**2.d0)+0.5d0*(bxc**2.d0+byc**2.d0+bzc**2.d0)
    q(i,9)=bxr
    q(i,10)=byr
    q(i,11)=bzr
enddo

end subroutine initBrioWuShockTube1D

subroutine bdryBrioWuShockTube1D(this,q)
class(grid)::this
double precision,dimension(1-this%nbuf:this%nMesh(1)+this%nbuf,this%nvar)::q
integer::nbuf,nx,i

nbuf=this%nbuf
nx=this%nMesh(1)

if(this%boundaryType .eq. 1) then !! zero gradient
  if(this%left_mpi .lt. 0) then
    do i=1,nbuf
       q(i-nbuf,:)=q(nbuf-i+1,:)
    enddo
  endif
  if(this%right_mpi .lt. 0) then
    do i=1,nbuf
       q(nx+i,:)=q(nx-i+1,:)
    enddo
  endif
endif
end subroutine bdryBrioWuSHockTube1D
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!! Woodward-Colella blast wave problem (1984, JCP, 54, 115)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine WCShockTube1D(gridID)
integer::gridID
type(grid)::g1
integer::ndim,nbuf,coordType,variable(8)
integer::nMesh(1),dims(1)
double precision::leftBdry(1),rightBdry(1)
logical::periods(1),reorder
double precision::t0,t1
integer::nstep

nstep=0
variable=0
ndim=1
nbuf=2
coordType=1
variable(1)=1
variable(2)=1
variable(8)=1

    nMesh(1)=1201
 leftBdry(1)=0.d0
rightBdry(1)=1.d0
dims(1)=nprocs
periods(1)=.false.
reorder=.false.

call g1%setTopologyMPI(ndim,dims,periods,reorder)
call g1%setMesh(nMesh,leftBdry,rightBdry,nbuf,coordType,gridID)
call g1%setVariable(variable) !! den,vx,ene
call g1%setMPIWindows()
call g1%setEoS(eosType=2) !! adiabatic
call g1%setadiGamma(gam=1.4d0) !! gamma
call g1%setCFL(CFL=0.9d0) !! Courant number
call g1%setTime(fstart=0,tend=0.038d0,dtout=0.038d0/200.d0)
call g1%setSlopeLimiter(limiterType=3) !! won't be correct  with limiterType=0
call g1%setSolverType(solverType=3) !! won't be correct with solverType=2
call g1%setBoundaryType(boundaryType=2)
call g1%initVariable()
call g1%exchangeBdryMPI(g1%q,g1%winq)
call g1%setBoundary(g1%q)
call g1%writeGrid()
g1%writeFlag=.false.

t0=MPI_WTIME()
do while (g1%t .lt. g1%tend)
   call g1%griddt()
   call g1%evolveGridRK2()
   nstep=nstep+1
   if(myid .eq. 0) then
     print *,nstep," testSuite.f03: gridID=",g1%gridID,'t=',g1%t,"dt=",g1%dt
   endif
   if(g1%writeFlag .eqv. .true.) then
      t1=MPI_WTIME()
      if(myid .eq. 0) then
        print *,"testSuiteMPI.f03: elapse time=",t1-t0,"secs"
      endif
      call g1%writeGrid()
      g1%writeFlag=.false.
      t0=MPI_WTIME()
   endif
enddo !! end do while

end subroutine WCShockTube1D

subroutine initWCShockTube1D(this,q)
class(grid)::this
double precision,dimension(1-this%nbuf:this%nMesh(1)+this%nbuf,this%nvar)::q
double precision,dimension(1-this%nbuf:this%nMesh(1)+this%nbuf)::xc
integer::nx,nvar,nbuf,gridID,eosType
integer::i

xc=this%xc(1)%coords
nx=this%nMesh(1)
   do i=1,nx
      if(xc(i) .lt. 0.1d0) then
         q(i,3)=2500.d0
      elseif(xc(i) .gt. 0.1d0 .and. xc(i) .lt. 0.9d0) then
         q(i,3)=0.025d0
      else
         q(i,3)=250.d0
      endif
         q(i,1)=1.d0
         q(i,2)=0.d0
   enddo
end subroutine initWCShockTube1D

subroutine bdryWCShockTube1D(this,q)
class(grid)::this
double precision,dimension(1-this%nbuf:this%nMesh(1)+this%nbuf,this%nvar)::q
integer::nbuf,nx,i

nbuf=this%nbuf
nx=this%nMesh(1)

if(this%boundaryType .eq. 2) then !! reflective boundary condition
  if(this%left_mpi .lt. 0) then
    do i=1,nbuf
       q(i-nbuf,:)=q(nbuf-i+1,:)
    enddo
    do i=1,nbuf
       q(i-nbuf,2)=-q(nbuf-i+1,2)
    enddo
  endif
  if(this%right_mpi .lt. 0) then
    do i=1,nbuf
       q(nx+i,:)=q(nx-i+1,:)
    enddo
    do i=1,nbuf
       q(nx+i,2)  =-q(nx-i+1,2)
    enddo
  endif
endif
end subroutine bdryWCShockTube1D

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!! Sod Shocktube problem (Sod, 1978, JCP, 27, 1)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine sodShockTube1D(gridID)
integer::gridID
type(grid)::g1
integer::ndim,nbuf,coordType,variable(8)
integer::nMesh(1),dims(1)
double precision::leftBdry(1),rightBdry(1)
logical::periods(1),reorder
double precision::t0,t1
integer::nstep

nstep=0
variable=0
ndim=1
nbuf=2
coordType=1
variable(1)=1
variable(2)=1
variable(8)=1

nMesh(1)=400
leftBdry(1)=-0.5d0
rightBdry(1)=0.5d0
dims(1)=nprocs
periods(1)=.false.
reorder=.false.

call g1%setTopologyMPI(ndim,dims,periods,reorder)
call g1%setMesh(nMesh,leftBdry,rightBdry,nbuf,coordType,gridID)
call g1%setVariable(variable) !! den,vx,ene
call g1%setMPIWindows()
call g1%setEoS(eosType=2) !! adiabatic
call g1%setadiGamma(gam=1.4d0) !! gamma
call g1%setCFL(CFL=0.6d0) !! Courant number
call g1%setTime(fstart=0,tend=0.2d0,dtout=0.2d0)
call g1%setSlopeLimiter(limiterType=1)
call g1%setSolverType(solverType=3)
call g1%setBoundaryType(boundaryType=1)
call g1%initVariable()
call g1%exchangeBdryMPI(g1%q,g1%winq)
call g1%setBoundary(g1%q)
call g1%writeGrid()
g1%writeFlag=.false.

t0=MPI_WTIME()
do while (g1%t .lt. g1%tend)
   call g1%griddt()
   call g1%evolveGridRK2()
   nstep=nstep+1
   if(myid .eq. 0) then
     !print *,nstep," testSuite.f03: gridID=",g1%gridID,'t=',g1%t,"dt=",g1%dt
   endif
   if(g1%writeFlag .eqv. .true.) then
      t1=MPI_WTIME()
      if(myid .eq. 0) then
        print *,"testSuiteMPI.f03: elapse time=",t1-t0,"secs"
      endif
      call g1%writeGrid()
      g1%writeFlag=.false.
      t0=MPI_WTIME()
   endif
enddo !! end do while
end subroutine sodShockTube1D

subroutine initSodShockTube1D(this,q)
class(grid)::this
double precision,dimension(1-this%nbuf:this%nMesh(1)+this%nbuf,this%nvar)::q
double precision,dimension(1-this%nbuf:this%nMesh(1)+this%nbuf)::xc
integer::nx,nvar,nbuf,gridID,eosType
integer::i

xc=this%xc(1)%coords
nx=this%nMesh(1)
   do i=1,nx
      if(xc(i) .lt. 0.d0) then
         q(i,1)=1.d0
         q(i,3)=2.5d0
      else
         q(i,1)=0.125d0
         q(i,3)=0.25d0
      endif
         q(i,2)=0.d0
   enddo
end subroutine initSodShockTube1D

subroutine bdrySodShockTube1D(this,q)
class(grid)::this
double precision,dimension(1-this%nbuf:this%nMesh(1)+this%nbuf,this%nvar)::q
integer::nbuf,nx,i

nbuf=this%nbuf
nx=this%nMesh(1)

if(this%boundaryType .eq. 1) then !! zero gradient
  if(this%left_mpi .lt. 0) then
    do i=1,nbuf
       q(i-nbuf,:)=q(nbuf-i+1,:)
    enddo
  endif
  if(this%right_mpi .lt. 0) then
    do i=1,nbuf
       q(nx+i,:)=q(nx-i+1,:)
    enddo
  endif
endif
end subroutine bdrySodShockTube1D

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!! isothermal shocktube problem
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine isoShockTube1D(gridID)

integer::gridID
type(grid)::g1
integer::ndim,nbuf,coordType,variable(8)
integer::nMesh(1),dims(1)
double precision::leftBdry(1),rightBdry(1)
logical::periods(1),reorder
double precision::t0,t1
integer::nstep

nstep=0
variable=0
ndim=1
nbuf=2
coordType=1
variable(1)=1
variable(2)=1

nMesh(1)=400
leftBdry(1)=-0.5d0
rightBdry(1)=0.5d0
dims(1)=nprocs
periods(1)=.false.
reorder=.false.

call g1%setTopologyMPI(ndim,dims,periods,reorder)
call g1%setMesh(nMesh,leftBdry,rightBdry,nbuf,coordType,gridID)
call g1%setVariable(variable) !! den,vx
call g1%setMPIWindows()
call g1%setEoS(eosType=1) !! isothermal
call g1%setSoundSpeed(snd=1.d0) !! isothermal sound speed
call g1%setCFL(CFL=0.6d0) !! Courant number
call g1%setTime(fstart=0,tend=0.2d0,dtout=0.2d0)
call g1%setSlopeLimiter(limiterType=1)
call g1%setSolverType(solverType=1)
call g1%setBoundaryType(boundaryType=1)
call g1%initVariable()
call g1%exchangeBdryMPI(g1%q,g1%winq)
call g1%setBoundary(g1%q)
call g1%writeGrid()
g1%writeFlag=.false.

t0=MPI_WTIME()
do while (g1%t .lt. g1%tend)
   call g1%griddt()
   call g1%evolveGridRK2()
   nstep=nstep+1
   if(myid .eq. 0) then
     print *,nstep," testSuite.f03: gridID=",g1%gridID,'t=',g1%t,"dt=",g1%dt
   endif
   if(g1%writeFlag .eqv. .true.) then
      t1=MPI_WTIME()
      if(myid .eq. 0) then
         print *,"testSuiteMPI.f03: elapse time=",t1-t0,"secs"
      endif
      call g1%writeGrid()
      g1%writeFlag=.false.
      t0=MPI_WTIME()
   endif
enddo !! end do while


end subroutine isoShockTube1D

subroutine initIsoShockTube1D(this,q)
class(grid)::this
double precision,dimension(1-this%nbuf:this%nMesh(1)+this%nbuf,this%nvar)::q
double precision,dimension(1-this%nbuf:this%nMesh(1)+this%nbuf)::xc
integer::nx,nvar,nbuf,gridID,eosType
integer::i

xc=this%xc(1)%coords
nx=this%nMesh(1)
   do i=1,nx
      if(xc(i) .lt. 0.d0) then
         q(i,1)=1.d0
      else
         q(i,1)=0.25d0
      endif
         q(i,2)=0.d0
   enddo
end subroutine initIsoShockTube1D

subroutine bdryIsoShockTube1D(this,q)
class(grid)::this
double precision,dimension(1-this%nbuf:this%nMesh(1)+this%nbuf,this%nvar)::q
integer::nbuf,nx,i

nbuf=this%nbuf
nx=this%nMesh(1)

if(this%boundaryType .eq. 1) then !! zero gradient
  if(this%left_mpi .lt. 0) then
    do i=1,nbuf
       q(i-nbuf,:)=q(nbuf-i+1,:)
    enddo
  endif
  if(this%right_mpi .lt. 0) then
    do i=1,nbuf
       q(nx+i,:)=q(nx-i+1,:)
    enddo
  endif
endif
end subroutine bdryIsoShockTube1D



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!! polytropic shocktube problem
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine polyShockTubeHD1D(gridID)

integer::gridID
type(grid)::g1
integer::ndim,nbuf,coordType,variable(8)
integer::nMesh(1),dims(1)
double precision::leftBdry(1),rightBdry(1)
logical::periods(1),reorder
double precision::t0,t1
integer::nstep

nstep=0
variable=0
ndim=1
nbuf=2
coordType=1
variable(1)=1
variable(2)=1

nMesh(1)=800
leftBdry(1)=-2.d0
rightBdry(1)=2.d0
dims(1)=nprocs
periods(1)=.false.
reorder=.false.
g1%polyK=0.5d0
g1%polyGamma=2.d0

call g1%setTopologyMPI(ndim,dims,periods,reorder)
call g1%setMesh(nMesh,leftBdry,rightBdry,nbuf,coordType,gridID)
call g1%setVariable(variable) !! den,vx
call g1%setMPIWindows()
call g1%setEoS(eosType=3) !! isothermal
!call g1%setSoundSpeed(snd=1.d0) !! isothermal sound speed
call g1%setCFL(CFL=0.6d0) !! Courant number
call g1%setTime(fstart=0,tend=0.5d0,dtout=0.5d0)
call g1%setSlopeLimiter(limiterType=1)
call g1%setSolverType(solverType=6)
call g1%setBoundaryType(boundaryType=1)
call g1%initVariable()
call g1%exchangeBdryMPI(g1%q,g1%winq)
call g1%setBoundary(g1%q)
call g1%writeGrid()
g1%writeFlag=.false.

t0=MPI_WTIME()
do while (g1%t .lt. g1%tend)
   call g1%griddt()
   call g1%evolveGridRK2()
   
   if(myid .eq. 0) then
     print *,nstep," testSuite.f03: gridID=",g1%gridID,'t=',g1%t,"dt=",g1%dt
   endif
   if(g1%writeFlag .eqv. .true.) then
      t1=MPI_WTIME()
      if(myid .eq. 0) then
         print *,"testSuiteMPI.f03: elapse time=",t1-t0,"secs"
      endif
      call g1%writeGrid()
      g1%writeFlag=.false.
      t0=MPI_WTIME()
   endif
   nstep=nstep+1
enddo !! end do while


end subroutine polyShockTubeHD1D

subroutine initPolyShockTubeHD1D(this,q)
class(grid)::this
double precision,dimension(1-this%nbuf:this%nMesh(1)+this%nbuf,this%nvar)::q
double precision,dimension(1-this%nbuf:this%nMesh(1)+this%nbuf)::xc
integer::nx,nvar,nbuf,gridID,eosType
integer::i

xc=this%xc(1)%coords
nx=this%nMesh(1)
   do i=1,nx
      if(xc(i) .lt. 0.d0) then
         q(i,1)=3.d0
      else
         q(i,1)=1.d0
      endif
         q(i,2)=0.d0
   enddo
end subroutine initPolyShockTubeHD1D

subroutine bdryPolyShockTubeHD1D(this,q)
class(grid)::this
double precision,dimension(1-this%nbuf:this%nMesh(1)+this%nbuf,this%nvar)::q
integer::nbuf,nx,i

nbuf=this%nbuf
nx=this%nMesh(1)

if(this%boundaryType .eq. 1) then !! zero gradient
  if(this%left_mpi .lt. 0) then
    do i=1,nbuf
       q(i-nbuf,:)=q(nbuf-i+1,:)
    enddo
  endif
  if(this%right_mpi .lt. 0) then
    do i=1,nbuf
       q(nx+i,:)=q(nx-i+1,:)
    enddo
  endif
endif
end subroutine bdryPolyShockTubeHD1D

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!! C-Shock test 3D
!!!!!! Tilley, Balsara, Meyer, 2012, New Astronomy, 17, 368
!!!!!! Toth, 1994, ApJ, 425, 171
!!!!!! Tilley & Balsara, 2008, MNRAS, 389, 1058
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine CShockTest3D(gridIDn,gridIDi)
implicit none
integer::gridIDn,gridIDi
type(grid)::gn,gi
integer::ndim,nbuf,coordType,nvariable(8),ivariable(8)
integer::nMesh(3),dims(3)
double precision::leftBdry(3),rightBdry(3)
logical::periods(3),reorder
double precision::t0,t1,dt,tend
integer::nstep,ierr

nstep=0
nvariable=0
ivariable=0
ndim=3
nbuf=2
coordType=1
nvariable(1)=1 !!! density
nvariable(2)=1 !!! momx
nvariable(3)=1 !!! momy
nvariable(4)=1 !!! momz
nvariable(8)=1 !!! ene

ivariable=1  !!! all turned on

nMesh(1)=4
nMesh(2)=512
nMesh(3)=512
leftBdry(3) =0.d0   
rightBdry(3)=0.02d0

leftBdry(1) =0.d0
rightBdry(1)=(rightBdry(3)-leftBdry(3))/dble(nMesh(3))*dble(nMesh(1))
leftBdry(2) =0.d0
rightBdry(2)=(rightBdry(3)-leftBdry(3))/dble(nMesh(3))*dble(nMesh(2))

dims(1)=1
dims(2)=1
dims(3)=nprocs
periods(1)=.true.
periods(2)=.true.
periods(3)=.false.
reorder=.true.
tend=0.01636d0
!tend=1.34d-5

call gn%enableAD(enable_ad=.true.)
call gi%enableAD(enable_ad=.true.)
call gn%setADparams(mu_ad=2.33d0,alpha_ad=7.31d4)
call gi%setADparams(mu_ad=30.d0 ,alpha_ad=7.31d4)
call gn%setTopologyMPI(ndim,dims,periods,reorder)
call gi%setTopologyMPI(ndim,dims,periods,reorder)
call gn%setMesh(nMesh,leftBdry,rightBdry,nbuf,coordType,gridIDn)
call gi%setMesh(nMesh,leftBdry,rightBdry,nbuf,coordType,gridIDi)
call gn%setVariable(nvariable)
call gi%setVariable(ivariable)
call gn%setMPIWindows()
call gi%setMPIWIndows()
call gn%setEOS(eosType=2)
call gi%setEOS(eosType=2)
call gn%setadiGamma(gam=1.0001d0)
call gi%setadiGamma(gam=5.d0/3.d0)
call gn%setCFL(CFL=0.6d0)
call gi%setCFL(CFL=0.6d0)
call gn%setTime(fstart=0,tend=tend,dtout=tend)
call gi%setTime(fstart=0,tend=tend,dtout=tend)
call gn%setSlopeLimiter(limiterType=3)
call gi%setSlopeLimiter(limiterType=3)
call gn%setSolverType(solverType=3)
call gi%setSolverType(solverType=4)
call gn%setBoundaryType(boundaryType=0)
call gi%setBoundaryType(boundaryType=0)
call gn%initVariable()
call gi%initVariable()
call gn%exchangeBdryMPI(gn%q,gn%winq)
call gi%exchangeBdryMPI(gi%q,gi%winq)
call gn%setBoundary(gn%q)
call gi%setBoundary(gi%q)
call gn%writeGrid()
call gi%writeGrid()
call gn%writeGrid_HD_vtk()
call gi%writeGrid_MHD_vtk()

gn%writeFlag=.false.
gi%writeFlag=.false.

t0=MPI_WTIME()
do while (gn%t .lt. gn%tend)
   call MPI_BARRIER(MPI_COMM_WORLD,ierr)
   call gn%griddt()
   call gi%griddt()
   call MPI_BARRIER(MPI_COMM_WORLD,ierr)
   dt=dmin1(gn%dt,gi%dt)
   if(dt .gt. gn%toutput-gn%t) then
     dt=gn%toutput-gn%t
     gn%writeFlag=.true.
     gi%writeFlag=.true.
     gn%toutput=gn%toutput+gn%dtout
     gi%toutput=gi%toutput+gi%dtout
     gn%fnum=gn%fnum+1
     gi%fnum=gi%fnum+1
   elseif (dt .gt. gn%tend-gn%t) then
     dt=gn%tend-gn%t
     gn%writeFlag=.true.
     gi%writeFlag=.true.
     gn%fnum=gn%fnum+1
     gi%fnum=gi%fnum+1
   endif
   gn%dt=dt
   gi%dt=dt
   call MPI_BARRIER(MPI_COMM_WORLD,ierr)
   call rk2AD_3D_HSHSMD(gn,gn%q,gn%q1,gn%q2,gi,gi%q,gi%q1,gi%q2)
   call MPI_BARRIER(MPI_COMM_WORLD,ierr)
   call gn%exchangeBdryMPI(gn%q,gn%winq)
   call gi%exchangeBdryMPI(gi%q,gi%winq)
   call gn%setBoundary(gn%q)
   call gi%setBoundary(gi%q)
   if(myid .eq. 0) then
     print *,"testSuite.f03: gridID=",gn%gridID,'t=',gn%t,"dt=",gn%dt
   endif
   if(gn%writeFlag .eqv. .true.) then
     if(gi%writeFlag .eqv. .true.) then
       t1=MPI_WTIME()
       if(myid .eq. 0) then
          print *,"testSuiteMPI.f03: elapse time=",t1-t0,"secs"
       endif
       call gn%writeGrid()
       call gi%writeGrid()
       call gn%writeGrid_HD_vtk()
       call gi%writeGrid_MHD_vtk()
       gn%writeFlag=.false.
       gi%writeFlag=.false.
       t0=MPI_WTIME()
     endif
   endif
enddo !! end do while

end subroutine CShockTest3D

subroutine initCShockTest3Dn(this,q)
implicit none
class(grid)::this
double precision,dimension(1-this%nbuf:this%nMesh(1)+this%nbuf,1-this%nbuf:this%nMesh(2)+this%nbuf, &
                           1-this%nbuf:this%nMesh(3)+this%nbuf,this%nvar)::q
integer::i,j,k,nx,ny,nz
double precision::vx,vy,vz,p,snd,den

nx=this%nMesh(1)
ny=this%nMesh(2)
nz=this%nMesh(3)
den=3.45d0 !! Msun/pc^3
vz=-2.2605d0 !! km/s
snd=0.344d0 !! km/s
p=snd**2*den/this%adiGamma
do k=1,nz
  do j=1,ny
    do i=1,nx
      q(i,j,k,1)=3.45d0
      q(i,j,k,2)=0.d0
      q(i,j,k,3)=0.d0
      q(i,j,k,4)=q(i,j,k,1)*vz
      q(i,j,k,5)=p/(this%adiGamma-1.d0)+0.5d0*den*(vz**2)
    enddo
  enddo 
enddo
end subroutine initCShockTest3Dn

subroutine initCShockTest3Di(this,q)
implicit none
class(grid)::this
double precision,dimension(1-this%nbuf:this%nMesh(1)+this%nbuf,1-this%nbuf:this%nMesh(2)+this%nbuf, &
                           1-this%nbuf:this%nMesh(3)+this%nbuf,this%nvar)::q
integer::i,j,k,nx,ny,nz
double precision::vx,vy,vz,p,snd,den,bxl,byl,bzl,bxr,byr,bzr
double precision::pi
double precision::bxc,byc,bzc
pi=4.d0*datan2(1.d0,1.d0)
nx=this%nMesh(1)
ny=this%nMesh(2)
nz=this%nMesh(3)
den=7.4d-3 !! Msun/pc^3
vx=0.d0
vy=0.d0
vz=-2.2605d0
bxl=-1.71d0
bzl=0.d0
byl=0.d0
bxr=-1.71d0
bzr=0.d0
byr=0.d0
snd=9.6d-2 !! km/s
p=snd**2*den/this%adiGamma
bxc=0.5d0*(bxl+bxr)
byc=0.5d0*(byl+byr)
bzc=0.5d0*(bzl+bzr)
do k=1,nz
  do j=1,ny
    do i=1,nx
      q(i,j,k,1)=den
      q(i,j,k,2)=den*vx
      q(i,j,k,3)=den*vy
      q(i,j,k,4)=den*vz
      q(i,j,k,5)=bxl
      q(i,j,k,6)=byl
      q(i,j,k,7)=bzl
      q(i,j,k,8)=p/(this%adiGamma-1.d0)+0.5d0*den*(vx**2+vy**2+vz**2)+0.5d0*(bxc**2+byc**2+bzc**2)
      q(i,j,k,9)=bxr
      q(i,j,k,10)=byr
      q(i,j,k,11)=bzr
    enddo
  enddo 
enddo
end subroutine initCShockTest3Di

subroutine bdryCShockTest3Dn(this,q)
implicit none
class(grid)::this
double precision,dimension(1-this%nbuf:this%nMesh(1)+this%nbuf,1-this%nbuf:this%nMesh(2)+this%nbuf, &
                           1-this%nbuf:this%nMesh(3)+this%nbuf,this%nvar)::q
integer::i,j,k,nx,ny,nz,nbuf
double precision::vx,vy,vz,p,snd,den

nx=this%nMesh(1)
ny=this%nMesh(2)
nz=this%nMesh(3)
nbuf=this%nbuf
den=3.45d0
!vx=-2.2605d0
!vy=-2.2605d0
vz=-2.2605d0
snd=0.344d0
p=snd**2*den/this%adiGamma

!if(this%left_mpi .lt. 0) then
!  do i=1,nbuf
!    q(i-nbuf,:,:,:)=q(nbuf-i+1,:,:,:)
!  enddo
!
!  do i=1,nbuf
!    q(i-nbuf,:,:,2)=-q(nbuf-i+1,:,:,2)
!  enddo
!endif

!if(this%right_mpi .lt. 0) then
!  do i=1,nbuf
!    q(nx+i,:,:,1)=den
!    q(nx+i,:,:,2)=den*vx
!    q(nx+i,:,:,3)=0.d0
!    q(nx+i,:,:,4)=0.d0
!    q(nx+i,:,:,5)=p/(this%adiGamma-1.d0)+0.5d0*den*(vx**2)
!  enddo
!endif

!if(this%up_mpi .lt. 0) then
!   do j=1,nbuf
!     q(:,j-nbuf,:,:)=q(:,nbuf-j+1,:,:)
!   enddo
!
!   do j=1,nbuf
!     q(:,j-nbuf,:,3)=-q(:,nbuf-j+1,:,3)
!   enddo
!endif

!if(this%down_mpi .lt. 0) then
!   do j=1,nbuf
!     q(:,ny+j,:,1)=den
!     q(:,ny+j,:,2)=0.d0
!     q(:,ny+j,:,3)=den*vy
!     q(:,ny+j,:,4)=0.d0
!     q(:,ny+j,:,5)=p/(this%adiGamma-1.d0)+0.5d0*den*(vy**2)
!   enddo
!endif

if(this%top_mpi .lt. 0) then
   do k=1,nbuf
     q(:,:,k-nbuf,:)=q(:,:,nbuf-k+1,:)
   enddo
   
   do k=1,nbuf
     q(:,:,k-nbuf,4)=-q(:,:,nbuf-k+1,4)
   enddo
endif

if(this%bottom_mpi .lt. 0) then
   do k=1,nbuf
     q(:,:,nz+k,1)=den
     q(:,:,nz+k,2)=0.d0
     q(:,:,nz+k,3)=0.d0
     q(:,:,nz+k,4)=den*vz
     q(:,:,nz+k,5)=p/(this%adiGamma-1.d0)+0.5*den*(vz**2)
   enddo
endif
end subroutine bdryCShockTest3Dn

subroutine bdryCShockTest3Di(this,q)
implicit none
class(grid)::this
double precision,dimension(1-this%nbuf:this%nMesh(1)+this%nbuf,1-this%nbuf:this%nMesh(2)+this%nbuf, &
                           1-this%nbuf:this%nMesh(3)+this%nbuf,this%nvar)::q
integer::i,j,k,nx,ny,nz,nbuf
double precision::vx,vy,vz,den,bxl,byl,bzl,bxr,byr,bzr,bxc,byc,bzc,snd,p,pi

pi=datan2(1.d0,1.d0)*4.d0
nx=this%nMesh(1)
ny=this%nMesh(2)
nz=this%nMesh(3)
nbuf=this%nbuf
den= 7.4d-3
vx = 0.d0
vy = 0.d0
vz =-2.2605d0

bxl=-1.71d0
bzl=0.d0
byl=0.d0

bxr=-1.71d0
bzr=0.d0
byr=0.d0

snd=9.6d-2
p=snd**2*den/this%adiGamma
bxc=0.5d0*(bxl+bxr)
byc=0.5d0*(byl+byr)
bzc=0.5d0*(bzl+bzr)

!if(this%left_mpi .lt. 0) then
!  do i=1,nbuf
!    q(i-nbuf,:,:,:)=q(nbuf-i+1,:,:,:)
!  enddo
!  do i=1,nbuf
!    q(i-nbuf,:,:,2)=-q(nbuf-i+1,:,:,2)
!  enddo
!endif

!if(this%right_mpi .lt. 0) then
!      do i=1,nbuf
!        q(nx+i,:,:,1)=den
!        q(nx+i,:,:,2)=den*vx
!        q(nx+i,:,:,3)=den*vy
!        q(nx+i,:,:,4)=den*vz
!        q(nx+i,:,:,5)=bxl
!        q(nx+i,:,:,6)=byl
!        q(nx+i,:,:,7)=bzl
!        q(nx+i,:,:,8)=p/(this%adiGamma-1.d0)+0.5d0*den*(vx**2+vy**2+vz**2)+0.5d0*(bxc**2+byc**2+bzc**2)
!        q(nx+i,:,:,9)=bxr
!        q(nx+i,:,:,10)=byr
!        q(nx+i,:,:,11)=bzr
!      enddo
!endif


!if(this%up_mpi .lt. 0) then
!   do j=1,nbuf
!     q(:,j-nbuf,:,:) =  q(:,nbuf-j+1,:,:)
!   enddo
!   do j=1,nbuf
!     q(:,j-nbuf,:,3) = -q(:,nbuf-j+1,:,3)
!   enddo
!endif

!if(this%down_mpi .lt. 0) then
!   do j=1,nbuf
!     q(:,ny+j,:,1) = den
!     q(:,ny+j,:,2) = den*vx
!     q(:,ny+j,:,3) = den*vy
!     q(:,ny+j,:,4) = den*vz
!     q(:,ny+j,:,5) = bxl
!     q(:,ny+j,:,6) = byl
!     q(:,ny+j,:,7) = bzl
!     q(:,ny+j,:,8) = p/(this%adiGamma-1.d0)+0.5d0*den*(vx**2+vy**2+vz**2)+0.5d0*(bxc**2+byc**2+bzc**2) 
!     q(:,ny+j,:,9) = bxr
!     q(:,ny+j,:,10) = byr
!     q(:,ny+j,:,11) = bzr    
!   enddo
!endif

if(this%top_mpi .lt. 0) then
   do k=1,nbuf
     q(:,:,k-nbuf,:)=q(:,:,nbuf-k+1,:)
   enddo
   do k=1,nbuf
     q(:,:,k-nbuf,4)=-q(:,:,nbuf-k+1,4)
   enddo
endif

if(this%bottom_mpi .lt. 0) then
   do k=1,nbuf
     q(:,:,nz+k,1)=den
     q(:,:,nz+k,2)=den*vx
     q(:,:,nz+k,3)=den*vy
     q(:,:,nz+k,4)=den*vz
     q(:,:,nz+k,5)=bxl
     q(:,:,nz+k,6)=byl
     q(:,:,nz+k,7)=bzl
     q(:,:,nz+k,8)=p/(this%adiGamma-1.d0)+0.5d0*den*(vx**2+vy**2+vz**2)+0.5d0*(bxc**2+byc**2+bzc**2)
     q(:,:,nz+k,9)=bxr
     q(:,:,nz+k,10)=byr
     q(:,:,nz+k,11)=bzr
   enddo
endif

end subroutine bdryCShockTest3Di

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!! Two-fluid Spherical Blast Wave
!!!!! Ref. Tilly & Balsara & Meyer, New Astronomy, 2012, 17, 368
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine sphericalBlastWaveAD2D(gridIDn,gridIDi)
implicit none

integer::gridIDn,gridIDi

type(grid)::gn,gi

integer::ndim,nbuf,coordType,nvariable(8),ivariable(8)
integer::nMesh(2),dims(2)
double precision::leftBdry(2),rightBdry(2)
logical::periods(2),reorder
double precision::t0,t1,dt,tend
integer::nstep,ierr

nstep=0
nvariable=0
ivariable=0
ndim=2
nbuf=2
coordType=1
nvariable(1)=1  ! density
nvariable(2)=1  ! x-momentum
nvariable(3)=1  ! y-momentum
nvariable(8)=1  ! energy 

ivariable=1

nMesh(1)=400
nMesh(2)=400
leftBdry(1)=0.d0
leftBdry(2)=0.d0
rightBdry(1)=1.d0
rightBdry(2)=1.d0
dims=(/2,2/)
call MPI_DIMS_CREATE(nprocs,ndim,dims,ierr)
periods(1)=.true.
periods(2)=.true.
reorder=.true.
tend = 0.006d0

call gn%enableAD(enable_ad=.true.)
call gi%enableAD(enable_ad=.true.)
call gn%setADparams(mu_ad=2.3d0,alpha_ad=3.5469d2)
call gi%setADParams(mu_ad=29.d0,alpha_ad=3.5469d2)
call gn%setTopologyMPI(ndim,dims,periods,reorder)
call gi%setTopologyMPI(ndim,dims,periods,reorder)

call gn%setMesh(nMesh,leftBdry,rightBdry,nbuf,coordType,gridIDn)
call gi%setMesh(nMesh,leftBdry,rightBdry,nbuf,coordType,gridIDi)
call gn%setVariable(nvariable)
call gi%setVariable(ivariable)
call gn%setMPIWindows()
call gi%setMPIWindows()
call gn%setEOS(eosType=2)
call gi%setEOS(eosType=2)
call gn%setadiGamma(gam=1.4d0)
call gi%setadiGamma(gam=5.d0/3.d0)
call gn%setCFL(CFL=0.2d0)
call gi%setCFL(CFL=0.2d0)
call gn%setTime(fstart=0,tend=tend,dtout=tend/1.d0)
call gi%setTime(fstart=0,tend=tend,dtout=tend/1.d0)
call gn%setSlopeLimiter(limiterType=1)
call gi%setSlopeLimiter(limiterType=1)
call gn%setSolverType(solverType=3)
call gi%setSolverType(solverType=4)
call gn%setBoundaryType(boundaryType=0)
call gi%setBoundaryType(boundaryType=0)
call gn%initVariable()
call gi%initVariable()
call gn%exchangeBdryMPI(gn%q,gn%winq)
call gi%exchangeBdryMPI(gi%q,gi%winq)
call gn%setBoundary(gn%q)
call gi%setBoundary(gi%q)
call gn%writeGrid()
call gi%writeGrid()
call gn%writeGrid_HD_vtk()
call gi%writeGrid_MHD_vtk()

t0=MPI_WTIME()
do while(gn%t .lt. gn%tend)
   call MPI_BARRIER(MPI_COMM_WORLD,ierr)
   call gn%griddt()
   call gi%griddt()
   call MPI_BARRIER(MPI_COMM_WORLD,ierr)
   dt=dmin1(gn%dt,gi%dt)
   if(dt .gt. gn%toutput-gn%t) then
     dt=gn%toutput-gn%t
     gn%writeFlag=.true.
     gi%writeFlag=.true.
     gn%toutput=gn%toutput+gn%dtout
     gi%toutput=gi%toutput+gi%dtout
     gn%fnum=gn%fnum+1
     gi%fnum=gi%fnum+1
   elseif(dt .gt. gn%tend-gn%t) then
     dt=gn%tend-gn%t
     gn%writeFlag=.true.
     gi%writeFlag=.true.
     gn%fnum=gn%fnum+1
     gi%fnum=gi%fnum+1
   endif   
   gn%dt=dt
   gi%dt=dt
   call MPI_BARRIER(MPI_COMM_WORLD,ierr)
   call rk2AD_2D_HSHSMD(gn,gn%q,gn%q1,gn%q2,gi,gi%q,gi%q1,gi%q2)
   call MPI_BARRIER(MPI_COMM_WORLD,ierr)
   call gn%exchangeBdryMPI(gn%q,gn%winq)
   call gi%exchangeBdryMPI(gi%q,gi%winq)
   call gn%setBoundary(gn%q)
   call gi%setBoundary(gi%q)
   
  if(myid .eq. 0) then
     print *,"testSuite.f03: gridID=",gn%gridID,'t=',gn%t,"dt=",gn%dt
   endif
   if(gn%writeFlag .eqv. .true.) then
     if(gi%writeFlag .eqv. .true.) then
       t1=MPI_WTIME()
       if(myid .eq. 0) then
          print *,"testSuiteMPI.f03: elapse time=",t1-t0,"secs"
       endif
       call gn%writeGrid()
       call gi%writeGrid()
       call gn%writeGrid_HD_vtk()
       call gi%writeGrid_MHD_vtk()
       gn%writeFlag=.false.
       gi%writeFlag=.false.
       t0=MPI_WTIME()
     endif
   endif

enddo !!!! end do while


end subroutine sphericalBlastWaveAD2D

subroutine initSphericalBlastWaveAD2Dn(this,q)
implicit none
class(grid)::this
double precision,dimension(1-this%nbuf:this%nMesh(1)+this%nbuf,1-this%nbuf:this%nMesh(2)+this%nbuf,this%nvar)::q
integer::nx,ny,i,j
double precision::xc,yc,rc,x,y,r,deni,deno,pi,po

nx=this%nMesh(1)
ny=this%nMesh(2)
xc=0.5d0
yc=0.5d0
rc=0.1d0
deni=10.d0
deno = 1.d0
pi = 4347.8d0
po = 4.3478d0

do j=1,ny
  do i=1,nx
    x=this%xc(1)%coords(i)
    y=this%xc(2)%coords(j)
    r=dsqrt((x-xc)**2+(y-yc)**2)
    if(r .le. rc) then
      q(i,j,1)=deni 
      q(i,j,4)=pi/(this%adiGamma-1.d0)
    else
      q(i,j,1)=deno
      q(i,j,4)=po/(this%adiGamma-1.d0)
    endif
    q(i,j,2)=0.d0
    q(i,j,3)=0.d0
  enddo
enddo
end subroutine initSphericalBlastWaveAD2Dn

subroutine initSphericalBlastWaveAD2Di(this,q)
implicit none
class(grid)::this
double precision,dimension(1-this%nbuf:this%nMesh(1)+this%nbuf,1-this%nbuf:this%nMesh(2)+this%nbuf,this%nvar)::q
integer::nx,ny,i,j
double precision::xc,yc,rc,x,y,r,deni,deno,pi,po,By0

pi = 4.d0*datan2(1.d0,1.d0)
nx=this%nMesh(1)
ny=this%nMesh(2)
xc=0.5d0
yc=0.5d0
rc=0.1d0
deni=0.0126d0
deno=0.00126d0
pi=0.43478d0
po=4.3478d-4
By0=dsqrt(8.d0*pi)/dsqrt(4.d0*pi)

do j=1,ny
  do i=1,nx
    x=this%xc(1)%coords(i)
    y=this%xc(2)%coords(j)
    r=dsqrt((x-xc)**2+(y-yc)**2)
    if(r .le. rc) then
      q(i,j,1)=deni
      q(i,j,8)=pi/(this%adiGamma-1.d0)+0.5d0*(By0**2)
    else
      q(i,j,1)=deno
      q(i,j,8)=po/(this%adiGamma-1.d0)+0.5d0*(By0**2)
    endif
    q(i,j,2)=0.d0  !!! x-momentum
    q(i,j,3)=0.d0  !!! y-momentum
    q(i,j,4)=0.d0  !!! z-momentum
    q(i,j,5)=0.d0  !!! bxl
    q(i,j,6)=By0   !!! byl
    q(i,j,7)=0.d0  !!! bzl
    q(i,j,9)=0.d0  !!! bxr
    q(i,j,10)=By0  !!! byr
    q(i,j,11)=0.d0 !!! bzr
  enddo
enddo
end subroutine initSphericalBlastWaveAD2Di

subroutine bdrySphericalBlastWaveAD2Dn(this,q)
implicit none
class(grid)::this
double precision,dimension(1-this%nbuf:this%nMesh(1)+this%nbuf,1-this%nbuf:this%nMesh(2)+this%nbuf,this%nvar)::q

!!!!! empty due to periodic boundary condition !!!!!
end subroutine bdrySphericalBlastWaveAD2Dn

subroutine bdrySphericalBlastWaveAD2Di(this,q)
implicit none
class(grid)::this
double precision,dimension(1-this%nbuf:this%nMesh(1)+this%nbuf,1-this%nbuf:this%nMesh(2)+this%nbuf,this%nvar)::q

!!!! empty due to periodic boundary condition !!!!!

end subroutine bdrySphericalBlastWaveAD2Di



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!! Testing 2D Driving turbulence module for 2D HD cases added by WGZeng
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine TestDrivingTurbulence2DHD(gridID)
    integer::gridID
    type(grid)::g1
    integer::ndim,nbuf,coordType,variable(8)
    integer::nMesh(2),dims(2)
    double precision::leftBdry(2),rightBdry(2)
    double precision::pi
    logical::periods(2),reorder
    double precision::t0,t1
    double precision:: dt_turb,t_accum_turb !!!zwg:: added for turbulence driving!!!!!
    integer::t_count_turb,DT_mode   !!!zwg:: added for turbulence driving!!!!!
    integer::nstep,ier

    nstep=0
    pi=4.d0*atan2(1.d0,1.d0)
    variable=0
    ndim=2
    nbuf=2
    coordType=1 !! Cartesian

    dt_turb=3.0E11	!0.d0	!5.0E10	!5e10s
    t_count_turb=1
    t_accum_turb=dt_turb

    variable(1)=1  !! rho
    variable(2)=1  !! vx
    variable(3)=1  !! vy
    variable(8)=1  !! energy, also required for isothermal simulation

    nMesh(1)=128
    nMesh(2)=128

    leftBdry(1)=0.d0
    rightBdry(1)=3.0857E17
    leftBdry(2)=0.d0
    rightBdry(2)=3.0857E17


    !!!! zwg: if use FFTW, then we should divide the domain along the last coordinate
    dims(1)=1  
    dims(2)=nprocs

    periods(1)=.true.
    periods(2)=.true.
    reorder=.true.

    call g1%setGridID(gridID=gridID)	
    call g1%enableDrivingTurbulence(DT_mode=0)
!!!!!!!!!!!enableSelfgravity!!!!!!!!!!
!!call g1%enableSelfgravity() 
!!call g1%setSgBdryType(sgBdryType=1) !! 0:isolated, 1:periodic
!!call g1%setGravConst(GravConst=4.3011d-3)
    call g1%setTopologyMPI(ndim,dims,periods,reorder)
    call g1%setMesh(nMesh,leftBdry,rightBdry,nbuf,coordType,gridID)
    call g1%setVariable(variable) 
!!!!!!!!!!! set MPI communication for selfgravity !!!!!
!!call g1%setSelfgravityKernel()
!!call g1%setSGMPIWindows()

    call g1%setMPIWindows()
    call g1%setEoS(eosType=2) !! adiabatic
    !call g1%setEoS(eosType=1) !! isothermal
    call g1%setadiGamma(gam=1.00001d0) !! ratio of heat capacity
    !call g1%setSoundSpeed(snd=186.d0) !! isothermal sound speed
    call g1%setCFL(CFL=0.4d0) !! courant number
    call g1%setTime(fstart=0,tend=3000000000000000.d0,dtout=1000000000000.d0)
    call g1%setSlopeLimiter(limiterType=3)
    call g1%setSolverType(solverType=3) !! 3=AdiHLLCHD
    call g1%setBoundaryType(boundaryType=3)
    call g1%initVariable()

    if(g1%enable_DT)then
      call g1%setDTenergyfaction(t_accum_turb, dt_turb)
      call g1%calcDrivingTurbulence_MD(g1%q)
    endif

    call g1%exchangeBdryMPI(g1%q,g1%winq)
    call g1%setBoundary(g1%q)

!!!!! calculate selfgravity for initial condition !!!!!
!!call g1%calcSelfgravity(g1%q) !!! here only neutral will be calculated
    call g1%writeGrid()
    call g1%writeGrid_HD_vtk()
    g1%writeFlag=.false.
    t0=MPI_WTIME()

    t_accum_turb=0.0
do while (g1%t .lt. g1%tend)

   call g1%griddt()

   if(g1%DT_mode==1)then  !!!!!!! zwg: added for continue driving turbulence!!!!!
     if((t_accum_turb .le. dt_turb) .and. (t_accum_turb+g1%dt .gt. dt_turb))  then
       call g1%setDTenergyfaction(t_accum_turb, dt_turb)
       call g1%calcDrivingTurbulence_MD(g1%q)
       call g1%exchangeBdryMPI(g1%q,g1%winq)
       call g1%setBoundary(g1%q)
       t_accum_turb=0
       t_count_turb=t_count_turb+1
     endif
     t_accum_turb=t_accum_turb+g1%dt
   endif !!!! end of if(g1%DT_mode==1)then

   call g1%evolveGridRK2()

   if(myid .eq. 0) then
     print *,"testSuite.f03: gridID=",g1%gridID,'t=',g1%t,"dt=",g1%dt
   endif
   if(g1%writeFlag .eqv. .true.) then
      t1=MPI_WTIME()
      if(myid .eq. 0) then
        print *,"testSuiteMPI.f03: elapse time=",t1-t0,"secs"
      endif
      call g1%writeGrid()
      call g1%writeGrid_HD_vtk()
      g1%writeFlag=.false.
      t0=MPI_WTIME()
   endif
enddo !! end do while

end subroutine TestDrivingTurbulence2DHD




subroutine initTurbulenceDriving2DHD(this,q)
use gridModule
class(grid)::this
integer:: Nx, Ny,i,j
double precision::rho0,u0,v0,p0,snd

double precision,dimension(1-this%nbuf:this%nMesh(1)+this%nbuf,1-this%nbuf:this%nMesh(2)+this%nbuf,this%nvar)::q
!!double precision,dimension(1:this%nMesh(1),1:this%nMesh(2))::Den_ini

Nx=this%nMesh(1)
Ny=this%nMesh(2)

if(this%enable_DT)then  !!! set DT parameters here!!!
  this%drivingWN_DT=2.0 ! Driving scale, i.e nx=128, scale = 2 -> 64 scale.
  this%Energy_DT=3.437068E6
  this%zeta_DT=1.0
  this%netmomx_DT=0.0    
  this%netmomy_DT=0.0
  !this%netmomz_DT=0.0  !!zwg for 3D case
endif

rho0=0.0001
u0=0.0
v0=0.0
snd=186.d0  !!!m/s
p0=snd**2*rho0/this%adiGamma

do i=1,Nx
  do j=1, Ny
    q(i,j,1)=rho0
    q(i,j,2)=rho0*u0
    q(i,j,3)=rho0*v0
    q(i,j,4)=p0/(this%adiGamma-1.d0)+0.5*rho0*(u0**2+v0**2)
   !q(i,j,4)=p_tmp/(this%adiGamma-1.d0)+0.5*rho_tmp*(u_tmp**2+v_tmp**2)
  enddo
enddo



end subroutine initTurbulenceDriving2DHD

subroutine bdryTurbulenceDriving2DHD(this,q)
class(grid)::this
double precision,dimension(1-this%nbuf:this%nMesh(1)+this%nbuf,1-this%nbuf:this%nMesh(2)+this%nbuf,this%nvar)::q
integer::nx,ny,nbuf,i,j

nbuf=this%nbuf
nx=this%nMesh(1)
ny=this%nMesh(2)

if(this%boundaryType .eq. 3) then !! periodic boundary condition
!!! do nothing for periodic boundary conditions
endif

end subroutine bdryTurbulenceDriving2DHD




!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!! Testing 3D Driving turbulence module for 3D HD cases added by WGZeng
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 
  subroutine TestDrivingTurbulence3DHD(gridID)
    integer::gridID
    type(grid)::g1
    integer::ndim,nbuf,coordType,variable(8)
    integer::nMesh(3),dims(3)
    double precision::leftBdry(3),rightBdry(3)
    double precision::pi
    logical::periods(3),reorder
    double precision::t0,t1
    double precision:: dt_turb,t_accum_turb !!!zwg:: added for turbulence driving!!!!!
    integer::t_count_turb,DT_mode  !!!zwg:: added for turbulence driving!!!!!
    integer::nstep,ier

    nstep=0
    pi=4.d0*atan2(1.d0,1.d0)
    variable=0
    ndim=3
    nbuf=2
    coordType=1 !! Cartesian

    !!!zwg:: added for turbulence driving!!!!!
    dt_turb=3.0E11	!0.d0	!5.0E10	!5e10s
    t_count_turb=1
    t_accum_turb=dt_turb

    variable(1)=1  !! rho
    variable(2)=1  !! vx
    variable(3)=1  !! vy
    variable(4)=1  !! vz
    variable(8)=1  !! energy, also required for isothermal simulation

    nMesh(1)=128
    nMesh(2)=128
    nMesh(3)=128

    leftBdry(1)=0.d0
    rightBdry(1)=3.0857E17
    leftBdry(2)=0.d0
    rightBdry(2)=3.0857E17
    leftBdry(3)=0.d0
    rightBdry(3)=3.0857E17


    !!!! zwg: if use FFTW, then we should divide the domain along the last coordinate
    dims(1)=1  
    dims(2)=1
    dims(3)=nprocs

    periods(1)=.true.
    periods(2)=.true.
    periods(2)=.true.
    reorder=.true.


    call g1%setGridID(gridID=gridID)	
    call g1%enableDrivingTurbulence(DT_mode=0)
!!!!!!!!!!!enableSelfgravity!!!!!!!!!!
!!call g1%enableSelfgravity() 
!!call g1%setSgBdryType(sgBdryType=1) !! 0:isolated, 1:periodic
!!call g1%setGravConst(GravConst=4.3011d-3)
    call g1%setTopologyMPI(ndim,dims,periods,reorder)
    call g1%setMesh(nMesh,leftBdry,rightBdry,nbuf,coordType,gridID)
    call g1%setVariable(variable) 
!!!!!!!!!!! set MPI communication for selfgravity !!!!!
!!call g1%setSelfgravityKernel()
!!call g1%setSGMPIWindows()

    call g1%setMPIWindows()
    call g1%setEoS(eosType=2) !! adiabatic
    !call g1%setEoS(eosType=1) !! isothermal
    call g1%setadiGamma(gam=1.00001d0) !! ratio of heat capacity
    !call g1%setSoundSpeed(snd=186.d0) !! isothermal sound speed
    call g1%setCFL(CFL=0.1d0) !! courant number
    call g1%setTime(fstart=0,tend=3000000000000000.d0,dtout=1000000000000.d0)

    !call g1%setRestart(fstart=5,tend=3000000000000000.d0,dtout=1000000000000.d0)
                                     !3000000000000000 1000000000000
    call g1%setSlopeLimiter(limiterType=3)
    call g1%setSolverType(solverType=3) !! 3=AdiHLLCHD
    call g1%setBoundaryType(boundaryType=3)


    !call g1%initVariable()
if(g1%isRestart==0)then
  call g1%initVariable()
elseif(g1%isRestart==1)then
 call g1%readGrid_HD_vtk()
  !call g1%readGrid_MHD_vtk()
endif

if(g1%isRestart==0)then
    if(g1%enable_DT)then
      call g1%setDTenergyfaction(t_accum_turb, dt_turb)
      call g1%calcDrivingTurbulence_MD(g1%q)
    endif
endif


    call g1%exchangeBdryMPI(g1%q,g1%winq)
    call g1%setBoundary(g1%q)

!!!!! calculate selfgravity for initial condition !!!!!
!!call g1%calcSelfgravity(g1%q) !!! here only neutral will be calculated

if(g1%isRestart==0)then
    call g1%writeGrid()
    call g1%writeGrid_HD_vtk()
endif


    g1%writeFlag=.false.
    t0=MPI_WTIME()

    t_accum_turb=0.0  !!! set t_accum_turb=0.0 if g1%isRestart==0

do while (g1%t .lt. g1%tend)
   nstep=nstep+1
   call g1%griddt()

   if(g1%DT_mode==1)then  !!!!!!! zwg: added for continue driving turbulence!!!!!
     if((t_accum_turb .le. dt_turb) .and. (t_accum_turb+g1%dt .gt. dt_turb))  then
       call g1%setDTenergyfaction(t_accum_turb, dt_turb)
       call g1%calcDrivingTurbulence_MD(g1%q)
       call g1%exchangeBdryMPI(g1%q,g1%winq)
       call g1%setBoundary(g1%q)
       t_accum_turb=0
       t_count_turb=t_count_turb+1
     endif
     t_accum_turb=t_accum_turb+g1%dt
   endif !!!! end of if(g1%DT_mode==1)then

   call g1%evolveGridRK2()

   if(myid .eq. 0) then
     print *,"testSuite.f03: gridID=",g1%gridID,'t=',g1%t,"dt=",g1%dt
   endif
   if(g1%writeFlag .eqv. .true.) then
      t1=MPI_WTIME()
      if(myid .eq. 0) then
        print *,"testSuiteMPI.f03: elapse time=",t1-t0,"secs"
      endif
      call g1%writeGrid()
      call g1%writeGrid_HD_vtk()
      g1%writeFlag=.false.
      t0=MPI_WTIME()
   endif
enddo !! end do while

end subroutine TestDrivingTurbulence3DHD
               

          

subroutine initTurbulenceDriving3DHD(this,q)
use gridModule
class(grid)::this
integer:: Nx,Ny,Nz,i,j,k
double precision::rho0,u0,v0,w0,p0,snd

!double precision,dimension(1-this%nbuf:this%nMesh(1)+this%nbuf,1-this%nbuf:this%nMesh(2)+this%nbuf,this%nvar)::q
double precision,dimension(1-this%nbuf:this%nMesh(1)+this%nbuf,1-this%nbuf:this%nMesh(2)+this%nbuf, &
                         & 1-this%nbuf:this%nMesh(3)+this%nbuf,this%nvar)::q

!!double precision,dimension(1:this%nMesh(1),1:this%nMesh(2))::Den_ini

Nx=this%nMesh(1)
Ny=this%nMesh(2)
Nz=this%nMesh(3)

if(this%enable_DT)then  !!! set DT parameters here!!!
  this%drivingWN_DT=8.0
  this%Energy_DT=10160.64*0.5
  this%zeta_DT=1.0
  this%netmomx_DT=0.0    
  this%netmomy_DT=0.0
  this%netmomz_DT=0.0  !!zwg for 3D case
endif

rho0=0.0001
u0=0.0
v0=0.0
w0=0.0
snd=186.d0  !!!m/s
p0=snd**2*rho0/this%adiGamma

  do i=1,Nx
    do j=1,Ny
      do k=1,Nz
      q(i,j,k,1)=rho0
      q(i,j,k,2)=rho0*u0
      q(i,j,k,3)=rho0*v0
      q(i,j,k,4)=rho0*w0
      q(i,j,k,5)=p0/(this%adiGamma-1.d0)+0.5*rho0*(u0**2+v0**2+w0**2)
      enddo
    enddo
  enddo



end subroutine initTurbulenceDriving3DHD

          !bdryTurbulenceDriving3DHD
subroutine bdryTurbulenceDriving3DHD(this,q)
class(grid)::this
double precision,dimension(1-this%nbuf:this%nMesh(1)+this%nbuf,1-this%nbuf:this%nMesh(2)+this%nbuf, &
                         & 1-this%nbuf:this%nMesh(3)+this%nbuf,this%nvar)::q
integer::nx,ny,nz,nbuf,i,j,k

nbuf=this%nbuf
nx=this%nMesh(1)
ny=this%nMesh(2)
nz=this%nMesh(3)

if(this%boundaryType .eq. 3) then !! periodic boundary condition
!!! do nothing for periodic boundary conditions
endif

end subroutine bdryTurbulenceDriving3DHD
              !bdryTurbulenceDriving3DHD







!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!! Testing 3D Driving turbulence module for 3D MHD cases added by WGZeng
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 
  subroutine TestDrivingTurbulence3DMHD(gridID)
    integer::gridID
    type(grid)::g1
    integer::ndim,nbuf,coordType,variable(8)
    integer::nMesh(3),dims(3)
    double precision::leftBdry(3),rightBdry(3)
    double precision::pi
    logical::periods(3),reorder
    double precision::t0,t1
    integer::nstep,ier

    nstep=0
    pi=4.d0*atan2(1.d0,1.d0)
    variable=0
    ndim=3
    nbuf=2
    coordType=1 !! Cartesian

    variable(1)=1  !! rho
    variable(2)=1  !! vx
    variable(3)=1  !! vy
    variable(4)=1  !! vz
    variable(5)=1  !! bx
    variable(6)=1  !! by
    variable(7)=1  !! bz
    variable(8)=1  !! energy, also required for isothermal simulation

    nMesh(1)=4
    nMesh(2)=2
    nMesh(3)=4

    leftBdry(1)=0.d0
    rightBdry(1)=2.d0*pi
    leftBdry(2)=0.d0
    rightBdry(2)=2.d0*pi
    leftBdry(3)=0.d0
    rightBdry(3)=2.d0*pi


    !!!! zwg: if use FFTW, then we should divide the domain along the last coordinate
    dims(1)=1  
    dims(2)=1
    dims(3)=nprocs

    periods(1)=.true.
    periods(2)=.true.
    periods(2)=.true.
    reorder=.true.

    !call g1%enableDrivingTurbulence()
    !write(*,*)"zwg :I am here!!!"
    !stop
    call g1%setTopologyMPI(ndim,dims,periods,reorder)
    call g1%setMesh(nMesh,leftBdry,rightBdry,nbuf,coordType,gridID)

    call g1%setVariable(variable) !! den,vx,vy,vz,bx,by,bz,ene

    !call g1%initVariable()
    call g1%init_For_FFTW()

    call g1%calcDrivingTurbulence(g1%q)

    


    !!call g1%enableSelfgravity()
    !!call g1%setSgBdryType(sgBdryType=1) ! 0=isolated, 1=periodic
    !!call g1%setGravConst(GravConst=1.d0)




  end subroutine TestDrivingTurbulence3DMHD



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!
!!!  2DHD Richtmyer-Meshkov Instability
!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine HDRichtmyerMeshkovInstability2D(gridID)

integer::gridID
type(grid)::g1
integer::ndim,nbuf,coordType,variable(8)
integer::nMesh(2),dims(2)
double precision::leftBdry(2),rightBdry(2)
logical::periods(2),reorder
double precision::t0,t1
integer::nstep,ierr

nstep=0
variable=0
ndim=2
nbuf=2
coordType=1
variable(1)=1
variable(2)=1
variable(3)=1
variable(8)=1

nMesh(1)=3072
nMesh(2)=256
leftBdry(1)=-2.0d0*3.1415926d0                !-0.05d0
leftBdry(2)=0.0d0                             !-0.5d0
rightBdry(1)=10.0d0*3.1415926d0                !0.05d0
rightBdry(2)=1.0d0*3.1415926d0                !0.5d0
dims=(/0,0/)
call MPI_DIMS_CREATE(nprocs,ndim,dims,ierr)
dims(1)=nprocs
dims(2)=1
periods(1)=.false.
periods(2)=.true.
reorder=.true.


call g1%setTopologyMPI(ndim,dims,periods,reorder)
call g1%setMesh(nMesh,leftBdry,rightBdry,nbuf,coordType,gridID)
call g1%setVariable(variable) !! den,vx
call g1%setMPIWindows()
call g1%setEoS(eosType=2) !! adiabatic
call g1%setadiGamma(gam=1.4d0) 
call g1%setCFL(CFL=0.8d0) !! Courant number
call g1%setTime(fstart=0,tend=50.0d0,dtout=1.0d0)
call g1%setSlopeLimiter(limiterType=1)
call g1%setSolverType(solverType=3)  !!!3=AdiHLLC,
call g1%setBoundaryType(boundaryType=3)
call g1%initVariable()
call g1%exchangeBdryMPI(g1%q,g1%winq)
call g1%setBoundary(g1%q)
call g1%writeGrid()
call g1%writeGrid_HD_vtk()
g1%writeFlag=.false.


t0=MPI_WTIME()
do while (g1%t .lt. g1%tend)
   call g1%griddt()
   call g1%evolveGridRK2()
   if(myid .eq. 0 ) then
     print *,"testSuite.f03: gridID=",g1%gridID,'t=',g1%t,"dt=",g1%dt
   endif
   if(g1%writeFlag .eqv. .true.) then
      t1=MPI_WTIME()
      if(myid .eq. 0) then
        print *,"testSuiteMPI.f03: elapse time=",t1-t0,"secs"
      endif
      call g1%writeGrid()
      call g1%writeGrid_HD_vtk()
      g1%writeFlag=.false.
      t0=MPI_WTIME()
   endif
enddo !! end do while


end subroutine HDRichtmyerMeshkovInstability2D




subroutine initHDRichtmyerMeshkovInstability2D(this,q)
class(grid)::this
double precision,dimension(1-this%nbuf:this%nMesh(1)+this%nbuf,1-this%nbuf:this%nMesh(2)+this%nbuf,this%nvar)::q
double precision,dimension(1-this%nbuf:this%nMesh(1)+this%nbuf)::xc
double precision,dimension(1-this%nbuf:this%nMesh(2)+this%nbuf)::yc
integer::nx,ny,nvar,nbuf,gridID,eosType
integer::i,j

double precision:: Ms=2.0d0
double precision:: p_tmp,rho_tmp,a_tmp,u_tmp,v_tmp
double precision:: rho_ratio=2.0,yitaY,L_rho=0.1

xc=this%xc(1)%coords
yc=this%xc(2)%coords
nx=this%nMesh(1)
ny=this%nMesh(2)

 do j=1,ny
   do i=1,nx
     if(xc(i)<=-0.42) then
       rho_tmp=1.4d0*8.0d0/3.0d0
       u_tmp=1.25
       v_tmp=0.0d0
       p_tmp=4.5d0
  
       q(i,j,1)=rho_tmp
       q(i,j,2)=rho_tmp*u_tmp
       q(i,j,3)=rho_tmp*v_tmp
       q(i,j,4)=p_tmp/(this%adiGamma-1.d0)+0.5*rho_tmp*(u_tmp**2+v_tmp**2)
       
     else     !!!if (xc(i)<=0.31415926d0*cos(2.0*yc(j))) then

       yitaY=0.31415926*( cos(2.0*yc(j))-1.0 )/2.0
       rho_tmp=(1.4d0+1.4*rho_ratio)/2.0-((1.4d0-1.4*rho_ratio)/2.0)*tanh( (xc(i)-yitaY)/L_rho  )
       u_tmp=0.0d0
       v_tmp=0.0d0
       p_tmp=1.0d0

       q(i,j,1)=rho_tmp
       q(i,j,2)=rho_tmp*u_tmp
       q(i,j,3)=rho_tmp*v_tmp
       q(i,j,4)=p_tmp/(this%adiGamma-1.d0)+0.5*rho_tmp*(u_tmp**2+v_tmp**2)



     endif

   enddo
 enddo

end subroutine initHDRichtmyerMeshkovInstability2D


subroutine initHDRichtmyerMeshkovInstability2D_old(this,q)
class(grid)::this
double precision,dimension(1-this%nbuf:this%nMesh(1)+this%nbuf,1-this%nbuf:this%nMesh(2)+this%nbuf,this%nvar)::q
double precision,dimension(1-this%nbuf:this%nMesh(1)+this%nbuf)::xc
double precision,dimension(1-this%nbuf:this%nMesh(2)+this%nbuf)::yc
integer::nx,ny,nvar,nbuf,gridID,eosType
integer::i,j

double precision:: Ms=2.0d0
double precision:: p_tmp,rho_tmp,a_tmp,u_tmp,v_tmp
double precision:: rho_ratio=2.0

xc=this%xc(1)%coords
yc=this%xc(2)%coords
nx=this%nMesh(1)
ny=this%nMesh(2)

 do j=1,ny
   do i=1,nx
     if(xc(i)<=-0.315) then
       rho_tmp=1.4d0*8.0d0/3.0d0
       u_tmp=1.25
       v_tmp=0.0d0
       p_tmp=4.5d0
  
       q(i,j,1)=rho_tmp
       q(i,j,2)=rho_tmp*u_tmp
       q(i,j,3)=rho_tmp*v_tmp
       q(i,j,4)=p_tmp/(this%adiGamma-1.d0)+0.5*rho_tmp*(u_tmp**2+v_tmp**2)
       
     elseif (xc(i)<=0.31415926d0*cos(2.0*yc(j))) then

       rho_tmp=1.4d0
       u_tmp=0.0d0
       v_tmp=0.0d0
       p_tmp=1.0d0

       q(i,j,1)=rho_tmp
       q(i,j,2)=rho_tmp*u_tmp
       q(i,j,3)=rho_tmp*v_tmp
       q(i,j,4)=p_tmp/(this%adiGamma-1.d0)+0.5*rho_tmp*(u_tmp**2+v_tmp**2)


     else
       rho_tmp=1.4d0*rho_ratio
       u_tmp=0.0d0
       v_tmp=0.0d0
       p_tmp=1.0d0

       q(i,j,1)=rho_tmp
       q(i,j,2)=rho_tmp*u_tmp
       q(i,j,3)=rho_tmp*v_tmp
       q(i,j,4)=p_tmp/(this%adiGamma-1.d0)+0.5*rho_tmp*(u_tmp**2+v_tmp**2)


     endif

   enddo
 enddo

end subroutine initHDRichtmyerMeshkovInstability2D_old



subroutine bdryHDRichtmyerMeshkovInstability2D(this,q)
class(grid)::this
double precision,dimension(1-this%nbuf:this%nMesh(1)+this%nbuf,1-this%nbuf:this%nMesh(2)+this%nbuf,this%nvar)::q
integer::nx,ny,nbuf,i,j

nbuf=this%nbuf
nx=this%nMesh(1)
ny=this%nMesh(2)

!if(this%boundaryType .eq. 1) then !! zero gradient
!!!!!! left boundary
  if(this%left_mpi .lt. 0) then
    do i=1,nbuf
       q(1-i,:,:)=q(1,:,:)
    enddo
  endif
!!!!!! right boundary
  if(this%right_mpi .lt. 0) then
    do i=1,nbuf
       q(nx+i,:,:)=q(nx,:,:)
    enddo
  endif
!endif


end subroutine bdryHDRichtmyerMeshkovInstability2D








!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!
!!!  2DMHD Richtmyer-Meshkov Instability 
!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine MHDRichtmyerMeshkovInstability2D(gridID)
integer::gridID
type(grid)::g1
integer::ndim,nbuf,coordType,variable(8)
integer::nMesh(2),dims(2)
double precision::leftBdry(2),rightBdry(2)
double precision::pi
logical::periods(2),reorder
double precision::t0,t1
integer::nstep,ierr

nstep=0
pi=4.d0*atan2(1.d0,1.d0)
variable=0
ndim=2
nbuf=2
coordType=1 !! Cartesian
variable(1)=1  !! rho
variable(2)=1  !! vx
variable(3)=1  !! vy
variable(4)=1  !! vz
variable(5)=1  !! bx
variable(6)=1  !! by
variable(7)=1  !! bz
variable(8)=1  !! energy


nMesh(1)=3072
nMesh(2)=256
leftBdry(1)=-2.0d0*3.1415926d0                !-0.05d0
leftBdry(2)=0.0d0                             !-0.5d0
rightBdry(1)=10.0d0*3.1415926d0                !0.05d0
rightBdry(2)=1.0d0*3.1415926d0                !0.5d0


!nMesh(1)=6400
!nMesh(2)=400
!leftBdry(1)=-3.0d0*3.1415926d0                !-0.05d0
!leftBdry(2)=0.0d0                             !-0.5d0
!rightBdry(1)=13.0d0*3.1415926d0                !0.05d0
!rightBdry(2)=1.0d0*3.1415926d0                !0.5d0
dims=(/0,0/)
call MPI_DIMS_CREATE(nprocs,ndim,dims,ierr)
dims(1)=nprocs
dims(2)=1
periods(1)=.false.
periods(2)=.true.
reorder=.true.


!call g1%setTopologyMPI(ndim,dims,periods,reorder)
!call g1%setMesh(nMesh,leftBdry,rightBdry,nbuf,coordType,gridID)
!call g1%setVariable(variable) !! den,vx,vy,vz,bx,by,bz,ene
!call g1%setMPIWindows()
!call g1%setEoS(eosType=2) !! adiabatic
!call g1%setadiGamma(gam=1.4d0) !! ratio of heat capacity
!call g1%setCFL(CFL=0.8d0) !! courant number
!call g1%setTime(fstart=0,tend=100.0d0,dtout=1.0d0)
!call g1%setSlopeLimiter(limiterType=1)
!call g1%setSolverType(solverType=5) !! 4=AdiHLLMHD, 5=AdiHLLDMHD
!call g1%setBoundaryType(boundaryType=3) !! 3=periodic
!call g1%initVariable()
!call g1%exchangeBdryMPI(g1%q,g1%winq)
!call g1%setBoundary(g1%q)
!call g1%writeGrid()
!call g1%writeGrid_MHD_vtk()
!g1%writeFlag=.false.



call g1%setTopologyMPI(ndim,dims,periods,reorder)
call g1%setMesh(nMesh,leftBdry,rightBdry,nbuf,coordType,gridID)
call g1%setVariable(variable) !! den,vx,vy,vz,bx,by,bz,ene
call g1%setMPIWindows()
call g1%setEoS(eosType=2) !! adiabatic
call g1%setadiGamma(gam=1.4d0) !! ratio of heat capacity
call g1%setCFL(CFL=0.8d0) !! courant number
call g1%setTime(fstart=0,tend=100.0d0,dtout=1.0d0)
call g1%setSlopeLimiter(limiterType=1)
call g1%setSolverType(solverType=5) !! 4=AdiHLLMHD, 5=AdiHLLDMHD
call g1%setBoundaryType(boundaryType=3) !! 3=periodic
call g1%initVariable()
call g1%exchangeBdryMPI(g1%q,g1%winq)
call g1%setBoundary(g1%q)
call g1%writeGrid()
call g1%writeGrid_MHD_vtk()
g1%writeFlag=.false.





t0=MPI_WTIME()
do while (g1%t .lt. g1%tend)
   call g1%griddt()
   call g1%evolveGridRK2()
   if(myid .eq. 0) then
     print *,"testSuite.f03: gridID=",g1%gridID,'t=',g1%t,"dt=",g1%dt
   endif
   if(g1%writeFlag .eqv. .true.) then
      t1=MPI_WTIME()
      if(myid .eq. 0) then
        print *,"testSuiteMPI.f03: elapse time=",t1-t0,"secs"
      endif
      call g1%writeGrid()
      call g1%writeGrid_MHD_vtk()
      g1%writeFlag=.false.
      t0=MPI_WTIME()
   endif
enddo !! end do while



end subroutine MHDRichtmyerMeshkovInstability2D





subroutine initMHDRichtmyerMeshkovInstability2D(this,q)
class(grid)::this
double precision,dimension(1-this%nbuf:this%nMesh(1)+this%nbuf,1-this%nbuf:this%nMesh(2)+this%nbuf,this%nvar)::q
double precision,dimension(1-this%nbuf:this%nMesh(1)+this%nbuf)::xc,xl,xr
double precision,dimension(1-this%nbuf:this%nMesh(2)+this%nbuf)::yc,yl,yr
integer::nx,ny,nvar,nbuf,gridID,eosType
integer::i,j
double precision::gam,rho,vx,vy,vz,bxl,byl,bzl,bxr,byr,bzr,p
double precision::bxc,byc,bzc
double precision::pi,b0
double precision::rc
double precision:: rho_ratio=2.0,yitaY,L_rho=0.1

pi=4.d0*atan2(1.d0,1.d0)

b0=0.1d0      !dsqrt(4.d0*pi/2.d0)

xc=this%xc(1)%coords
yc=this%xc(2)%coords
xl=this%xl(1)%coords
yl=this%xl(2)%coords
xr=this%xr(1)%coords
yr=this%xr(2)%coords
nx=this%nMesh(1)
ny=this%nMesh(2)
gam=this%adiGamma





do j=1,ny
  do i=1,nx

     if(xc(i)<=-0.42) then
       rho=1.4d0*8.0d0/3.0d0
       vx=1.25
       vy=0.0d0
       vz=0.0d0
        p=4.5d0

   bxl=b0
   byl=0.0d0
   bzl=0.0d0
   bxr=bxl
   byr=byl
   bzr=bzl


    bxc=0.5d0*(bxl+bxr)
    byc=0.5d0*(byl+byr)
    bzc=0.5d0*(bzl+bzr)



    q(i,j,1)=rho
    q(i,j,2)=rho*vx
    q(i,j,3)=rho*vy
    q(i,j,4)=rho*vz
    q(i,j,5)=bxl
    q(i,j,6)=byl
    q(i,j,7)=bzl
    q(i,j,8)=p/(gam-1.d0)+0.5d0*rho*(vx**2.d0+vy**2.d0+vz**2.d0)+0.5d0*(bxc**2.d0+byc**2.d0+bzc**2.d0)
    q(i,j,9)=bxr
    q(i,j,10)=byr
    q(i,j,11)=bzr


   else   !!!if (xc(i)<=0.31415926d0*cos(2.0*yc(j))) then

       yitaY=0.31415926*( cos(2.0*yc(j))-1.0 )/2.0
       rho=(1.4d0+1.4*rho_ratio)/2.0-((1.4d0-1.4*rho_ratio)/2.0)*tanh( (xc(i)-yitaY)/L_rho  )
       vx=0.0d0
       vy=0.0d0
       vz=0.0d0
        p=1.0d0

   bxl=b0
   byl=0.0d0
   bzl=0.0d0
   bxr=bxl
   byr=byl
   bzr=bzl


    bxc=0.5d0*(bxl+bxr)
    byc=0.5d0*(byl+byr)
    bzc=0.5d0*(bzl+bzr)



    q(i,j,1)=rho
    q(i,j,2)=rho*vx
    q(i,j,3)=rho*vy
    q(i,j,4)=rho*vz
    q(i,j,5)=bxl
    q(i,j,6)=byl
    q(i,j,7)=bzl
    q(i,j,8)=p/(gam-1.d0)+0.5d0*rho*(vx**2.d0+vy**2.d0+vz**2.d0)+0.5d0*(bxc**2.d0+byc**2.d0+bzc**2.d0)
    q(i,j,9)=bxr
    q(i,j,10)=byr
    q(i,j,11)=bzr

   endif



  enddo
enddo

end subroutine initMHDRichtmyerMeshkovInstability2D





subroutine initMHDRichtmyerMeshkovInstability2D_old(this,q)
class(grid)::this
double precision,dimension(1-this%nbuf:this%nMesh(1)+this%nbuf,1-this%nbuf:this%nMesh(2)+this%nbuf,this%nvar)::q
double precision,dimension(1-this%nbuf:this%nMesh(1)+this%nbuf)::xc,xl,xr
double precision,dimension(1-this%nbuf:this%nMesh(2)+this%nbuf)::yc,yl,yr
integer::nx,ny,nvar,nbuf,gridID,eosType
integer::i,j
double precision::gam,rho,vx,vy,vz,bxl,byl,bzl,bxr,byr,bzr,p
double precision::bxc,byc,bzc
double precision::pi,b0
double precision::rc
double precision:: rho_ratio=2.0

pi=4.d0*atan2(1.d0,1.d0)

b0=0.1d0      !dsqrt(4.d0*pi/2.d0)

xc=this%xc(1)%coords
yc=this%xc(2)%coords
xl=this%xl(1)%coords
yl=this%xl(2)%coords
xr=this%xr(1)%coords
yr=this%xr(2)%coords
nx=this%nMesh(1)
ny=this%nMesh(2)
gam=this%adiGamma





do j=1,ny
  do i=1,nx

     if(xc(i)<=-0.315) then
       rho=1.4d0*8.0d0/3.0d0
       vx=1.25
       vy=0.0d0
       vz=0.0d0
        p=4.5d0

   bxl=b0
   byl=0.0d0
   bzl=0.0d0
   bxr=bxl
   byr=byl
   bzr=bzl


    bxc=0.5d0*(bxl+bxr)
    byc=0.5d0*(byl+byr)
    bzc=0.5d0*(bzl+bzr)



    q(i,j,1)=rho
    q(i,j,2)=rho*vx
    q(i,j,3)=rho*vy
    q(i,j,4)=rho*vz
    q(i,j,5)=bxl
    q(i,j,6)=byl
    q(i,j,7)=bzl
    q(i,j,8)=p/(gam-1.d0)+0.5d0*rho*(vx**2.d0+vy**2.d0+vz**2.d0)+0.5d0*(bxc**2.d0+byc**2.d0+bzc**2.d0)
    q(i,j,9)=bxr
    q(i,j,10)=byr
    q(i,j,11)=bzr


   elseif (xc(i)<=0.31415926d0*cos(2.0*yc(j))) then

       rho=1.4d0
       vx=0.0d0
       vy=0.0d0
       vz=0.0d0
        p=1.0d0

   bxl=b0
   byl=0.0d0
   bzl=0.0d0
   bxr=bxl
   byr=byl
   bzr=bzl


    bxc=0.5d0*(bxl+bxr)
    byc=0.5d0*(byl+byr)
    bzc=0.5d0*(bzl+bzr)



    q(i,j,1)=rho
    q(i,j,2)=rho*vx
    q(i,j,3)=rho*vy
    q(i,j,4)=rho*vz
    q(i,j,5)=bxl
    q(i,j,6)=byl
    q(i,j,7)=bzl
    q(i,j,8)=p/(gam-1.d0)+0.5d0*rho*(vx**2.d0+vy**2.d0+vz**2.d0)+0.5d0*(bxc**2.d0+byc**2.d0+bzc**2.d0)
    q(i,j,9)=bxr
    q(i,j,10)=byr
    q(i,j,11)=bzr




   else

   
       rho=1.4d0*rho_ratio
       vx=0.0d0
       vy=0.0d0
       vz=0.0d0
        p=1.0d0

   bxl=b0
   byl=0.0d0
   bzl=0.0d0
   bxr=bxl
   byr=byl
   bzr=bzl


    bxc=0.5d0*(bxl+bxr)
    byc=0.5d0*(byl+byr)
    bzc=0.5d0*(bzl+bzr)



    q(i,j,1)=rho
    q(i,j,2)=rho*vx
    q(i,j,3)=rho*vy
    q(i,j,4)=rho*vz
    q(i,j,5)=bxl
    q(i,j,6)=byl
    q(i,j,7)=bzl
    q(i,j,8)=p/(gam-1.d0)+0.5d0*rho*(vx**2.d0+vy**2.d0+vz**2.d0)+0.5d0*(bxc**2.d0+byc**2.d0+bzc**2.d0)
    q(i,j,9)=bxr
    q(i,j,10)=byr
    q(i,j,11)=bzr






   endif



  enddo
enddo

end subroutine initMHDRichtmyerMeshkovInstability2D_old

subroutine bdryMHDRichtmyerMeshkovInstability2D(this,q)
class(grid)::this
double precision,dimension(1-this%nbuf:this%nMesh(1)+this%nbuf,1-this%nbuf:this%nMesh(2)+this%nbuf,this%nvar)::q
integer::nx,ny,nbuf,i,j

nbuf=this%nbuf
nx=this%nMesh(1)
ny=this%nMesh(2)

!if(this%boundaryType .eq. 3) then !! periodic boundary condition
!endif

!!!!!! left boundary
  if(this%left_mpi .lt. 0) then
    do i=1,nbuf
       q(1-i,:,:)=q(1,:,:)
    enddo
  endif
!!!!!! right boundary
  if(this%right_mpi .lt. 0) then
    do i=1,nbuf
       q(nx+i,:,:)=q(nx,:,:)
    enddo
  endif

end subroutine bdryMHDRichtmyerMeshkovInstability2D






!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!
!!!  3DHD Richtmyer-Meshkov Instability
!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine HDRichtmyerMeshkovInstability3D(gridID)

integer::gridID
type(grid)::g1
integer::ndim,nbuf,coordType,variable(8)
integer::nMesh(3),dims(3)
double precision::leftBdry(3),rightBdry(3)
logical::periods(3),reorder
double precision::t0,t1
integer::nstep,ierr

nstep=0
variable=0
ndim=3
nbuf=2
coordType=1
variable(1)=1
variable(2)=1
variable(3)=1
variable(4)=1
variable(8)=1

nMesh(1)=3200
nMesh(2)=200
nMesh(3)=4

leftBdry(1)=-3.0d0*3.1415926d0                !-0.05d0
leftBdry(2)=0.0d0 
leftBdry(3)=0.0d0                            !-0.5d0
rightBdry(1)=13.0d0*3.1415926d0                !0.05d0
rightBdry(2)=1.0d0*3.1415926d0                !0.5d0
rightBdry(3)=1.0d0
dims=(/0,0,0/)
call MPI_DIMS_CREATE(nprocs,ndim,dims,ierr)
dims(1)=nprocs
dims(2)=1
dims(3)=1
periods(1)=.false.
periods(2)=.true.
periods(3)=.true.

reorder=.true.



call g1%setTopologyMPI(ndim,dims,periods,reorder)
call g1%setMesh(nMesh,leftBdry,rightBdry,nbuf,coordType,gridID)
call g1%setVariable(variable) !! den,vx
call g1%setMPIWindows()
call g1%setEoS(eosType=2) !! adiabatic
call g1%setadiGamma(gam=1.4d0) 
call g1%setCFL(CFL=0.8d0) !! Courant number
call g1%setTime(fstart=0,tend=100.0d0,dtout=1.0d0)

!call g1%setRestart(fstart=1,tend=100.0d0,dtout=1.0d0)

call g1%setSlopeLimiter(limiterType=1)
call g1%setSolverType(solverType=3)  !!!3=AdiHLLC,
call g1%setBoundaryType(boundaryType=3)

if(g1%isRestart==0)then
call g1%initVariable()
elseif(g1%isRestart==1)then
  call g1%readGrid_HD_vtk()
  !call g1%readGrid_MHD_vtk()
endif

call g1%exchangeBdryMPI(g1%q,g1%winq)
call g1%setBoundary(g1%q)

if(g1%isRestart==0)then
call g1%writeGrid()
call g1%writeGrid_HD_vtk()
endif

g1%writeFlag=.false.

t0=MPI_WTIME()
do while (g1%t .lt. g1%tend)
   call g1%griddt()
   call g1%evolveGridRK2()
   if(myid .eq. 0 ) then
     print *,"testSuite.f03: gridID=",g1%gridID,'t=',g1%t,"dt=",g1%dt
   endif
   if(g1%writeFlag .eqv. .true.) then
      t1=MPI_WTIME()
      if(myid .eq. 0) then
        print *,"testSuiteMPI.f03: elapse time=",t1-t0,"secs"
      endif
      call g1%writeGrid()
      call g1%writeGrid_HD_vtk()
      g1%writeFlag=.false.
      t0=MPI_WTIME()
   endif
enddo !! end do while


end subroutine HDRichtmyerMeshkovInstability3D


subroutine initHDRichtmyerMeshkovInstability3D(this,q)
class(grid)::this
double precision,dimension(1-this%nbuf:this%nMesh(1)+this%nbuf,1-this%nbuf:this%nMesh(2)+this%nbuf,&
& 1-this%nbuf:this%nMesh(3)+this%nbuf,this%nvar)::q
double precision,dimension(1-this%nbuf:this%nMesh(1)+this%nbuf)::xc
double precision,dimension(1-this%nbuf:this%nMesh(2)+this%nbuf)::yc
double precision,dimension(1-this%nbuf:this%nMesh(3)+this%nbuf)::zc
integer::nx,ny,nz,nvar,nbuf,gridID,eosType
integer::i,j,k

double precision:: Ms=2.0d0
double precision:: p_tmp,rho_tmp,a_tmp,u_tmp,v_tmp,w_tmp
double precision:: rho_ratio=2.0

xc=this%xc(1)%coords
yc=this%xc(2)%coords
zc=this%xc(3)%coords
nx=this%nMesh(1)
ny=this%nMesh(2)
nz=this%nMesh(3)

do k=1,nz
 do j=1,ny
   do i=1,nx
     if(xc(i)<=-0.315) then
       rho_tmp=1.4d0*8.0d0/3.0d0
       u_tmp=1.25
       v_tmp=0.0d0
       w_tmp=0.0d0
       p_tmp=4.5d0
  
       q(i,j,k,1)=rho_tmp
       q(i,j,k,2)=rho_tmp*u_tmp
       q(i,j,k,3)=rho_tmp*v_tmp
       q(i,j,k,4)=rho_tmp*w_tmp
       q(i,j,k,5)=p_tmp/(this%adiGamma-1.d0)+0.5*rho_tmp*(u_tmp**2+v_tmp**2+w_tmp**2)
       
     elseif (xc(i)<=0.31415926d0*cos(2.0*yc(j))) then

       rho_tmp=1.4d0
       u_tmp=0.0d0
       v_tmp=0.0d0
       p_tmp=1.0d0

       q(i,j,k,1)=rho_tmp
       q(i,j,k,2)=rho_tmp*u_tmp
       q(i,j,k,3)=rho_tmp*v_tmp
       q(i,j,k,4)=rho_tmp*w_tmp
       q(i,j,k,5)=p_tmp/(this%adiGamma-1.d0)+0.5*rho_tmp*(u_tmp**2+v_tmp**2+w_tmp**2)


     else
       rho_tmp=1.4d0*rho_ratio
       u_tmp=0.0d0
       v_tmp=0.0d0
       p_tmp=1.0d0

       q(i,j,k,1)=rho_tmp
       q(i,j,k,2)=rho_tmp*u_tmp
       q(i,j,k,3)=rho_tmp*v_tmp
       q(i,j,k,4)=rho_tmp*w_tmp
       q(i,j,k,5)=p_tmp/(this%adiGamma-1.d0)+0.5*rho_tmp*(u_tmp**2+v_tmp**2+w_tmp**2)

     endif

   enddo
 enddo
enddo

end subroutine initHDRichtmyerMeshkovInstability3D



subroutine bdryHDRichtmyerMeshkovInstability3D(this,q)
class(grid)::this
double precision,dimension(1-this%nbuf:this%nMesh(1)+this%nbuf,1-this%nbuf:this%nMesh(2)+this%nbuf,&
& 1-this%nbuf:this%nMesh(3)+this%nbuf,this%nvar)::q
integer::nx,ny,nz,nbuf,i,j,k

nbuf=this%nbuf
nx=this%nMesh(1)
ny=this%nMesh(2)
nz=this%nMesh(3)


!if(this%boundaryType .eq. 1) then !! zero gradient
!!!!!! left boundary
  if(this%left_mpi .lt. 0) then
    do i=1,nbuf
       q(1-i,:,:,:)=q(1,:,:,:)
    enddo
  endif
!!!!!! right boundary
  if(this%right_mpi .lt. 0) then
    do i=1,nbuf
       q(nx+i,:,:,:)=q(nx,:,:,:)
    enddo
  endif
!endif


end subroutine bdryHDRichtmyerMeshkovInstability3D




!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!
!!!  3DMHD Richtmyer-Meshkov Instability
!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine MHDRichtmyerMeshkovInstability3D(gridID)
integer::gridID
type(grid)::g1
integer::ndim,nbuf,coordType,variable(8)
integer::nMesh(3),dims(3)
double precision::leftBdry(3),rightBdry(3)
logical::periods(3),reorder
double precision::t0,t1, dt,tend
integer::nstep,ierr
double precision::pi



nstep=0
variable=0
ndim=3
nbuf=2
coordType=1 !! Cartesian
variable(1)=1  !! rho
variable(2)=1  !! vx
variable(3)=1  !! vy
variable(4)=1  !! vz
variable(5)=1  !! bx
variable(6)=1  !! by
variable(7)=1  !! bz
variable(8)=1  !! energy

pi=4.d0*atan2(1.d0,1.d0)


nMesh(1)=3200
nMesh(2)=200
nMesh(3)=4
leftBdry(1)=-3.0d0*3.1415926d0                !-0.05d0
leftBdry(2)=0.0d0                             !-0.5d0
leftBdry(3)=0.0d0
rightBdry(1)=13.0d0*3.1415926d0                !0.05d0
rightBdry(2)=1.0d0*3.1415926d0                !0.5d0
rightBdry(3)=1.0d0
dims=(/0,0,0/)
call MPI_DIMS_CREATE(nprocs,ndim,dims,ierr)
dims(1)=nprocs
dims(2)=1
dims(3)=1
periods(1)=.false.
periods(2)=.true.
periods(3)=.true.
reorder=.true.

call g1%setGridID(gridID)
call g1%setTopologyMPI(ndim,dims,periods,reorder)
call g1%setMesh(nMesh,leftBdry,rightBdry,nbuf,coordType,gridID)
call g1%setVariable(variable) !! den,vx,vy,vz,bx,by,bz,ene
call g1%setMPIWindows()

!call g1%setEoS(eosType=1) !!1=isothermal, 2=adiabatic, 3=polytropic
!call g1%setSoundSpeed(snd=1.0d0) !! sound speed

call g1%setEoS(eosType=2) !! adiabatic
call g1%setadiGamma(gam=1.4d0) !! ratio of heat capacity

call g1%setCFL(CFL=0.8d0) !! courant number
call g1%setTime(fstart=0,tend=tend,dtout=1.0d0)
call g1%setSlopeLimiter(limiterType=3)
call g1%setSolverType(solverType=5) !! 4=AdiHLLMHD, 5=AdiHLLDMHD
call g1%setBoundaryType(boundaryType=3) !! 3=periodic
call g1%initVariable()
call g1%setBoundary(g1%q)
call g1%exchangeBdryMPI(g1%q,g1%winq)
call g1%writeGrid()
call g1%writeGrid_MHD_vtk()


g1%writeFlag=.false.
t0=MPI_WTIME()
do while (g1%t .lt. g1%tend)
   call g1%griddt()
   dt=g1%dt


   if(dt .gt. g1%toutput-g1%dt)then
     dt=g1%toutput-g1%dt
     g1%writeFlag=.true.
     g1%toutput=g1%toutput+g1%dtout
     g1%fnum=g1%fnum+1
   elseif (dt .gt. g1%tend-g1%t) then
     dt=g1%tend-g1%t
     g1%writeFlag=.true.
     g1%fnum=g1%fnum+1
   endif
   g1%dt=dt




   call rk2MHD_3D(g1,g1%q,g1%q1,g1%q2)
   call g1%exchangeBdryMPI(g1%q,g1%winq)
   call g1%setBoundary(g1%q)


   if(myid .eq. 0) then
     print *,"testSuite.f03: gridID=",g1%gridID,'t=',g1%t,"dt=",g1%dt
   endif
   if(g1%writeFlag .eqv. .true.) then
      t1=MPI_WTIME()
      if(myid .eq. 0) then
        print *,"testSuiteMPI.f03: elapse time=",t1-t0,"secs"
      endif
      call g1%writeGrid()
      call g1%writeGrid_MHD_vtk()
      g1%writeFlag=.false.
      t0=MPI_WTIME()
   endif
enddo !! end do while



end subroutine MHDRichtmyerMeshkovInstability3D





subroutine initMHDRichtmyerMeshkovInstability3D(this,q)
class(grid)::this
double precision,dimension(1-this%nbuf:this%nMesh(1)+this%nbuf,1-this%nbuf:this%nMesh(2)+this%nbuf, &
                         & 1-this%nbuf:this%nMesh(3)+this%nbuf,this%nvar)::q
double precision,dimension(1-this%nbuf:this%nMesh(1)+this%nbuf)::xc,xl,xr
double precision,dimension(1-this%nbuf:this%nMesh(2)+this%nbuf)::yc,yl,yr
double precision,dimension(1-this%nbuf:this%nMesh(3)+this%nbuf)::zc,zl,zr

integer::nx,ny,nz,nvar,nbuf,gridID,eosType
integer::i,j,k
double precision::gam,rho,vx,vy,vz,bxl,byl,bzl,bxr,byr,bzr,p
double precision::bxc,byc,bzc
double precision::pi,b0
double precision::rc
double precision:: rho_ratio=2.0

pi=4.d0*atan2(1.d0,1.d0)

b0=0.1d0      !dsqrt(4.d0*pi/2.d0)


xc=this%xc(1)%coords
yc=this%xc(2)%coords
zc=this%xc(3)%coords
xl=this%xl(1)%coords
yl=this%xl(2)%coords
zl=this%xl(3)%coords
xr=this%xr(1)%coords
yr=this%xr(2)%coords
zr=this%xr(3)%coords
nx=this%nMesh(1)
ny=this%nMesh(2)
nz=this%nMesh(3)
gam=this%adiGamma




do k=1,nz
do j=1,ny
  do i=1,nx

     if(xc(i)<=-0.315) then
       rho=1.4d0*8.0d0/3.0d0
       vx=1.25
       vy=0.0d0
       vz=0.0d0
        p=4.5d0

   bxl=b0
   byl=0.0d0
   bzl=0.0d0
   bxr=bxl
   byr=byl
   bzr=bzl


    bxc=0.5d0*(bxl+bxr)
    byc=0.5d0*(byl+byr)
    bzc=0.5d0*(bzl+bzr)



    q(i,j,k,1)=rho
    q(i,j,k,2)=rho*vx
    q(i,j,k,3)=rho*vy
    q(i,j,k,4)=rho*vz
    q(i,j,k,5)=bxl
    q(i,j,k,6)=byl
    q(i,j,k,7)=bzl
    !q(i,j,k,8)=0.0d0
    q(i,j,k,8)=p/(gam-1.d0)+0.5d0*rho*(vx**2.d0+vy**2.d0+vz**2.d0)+0.5d0*(bxc**2.d0+byc**2.d0+bzc**2.d0)
    q(i,j,k,9)=bxr
    q(i,j,k,10)=byr
    q(i,j,k,11)=bzr


    elseif (xc(i)<=0.31415926d0*cos(2.0*yc(j))) then
    !elseif (xc(i)<=cos(2.0*yc(j))) then

       rho=1.4d0
       vx=0.0d0
       vy=0.0d0
       vz=0.0d0
        p=1.0d0

   bxl=b0
   byl=0.0d0
   bzl=0.0d0
   bxr=bxl
   byr=byl
   bzr=bzl


    bxc=0.5d0*(bxl+bxr)
    byc=0.5d0*(byl+byr)
    bzc=0.5d0*(bzl+bzr)



    q(i,j,k,1)=rho
    q(i,j,k,2)=rho*vx
    q(i,j,k,3)=rho*vy
    q(i,j,k,4)=rho*vz
    q(i,j,k,5)=bxl
    q(i,j,k,6)=byl
    q(i,j,k,7)=bzl
    !q(i,j,k,8)=0.0d0
    q(i,j,k,8)=p/(gam-1.d0)+0.5d0*rho*(vx**2.d0+vy**2.d0+vz**2.d0)+0.5d0*(bxc**2.d0+byc**2.d0+bzc**2.d0)
    q(i,j,k,9)=bxr
    q(i,j,k,10)=byr
    q(i,j,k,11)=bzr




   else

   
       rho=1.4d0*rho_ratio
       vx=0.0d0
       vy=0.0d0
       vz=0.0d0
        p=1.0d0

   bxl=b0
   byl=0.0d0
   bzl=0.0d0
   bxr=bxl
   byr=byl
   bzr=bzl


    bxc=0.5d0*(bxl+bxr)
    byc=0.5d0*(byl+byr)
    bzc=0.5d0*(bzl+bzr)



    q(i,j,k,1)=rho
    q(i,j,k,2)=rho*vx
    q(i,j,k,3)=rho*vy
    q(i,j,k,4)=rho*vz
    q(i,j,k,5)=bxl
    q(i,j,k,6)=byl
    q(i,j,k,7)=bzl
    !q(i,j,k,8)=0.0d0
    q(i,j,k,8)=p/(gam-1.d0)+0.5d0*rho*(vx**2.d0+vy**2.d0+vz**2.d0)+0.5d0*(bxc**2.d0+byc**2.d0+bzc**2.d0)
    q(i,j,k,9)=bxr
    q(i,j,k,10)=byr
    q(i,j,k,11)=bzr



   endif



  enddo
enddo
enddo

end subroutine initMHDRichtmyerMeshkovInstability3D





subroutine bdryMHDRichtmyerMeshkovInstability3D(this,q)
class(grid)::this
double precision,dimension(1-this%nbuf:this%nMesh(1)+this%nbuf,1-this%nbuf:this%nMesh(2)+this%nbuf, &
                         & 1-this%nbuf:this%nMesh(3)+this%nbuf,this%nvar)::q
integer::nx,ny,nz,nbuf,i,j,k

nbuf=this%nbuf
nx=this%nMesh(1)
ny=this%nMesh(2)
ny=this%nMesh(3)

!if(this%boundaryType .eq. 3) then !! periodic boundary condition
!endif

!!!!!! left boundary
  if(this%left_mpi .lt. 0) then
    do i=1,nbuf
       q(1-i,:,:,:)=q(1,:,:,:)
    enddo
  endif
!!!!!! right boundary
  if(this%right_mpi .lt. 0) then
    do i=1,nbuf
       q(nx+i,:,:,:)=q(nx,:,:,:)
    enddo
  endif

end subroutine bdryMHDRichtmyerMeshkovInstability3D





!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!! Adi MHD Shocktube 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine AdiShockTubeMHD3D(gridID)
integer::gridID
type(grid)::g1
integer::ndim,nbuf,coordType,variable(8)
integer::nMesh(3),dims(3)
double precision::leftBdry(3),rightBdry(3)
logical::periods(3),reorder
double precision::t0,t1,dt,tend
integer::nstep,ierr

nstep=0
variable=0
ndim=3
nbuf=2
coordType=1 !! Cartesian
variable(1)=1  !! rho
variable(2)=1  !! vx
variable(3)=1  !! vy
variable(4)=1  !! vz
variable(5)=1  !! bx
variable(6)=1  !! by
variable(7)=1  !! bz
variable(8)=1  !! energy

nMesh(1)=3200
nMesh(2)=200
nMesh(3)=4
leftBdry(1)=-3.0d0*3.1415926d0                !-0.05d0
leftBdry(2)=0.0d0                             !-0.5d0
leftBdry(3)=0.0d0
rightBdry(1)=13.0d0*3.1415926d0                !0.05d0
rightBdry(2)=1.0d0*3.1415926d0                !0.5d0
rightBdry(3)=1.0d0
dims=(/0,0,0/)
call MPI_DIMS_CREATE(nprocs,ndim,dims,ierr)
dims(1)=nprocs
dims(2)=1
dims(3)=1

periods(1)=.false.
periods(2)=.true.
periods(3)=.true.
reorder=.true.

tend=20.0

call g1%setTopologyMPI(ndim,dims,periods,reorder)
call g1%setMesh(nMesh,leftBdry,rightBdry,nbuf,coordType,gridID)
call g1%setVariable(variable) !! den,vx,vy,vz,bx,by,bz,ene
call g1%setMPIWindows()

!call g1%setEoS(eosType=1) !!1=isothermal, 2=adiabatic, 3=polytropic
!call g1%setSoundSpeed(snd=sqrt(1.0d0/1.4d0)) !! sound speed

call g1%setEoS(eosType=2) !!1=isothermal, 2=adiabatic, 3=polytropic
call g1%setadiGamma(gam=1.4d0) !! ratio of heat capacity

call g1%setCFL(CFL=0.8d0) !! courant number

call g1%setTime(fstart=0,tend=tend,dtout=0.5d0)
call g1%setRestart(fstart=1,tend=tend,dtout=0.5d0)

call g1%setSlopeLimiter(limiterType=1)
call g1%setSolverType(solverType=5) !! 4=AdiHLLMHD, 5=AdiHLLDMHD
call g1%setBoundaryType(boundaryType=1)

if(g1%isRestart==0)then
  call g1%initVariable()
elseif(g1%isRestart==1)then
 !call g1%readGrid_HD_vtk()
  call g1%readGrid_MHD_vtk()
endif

call g1%setBoundary(g1%q)
call g1%exchangeBdryMPI(g1%q,g1%winq)

if(g1%isRestart==0)then
call g1%writeGrid()
call g1%writeGrid_MHD_vtk()
endif

g1%writeFlag=.false.

t0=MPI_WTIME()
do while (g1%t .lt. g1%tend)
   nstep=nstep+1
   call g1%griddt()

   dt=g1%dt


   if(dt .gt. g1%toutput-g1%dt)then
     dt=g1%toutput-g1%dt
     g1%writeFlag=.true.
     g1%toutput=g1%toutput+g1%dtout
     g1%fnum=g1%fnum+1
   elseif (dt .gt. g1%tend-g1%t) then
     dt=g1%tend-g1%t
     g1%writeFlag=.true.
     g1%fnum=g1%fnum+1
   endif
   g1%dt=dt


   call rk2MHD_3D(g1,g1%q,g1%q1,g1%q2)

   call g1%exchangeBdryMPI(g1%q,g1%winq)
   call g1%setBoundary(g1%q)

   if(myid .eq. 0 ) then
     print *,"testSuite.f03: gridID=",g1%gridID,'t=',g1%t,"dt=",g1%dt
   endif
   if(g1%writeFlag .eqv. .true.) then
      t1=MPI_WTIME()
      if(myid .eq. 0) then
        print *,"testSuiteMPI.f03: elapse time=",t1-t0,"secs"
      endif
      call g1%writeGrid()
      call g1%writeGrid_MHD_vtk()
      g1%writeFlag=.false.
      t0=MPI_WTIME()
   endif

enddo !! end do while

end subroutine AdiShockTubeMHD3D

subroutine initAdiShockTubeMHD3D(this,q)
class(grid)::this
double precision,dimension(1-this%nbuf:this%nMesh(1)+this%nbuf,1-this%nbuf:this%nMesh(2)+this%nbuf,&
                           1-this%nbuf:this%nMesh(3)+this%nbuf,this%nvar)::q
double precision,dimension(1-this%nbuf:this%nMesh(1)+this%nbuf)::xc
double precision,dimension(1-this%nbuf:this%nMesh(2)+this%nbuf)::yc
double precision,dimension(1-this%nbuf:this%nMesh(3)+this%nbuf)::zc

integer::nx,ny,nz,nvar,nbuf,gridID,eosType
integer::i,j,k
   double precision::gam,rho,vx,vy,vz,bxl,byl,bzl,bxr,byr,bzr,p
   double precision::bxc,byc,bzc,pi,snd,rho_ratio

xc=this%xc(1)%coords
yc=this%xc(2)%coords
zc=this%xc(3)%coords
nx=this%nMesh(1)
ny=this%nMesh(2)
nz=this%nMesh(3)
snd=this%snd
pi=4.d0*datan2(1.d0,1.d0)
gam=this%adiGamma
rho_ratio=2.0

do k=1,nz
  do j=1,ny
    do i=1,nx
      if(xc(i) .lt. -0.315d0) then
         rho=1.4d0*8.0d0/3.0d0
         vx=1.25d0
         vy=0.d0
         vz=0.0d0
         p=4.5d0
         bxl=0.1d0
         bxr=bxl
         byl=0.0d0
         byr=0.0d0
         bzl=0.0d0
         bzr=0.0d0
      elseif (xc(i)<=0.31415926d0*cos(2.0*yc(j))) then
          rho=1.4d0
          vx=0.d0
          vy=0.d0
          vz=0.d0
          p=1.0d0

         !rho=1.4d0*8.0d0/3.0d0
         !vx=1.25d0
         !vy=0.d0
         !vz=0.0d0
         !p=4.5d0

          bxl=0.1d0
          bxr=bxl
          byl=0.0d0
          byr=0.0d0
          bzl=0.0d0
          bzr=0.0d0

      else

          rho=1.4d0*rho_ratio
          vx=0.d0
          vy=0.d0
          vz=0.d0
          p=1.0d0

         !rho=1.4d0*8.0d0/3.0d0
         !vx=1.25d0
         !vy=0.d0
         !vz=0.0d0
         !p=4.5d0

          bxl=0.1d0
          bxr=bxl
          byl=0.0d0
          byr=0.0d0
          bzl=0.0d0
          bzr=0.0d0

      endif

      bxc=0.5d0*(bxl+bxr)
      byc=0.5d0*(byl+byr)
      bzc=0.5d0*(bzl+bzr)

      q(i,j,k,1)=rho
      q(i,j,k,2)=rho*vx
      q(i,j,k,3)=rho*vy
      q(i,j,k,4)=rho*vz
      q(i,j,k,5)=bxl
      q(i,j,k,6)=byl
      q(i,j,k,7)=bzl
      !q(i,j,k,8)=0.d0
      q(i,j,k,8)=p/(gam-1.d0)+0.5d0*rho*(vx**2.d0+vy**2.d0+vz**2.d0)+0.5d0*(bxc**2.d0+byc**2.d0+bzc**2.d0)
      q(i,j,k,9)=bxr
      q(i,j,k,10)=byr
      q(i,j,k,11)=bzr
    enddo
  enddo
enddo

end subroutine initAdiShockTubeMHD3D

subroutine bdryAdiShockTubeMHD3D(this,q)
class(grid)::this
double precision,dimension(1-this%nbuf:this%nMesh(1)+this%nbuf,1-this%nbuf:this%nMesh(2)+this%nbuf,&
                           1-this%nbuf:this%nMesh(3)+this%nbuf,this%nvar)::q
integer::nx,ny,nz,nbuf,i,j,k

nbuf=this%nbuf
nx=this%nMesh(1)
ny=this%nMesh(2)
nz=this%nMesh(3)

!if(this%boundaryType .eq. 1) then !! zero gradient

!!!!!! left boundary
  if(this%left_mpi .lt. 0) then
    do i=1,nbuf
       q(1-i,:,:,:)=q(1,:,:,:)
    enddo
  endif
!!!!!! right boundary
  if(this%right_mpi .lt. 0) then
    do i=1,nbuf
       q(nx+i,:,:,:)=q(nx,:,:,:)
    enddo
  endif


!!!!!! up boundary
  !if(this%up_mpi .lt. 0) then
    !do j=1,nbuf
       !q(:,1-j,:,:)=q(:,1,:,:)
    !enddo
  !endif
!!!!!! down boundary
  !if(this%down_mpi .lt. 0) then
    !do j=1,nbuf
       !q(:,ny+j,:,:)=q(:,ny,:,:)
    !enddo
  !endif
!!!! top boundary
   !if(this%top_mpi .lt. 0) then
     !do k=1,nbuf
       !q(:,:,1-k,:)=q(:,:,1,:)
     !enddo
   !endif
!!!! bottom boundary
   !if(this%bottom_mpi .lt. 0) then
     !do k=1,nbuf
       !q(:,:,nz+k,:)=q(:,:,nz,:)
     !enddo
   !endif
!!!!

!endif

end subroutine bdryAdiShockTubeMHD3D

end module testSuiteMPI
