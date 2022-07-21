module gridModule
use mpi
use, intrinsic::iso_c_binding
implicit none

include 'fftw3-mpi.f03'

integer::nprocs,myid

type coord1D
   double precision,dimension(:),allocatable::coords
end type coord1D

type grid
!! solverType (1=exactHD,2=HLLHD,3=HLLCHD,4=HLLMHD,5=HLLDMHD,6=polyHLL)
!! limiterType (0=zero, 1=van Leer, 2=fslop, 3=minmod)
!! eosType (1=isothermal, 2=adiabatic, 3=polytropic)
!! boundaryType (0=user defined, 1=zero gradient, 2=reflective,3=periodic)
!! sgBdryType (0=isolated [default], 1=periodic)
!! ***NOTE**** not all the combinations of solverType and eosType are meaningful and available, the following has been (to be) implemented
!! exactHD => isothermal, (adiabatic,polytropic)
!! HLLHD => isothermal, adiabatic, polytropic
!! HLLCHD => adiabatic
!! HLLDMHD => adiabatic, isothermal
!! coordType (1=Cartesian, 2=cylindrical log in r, 3=cylindrical uniform in r)

   integer,dimension(:),allocatable::nMesh,mpiCoord,nMesh_global

   integer::nbuf,coordType,gridID,nvar,ndim,fstart,fnum,limiterType,solverType,eosType,boundaryType,variable(8)
  
   !!!!zwg added for spactial dimension index
   integer:: Sp_D
  

   double precision,dimension(:),allocatable::leftBdry,rightBdry,q,q1,q2

   double precision,dimension(:),allocatable::sgfx,sgfy,sgfz,sgDensityBuffer

   complex(C_DOUBLE_COMPLEX),dimension(:),pointer::sgfxKernel,sgfyKernel,sgfzKernel,sgDenCmplx

   complex(C_DOUBLE_COMPLEX),dimension(:),pointer::sgfxCmplx,sgfyCmplx,sgfzCmplx  

   double precision,dimension(:),allocatable::databuf1,databuf2,leftbuf,rightbuf,frontbuf,backbuf

   double precision,dimension(:),allocatable::leftBdry_global,rightBdry_global

   type(coord1D),dimension(:),allocatable::xl,xr,xc,dx

   double precision::t,tend,CFL,dtout,dt,snd,toutput,adiGamma,GravConst,polyGamma,polyK

   logical:: writeFlag,changeSolver,enable_sg=.false.,enable_ad=.false.,neg_pressure

   integer::vu_mpi,left_mpi,right_mpi,up_mpi,down_mpi,top_mpi,bottom_mpi

   integer,dimension(:),allocatable::dims_mpi

   integer::winq,winq1,winq2,winbuf1,winbuf2,winsg

   integer,dimension(:),allocatable::nx_list,ny_list,nz_list
   integer::gnx_list_R,gny_list_R,gnz_list_R

   type(C_PTR)::sgPlanDen,sgPlanFx,sgPlanFy,sgPlanFz

   double precision::mu_ad,alpha_ad

   integer::sgBdryType

   !!! zwg: the following avariables are added for driving turbulence by zwg

   logical:: enable_DT=.false.
   integer:: DT_mode
   integer:: enable_FFTW=0
   integer:: isRestart=0

   double precision,dimension(:),allocatable::Den_ini  !! zwg added for checking R2C and C2R interfaces of FFTW
   real(C_DOUBLE), dimension(:),pointer:: Den_R2C_in   !! zwg added for checking R2C and C2R interfaces of FFTW
   complex(C_DOUBLE_COMPLEX),dimension(:),pointer::Den_R2C_out  !! zwg added for checking R2C and C2R interfaces of FFTW
   real(C_DOUBLE), dimension(:),pointer:: Den_C2R_out   !! zwg added for checking R2C and C2R interfaces of FFTW
   type(C_PTR)::Den_Plan_R2C,Den_Plan_C2R  !! zwg added for checking R2C and C2R interfaces of FFTW

   complex(C_DOUBLE_COMPLEX),dimension(:),pointer::DTVx_C2R_in,DTVy_C2R_in,DTVz_C2R_in
   real(C_DOUBLE), dimension(:),pointer:: DTVx_C2R_out,DTVy_C2R_out,DTVz_C2R_out

   type(C_PTR)::DTVx_Plan_C2R,DTVy_Plan_C2R,DTVz_Plan_C2R

   double precision::DTenergyfaction
   double precision,dimension(:),allocatable::DTVx,DTVy,DTVz 
   double precision::drivingWN_DT,energy_DT,zeta_DT  !!! zwg : this variables are setted in subroutine init2d_for_FFTW of file init2d.f03 or init3D.f03
   double precision::netmomx_DT,netmomy_DT,netmomz_DT


   contains
      procedure::enableAD
      procedure::setADparams
      procedure::setGravConst
      procedure::setSelfgravityKernel
      procedure::enableSelfgravity
      procedure::setTopologyMPI
      procedure::setGridID
      procedure::setMesh
      procedure::setVariable
      procedure::setMPIWindows
      procedure::setSGMPIWindows
      procedure::setEoS 
      procedure::setSoundSpeed     
      procedure::setCFL
      procedure::setTime
      procedure::setRestart
      procedure::setSlopeLimiter
      procedure::setSolverType
      procedure::setBoundaryType
      procedure::initVariable
      procedure::exchangeBdryMPI
      procedure::setBoundary
      procedure::calcSelfgravity
      procedure::griddt
      procedure::writeGrid 
      procedure::evolveGridRK2
      procedure::setAdiGamma
      procedure::setSgBdryType

      !!! zwg added for driving turbulence
      procedure::enableDrivingTurbulence
      procedure::setDTenergyfaction
      procedure::calcDrivingTurbulence
      procedure::calcDrivingTurbulence_MD 
      procedure::init_For_FFTW 

      !!! subroutines added by zwg
      procedure:: writeGrid_HD_vtk
      procedure:: writeGrid_MHD_vtk
      procedure:: writeGrid_tecplot

      procedure:: readGrid_HD_vtk
      procedure:: readGrid_MHD_vtk
end type grid

contains

subroutine enableAD(this,enable_ad)
class(grid)::this
logical::enable_ad
this%enable_ad=enable_ad
end subroutine enableAD

subroutine setADparams(this,mu_ad,alpha_ad)
class(grid)::this
double precision::mu_ad,alpha_ad
this%mu_ad=mu_ad
this%alpha_ad=alpha_ad
end subroutine setADparams

subroutine setGravConst(this,GravConst)
class(grid)::this
double precision::GravConst
this%GravConst=GravConst
end subroutine setGravConst

subroutine setSgBdryType(this,sgBdryType)
class(grid)::this
integer::sgBdryType
this%sgBdryType=sgBdryType  !! 0=isolated [default], 1=periodic
end subroutine

subroutine calcSelfgravity(this,q)
class(grid)::this
double precision,dimension(:)::q
integer::ndim
ndim=this%ndim

if(ndim .eq. 2) then

 if(this%sgBdryType .eq. 0 ) then  !!! zwg added for poisson solver with 2D periodic BCs
  call calcSG2D(this,q,this%sgDensityBuffer,this%sgDenCmplx,this%sgfxKernel,this%sgfyKernel, &
       this%sgfxCmplx,this%sgfyCmplx,this%sgfx,this%sgfy)
 elseif(this%sgBdryType .eq. 1) then  !!zwg added for poisson solver with 2D periodic BCs
  call calcSG2Dperiodic(this,q,this%sgDenCmplx,this%sgfxCmplx,this%sgfyCmplx,this%sgfx,this%sgfy)  !!!zwg added for poisson solver with 2D periodic BCs
 endif   !!!zwg added for poisson solver with 2D periodic BCs

elseif(ndim .eq. 3) then
  if(this%sgBdryType .eq. 0 ) then
    call calcSG3D(this,q,this%sgDensityBuffer,this%sgDenCmplx,this%sgfxKernel,this%sgfyKernel, &
    this%sgfzKernel,this%sgfxCmplx,this%sgfyCmplx,this%sgfzCmplx,this%sgfx,this%sgfy,this%sgfz)
  elseif(this%sgBdryType .eq. 1) then
    call calcSG3Dperiodic(this,q,this%sgDenCmplx,this%sgfxCmplx,this%sgfyCmplx,this%sgfzCmplx, &
                          this%sgfx,this%sgfy,this%sgfz)
  endif
else
  print *,"gridModule.f03: selfgravity for",ndim,"D problems has not yet implemented...."
  stop
endif
end subroutine calcSelfgravity



subroutine calcDrivingTurbulence(this,q)

class(grid)::this
double precision,dimension(:)::q
integer::ndim
ndim=this%ndim


if(ndim .eq. 2) then
call calcDT2D(this,this%Den_ini,this%Den_R2C_in,this%Den_R2C_out,this%Den_C2R_out, & 
            & this%DTVx_C2R_in,this%DTVx_C2R_out,this%DTVy_C2R_in,this%DTVy_C2R_out,this%DTVx,this%DTVy)


elseif(ndim .eq. 3) then

call calcDT3D(this,this%Den_ini,this%Den_R2C_in,this%Den_R2C_out,this%Den_C2R_out, this%DTVx_C2R_in,this%DTVx_C2R_out, &
            & this%DTVy_C2R_in,this%DTVy_C2R_out,this%DTVz_C2R_in,this%DTVz_C2R_out,this%DTVx,this%DTVy,this%DTVz)

endif

end subroutine calcDrivingTurbulence



subroutine calcDrivingTurbulence_MD(this,q)

class(grid)::this
double precision,dimension(:)::q
integer::ndim
ndim=this%ndim


if(ndim .eq. 2) then
  call calcDT2D_MD(this,this%DTVx_C2R_in,this%DTVx_C2R_out, &
                      & this%DTVy_C2R_in,this%DTVy_C2R_out, &
                      & this%DTVx,this%DTVy,this%q)
elseif(ndim .eq. 3) then

  call calcDT3D_MD(this,this%DTVx_C2R_in,this%DTVx_C2R_out, &
                      & this%DTVy_C2R_in,this%DTVy_C2R_out, &
                      & this%DTVz_C2R_in,this%DTVz_C2R_out, &
                      & this%DTVx,this%DTVy,this%DTVz,this%q)
endif

end subroutine calcDrivingTurbulence_MD




subroutine setSGMPIWindows(this)
class(grid)::this
if(this%enable_sg) then
  if(this%ndim .eq. 2) then
    call initSGWindows2D(this,this%sgDensityBuffer)
  endif
  if(this%ndim .eq. 3) then
    call initSGWindows3D(this,this%sgDensityBuffer)
  endif
else
  print *,"gridModule.f03: selfgravity is not enabled..."
endif
end subroutine setSGMPIWindows

subroutine setSelfgravityKernel(this)
class(grid)::this
if(this%ndim .eq. 2)then
  call sgkernel2d(this,this%sgfxKernel,this%sgfyKernel)
endif
if(this%ndim .eq. 3)then
  if(this%sgBdryType .eq. 0) then
    call sgkernel3d(this,this%sgfxKernel,this%sgfyKernel,this%sgfzKernel)
  endif
endif
end subroutine setSelfgravityKernel

subroutine enableSelfgravity(this)
class(grid)::this
this%enable_sg=.true.

this%sgBdryType=0  !! isolated boundary condition by default
print *,"gridModule.f03: selfgravity has been enabled ....."

if(this%enable_FFTW==0) then
call fftw_mpi_init()
!!! zwg added for no more calling fftw_mpi_init() in subroutine enableDrivingTurbulence
this%enable_FFTW=1     
endif  

end subroutine enableSelfgravity


!!!!! this subroutine is added by zwg to drive turbulence for Scorpio
subroutine enableDrivingTurbulence(this,DT_mode)

class(grid)::this
integer:: DT_mode
this%enable_DT=.true.
this%DT_mode=DT_mode
!!! for drive turbulence case, the boundaryType should be periodic: boundaryType=3
!write(*,*)"zwg: this%enable_FFTW before setting",this%enable_FFTW

if(this%enable_FFTW==0) then
call fftw_mpi_init()
this%enable_FFTW=1
endif

!write(*,*)"zwg: this%enable_FFTW after setting",this%enable_FFTW

end subroutine enableDrivingTurbulence


!!subroutine setDTenergyfaction(this,t_accum_turb,dt_turb)
subroutine setDTenergyfaction(this,t_accum_turb,dt_turb)
class(grid)::this
  double precision::t_accum_turb, dt_turb

  this%DTenergyfaction=t_accum_turb/dt_turb

end subroutine setDTenergyfaction


subroutine exchangeBdryMPI(this,q,win)
class(grid)::this
double precision,dimension(:)::q
integer::win,ndim
ndim=this%ndim

!print *,"gridModule.f03: myid=",myid,"ndim=",ndim,"win=",win
if(ndim .eq. 1) then
   call exchgBdryMPI1D(this,q,win)
endif
if(ndim .eq. 2) then
   call exchgBdryMPI2D(this,q,win,this%winbuf1,this%databuf1,this%leftbuf,this%rightbuf)
endif
if(ndim .eq. 3) then
   call exchgBdryMPI3D(this,q,win,this%winbuf1,this%winbuf2,this%databuf1,this%databuf2,&
&                      this%leftbuf,this%rightbuf,this%frontbuf,this%backbuf)
endif

end subroutine exchangeBdryMPI

subroutine setMPIWindows(this)
class(grid)::this
integer::ndim

ndim=this%ndim
if(ndim .eq. 1) then
   call initMPIWindows1D(this,this%q,this%q1,this%q2)
endif
if(ndim .eq. 2) then
   call initMPIWindows2D(this,this%q,this%q1,this%q2,this%databuf1)
endif
if(ndim .eq. 3) then
   call initMPIWindows3D(this,this%q,this%q1,this%q2,this%databuf1,this%databuf2)
endif
end subroutine setMPIWindows



!!!!!!!!  zwg set mpi process for each blocks
!!!!!!!!
subroutine setTopologyMPI(this,ndim,dims,periods,reorder)

class(grid)::this

integer::ndim

integer::dims(ndim),vu,ierr

logical::periods(ndim),reorder

integer::mpiCoord(ndim),i,np

!!!! to get the total number of processes!!!!

np=dims(1)  !!! zwg : np is just a local index!!!!
do i=2,ndim !!! zwg :should be noted that here i should be starts with 2 and end with ndim
 np=np*dims(i) 
enddo  !!

allocate(this%dims_mpi(ndim))
this%dims_mpi=dims
!!! use all the arry for this%dims_mpi

call MPI_CART_CREATE(MPI_COMM_WORLD,ndim,dims,periods,reorder,vu,ierr) 

!! get a new communicator vu for the newly created MPI topology
!!! zwg : not MPI_COMM_WORLD anymore!!!

this%vu_mpi=vu
!!! zwg :use this new communicator vn for this%vn_mpi

call MPI_COMM_RANK(vu,myid,ierr)

!!! zwg : get the current id for present core 


if(np .eq. nprocs) then
   if(myid .eq. 0) then
     print *,"gridModule.f03: The MPI topology is self-consistent..."
   endif
else
   if(myid .eq. 0) then
     print *,"gridModule.f03: The MPI topolgy is not self-consistent..."
   endif
   stop
endif


!!!!! zwg :: the following routines is not well understood!!!
allocate(this%mpiCoord(ndim))

call MPI_CART_COORDS(vu,myid,ndim,mpiCoord,ierr)

!!!zwg: Determines process coords in cartesian topology given rank in group
!!!zwg : the exat mean of mpiCoord still not clear for me currently!!!!

this%mpiCoord=mpiCoord
!print *,"gridModule.f03: myid=",myid, "Coord=",mpiCoord  !!!, this%mpiCoord

 !stop



!!!  zwg : to get the process index for neighbors cores in x diraction
  call MPI_CART_SHIFT(this%vu_mpi,0,1,this%left_mpi,this%right_mpi,ierr)



!!!  zwg : to get the process index for neighbors cores in Y diraction
if(ndim .gt. 1) then
  call MPI_CART_SHIFT(this%vu_mpi,1,1,this%up_mpi,this%down_mpi,ierr)

endif

!!!  zwg : to get the process index for neighbors cores in z diraction
if(ndim .gt. 2) then
  call MPI_CART_SHIFT(this%vu_mpi,2,1,this%top_mpi,this%bottom_mpi,ierr)

endif

!!! zwg try to print the results, 
!!! the results shows that it seems that the return values of subroutine MPI_CART_SHIFT should plus 1 for each process ID

!if(myid==0) then
  !print *,"gridModule.f03: myid=",myid,"left=",this%left_mpi,"right=",this%right_mpi
  !print *,"gridModule.f03: myid=",myid,"up=",this%up_mpi,"down=",this%down_mpi
  !print *,"gridModule.f03: myid",myid,"top=",this%top_mpi,"botton=",this%bottom_mpi
!endif

 !!stop  !!! zwg :: added for debugging
end subroutine setTopologyMPI



subroutine setMPI(np)
integer::np
nprocs=np
end subroutine setMPI



subroutine setGridID(this,gridID)
class(grid)::this
integer::gridID
this%gridID=gridID
end subroutine setGridID



subroutine setAdiGamma(this,gam)
class(grid)::this
double precision::gam
this%adiGamma=gam
end subroutine setAdiGamma



subroutine setBoundaryType(this,boundaryType)
class(grid)::this
integer::boundaryType
this%boundaryType=boundaryType
end subroutine setBoundaryType



subroutine evolveGridRK2(this)
class(grid)::this
integer::ndim
ndim=this%ndim
if(ndim .eq. 1) then
   call rk2_1D(this,this%q,this%q1,this%q2)
endif
if(ndim .eq. 2) then
   call rk2_2D(this,this%q,this%q1,this%q2)
endif
if(ndim .eq. 3) then
   call rk2_3D(this,this%q,this%q1,this%q2)
endif
end subroutine  



!!!!output data using vtk formate added by zwg
!!!! just for 2D and 3D cases

subroutine writeGrid_MHD_vtk(this)  
class(grid)::this
integer:: ndim

ndim=this%ndim


 if(ndim==2)then
  !call output2d_MHD_vtk(this,this%q)  !!! not use anymore

   call output2d_MHD_pvts(this,this%q)
  !call output2d_MHD_pvtu(this,this%q)

 elseif(ndim==3)then
   call output3d_MHD_pvts(this,this%q)
  !call output3d_MHD_pvtu(this,this%q)

 endif !! end of if(ndim==2)then

end subroutine writeGrid_MHD_vtk


subroutine writeGrid_HD_vtk(this)  
class(grid)::this
integer:: ndim

ndim=this%ndim


 if(ndim==2)then
  !call output2d_MHD_vtk(this,this%q)   !! not use anymore

   call output2d_HD_pvts(this,this%q)   !! Structured  output, already implemented
  !call output2d_HD_pvtu(this,this%q)   !! Unstructured output, not implemented yet

 elseif(ndim==3)then
   call output3d_HD_pvts(this,this%q)
  !call output3d_HD_pvtu(this,this%q)

 endif !! end of if(ndim==2)then



end subroutine writeGrid_HD_vtk



subroutine readGrid_MHD_vtk(this)  
class(grid)::this
integer:: ndim

ndim=this%ndim


 if(ndim==2)then
  
   !call input2d_MHD_pvts(this,this%q)
  
 elseif(ndim==3)then

   call input3d_MHD_pvts(this,this%q)

 endif !! end of if(ndim==2)then

end subroutine readGrid_MHD_vtk


subroutine readGrid_HD_vtk(this)  
class(grid)::this
integer:: ndim

ndim=this%ndim


 if(ndim==2)then
  

   !call input2d_HD_pvts(this,this%q)   !! Structured  output, already implemented
  

 elseif(ndim==3)then
   call input3d_HD_pvts(this,this%q)
  

 endif !! end of if(ndim==2)then



end subroutine readGrid_HD_vtk





!!!!output data using tecplot formate added by zwg
!!!! just for 2D and 3D cases
subroutine writeGrid_tecplot(this)  
class(grid)::this



end subroutine writeGrid_tecplot



subroutine writeGrid(this)
class(grid)::this
integer::ntotal,fnum,gridID,nx,ny,nz,nbuf,nvar,ndim,i

character(len=13)::flnm
character(len=4)::gridnum
character(len=4)::filenum

character(len=20)::dsetname
character(len=5)::gridname
character(len=2)::varnum  !!! zwg no use!!!
double precision:: variable(8)

integer::xstart,xend
integer::ystart,yend
integer::zstart,zend

integer::varindex

fnum=this%fnum

gridID=this%gridID

ndim=this%ndim

nvar=this%nvar
nbuf=this%nbuf
nx=this%nMesh(1)
variable=dble(this%variable)
if(ndim .gt. 1) then
   ny=this%nMesh(2)
endif
if(ndim .gt. 2) then
   nz=this%nMesh(3)
endif


!!! zwg :: it seems that ntotal is useless
if(ndim .eq. 1) then
  ntotal=(nx+2*nbuf)*nvar
elseif(ndim .eq. 2) then
  ntotal=(nx+2*nbuf)*(ny+2*nbuf)*nvar
elseif(ndim .eq. 3) then
  ntotal=(nx+2*nbuf)*(ny+2*nbuf)*(nz+2*nbuf)*nvar
endif

write(gridnum,'(I4.4)') gridID
write(filenum,'(I4.4)') fnum
flnm='g'//gridnum//'_'//filenum//'.h5'
gridname='g'//gridnum

if(myid .eq. 0)then
  print *,"gridModule.f03: writing data to ",flnm
endif

dsetname='nMesh'
call output1d(this,0,dble(this%nMesh_global),this%ndim,this%ndim,0,0,flnm,dsetname,1)
dsetname='nbuf'
call output1d(this,0,dble(this%nbuf),1,1,0,0,flnm,dsetname,0)
dsetname='coordType'
call output1d(this,0,dble(this%coordType),1,1,0,0,flnm,dsetname,0)
dsetname='gridID'
call output1d(this,0,dble(this%gridID),1,1,0,0,flnm,dsetname,0)
dsetname='leftBdry'
call output1d(this,0,this%leftBdry_global,this%ndim,this%ndim,0,0,flnm,dsetname,0)
dsetname='rightBdry'
call output1d(this,0,this%rightBdry_global,this%ndim,this%ndim,0,0,flnm,dsetname,0)
dsetname='variable'
call output1d(this,0,dble(this%variable),8,8,0,0,flnm,dsetname,0)
dsetname='nvar'
call output1d(this,0,dble(this%nvar),1,1,0,0,flnm,dsetname,0)
dsetname='ndim'
call output1d(this,0,dble(this%ndim),1,1,0,0,flnm,dsetname,0)
dsetname='t'
call output1d(this,0,this%t,1,1,0,0,flnm,dsetname,0)
dsetname='CFL'
call output1d(this,0,this%CFL,1,1,0,0,flnm,dsetname,0)
dsetname='snd'
call output1d(this,0,this%snd,1,1,0,0,flnm,dsetname,0)
dsetname='limiterType'
call output1d(this,0,dble(this%limiterType),1,1,0,0,flnm,dsetname,0)
dsetname='solverType'
call output1d(this,0,dble(this%solverType),1,1,0,0,flnm,dsetname,0)
dsetname='eosType'
call output1d(this,0,dble(this%eosType),1,1,0,0,flnm,dsetname,0)
dsetname='adiGamma'
call output1d(this,0,this%adiGamma,1,1,0,0,flnm,dsetname,0)
dsetname='boundaryType'
call output1d(this,0,dble(this%boundaryTYpe),1,1,0,0,flnm,dsetname,0)

dsetname='xl1'
call output1d(this,1,this%xl(1)%coords,this%nMesh_global(1),this%nMesh(1),this%nbuf,2,flnm,dsetname,0)
dsetname='xr1'
call output1d(this,1,this%xr(1)%coords,this%nMesh_global(1),this%nMesh(1),this%nbuf,2,flnm,dsetname,0)
dsetname='xc1'
call output1d(this,1,this%xc(1)%coords,this%nMesh_global(1),this%nMesh(1),this%nbuf,2,flnm,dsetname,0)
dsetname='dx1'
call output1d(this,1,this%dx(1)%coords,this%nMesh_global(1),this%nMesh(1),this%nbuf,2,flnm,dsetname,0)


if(this%ndim .gt. 1) then
dsetname='xl2'
call output1d(this,2,this%xl(2)%coords,this%nMesh_global(2),this%nMesh(2),this%nbuf,2,flnm,dsetname,0)
dsetname='xr2'
call output1d(this,2,this%xr(2)%coords,this%nMesh_global(2),this%nMesh(2),this%nbuf,2,flnm,dsetname,0)
dsetname='xc2'
call output1d(this,2,this%xc(2)%coords,this%nMesh_global(2),this%nMesh(2),this%nbuf,2,flnm,dsetname,0)
dsetname='dx2'
call output1d(this,2,this%dx(2)%coords,this%nMesh_global(2),this%nMesh(2),this%nbuf,2,flnm,dsetname,0)
endif

if(this%ndim .gt. 2) then
dsetname='xl3'
call output1d(this,3,this%xl(3)%coords,this%nMesh_global(3),this%nMesh(3),this%nbuf,2,flnm,dsetname,0)
dsetname='xr3'
call output1d(this,3,this%xr(3)%coords,this%nMesh_global(3),this%nMesh(3),this%nbuf,2,flnm,dsetname,0)
dsetname='xc3'
call output1d(this,3,this%xc(3)%coords,this%nMesh_global(3),this%nMesh(3),this%nbuf,2,flnm,dsetname,0)
dsetname='dx3'
call output1d(this,3,this%dx(3)%coords,this%nMesh_global(3),this%nMesh(3),this%nbuf,2,flnm,dsetname,0)
endif

if(this%ndim .eq. 1) then
   varindex=0
   if(variable(1) .eq. 1) then
     varindex=varindex+1
     dsetname='den'
     xstart=1+(varindex-1)*(nx+2*nbuf)
     xend=xstart+nx+2*nbuf-1
     call output1d(this,1,this%q(xstart:xend),this%nMesh_global(1),this%nMesh(1),this%nbuf,2,flnm,dsetname,0)
   endif
   
   if(variable(2) .eq. 1)then
     varindex=varindex+1
     dsetname='momx'
     xstart=1+(varindex-1)*(nx+2*nbuf)
     xend=xstart+nx+2*nbuf-1
     call output1d(this,1,this%q(xstart:xend),this%nMesh_global(1),this%nMesh(1),this%nbuf,2,flnm,dsetname,0)
   endif

   if(variable(3) .eq. 1)then
     varindex=varindex+1
     dsetname='momy'
     xstart=1+(varindex-1)*(nx+2*nbuf)
     xend=xstart+nx+2*nbuf-1
     call output1d(this,1,this%q(xstart:xend),this%nMesh_global(1),this%nMesh(1),this%nbuf,2,flnm,dsetname,0)
   endif

   if(variable(4) .eq. 1)then
     varindex=varindex+1
     dsetname='momz'
     xstart=1+(varindex-1)*(nx+2*nbuf)
     xend=xstart+nx+2*nbuf-1
     call output1d(this,1,this%q(xstart:xend),this%nMesh_global(1),this%nMesh(1),this%nbuf,2,flnm,dsetname,0)
   endif

   if(variable(5) .eq. 1)then
     varindex=varindex+1
     dsetname='bxl'
     xstart=1+(varindex-1)*(nx+2*nbuf)
     xend=xstart+nx+2*nbuf-1
     call output1d(this,1,this%q(xstart:xend),this%nMesh_global(1),this%nMesh(1),this%nbuf,2,flnm,dsetname,0)
   endif

   if(variable(6) .eq. 1)then
     varindex=varindex+1
     dsetname='byl'
     xstart=1+(varindex-1)*(nx+2*nbuf)
     xend=xstart+nx+2*nbuf-1
     call output1d(this,1,this%q(xstart:xend),this%nMesh_global(1),this%nMesh(1),this%nbuf,2,flnm,dsetname,0)
   endif

   if(variable(7) .eq. 1)then
     varindex=varindex+1
     dsetname='bzl'
     xstart=1+(varindex-1)*(nx+2*nbuf)
     xend=xstart+nx+2*nbuf-1
     call output1d(this,1,this%q(xstart:xend),this%nMesh_global(1),this%nMesh(1),this%nbuf,2,flnm,dsetname,0)
   endif

   if(variable(8) .eq. 1)then
     varindex=varindex+1
     dsetname='ene'
     xstart=1+(varindex-1)*(nx+2*nbuf)
     xend=xstart+nx+2*nbuf-1
     call output1d(this,1,this%q(xstart:xend),this%nMesh_global(1),this%nMesh(1),this%nbuf,2,flnm,dsetname,0)
   endif

   if(variable(5) .eq. 1)then
     varindex=varindex+1
     dsetname='bxr'
     xstart=1+(varindex-1)*(nx+2*nbuf)
     xend=xstart+nx+2*nbuf-1
     call output1d(this,1,this%q(xstart:xend),this%nMesh_global(1),this%nMesh(1),this%nbuf,2,flnm,dsetname,0)
   endif

   if(variable(6) .eq. 1)then
     varindex=varindex+1
     dsetname='byr'
     xstart=1+(varindex-1)*(nx+2*nbuf)
     xend=xstart+nx+2*nbuf-1
     call output1d(this,1,this%q(xstart:xend),this%nMesh_global(1),this%nMesh(1),this%nbuf,2,flnm,dsetname,0)
   endif

   if(variable(7) .eq. 1)then
     varindex=varindex+1
     dsetname='bzr'
     xstart=1+(varindex-1)*(nx+2*nbuf)
     xend=xstart+nx+2*nbuf-1
     call output1d(this,1,this%q(xstart:xend),this%nMesh_global(1),this%nMesh(1),this%nbuf,2,flnm,dsetname,0)
   endif
endif  !!! zwg : endof if(this%ndim .eq. 1) then

if(this%ndim .eq. 2) then
   varindex=0
   if(variable(1) .eq. 1) then
     varindex=varindex+1
     dsetname='den'
     xstart=1+(varindex-1)*(nx+2*nbuf)*(ny+2*nbuf)
     xend  =varindex*(nx+2*nbuf)*(ny+2*nbuf)
     call output2d(this,this%q(xstart:xend),this%nMesh_global(1),this%nMesh_global(2),this%nMesh(1)& 
     ,this%nMesh(2),this%nbuf,this%nbuf,2,2,flnm,dsetname,0) 

         !output2d(this,q_loc,nx_global,ny_global,nx_loc,ny_loc,inputbufx,inputbufy,bufx,bufy,flnm,dsetname,flag) 

   endif

   if(variable(2) .eq. 1) then
     varindex=varindex+1
     dsetname='momx'
     xstart=1+(varindex-1)*(nx+2*nbuf)*(ny+2*nbuf)
     xend  =varindex*(nx+2*nbuf)*(ny+2*nbuf)
     call output2d(this,this%q(xstart:xend),this%nMesh_global(1),this%nMesh_global(2),this%nMesh(1)& 
     ,this%nMesh(2),this%nbuf,this%nbuf,2,2,flnm,dsetname,0)  
   endif

   if(variable(3) .eq. 1) then
     varindex=varindex+1
     dsetname='momy'
     xstart=1+(varindex-1)*(nx+2*nbuf)*(ny+2*nbuf)
     xend  =varindex*(nx+2*nbuf)*(ny+2*nbuf)
     call output2d(this,this%q(xstart:xend),this%nMesh_global(1),this%nMesh_global(2),this%nMesh(1)& 
     ,this%nMesh(2),this%nbuf,this%nbuf,2,2,flnm,dsetname,0)  
   endif

   if(variable(4) .eq. 1) then
     varindex=varindex+1
     dsetname='momz'
     xstart=1+(varindex-1)*(nx+2*nbuf)*(ny+2*nbuf)
     xend  =varindex*(nx+2*nbuf)*(ny+2*nbuf)
     call output2d(this,this%q(xstart:xend),this%nMesh_global(1),this%nMesh_global(2),this%nMesh(1)& 
     ,this%nMesh(2),this%nbuf,this%nbuf,2,2,flnm,dsetname,0)  
   endif

   if(variable(5) .eq. 1) then
     varindex=varindex+1
     dsetname='bxl'
     xstart=1+(varindex-1)*(nx+2*nbuf)*(ny+2*nbuf)
     xend  =varindex*(nx+2*nbuf)*(ny+2*nbuf)
     call output2d(this,this%q(xstart:xend),this%nMesh_global(1),this%nMesh_global(2),this%nMesh(1)& 
     ,this%nMesh(2),this%nbuf,this%nbuf,2,2,flnm,dsetname,0)  
   endif

   if(variable(6) .eq. 1) then
     varindex=varindex+1
     dsetname='byl'
     xstart=1+(varindex-1)*(nx+2*nbuf)*(ny+2*nbuf)
     xend  =varindex*(nx+2*nbuf)*(ny+2*nbuf)
     call output2d(this,this%q(xstart:xend),this%nMesh_global(1),this%nMesh_global(2),this%nMesh(1)& 
     ,this%nMesh(2),this%nbuf,this%nbuf,2,2,flnm,dsetname,0)  
   endif

   if(variable(7) .eq. 1) then
     varindex=varindex+1
     dsetname='bzl'
     xstart=1+(varindex-1)*(nx+2*nbuf)*(ny+2*nbuf)
     xend  =varindex*(nx+2*nbuf)*(ny+2*nbuf)
     call output2d(this,this%q(xstart:xend),this%nMesh_global(1),this%nMesh_global(2),this%nMesh(1)& 
     ,this%nMesh(2),this%nbuf,this%nbuf,2,2,flnm,dsetname,0)  
   endif

   if(variable(8) .eq. 1) then
     varindex=varindex+1
     dsetname='ene'
     xstart=1+(varindex-1)*(nx+2*nbuf)*(ny+2*nbuf)
     xend  =varindex*(nx+2*nbuf)*(ny+2*nbuf)
     call output2d(this,this%q(xstart:xend),this%nMesh_global(1),this%nMesh_global(2),this%nMesh(1)& 
     ,this%nMesh(2),this%nbuf,this%nbuf,2,2,flnm,dsetname,0)  
   endif

   if(variable(5) .eq. 1) then
     varindex=varindex+1
     dsetname='bxr'
     xstart=1+(varindex-1)*(nx+2*nbuf)*(ny+2*nbuf)
     xend  =varindex*(nx+2*nbuf)*(ny+2*nbuf)
     call output2d(this,this%q(xstart:xend),this%nMesh_global(1),this%nMesh_global(2),this%nMesh(1)& 
     ,this%nMesh(2),this%nbuf,this%nbuf,2,2,flnm,dsetname,0)  
   endif

   if(variable(6) .eq. 1) then
     varindex=varindex+1
     dsetname='byr'
     xstart=1+(varindex-1)*(nx+2*nbuf)*(ny+2*nbuf)
     xend  =varindex*(nx+2*nbuf)*(ny+2*nbuf)
     call output2d(this,this%q(xstart:xend),this%nMesh_global(1),this%nMesh_global(2),this%nMesh(1)& 
     ,this%nMesh(2),this%nbuf,this%nbuf,2,2,flnm,dsetname,0)  
   endif

   if(variable(7) .eq. 1) then
     varindex=varindex+1
     dsetname='bzr'
     xstart=1+(varindex-1)*(nx+2*nbuf)*(ny+2*nbuf)
     xend  =varindex*(nx+2*nbuf)*(ny+2*nbuf)
     call output2d(this,this%q(xstart:xend),this%nMesh_global(1),this%nMesh_global(2),this%nMesh(1)& 
     ,this%nMesh(2),this%nbuf,this%nbuf,2,2,flnm,dsetname,0)  
   endif

endif !! zwg : if(this%ndim .eq. 2) then

if(this%ndim .eq. 3) then
   varindex=0
   if(variable(1) .eq. 1) then
     varindex=varindex+1
     dsetname='den'
     xstart=1+(varindex-1)*(nx+2*nbuf)*(ny+2*nbuf)*(nz+2*nbuf)
     xend  =varindex*(nx+2*nbuf)*(ny+2*nbuf)*(nz+2*nbuf)
     call output3d(this,this%q(xstart:xend),this%nMesh_global(1),this%nMesh_global(2),this%nMesh_global(3),this%nMesh(1) &
     ,this%nMesh(2),this%nMesh(3),this%nbuf,this%nbuf,this%nbuf,2,2,2,flnm,dsetname,0)
   endif

   if(variable(2) .eq. 1) then
     varindex=varindex+1
     dsetname='momx'
     xstart=1+(varindex-1)*(nx+2*nbuf)*(ny+2*nbuf)*(nz+2*nbuf)
     xend  =varindex*(nx+2*nbuf)*(ny+2*nbuf)*(nz+2*nbuf)
     call output3d(this,this%q(xstart:xend),this%nMesh_global(1),this%nMesh_global(2),this%nMesh_global(3),this%nMesh(1) &
     ,this%nMesh(2),this%nMesh(3),this%nbuf,this%nbuf,this%nbuf,2,2,2,flnm,dsetname,0)
   endif   

   if(variable(3) .eq. 1) then
     varindex=varindex+1
     dsetname='momy'
     xstart=1+(varindex-1)*(nx+2*nbuf)*(ny+2*nbuf)*(nz+2*nbuf)
     xend  =varindex*(nx+2*nbuf)*(ny+2*nbuf)*(nz+2*nbuf)
     call output3d(this,this%q(xstart:xend),this%nMesh_global(1),this%nMesh_global(2),this%nMesh_global(3),this%nMesh(1) &
     ,this%nMesh(2),this%nMesh(3),this%nbuf,this%nbuf,this%nbuf,2,2,2,flnm,dsetname,0)
   endif

   if(variable(4) .eq. 1) then
     varindex=varindex+1
     dsetname='momz'
     xstart=1+(varindex-1)*(nx+2*nbuf)*(ny+2*nbuf)*(nz+2*nbuf)
     xend  =varindex*(nx+2*nbuf)*(ny+2*nbuf)*(nz+2*nbuf)
     call output3d(this,this%q(xstart:xend),this%nMesh_global(1),this%nMesh_global(2),this%nMesh_global(3),this%nMesh(1) &
     ,this%nMesh(2),this%nMesh(3),this%nbuf,this%nbuf,this%nbuf,2,2,2,flnm,dsetname,0)
   endif

   if(variable(5) .eq. 1) then
     varindex=varindex+1
     dsetname='bxl'
     xstart=1+(varindex-1)*(nx+2*nbuf)*(ny+2*nbuf)*(nz+2*nbuf)
     xend  =varindex*(nx+2*nbuf)*(ny+2*nbuf)*(nz+2*nbuf)
     call output3d(this,this%q(xstart:xend),this%nMesh_global(1),this%nMesh_global(2),this%nMesh_global(3),this%nMesh(1) &
     ,this%nMesh(2),this%nMesh(3),this%nbuf,this%nbuf,this%nbuf,2,2,2,flnm,dsetname,0)
   endif

   if(variable(6) .eq. 1) then
     varindex=varindex+1
     dsetname='byl'
     xstart=1+(varindex-1)*(nx+2*nbuf)*(ny+2*nbuf)*(nz+2*nbuf)
     xend  =varindex*(nx+2*nbuf)*(ny+2*nbuf)*(nz+2*nbuf)
     call output3d(this,this%q(xstart:xend),this%nMesh_global(1),this%nMesh_global(2),this%nMesh_global(3),this%nMesh(1) &
     ,this%nMesh(2),this%nMesh(3),this%nbuf,this%nbuf,this%nbuf,2,2,2,flnm,dsetname,0)
   endif

   if(variable(7) .eq. 1) then
     varindex=varindex+1
     dsetname='bzl'
     xstart=1+(varindex-1)*(nx+2*nbuf)*(ny+2*nbuf)*(nz+2*nbuf)
     xend  =varindex*(nx+2*nbuf)*(ny+2*nbuf)*(nz+2*nbuf)
     call output3d(this,this%q(xstart:xend),this%nMesh_global(1),this%nMesh_global(2),this%nMesh_global(3),this%nMesh(1) &
     ,this%nMesh(2),this%nMesh(3),this%nbuf,this%nbuf,this%nbuf,2,2,2,flnm,dsetname,0)
   endif

   if(variable(8) .eq. 1) then
     varindex=varindex+1
     dsetname='ene'
     xstart=1+(varindex-1)*(nx+2*nbuf)*(ny+2*nbuf)*(nz+2*nbuf)
     xend  =varindex*(nx+2*nbuf)*(ny+2*nbuf)*(nz+2*nbuf)
     call output3d(this,this%q(xstart:xend),this%nMesh_global(1),this%nMesh_global(2),this%nMesh_global(3),this%nMesh(1) &
     ,this%nMesh(2),this%nMesh(3),this%nbuf,this%nbuf,this%nbuf,2,2,2,flnm,dsetname,0)
   endif


   if(variable(5) .eq. 1) then
     varindex=varindex+1
     dsetname='bxr'
     xstart=1+(varindex-1)*(nx+2*nbuf)*(ny+2*nbuf)*(nz+2*nbuf)
     xend  =varindex*(nx+2*nbuf)*(ny+2*nbuf)*(nz+2*nbuf)
     call output3d(this,this%q(xstart:xend),this%nMesh_global(1),this%nMesh_global(2),this%nMesh_global(3),this%nMesh(1) &
     ,this%nMesh(2),this%nMesh(3),this%nbuf,this%nbuf,this%nbuf,2,2,2,flnm,dsetname,0)
   endif

   if(variable(6) .eq. 1) then
     varindex=varindex+1
     dsetname='byr'
     xstart=1+(varindex-1)*(nx+2*nbuf)*(ny+2*nbuf)*(nz+2*nbuf)
     xend  =varindex*(nx+2*nbuf)*(ny+2*nbuf)*(nz+2*nbuf)
     call output3d(this,this%q(xstart:xend),this%nMesh_global(1),this%nMesh_global(2),this%nMesh_global(3),this%nMesh(1) &
     ,this%nMesh(2),this%nMesh(3),this%nbuf,this%nbuf,this%nbuf,2,2,2,flnm,dsetname,0)
   endif

   if(variable(7) .eq. 1) then
     varindex=varindex+1
     dsetname='bzr'
     xstart=1+(varindex-1)*(nx+2*nbuf)*(ny+2*nbuf)*(nz+2*nbuf)
     xend  =varindex*(nx+2*nbuf)*(ny+2*nbuf)*(nz+2*nbuf)
     call output3d(this,this%q(xstart:xend),this%nMesh_global(1),this%nMesh_global(2),this%nMesh_global(3),this%nMesh(1) &
     ,this%nMesh(2),this%nMesh(3),this%nbuf,this%nbuf,this%nbuf,2,2,2,flnm,dsetname,0)
   endif

endif !! zwg : if(this%ndim .eq. 2) then

if(this%enable_sg) then
  if(ndim .eq. 2) then
    dsetname='sgfx'
    call output2d(this,this%sgfx,this%nMesh_global(1),this%nMesh_global(2),this%nMesh(1),&
                this%nMesh(2), 0,0,0,0,flnm,dsetname,0)
    dsetname='sgfy'
    call output2d(this,this%sgfy,this%nMesh_global(1),this%nMesh_global(2),this%nMesh(1),&
                this%nMesh(2), 0,0,0,0,flnm,dsetname,0)
    dsetname='GravConst'
    call output1d(this,0,this%GravConst,1,1,0,0,flnm,dsetname,0)
  endif
  if(ndim .eq. 3) then
    dsetname='sgfx'
    call output3d(this,this%sgfx,this%nMesh_global(1),this%nMesh_global(2),this%nMesh_global(3),this%nMesh(1),&
                this%nMesh(2), this%nMesh(3),0,0,0,0,0,0,flnm,dsetname,0)
    dsetname='sgfy'
    call output3d(this,this%sgfy,this%nMesh_global(1),this%nMesh_global(2),this%nMesh_global(3),this%nMesh(1),&
                this%nMesh(2), this%nMesh(3),0,0,0,0,0,0,flnm,dsetname,0)
    dsetname='sgfz'
    call output3d(this,this%sgfz,this%nMesh_global(1),this%nMesh_global(2),this%nMesh_global(3),this%nMesh(1),&
                this%nMesh(2), this%nMesh(3),0,0,0,0,0,0,flnm,dsetname,0)
    dsetname='GravConst'
    call output1d(this,0,this%GravConst,1,1,0,0,flnm,dsetname,0)
  endif
endif
end subroutine writeGrid

subroutine griddt(this)
class(grid)::this
integer::ndim
ndim=this%ndim

if(ndim .eq. 1) then
    call dt1D(this,this%q)
endif
if(ndim .eq. 2) then
    call dt2D(this,this%q)
endif
if(ndim .eq. 3) then
    call dt3D(this,this%q)
endif
end subroutine griddt

subroutine setBoundary(this,q)
class(grid)::this
double precision,dimension(:)::q
integer::ndim
ndim=this%ndim
if(ndim .eq. 1) then
     call setBdry1D(this,q)
endif
if(ndim .eq. 2) then
   call setBdry2D(this,q)
endif
if(ndim .eq. 3) then
   call setBdry3D(this,q)
endif
end subroutine setBoundary

subroutine setSolverType(this,solverType)
class(grid)::this
integer::solverType
this%solverType=solverType
this%changeSolver=.false.
this%neg_pressure=.false.
end subroutine setSolverType

subroutine setSlopeLimiter(this,limiterType)
class(grid)::this
integer::limiterType
this%limiterType=limiterType
end subroutine setSlopeLimiter




subroutine setTime(this,fstart,tend,dtout)
class(grid)::this
integer::fstart,gridID,nbuf,coordType,ndim
double precision::valuetmp,variabletmp(8)
integer::variable(8),ntotal,nx,ny,nz,nvar
double precision::tend,dtout
double precision,dimension(:),allocatable::nMeshdble,leftBdry,rightBdry
integer,dimension(:),allocatable::nMesh
integer::varindex
integer::xstart,xend
integer::ystart,yend
integer::zstart,zend
character(len=13)::flnm
character(len=4)::gridnum
character(len=4)::filenum
character(len=20)::dsetname
character(len=5)::gridname

this%fstart=fstart
this%tend=tend
this%dtout=dtout
this%toutput=dtout

if(fstart .eq. 0) then
   this%t=0.d0
   this%fnum=0
   this%writeFlag=.false.
else
   gridID=this%gridID
   this%fnum=fstart
   this%writeFlag=.false.
   this%changeSolver=.false.
   write(gridnum,'(I4.4)') gridID
   write(filenum,'(I4.4)') fstart
   flnm='g'//gridnum//'_'//filenum//'.h5'
   gridname='g'//gridnum
   print *,"gridModule.f03: reading data from ",flnm

   dsetname='ndim'
   call read1d(this,0,valuetmp,1,1,0,0,flnm,dsetname)   
   ndim=int(valuetmp)

   dsetname='nbuf'
   call read1d(this,0,valuetmp,1,1,0,0,flnm,dsetname)
   nbuf=int(valuetmp)

   dsetname='coordType'
   call read1d(this,0,valuetmp,1,1,0,0,flnm,dsetname)
   coordType=int(valuetmp)

   allocate(nMesh(ndim),nMeshdble(ndim),leftBdry(ndim),rightBdry(ndim))

   dsetname='nMesh'
   call read1d(this,0,nMeshdble,ndim,ndim,0,0,flnm,dsetname)
   nMesh=int(nMeshdble)

   dsetname='leftBdry'
   call read1d(this,0,leftBdry,ndim,ndim,0,0,flnm,dsetname)

   dsetname='rightBdry'
   call read1d(this,0,rightBdry,ndim,ndim,0,0,flnm,dsetname)

   call this%setMesh(nMesh,leftBdry,rightBdry,nbuf,coordType,gridID)
   print *,"myid=",myid,"gridID=",this%gridID
   print *,"myid=",myid,"ndim=",this%ndim
   print *,"myid=",myid,"nbuf=",this%nbuf
   print *,"myid=",myid,"coordType=",this%coordType
   print *,"myid=",myid,"nMesh=",this%nMesh
   print *,"myid=",myid,"leftBdry=",this%leftBdry
   print *,"myid=",myid,"rightBdry=",this%rightBdry

   nx=this%nMesh(1)
   if(ndim .gt. 1) then
     ny=this%nMesh(2)
   endif
   if(ndim .gt. 2) then
     nz=this%nMesh(3)
   endif

   dsetname='variable'
   call read1d(this,0,variabletmp,8,8,0,0,flnm,dsetname)
   variable = int(variabletmp)
   call this%setVariable(variable)
   print *,"myid=",myid,"variable=",this%variable
   print *,"myid=",myid,"nvar=",this%nvar
   nvar=this%nvar

   call this%setMPIWindows()

   dsetname='t'
   call read1d(this,0,valuetmp,1,1,0,0,flnm,dsetname)
   this%t=valuetmp
   print *,"myid=",myid,"t=",this%t
   this%toutput=this%t+dtout

   dsetname='CFL'
   call read1d(this,0,valuetmp,1,1,0,0,flnm,dsetname)
   this%CFL=valuetmp
   print *,"myid=",myid,"CFL=",this%CFL
  
   dsetname='snd'
   call read1d(this,0,valuetmp,1,1,0,0,flnm,dsetname)
   this%snd=valuetmp
   print *,"myid=",myid,"snd=", this%snd

   dsetname='limiterType'
   call read1d(this,0,valuetmp,1,1,0,0,flnm,dsetname)
   this%limiterType=int(valuetmp)
   print *,"myid=",myid,"limiterTYpe=",this%limiterType

   
   dsetname='solverType'
   call read1d(this,0,valuetmp,1,1,0,0,flnm,dsetname)
   this%solverType=int(valuetmp)
   print *,"myid=",myid,"solverType=",this%solverType

   dsetname='eosType'
   call read1d(this,0,valuetmp,1,1,0,0,flnm,dsetname)
   this%eosType=int(valuetmp)
   print *,"myid=",myid,"eosType=",this%eosType

   dsetname='boundaryType'
   call read1d(this,0,valuetmp,1,1,0,0,flnm,dsetname)
   this%boundaryType=int(valuetmp)
   print *,"myid=",myid,"boundaryType=",this%boundaryType

   dsetname='adiGamma'
   call read1d(this,0,valuetmp,1,1,0,0,flnm,dsetname)
   this%adiGamma=valuetmp
   print *,"myid=",myid,"gamma=",this%adiGamma

   if(ndim .eq. 1) then
     print *,"restart is not available for 1D problems....please run with a fresh start.."
     stop
   endif

   if(ndim .eq. 2) then
     varindex=0
     if(variable(1) .eq. 1) then
       varindex=varindex+1
       dsetname='den'
       xstart=1+(varindex-1)*(nx+2*nbuf)*(ny+2*nbuf)
       xend  =varindex*(nx+2*nbuf)*(ny+2*nbuf)
       call read2d(this,this%q(xstart:xend),this%nMesh_global(1),this%nMesh_global(2),this%nMesh(1)&
       ,this%nMesh(2),this%nbuf,this%nbuf,2,2,flnm,dsetname)
     endif

     if(variable(2) .eq. 1) then
       varindex=varindex+1
       dsetname='momx'
       xstart=1+(varindex-1)*(nx+2*nbuf)*(ny+2*nbuf)
       xend  =varindex*(nx+2*nbuf)*(ny+2*nbuf)
       call read2d(this,this%q(xstart:xend),this%nMesh_global(1),this%nMesh_global(2),this%nMesh(1)&
       ,this%nMesh(2),this%nbuf,this%nbuf,2,2,flnm,dsetname)
     endif

     if(variable(3) .eq. 1) then
       varindex=varindex+1
       dsetname='momy'
       xstart=1+(varindex-1)*(nx+2*nbuf)*(ny+2*nbuf)
       xend  =varindex*(nx+2*nbuf)*(ny+2*nbuf)
       call read2d(this,this%q(xstart:xend),this%nMesh_global(1),this%nMesh_global(2),this%nMesh(1)&
       ,this%nMesh(2),this%nbuf,this%nbuf,2,2,flnm,dsetname)
     endif

     if(variable(4) .eq. 1) then
       varindex=varindex+1
       dsetname='momz'
       xstart=1+(varindex-1)*(nx+2*nbuf)*(ny+2*nbuf)
       xend  =varindex*(nx+2*nbuf)*(ny+2*nbuf)
       call read2d(this,this%q(xstart:xend),this%nMesh_global(1),this%nMesh_global(2),this%nMesh(1)&
       ,this%nMesh(2),this%nbuf,this%nbuf,2,2,flnm,dsetname)
     endif

     if(variable(5) .eq. 1) then
       varindex=varindex+1
       dsetname='bxl'
       xstart=1+(varindex-1)*(nx+2*nbuf)*(ny+2*nbuf)
       xend  =varindex*(nx+2*nbuf)*(ny+2*nbuf)
       call read2d(this,this%q(xstart:xend),this%nMesh_global(1),this%nMesh_global(2),this%nMesh(1)&
       ,this%nMesh(2),this%nbuf,this%nbuf,2,2,flnm,dsetname)
     endif

     if(variable(6) .eq. 1) then
       varindex=varindex+1
       dsetname='byl'
       xstart=1+(varindex-1)*(nx+2*nbuf)*(ny+2*nbuf)
       xend  =varindex*(nx+2*nbuf)*(ny+2*nbuf)
       call read2d(this,this%q(xstart:xend),this%nMesh_global(1),this%nMesh_global(2),this%nMesh(1)&
       ,this%nMesh(2),this%nbuf,this%nbuf,2,2,flnm,dsetname)
     endif

     if(variable(7) .eq. 1) then
       varindex=varindex+1
       dsetname='bzl'
       xstart=1+(varindex-1)*(nx+2*nbuf)*(ny+2*nbuf)
       xend  =varindex*(nx+2*nbuf)*(ny+2*nbuf)
       call read2d(this,this%q(xstart:xend),this%nMesh_global(1),this%nMesh_global(2),this%nMesh(1)&
       ,this%nMesh(2),this%nbuf,this%nbuf,2,2,flnm,dsetname)
     endif
 
     if(variable(8) .eq. 1) then
       varindex=varindex+1
       dsetname='ene'
       xstart=1+(varindex-1)*(nx+2*nbuf)*(ny+2*nbuf)
       xend  =varindex*(nx+2*nbuf)*(ny+2*nbuf)
       call read2d(this,this%q(xstart:xend),this%nMesh_global(1),this%nMesh_global(2),this%nMesh(1)&
       ,this%nMesh(2),this%nbuf,this%nbuf,2,2,flnm,dsetname)
     endif

     if(variable(5) .eq. 1) then
       varindex=varindex+1
       dsetname='bxr'
       xstart=1+(varindex-1)*(nx+2*nbuf)*(ny+2*nbuf)
       xend  =varindex*(nx+2*nbuf)*(ny+2*nbuf)
       call read2d(this,this%q(xstart:xend),this%nMesh_global(1),this%nMesh_global(2),this%nMesh(1)&
       ,this%nMesh(2),this%nbuf,this%nbuf,2,2,flnm,dsetname)
     endif

     if(variable(6) .eq. 1) then
       varindex=varindex+1
       dsetname='byr'
       xstart=1+(varindex-1)*(nx+2*nbuf)*(ny+2*nbuf)
       xend  =varindex*(nx+2*nbuf)*(ny+2*nbuf)
       call read2d(this,this%q(xstart:xend),this%nMesh_global(1),this%nMesh_global(2),this%nMesh(1)&
       ,this%nMesh(2),this%nbuf,this%nbuf,2,2,flnm,dsetname)
     endif

     if(variable(7) .eq. 1) then
       varindex=varindex+1
       dsetname='bzr'
       xstart=1+(varindex-1)*(nx+2*nbuf)*(ny+2*nbuf)
       xend  =varindex*(nx+2*nbuf)*(ny+2*nbuf)
       call read2d(this,this%q(xstart:xend),this%nMesh_global(1),this%nMesh_global(2),this%nMesh(1)&
       ,this%nMesh(2),this%nbuf,this%nbuf,2,2,flnm,dsetname)
     endif
   endif !! if ndim .eq. 2


   if(ndim .eq. 3) then
     varindex=0
     if(variable(1) .eq. 1) then
       varindex=varindex+1
       dsetname='den'
       xstart=1+(varindex-1)*(nx+2*nbuf)*(ny+2*nbuf)*(nz+2*nbuf)
       xend  =varindex*(nx+2*nbuf)*(ny+2*nbuf)*(nz+2*nbuf)
       call read3d(this,this%q(xstart:xend),this%nMesh_global(1),this%nMesh_global(2),this%nMesh_global(3),this%nMesh(1)&
       ,this%nMesh(2),this%nMesh(3),this%nbuf,this%nbuf,this%nbuf,2,2,2,flnm,dsetname)
     endif

     if(variable(2) .eq. 1) then
       varindex=varindex+1
       dsetname='momx'
       xstart=1+(varindex-1)*(nx+2*nbuf)*(ny+2*nbuf)*(nz+2*nbuf)
       xend  =varindex*(nx+2*nbuf)*(ny+2*nbuf)*(nz+2*nbuf)
       call read3d(this,this%q(xstart:xend),this%nMesh_global(1),this%nMesh_global(2),this%nMesh_global(3),this%nMesh(1)&
       ,this%nMesh(2),this%nMesh(3),this%nbuf,this%nbuf,this%nbuf,2,2,2,flnm,dsetname)
     endif

     if(variable(3) .eq. 1) then
       varindex=varindex+1
       dsetname='momy'
       xstart=1+(varindex-1)*(nx+2*nbuf)*(ny+2*nbuf)*(nz+2*nbuf)
       xend  =varindex*(nx+2*nbuf)*(ny+2*nbuf)*(nz+2*nbuf)
       call read3d(this,this%q(xstart:xend),this%nMesh_global(1),this%nMesh_global(2),this%nMesh_global(3),this%nMesh(1)&
       ,this%nMesh(2),this%nMesh(3),this%nbuf,this%nbuf,this%nbuf,2,2,2,flnm,dsetname)
     endif

     if(variable(4) .eq. 1) then
       varindex=varindex+1
       dsetname='momz'
       xstart=1+(varindex-1)*(nx+2*nbuf)*(ny+2*nbuf)*(nz+2*nbuf)
       xend  =varindex*(nx+2*nbuf)*(ny+2*nbuf)*(nz+2*nbuf)
       call read3d(this,this%q(xstart:xend),this%nMesh_global(1),this%nMesh_global(2),this%nMesh_global(3),this%nMesh(1)&
       ,this%nMesh(2),this%nMesh(3),this%nbuf,this%nbuf,this%nbuf,2,2,2,flnm,dsetname)
     endif

     if(variable(5) .eq. 1) then
       varindex=varindex+1
       dsetname='bxl'
       xstart=1+(varindex-1)*(nx+2*nbuf)*(ny+2*nbuf)*(nz+2*nbuf)
       xend  =varindex*(nx+2*nbuf)*(ny+2*nbuf)*(nz+2*nbuf)
       call read3d(this,this%q(xstart:xend),this%nMesh_global(1),this%nMesh_global(2),this%nMesh_global(3),this%nMesh(1)&
       ,this%nMesh(2),this%nMesh(3),this%nbuf,this%nbuf,this%nbuf,2,2,2,flnm,dsetname)
     endif

     if(variable(6) .eq. 1) then
       varindex=varindex+1
       dsetname='byl'
       xstart=1+(varindex-1)*(nx+2*nbuf)*(ny+2*nbuf)*(nz+2*nbuf)
       xend  =varindex*(nx+2*nbuf)*(ny+2*nbuf)*(nz+2*nbuf)
       call read3d(this,this%q(xstart:xend),this%nMesh_global(1),this%nMesh_global(2),this%nMesh_global(3),this%nMesh(1)&
       ,this%nMesh(2),this%nMesh(3),this%nbuf,this%nbuf,this%nbuf,2,2,2,flnm,dsetname)
     endif

     if(variable(7) .eq. 1) then
       varindex=varindex+1
       dsetname='bzl'
       xstart=1+(varindex-1)*(nx+2*nbuf)*(ny+2*nbuf)*(nz+2*nbuf)
       xend  =varindex*(nx+2*nbuf)*(ny+2*nbuf)*(nz+2*nbuf)
       call read3d(this,this%q(xstart:xend),this%nMesh_global(1),this%nMesh_global(2),this%nMesh_global(3),this%nMesh(1)&
       ,this%nMesh(2),this%nMesh(3),this%nbuf,this%nbuf,this%nbuf,2,2,2,flnm,dsetname)
     endif
 
     if(variable(8) .eq. 1) then
       varindex=varindex+1
       dsetname='ene'
       xstart=1+(varindex-1)*(nx+2*nbuf)*(ny+2*nbuf)*(nz+2*nbuf)
       xend  =varindex*(nx+2*nbuf)*(ny+2*nbuf)*(nz+2*nbuf)
       call read3d(this,this%q(xstart:xend),this%nMesh_global(1),this%nMesh_global(2),this%nMesh_global(3),this%nMesh(1)&
       ,this%nMesh(2),this%nMesh(3),this%nbuf,this%nbuf,this%nbuf,2,2,2,flnm,dsetname)
     endif

     if(variable(5) .eq. 1) then
       varindex=varindex+1
       dsetname='bxr'
       xstart=1+(varindex-1)*(nx+2*nbuf)*(ny+2*nbuf)*(nz+2*nbuf)
       xend  =varindex*(nx+2*nbuf)*(ny+2*nbuf)*(nz+2*nbuf)
       call read3d(this,this%q(xstart:xend),this%nMesh_global(1),this%nMesh_global(2),this%nMesh_global(3),this%nMesh(1)&
       ,this%nMesh(2),this%nMesh(3),this%nbuf,this%nbuf,this%nbuf,2,2,2,flnm,dsetname)
     endif

     if(variable(6) .eq. 1) then
       varindex=varindex+1
       dsetname='byr'
       xstart=1+(varindex-1)*(nx+2*nbuf)*(ny+2*nbuf)*(nz+2*nbuf)
       xend  =varindex*(nx+2*nbuf)*(ny+2*nbuf)*(nz+2*nbuf)
       call read3d(this,this%q(xstart:xend),this%nMesh_global(1),this%nMesh_global(2),this%nMesh_global(3),this%nMesh(1)&
       ,this%nMesh(2),this%nMesh(3),this%nbuf,this%nbuf,this%nbuf,2,2,2,flnm,dsetname)
     endif

     if(variable(7) .eq. 1) then
       varindex=varindex+1
       dsetname='bzr'
       xstart=1+(varindex-1)*(nx+2*nbuf)*(ny+2*nbuf)*(nz+2*nbuf)
       xend  =varindex*(nx+2*nbuf)*(ny+2*nbuf)*(nz+2*nbuf)
       call read3d(this,this%q(xstart:xend),this%nMesh_global(1),this%nMesh_global(2),this%nMesh_global(3),this%nMesh(1)&
       ,this%nMesh(2),this%nMesh(3),this%nbuf,this%nbuf,this%nbuf,2,2,2,flnm,dsetname)
     endif
   endif !! if ndim .eq. 3

   call this%setBoundary(this%q)
   call this%exchangeBdryMPI(this%q,this%winq)
  !call this%setBoundary(this%q)
   deallocate(nMesh,nMeshdble,leftBdry,rightBdry)

   !zwg! call this%setBoundary(this%q)
   !zwg! call this%exchangeBdryMPI(this%q,this%winq)
   !zwg! if (this%enable_sg .eqv. .true.) then
     !zwg! call this%calcSelfgravity(this%q)
   !zwg! endif
   !zwg! deallocate(nMesh,nMeshdble,leftBdry,rightBdry)


endif
end subroutine setTime



subroutine setRestart(this,fstart,tend,dtout)
class(grid)::this
integer::fstart,gridID,nbuf,coordType,ndim
double precision::valuetmp,variabletmp(8)
integer::variable(8),ntotal,nx,ny,nz,nvar
double precision::tend,dtout
double precision,dimension(:),allocatable::nMeshdble,leftBdry,rightBdry
integer,dimension(:),allocatable::nMesh
integer::varindex
integer::xstart,xend
integer::ystart,yend
integer::zstart,zend
character(len=13)::flnm
character(len=4)::gridnum
character(len=4)::filenum
character(len=20)::dsetname
character(len=5)::gridname

this%fstart=fstart
this%tend=tend
this%dtout=dtout
this%toutput=dtout

if(fstart .eq. 0) then
   this%t=0.d0
   this%fnum=0
   this%writeFlag=.false.
else

   gridID=this%gridID
   this%fnum=fstart
   this%t=fstart*dtout
   this%toutput=this%t+dtout
   this%isRestart=1






   !this%writeFlag=.false.


   !this%changeSolver=.false.
   !write(gridnum,'(I4.4)') gridID
   !write(filenum,'(I4.4)') fstart
   !flnm='g'//gridnum//'_'//filenum//'.h5'
   !gridname='g'//gridnum
   !print *,"gridModule.f03: reading data from ",flnm

endif
end subroutine setRestart



subroutine setCFL(this,CFL)
class(grid)::this
double precision::CFL
this%CFL=CFL
end subroutine setCFL

subroutine initVariable(this)
class(grid)::this
integer::ndim

ndim=this%ndim
if(ndim .eq. 1) then
   call init1d(this,this%q)
endif
if(ndim .eq. 2) then
   call init2d(this,this%q)
endif
if(ndim .eq. 3) then
   call init3d(this,this%q)
endif
end subroutine initVariable


subroutine init_For_FFTW(this)
class(grid)::this
integer::ndim

ndim=this%ndim

if(ndim .eq. 2) then
   call init2d_for_FFTW(this,this%q,this%Den_ini)
endif

if(ndim .eq. 3) then
   call init3d_for_FFTW(this,this%q,this%Den_ini)
endif

end subroutine init_For_FFTW



subroutine setSoundSpeed(this,snd)
class(grid)::this
double precision::snd
this%snd=snd
end subroutine setSoundSpeed

subroutine setEoS(this,eosType)
class(grid)::this
integer::eosType
this%eosType=eosType
end subroutine setEoS




subroutine setVariable(this,variable)
class(grid)::this
integer,dimension(1:8)::variable
integer::den,vx,vy,vz,bx,by,bz,ene
integer::ndim,nbuf,nx,ny,nz,nvar
type(C_PTR)::cdata, cdata1
integer(C_INTPTR_T)::nx_ker,ny_ker,nz_ker

this%variable=variable
den=variable(1)
vx =variable(2)
vy =variable(3)
vz =variable(4)
bx =variable(5)
by =variable(6)
bz =variable(7)
ene=variable(8)

nbuf=this%nbuf
ndim=this%ndim

nvar=den+vx+vy+vz+(bx+by+bz)*2+ene
!!! zwg: why (bx+by+bz)*2  ???

this%nvar=nvar


!!!!  zwg: nx,ny,nz are local variables,
!!!!  zwg to compute nx ny nz here!!!
!!!!  zwg: ny ,nx,nz are all included the buffer cells

nx=this%nMesh(1)+2*nbuf

if(ndim .gt. 1) then
   ny=this%nMesh(2)+2*nbuf
endif

if(ndim .gt. 2) then
   nz=this%nMesh(3)+2*nbuf
endif

if(ndim .eq. 1) then
                        !print *,"gridModule.f03: nx*nvar,nx,nvar=",nx*nvar,nx,nvar
                        !print *,"gridModule.f03: variable=",variable

   allocate(this%q(nx*nvar))
   allocate(this%q1(nx*nvar))
   allocate(this%q2(nx*nvar))

   if(this%enable_sg) then
     print *,"gridModule:f03: selfgravity has not yet implemented for 1D problems..."
     stop
   endif

elseif(ndim .eq. 2) then
   allocate(this%q(nx*ny*nvar))
   allocate(this%q1(nx*ny*nvar))
   allocate(this%q2(nx*ny*nvar))

   allocate(this%databuf1(nbuf*ny))

   allocate(this%leftbuf(nbuf*ny))

   allocate(this%rightbuf(nbuf*ny))


   if(this%enable_sg) then

     if(this%sgBdryTYpe .eq. 0) then
       allocate(this%sgfx(this%nMesh(1)*this%nMesh(2)))
       allocate(this%sgfy(this%nMesh(1)*this%nMesh(2)))

       nx_ker=2*this%nMesh(1)
       ny_ker=2*this%nMesh(2)

       allocate(this%sgDensityBuffer(this%nMesh(1)*2*this%nMesh(2)))

       cdata=fftw_alloc_complex(nx_ker*ny_ker)
       call c_f_pointer(cdata,this%sgfxKernel,[nx_ker*ny_ker])

       cdata=fftw_alloc_complex(nx_ker*ny_ker)
       call c_f_pointer(cdata,this%sgfyKernel,[nx_ker*ny_ker])

       cdata=fftw_alloc_complex(nx_ker*ny_ker)
       call c_f_pointer(cdata,this%sgDenCmplx,[nx_ker*ny_ker])

       cdata=fftw_alloc_complex(nx_ker*ny_ker)
       call c_f_pointer(cdata,this%sgfxCmplx,[nx_ker*ny_ker])
       
       cdata=fftw_alloc_complex(nx_ker*ny_ker)
       call c_f_pointer(cdata,this%sgfyCmplx,[nx_ker*ny_ker])   
 
       call sgPlan2D(this,this%sgDenCmplx,this%sgfxCmplx,this%sgfyCmplx)
     endif  !! if this%sgBdryType .eq. 0 (isolated)

     if(this%sgBdryType .eq. 1) then
        !print *,"gridModule.f03: poisson solver for 2D periodic has not yet implemented!!"  !!!!zwg :: 
        !stop


       !!! The following code is implimented by zwg for 2D poisson solver with periodic BCs, and should be checked before using

       allocate(this%sgfx(this%nMesh(1)*this%nMesh(2)))
       allocate(this%sgfy(this%nMesh(1)*this%nMesh(2)))
       allocate(this%sgDensityBuffer(this%nMesh(1)*this%nMesh(2)))

       nx_ker=this%nMesh(1)
       ny_ker=this%nMesh(2)

       cdata=fftw_alloc_complex(nx_ker*ny_ker)
       call c_f_pointer(cdata,this%sgDenCmplx,[nx_ker*ny_ker])

       cdata=fftw_alloc_complex(nx_ker*ny_ker)
       call c_f_pointer(cdata,this%sgfxCmplx,[nx_ker*ny_ker])

       cdata=fftw_alloc_complex(nx_ker*ny_ker)
       call c_f_pointer(cdata,this%sgfyCmplx,[nx_ker*ny_ker])

       call sgPlan2D(this,this%sgDenCmplx,this%sgfxCmplx,this%sgfyCmplx)

       !!! The above code is implimented by zwg for 2D poisson solver with periodic BCs, and should be checked before using


     endif !! if this%sgBdryType .eq. 1 (periodic)

   endif  !! if this%enable_sg




    if(this%enable_DT) then  !!!zwg:  DT is short for Driving Turbulence (this is for 2D problems implimentation)

       allocate(this%Den_ini(this%nMesh(1)*this%nMesh(2)))
       allocate(this%DTVx(this%nMesh(1)*this%nMesh(2)))
       allocate(this%DTVy(this%nMesh(1)*this%nMesh(2)))
       !allocate(this%DTVz(this%nMesh(1)*this%nMesh(2)))

       nx_ker=this%nMesh(1)
       ny_ker=this%nMesh(2)
       
        cdata=fftw_alloc_complex(2*(nx_ker/2+1)*ny_ker)
        call c_f_pointer(cdata,this%Den_R2C_in, [2*(nx_ker/2+1)*ny_ker])

        cdata=fftw_alloc_complex((nx_ker/2+1)*ny_ker)
        call c_f_pointer(cdata,this%Den_R2C_out,[(nx_ker/2+1)*ny_ker] )

        cdata=fftw_alloc_complex(2*(nx_ker/2+1)*ny_ker)
        call c_f_pointer(cdata,this%Den_C2R_out,[2*(nx_ker/2+1)*ny_ker])
        

        cdata=fftw_alloc_complex((nx_ker/2+1)*ny_ker)
        call c_f_pointer(cdata,this%DTVx_C2R_in,[(nx_ker/2+1)*ny_ker] )
  
        cdata=fftw_alloc_complex(2*(nx_ker/2+1)*ny_ker)
        call c_f_pointer(cdata,this%DTVx_C2R_out,[2*(nx_ker/2+1)*ny_ker])

        cdata=fftw_alloc_complex((nx_ker/2+1)*ny_ker)
        call c_f_pointer(cdata,this%DTVy_C2R_in,[(nx_ker/2+1)*ny_ker] )
  
        cdata=fftw_alloc_complex(2*(nx_ker/2+1)*ny_ker)
        call c_f_pointer(cdata,this%DTVy_C2R_out,[2*(nx_ker/2+1)*ny_ker])

        !cdata=fftw_alloc_complex((nx_ker/2+1)*ny_ker)
        !call c_f_pointer(cdata,this%DTVz_C2R_in,[(nx_ker/2+1)*ny_ker] )
  
        !cdata=fftw_alloc_complex(2*(nx_ker/2+1)*ny_ker)
        !call c_f_pointer(cdata,this%DTVz_C2R_out,[2*(nx_ker/2+1)*ny_ker])


        call DTPlan2D(this,this%Den_R2C_in,this%Den_R2C_out,this%Den_C2R_out, &
        & this%DTVx_C2R_in,this%DTVx_C2R_out,this%DTVy_C2R_in,this%DTVy_C2R_out)

        !call DTPlan3D(this,this%Den_R2C_in,this%Den_R2C_out,this%Den_C2R_out, &
        !& this%DTVx_C2R_in,this%DTVx_C2R_out,this%DTVy_C2R_in,this%DTVy_C2R_out,this%DTVz_C2R_in,this%DTVz_C2R_out)


        !call Plan2D_RC(this,this%Den_R2C_in,this%Den_R2C_out,this%Den_C2R_out)
        !call Plan2D_RC(this,this%Den_R2C_in,this%Den_R2C_out,this%Den_C2R_out)
        !complex(C_DOUBLE_COMPLEX),dimension(:),pointer::DTVx_C2R_in,DTVy_C2R_in,DTVz_C2R_in
        !real(C_DOUBLE), dimension(:),pointer:: DTVx_C2R_out,DTVy_C2R_out,DTVz_C2R_out

        !DTVx,DTVy,DTVz Den_ini

    endif

elseif(ndim .eq. 3) then
   allocate(this%q(nx*ny*nz*nvar))
   allocate(this%q1(nx*ny*nz*nvar))
   allocate(this%q2(nx*ny*nz*nvar))

   allocate(this%databuf1(nbuf*ny*nz))
   allocate(this%databuf2(nbuf*nx*nz))

   allocate(this%leftbuf (nbuf*ny*nz))
   allocate(this%rightbuf(nbuf*ny*nz))

   allocate(this%frontbuf(nbuf*nx*nz))
   allocate(this%backbuf (nbuf*nx*nz))
 
   if(this%enable_sg) then
     if(this%sgBdryType .eq. 0) then
       allocate(this%sgfx(this%nMesh(1)*this%nMesh(2)*this%nMesh(3)))
       allocate(this%sgfy(this%nMesh(1)*this%nMesh(2)*this%nMesh(3)))
       allocate(this%sgfz(this%nMesh(1)*this%nMesh(2)*this%nMesh(3)))
     
       nx_ker=2*this%nMesh(1)
       ny_ker=2*this%nMesh(2)
       nz_ker=2*this%nMesh(3)

       allocate(this%sgDensityBuffer(this%nMesh(1)*this%nMesh(2)*2*this%nMesh(3))) 
   
       cdata=fftw_alloc_complex(nx_ker*ny_ker*nz_ker)
       call c_f_pointer(cdata,this%sgfxKernel,[nx_ker*ny_ker*nz_ker])

       cdata=fftw_alloc_complex(nx_ker*ny_ker*nz_ker)
       call c_f_pointer(cdata,this%sgfyKernel,[nx_ker*ny_ker*nz_ker])

       cdata=fftw_alloc_complex(nx_ker*ny_ker*nz_ker)
       call c_f_pointer(cdata,this%sgfzKernel,[nx_ker*ny_ker*nz_ker])

       cdata=fftw_alloc_complex(nx_ker*ny_ker*nz_ker)
       call c_f_pointer(cdata,this%sgDenCmplx,[nx_ker*ny_ker*nz_ker])
     
       cdata=fftw_alloc_complex(nx_ker*ny_ker*nz_ker)
       call c_f_pointer(cdata,this%sgfxCmplx,[nx_ker*ny_ker*nz_ker])

       cdata=fftw_alloc_complex(nx_ker*ny_ker*nz_ker)
       call c_f_pointer(cdata,this%sgfyCmplx,[nx_ker*ny_ker*nz_ker])

       cdata=fftw_alloc_complex(nx_ker*ny_ker*nz_ker)
       call c_f_pointer(cdata,this%sgfzCmplx,[nx_ker*ny_ker*nz_ker])

       call sgPlan3D(this,this%sgDenCmplx,this%sgfxCmplx,this%sgfyCmplx,this%sgfzCmplx)
     endif !! if this%sgBdryType .eq. 0  (isolated)

     if(this%sgBdryType .eq. 1) then
       allocate(this%sgfx(this%nMesh(1)*this%nMesh(2)*this%nMesh(3)))
       allocate(this%sgfy(this%nMesh(1)*this%nMesh(2)*this%nMesh(3)))
       allocate(this%sgfz(this%nMesh(1)*this%nMesh(2)*this%nMesh(3)))
       allocate(this%sgDensityBuffer(this%nMesh(1)*this%nMesh(2)*this%nMesh(3)))

       nx_ker=this%nMesh(1)
       ny_ker=this%nMesh(2)
       nz_ker=this%nMesh(3)
   
       cdata=fftw_alloc_complex(nx_ker*ny_ker*nz_ker)
       call c_f_pointer(cdata,this%sgDenCmplx,[nx_ker*ny_ker*nz_ker])

       cdata=fftw_alloc_complex(nx_ker*ny_ker*nz_ker)
       call c_f_pointer(cdata,this%sgfxCmplx,[nx_ker*ny_ker*nz_ker])

       cdata=fftw_alloc_complex(nx_ker*ny_ker*nz_ker)
       call c_f_pointer(cdata,this%sgfyCmplx,[nx_ker*ny_ker*nz_ker])

       cdata=fftw_alloc_complex(nx_ker*ny_ker*nz_ker)
       call c_f_pointer(cdata,this%sgfzCmplx,[nx_ker*ny_ker*nz_ker])

       call sgPlan3D(this,this%sgDenCmplx,this%sgfxCmplx,this%sgfyCmplx,this%sgfzCmplx)            
 
     endif !! zwg: if this%sgBdryType .eq. 1 (periodic)

   endif  !! zwg :if this%enable_sg




    if(this%enable_DT) then  !!!zwg:  DT is short for Driving Turbulence (this is for 3D problems implimentation)

       allocate(this%Den_ini(this%nMesh(1)*this%nMesh(2)*this%nMesh(3)))
       allocate(this%DTVx(this%nMesh(1)*this%nMesh(2)*this%nMesh(3)))
       allocate(this%DTVy(this%nMesh(1)*this%nMesh(2)*this%nMesh(3)))
       allocate(this%DTVz(this%nMesh(1)*this%nMesh(2)*this%nMesh(3)))

       nx_ker=this%nMesh(1)
       ny_ker=this%nMesh(2)
       nz_ker=this%nMesh(3)
       
        cdata=fftw_alloc_complex(2*(nx_ker/2+1)*ny_ker*nz_ker)
        call c_f_pointer(cdata,this%Den_R2C_in, [2*(nx_ker/2+1)*ny_ker*nz_ker])

        cdata=fftw_alloc_complex((nx_ker/2+1)*ny_ker*nz_ker)
        call c_f_pointer(cdata,this%Den_R2C_out,[(nx_ker/2+1)*ny_ker*nz_ker] )

        cdata=fftw_alloc_complex(2*(nx_ker/2+1)*ny_ker*nz_ker)
        call c_f_pointer(cdata,this%Den_C2R_out,[2*(nx_ker/2+1)*ny_ker*nz_ker])
        

        cdata=fftw_alloc_complex((nx_ker/2+1)*ny_ker*nz_ker)
        call c_f_pointer(cdata,this%DTVx_C2R_in,[(nx_ker/2+1)*ny_ker*nz_ker] )
  
        cdata=fftw_alloc_complex(2*(nx_ker/2+1)*ny_ker*nz_ker)
        call c_f_pointer(cdata,this%DTVx_C2R_out,[2*(nx_ker/2+1)*ny_ker*nz_ker])

        cdata=fftw_alloc_complex((nx_ker/2+1)*ny_ker*nz_ker)
        call c_f_pointer(cdata,this%DTVy_C2R_in,[(nx_ker/2+1)*ny_ker*nz_ker] )
  
        cdata=fftw_alloc_complex(2*(nx_ker/2+1)*ny_ker*nz_ker)
        call c_f_pointer(cdata,this%DTVy_C2R_out,[2*(nx_ker/2+1)*ny_ker*nz_ker])

        cdata=fftw_alloc_complex((nx_ker/2+1)*ny_ker*nz_ker)
        call c_f_pointer(cdata,this%DTVz_C2R_in,[(nx_ker/2+1)*ny_ker*nz_ker] )
  
        cdata=fftw_alloc_complex(2*(nx_ker/2+1)*ny_ker*nz_ker)
        call c_f_pointer(cdata,this%DTVz_C2R_out,[2*(nx_ker/2+1)*ny_ker*nz_ker])


        !call DTPlan2D(this,this%Den_R2C_in,this%Den_R2C_out,this%Den_C2R_out, &
        !& this%DTVx_C2R_in,this%DTVx_C2R_out,this%DTVy_C2R_in,this%DTVy_C2R_out)

        call DTPlan3D(this,this%Den_R2C_in,this%Den_R2C_out,this%Den_C2R_out, &
        & this%DTVx_C2R_in,this%DTVx_C2R_out,this%DTVy_C2R_in,this%DTVy_C2R_out,this%DTVz_C2R_in,this%DTVz_C2R_out)


        !call Plan2D_RC(this,this%Den_R2C_in,this%Den_R2C_out,this%Den_C2R_out)
        !call Plan2D_RC(this,this%Den_R2C_in,this%Den_R2C_out,this%Den_C2R_out)
        !complex(C_DOUBLE_COMPLEX),dimension(:),pointer::DTVx_C2R_in,DTVy_C2R_in,DTVz_C2R_in
        !real(C_DOUBLE), dimension(:),pointer:: DTVx_C2R_out,DTVy_C2R_out,DTVz_C2R_out

        !DTVx,DTVy,DTVz Den_ini

    endif






   

endif  !!!!! zwg : end of if(ndim .eq. 1) then

end subroutine setVariable


subroutine setMesh(this,nMesh_global,leftBdry_global,rightBdry_global,nbuf,coordType,gridID)

class(grid)::this

integer,dimension(:)::nMesh_global
double precision,dimension(:)::leftBdry_global,rightBdry_global

integer:: nbuf,coordType,gridID,ndim
integer:: remainderx,remaindery,remainderz,i

integer,dimension(:),allocatable::nx_list,ny_list,nz_list
integer::gnx_list_R,gny_list_R,gnz_list_R
integer::nprocx,nprocy,nprocz
integer::nxbase,nybase,nzbase
integer::coordx,coordy,coordz
integer::acc
double precision::dx,dy,dz
double precision::dlogr,logri,logro

integer::modtot

nprocx=0
nprocy=0
nprocz=0

nxbase=0
nybase=0
nzbase=0

ndim=size(nMesh_global)


!!! be caraful, this moduel is used just for the self-gravity case!!!!
if(this%enable_sg) then
  modtot=0
  !!! zwg :: if total number of mesh points can be divided by nprocs 
  !!! zwg :: then, modtot always be 0
  do i=1,ndim
    modtot=modtot+mod(nMesh_global(i),nprocs)
  enddo
 
  if(modtot .ne. 0) then
    print *,"gridModule.f03: mesh nubmer in all direction should be divisible by the number of processors=",nprocs
    !!! zwg :: modtot should be zero!!!  not flexiable!!!!
    stop
  endif
endif  !! zwg: end of if(this%enable_sg) then


!!!! zwg added for check mesh for driving turbulence case!!!
if(this%enable_DT) then

   modtot=0
   if(ndim==2)then
     modtot=mod(nMesh_global(2),nprocs)
   else if(ndim==3)then
     modtot=mod(nMesh_global(3),nprocs)
   endif !! if(ndim==2)then


  if(modtot .ne. 0) then
    !print *,""
    print *,"gridModule.f03: for DT case, mesh nubmer in last direction should be divisible by the number of processors=",nprocs
    stop
  endif


endif !!!  end of if(this%enable_DT) 





allocate(this%nMesh(ndim),this%leftBdry(ndim),this%rightBdry(ndim))  !!!! zwg: Local grid information 

allocate(this%nMesh_global(ndim),this%leftBdry_global(ndim),this%rightBdry_global(ndim)) !!!! zwg:  global grid information!!!!

!!! zwg: I do not have ideas about xl, xr, and xc currently.

allocate(this%xl(ndim),this%xr(ndim),this%xc(ndim),this%dx(ndim))

this%nbuf=nbuf

this%nMesh_global=nMesh_global

this%leftBdry_global=leftBdry_global
this%rightBdry_global=rightBdry_global
this%coordType=coordType

this%gridID=gridID
this%ndim=ndim


!!!! zwg: actually for each process,namely, the local grid information
!!!! calculate nMesh, leftBdry, rightBdry
 
nprocx=this%dims_mpi(1)  !!! zwg also dims(1) got by subroutine of MPI_DIMS_CREATE(nprocs,ndim,dims,ierr) in file testSuiteMPI.f03
!! nprocx: the processes numbers in x Direction

allocate(nx_list(0:nprocx-1))
!!! the index of processes in x direction

allocate(this%nx_list(0:nprocx-1))

nxbase=nMesh_global(1)

remainderx=mod(nxbase,nprocx)  !!! zwg :the remainder mesh points in x direction



!!! zwg: nx_list is the number of mesh points for each processes in x direction!!! 
do i=0,nprocx-1
  if(i .gt. (nprocx-remainderx-1)) then  !!! zwg if the 
    nx_list(i)=(nxbase-remainderx)/nprocx+1
  else
    nx_list(i)=(nxbase-remainderx)/nprocx
  endif
enddo

this%nx_list=nx_list

this%nMesh(1)=nx_list(this%mpiCoord(1))

if (coordType .eq. 1 .or. coordType .eq. 3) then !!! Cartesian coordinates

!!! in the case of polar coordinates (uniform mesh in r) x=>r, y=>phi, z=>z
dx=(this%rightBdry_global(1)-this%leftBdry_global(1))/dble(nxbase)  
!!! zwg : dble function of fortran to transfor a data to double data type


!!! each processe does the following loop , but their max loop index this%mpiCoord is locally dependented, 
!!! which causes a different acc !! hahaha (^_^)

acc=0
do i=0,this%mpiCoord(1)
  acc=acc+nx_list(i)
enddo
gnx_list_R=acc
!integer::gnx_list_R,gny_list_R,gnz_list_R

!!! zwg: actually acc is the global grid index for each block 



this%leftBdry(1) =this%leftBdry_global(1)+(acc-this%nMesh(1))*dx
this%rightBdry(1)=this%leftBdry_global(1)+acc*dx
endif !! zwg  end of if (coordType .eq. 1 .or. coordType .eq. 3) then !!! Cartesian coordinates

if (coordType .eq. 2) then !!! Cylindrical coordinates with logrithmic r
!!!! x=>r, y=>phi, z=>z
logri=dlog( this%leftBdry_global(1))
logro=dlog(this%rightBdry_global(1))
dlogr=(logro-logri)/dble(nxbase)
acc=0
do i=0,this%mpiCoord(1)
  acc=acc+nx_list(i)
enddo
this%leftBdry(1)=dexp(logri+(acc-this%nMesh(1))*dlogr)
this%rightBdry(1)=dexp(logri+acc*dlogr)
endif   !!! zwg : if (coordType .eq. 2) then !!! Cylindrical coordinates with logrithmic r


if(ndim .gt. 1) then
  nprocy=this%dims_mpi(2)
  allocate(ny_list(0:nprocy-1))
  allocate(this%ny_list(0:nprocy-1))
  nybase=nMesh_global(2)
  remaindery=mod(nybase,nprocy)
  do i=0,nprocy-1
    if(i .gt. (nprocy-remaindery-1)) then
      ny_list(i)=(nybase-remaindery)/nprocy+1
    else
      ny_list(i)=(nybase-remaindery)/nprocy
    endif
  enddo
  this%ny_list=ny_list
  this%nMesh(2)=ny_list(this%mpiCoord(2))
  dy=(this%rightBdry_global(2)-this%leftBdry_global(2))/dble(nybase)
  acc=0
  do i=0,this%mpiCoord(2)
    acc=acc+ny_list(i)
  enddo
  gny_list_R=acc
  this%leftBdry(2)=this%leftBdry_global(2)+(acc-this%nMesh(2))*dy
  this%rightBdry(2)=this%leftBdry_global(2)+acc*dy
endif


if(ndim .gt. 2) then
  nprocz=this%dims_mpi(3)
  allocate(nz_list(0:nprocz-1))
  allocate(this%nz_list(0:nprocz-1))
  nzbase=nMesh_global(3)
  remainderz=mod(nzbase,nprocz)
  
  do i=0,nprocz-1
    if(i .gt. (nprocz-remainderz-1)) then
      nz_list(i)=(nzbase-remainderz)/nprocz+1
    else
      nz_list(i)=(nzbase-remainderz)/nprocz
    endif
  enddo
  this%nz_list=nz_list
  this%nMesh(3)=nz_list(this%mpiCoord(3))
  dz=(this%rightBdry_global(3)-this%leftBdry_global(3))/dble(nzbase)
  acc=0
  do i=0,this%mpiCoord(3)
    acc=acc+nz_list(i)
  enddo
  gnz_list_R=acc

  this%leftBdry(3)=this%leftBdry_global(3)+(acc-this%nMesh(3))*dz
  this%rightBdry(3)=this%leftBdry_global(3)+acc*dz
endif

this%gnx_list_R=gnx_list_R
this%gny_list_R=gny_list_R
this%gnz_list_R=gnz_list_R

!write(*,*)"gnx_list_R=",gnx_list_R,"gny_list_R=",gny_list_R,"gnz_list_R=",gnz_list_R

!print *,"gridModule.f03: myid=",myid,"(nprocx,nprocy,nprocz)=",nprocx,nprocy,nprocz
!print *,"gridModule.f03: myid=",myid,"nMesh=",this%nMesh
!print *,"gridModule.f03: myid=",myid,nz_list
!print *,"gridModule.f03: myid=",myid,"leftBdry=",this%leftBdry,"rightBdry=",this%rightBdry

do i=1,ndim
   allocate(this%xl(i)%coords(1-nbuf:this%nMesh(i)+nbuf))
   allocate(this%xr(i)%coords(1-nbuf:this%nMesh(i)+nbuf))
   allocate(this%xc(i)%coords(1-nbuf:this%nMesh(i)+nbuf))
   allocate(this%dx(i)%coords(1-nbuf:this%nMesh(i)+nbuf))
enddo

call setCoordinates(this)

end subroutine setMesh





end module gridModule
