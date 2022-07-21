subroutine initSGWindows3D(this,sgDensityBuffer)
use gridModule
implicit none
class(grid)::this
double precision,dimension(this%nMesh(1),this%nMesh(2),2*this%nMesh(3))::sgDensityBuffer
integer::sizedouble,ierr,datasize
integer::nx,ny,nz

nx=this%nMesh(1)
ny=this%nMesh(2)
nz=this%nMesh(3)

if(this%sgBdryType .eq. 0) then
  call MPI_SIZEOF(sgDensityBuffer(1,1,1),sizedouble,ierr)
  datasize=(nx)*(ny)*(2*nz)*sizedouble
  call MPI_WIN_CREATE(sgDensityBuffer,datasize,sizedouble,MPI_INFO_NULL,MPI_COMM_WORLD,this%winsg,ierr)
elseif(this%sgBdryType .eq. 1) then
  if(myid .eq. 0) then
    print *,"initSGWindows3D.f03: no need to initialize MPI windows for perioic selfgravity..."
  endif
endif
end subroutine initSGWindows3D
