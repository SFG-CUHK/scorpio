subroutine initSGWindows2D(this,sgDensityBuffer)
use gridModule
implicit none
class(grid)::this
double precision,dimension(this%nMesh(1),2*this%nMesh(2))::sgDensityBuffer
integer::sizedouble,ierr,datasize
integer::nx,ny

nx=this%nMesh(1)
ny=this%nMesh(2)

call MPI_SIZEOF(sgDensityBuffer(1,1),sizedouble,ierr)
datasize=(nx)*(2*ny)*sizedouble
call MPI_WIN_CREATE(sgDensityBuffer,datasize,sizedouble,MPI_INFO_NULL,MPI_COMM_WORLD,this%winsg,ierr)
end subroutine initSGWindows2D
