subroutine initMPIWindows1D(this,q,q1,q2)
use gridModule
implicit none
class(grid)::this
double precision,dimension(1-this%nbuf:this%nMesh(1)+this%nbuf,this%nvar)::q,q1,q2
integer::nx,nbuf,nvar,sizedouble,ierr
integer(kind=MPI_ADDRESS_KIND)::datasize

nx=this%nMesh(1)
nbuf=this%nbuf
nvar=this%nvar
call MPI_SIZEOF(q(1,1),sizedouble,ierr)

datasize=(nx+2*nbuf)*nvar*sizedouble
call MPI_WIN_CREATE(q,datasize,sizedouble,MPI_INFO_NULL,MPI_COMM_WORLD,this%winq,ierr)
call MPI_WIN_CREATE(q1,datasize,sizedouble,MPI_INFO_NULL,MPI_COMM_WORLD,this%winq1,ierr)
call MPI_WIN_CREATE(q2,datasize,sizedouble,MPI_INFO_NULL,MPI_COMM_WORLD,this%winq2,ierr)
end subroutine initMPIWindows1D
