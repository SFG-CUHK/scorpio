subroutine initMPIWindows3D(this,q,q1,q2,databuf1,databuf2)
use gridModule
implicit none

class(grid)::this

double precision,dimension(1-this%nbuf:this%nMesh(1)+this%nbuf,1-this%nbuf:this%nMesh(2)+this%nbuf &
 & ,1-this%nbuf:this%nMesh(3)+this%nbuf,this%nvar)::q,q1,q2

double precision,dimension(1-this%nbuf:0,1-this%nbuf:this%nMesh(2)+this%nbuf,1-this%nbuf:this%nMesh(3)+this%nbuf)::databuf1

double precision,dimension(1-this%nbuf:this%nMesh(1)+this%nbuf,1-this%nbuf:0,1-this%nbuf:this%nMesh(3)+this%nbuf)::databuf2

integer::nx,ny,nz,nbuf,nvar,sizedouble,ierr

integer(kind=MPI_ADDRESS_KIND)::datasize

!!! databuf1 is used to exchange data which are not contiguously distributed in memory
!!! For 2D cases, data are contiguously distributed in y-direction, while scattered in x-direction

nx=this%nMesh(1)
ny=this%nMesh(2)
nz=this%nMesh(3)
nbuf=this%nbuf
nvar=this%nvar
databuf1=0.d0
databuf2=0.d0

call MPI_SIZEOF(q(1,1,1,1),sizedouble,ierr)

datasize=(nx+2*nbuf)*(ny+2*nbuf)*(nz+2*nbuf)*nvar*sizedouble

call MPI_WIN_CREATE(q ,datasize,sizedouble,MPI_INFO_NULL,MPI_COMM_WORLD,this%winq,ierr)
call MPI_WIN_CREATE(q1,datasize,sizedouble,MPI_INFO_NULL,MPI_COMM_WORLD,this%winq1,ierr)
call MPI_WIN_CREATE(q2,datasize,sizedouble,MPI_INFO_NULL,MPI_COMM_WORLD,this%winq2,ierr)

call MPI_SIZEOF(databuf1(0,0,0),sizedouble,ierr)
datasize=nbuf*(ny+2*nbuf)*(nz+2*nbuf)*sizedouble
call MPI_WIN_CREATE(databuf1,datasize,sizedouble,MPI_INFO_NULL,MPI_COMM_WORLD,this%winbuf1,ierr)
datasize=nbuf*(nx+2*nbuf)*(nz+2*nbuf)*sizedouble
call MPI_WIN_CREATE(databuf2,datasize,sizedouble,MPI_INFO_NULL,MPI_COMM_WORLD,this%winbuf2,ierr)
end subroutine initMPIWindows3D
