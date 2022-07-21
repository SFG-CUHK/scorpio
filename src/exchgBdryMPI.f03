subroutine exchgBdryMPI1D(this,q,win)
use gridModule
use mpi
implicit none
class(grid)::this
double precision,dimension(1-this%nbuf:this%nMesh(1)+this%nbuf,this%nvar)::q
integer::win,nx,nvar,nbuf,ierr,i,left,right
integer (kind=MPI_ADDRESS_KIND)::offset

nx=this%nMesh(1)
nvar=this%nvar
nbuf=this%nbuf
left=this%left_mpi
right=this%right_mpi

call MPI_WIN_FENCE(0,win,ierr)
  if(right .ge. 0) then
    do i=1,nvar
      offset=(i-1)*(this%nx_list(right)+2*nbuf)+nbuf
      call MPI_GET(q(nx+1,i),nbuf,MPI_DOUBLE,right,offset,nbuf,MPI_DOUBLE,win,ierr)
      offset=(i-1)*(this%nx_list(right)+2*nbuf)
      call MPI_PUT(q(nx-1,i),nbuf,MPI_DOUBLE,right,offset,nbuf,MPI_DOUBLE,win,ierr)
    enddo
  endif
call MPI_WIN_FENCE(0,win,ierr)

!print *,"exchgBdryMPI.f03: myid=",myid,"density:",q(:,1)

end subroutine exchgBdryMPI1D
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine exchgBdryMPI2D(this,q,win,winbuf1,databuf1,leftbuf,rightbuf)
use gridModule
use mpi
implicit none
class(grid)::this
double precision,dimension(1-this%nbuf:this%nMesh(1)+this%nbuf,1-this%nbuf:this%nMesh(2)+this%nbuf,this%nvar)::q
double precision,dimension(this%nbuf,this%nMesh(2)+2*this%nbuf)::databuf1
double precision,dimension(1-this%nbuf:0,1-this%nbuf:this%nMesh(2)+this%nbuf)::leftbuf
double precision,dimension(this%nMesh(1)+1:this%nMesh(1)+this%nbuf,1-this%nbuf:this%nMesh(2)+this%nbuf)::rightbuf
integer::win,winbuf1,nx,ny,nvar,nbuf,ierr,i,j,left,right,up,down,coords(2),coordx,coordy
integer(kind=MPI_ADDRESS_KIND)::offset
integer::datasize

nx=this%nMesh(1)
ny=this%nMesh(2)
nbuf=this%nbuf
nvar=this%nvar
left=this%left_mpi
right=this%right_mpi
up=this%up_mpi
down=this%down_mpi


call MPI_WIN_FENCE(0,win,ierr)
  if(down .ge. 0) then
     call MPI_CART_COORDS(this%vu_mpi,down,2,coords,ierr)
     coordx=coords(1)
     coordy=coords(2)
     datasize=nbuf*(nx+2*nbuf)
     do i=1,nvar
       offset=(i-1)*(this%nx_list(coordx)+2*nbuf)*(this%ny_list(coordy)+2*nbuf)+nbuf*(this%nx_list(coordx)+2*nbuf)
       call MPI_GET(q(1-nbuf,ny+1,i),datasize,MPI_DOUBLE,down,offset,datasize,MPI_DOUBLE,win,ierr)
       offset=offset-nbuf*(this%nx_list(coordx)+2*nbuf)
       call MPI_PUT(q(1-nbuf,ny-1,i),datasize,MPI_DOUBLE,down,offset,datasize,MPI_DOUBLE,win,ierr)
     enddo
  endif
call MPI_WIN_FENCE(0,win,ierr)

do i=1,nvar
  rightbuf=q(nx-nbuf+1:nx,1-nbuf:ny+nbuf,i)
  offset=0
  call MPI_WIN_FENCE(0,winbuf1,ierr)
    if (right .ge. 0) then
      call MPI_PUT(rightbuf,nbuf*(ny+2*nbuf),MPI_DOUBLE,right,offset,nbuf*(ny+2*nbuf),MPI_DOUBLE,winbuf1,ierr)
    endif
  call MPI_WIN_FENCE(0,winbuf1,ierr)
  if(left .ge. 0) then
    q(1-nbuf:0,1-nbuf:ny+nbuf,i)=databuf1
  endif
enddo

do i=1,nvar
  leftbuf=q(1:nbuf,1-nbuf:ny+nbuf,i)
  offset=0
  call MPI_WIN_FENCE(0,winbuf1,ierr)
    if(left .ge. 0) then
      call MPI_PUT(leftbuf,nbuf*(ny+2*nbuf),MPI_DOUBLE,left,offset,nbuf*(ny+2*nbuf),MPI_DOUBLE,winbuf1,ierr)
    endif
  call MPI_WIN_FENCE(0,winbuf1,ierr)
  if(right .ge. 0) then
    q(nx+1:nx+nbuf,1-nbuf:ny+nbuf,i)=databuf1
  endif
enddo

!print *,"exchgBdry2D.f03: myid=",myid,"den="
!do j=1-nbuf,ny+nbuf
!  print *,q(:,j,1)
!enddo

end subroutine exchgBdryMPI2D

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine exchgBdryMPI3D(this,q,win,winbuf1,winbuf2,databuf1,databuf2,leftbuf,rightbuf,frontbuf,backbuf)
use gridModule
use mpi
implicit none
class(grid)::this
double precision,dimension(1-this%nbuf:this%nMesh(1)+this%nbuf,1-this%nbuf:this%nMesh(2)+this%nbuf,&
&                          1-this%nbuf:this%nMesh(3)+this%nbuf,this%nvar)::q
double precision,dimension(this%nbuf,this%nMesh(2)+2*this%nbuf,this%nMesh(3)+2*this%nbuf)::databuf1
double precision,dimension(this%nMesh(1)+2*this%nbuf,this%nbuf,this%nMesh(3)+2*this%nbuf)::databuf2
double precision,dimension(1-this%nbuf:0,1-this%nbuf:this%nMesh(2)+this%nbuf,1-this%nbuf:this%nMesh(3)+this%nbuf)::leftbuf
double precision,dimension(1-this%nbuf:this%nMesh(1)+this%nbuf,1-this%nbuf:0,1-this%nbuf:this%nMesh(3)+this%nbuf)::frontbuf
double precision,dimension(this%nMesh(1)+1:this%nMesh(1)+this%nbuf,1-this%nbuf:this%nMesh(2)+this%nbuf,&
&                          1-this%nbuf:this%nMesh(3)+this%nbuf)::rightbuf
double precision,dimension(1-this%nbuf:this%nMesh(1)+this%nbuf,this%nMesh(2)+1:this%nMesh(2)+this%nbuf,&
&                          1-this%nbuf:this%nMesh(3)+this%nbuf)::backbuf
integer::win,winbuf1,winbuf2,nx,ny,nz,nvar,nbuf,ierr,i,j,k,left,right,up,down,top,bottom,coords(3),coordx,coordy,coordz,m
integer(kind=MPI_ADDRESS_KIND)::offset
integer::datasize

nx=this%nMesh(1)
ny=this%nMesh(2)
nz=this%nMesh(3)
nbuf=this%nbuf
nvar=this%nvar
left=this%left_mpi
right=this%right_mpi
up=this%up_mpi
down=this%down_mpi
top=this%top_mpi
bottom=this%bottom_mpi

!!!!!! top & bottom boundaries !!!!!!
call MPI_WIN_FENCE(0,win,ierr)
  if(bottom .ge. 0) then
     call MPI_CART_COORDS(this%vu_mpi,bottom,3,coords,ierr)
     coordx=coords(1)
     coordy=coords(2)
     coordz=coords(3)
     datasize=nbuf*(nx+2*nbuf)*(ny+2*nbuf)
     do i=1,nvar
       offset=(i-1)*(this%nx_list(coordx)+2*nbuf)*(this%ny_list(coordy)+2*nbuf)*(this%nz_list(coordz)+2*nbuf)&
&                  +nbuf*(this%nx_list(coordx)+2*nbuf)*(this%ny_list(coordy)+2*nbuf)
       call MPI_GET(q(1-nbuf,1-nbuf,nz+1,i),datasize,MPI_DOUBLE,bottom,offset,datasize,MPI_DOUBLE,win,ierr)
       offset=offset-nbuf*(this%nx_list(coordx)+2*nbuf)*(this%ny_list(coordy)+2*nbuf)
       call MPI_PUT(q(1-nbuf,1-nbuf,nz-1,i),datasize,MPI_DOUBLE,bottom,offset,datasize,MPI_DOUBLE,win,ierr)
     enddo
     
  endif
call MPI_WIN_FENCE(0,win,ierr)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!! left & right boundaries !!!!!!

do i=1,nvar
  rightbuf=q(nx-nbuf+1:nx,1-nbuf:ny+nbuf,1-nbuf:nz+nbuf,i)
  offset=0
  call MPI_WIN_FENCE(0,winbuf1,ierr)
   if (right .ge. 0) then
    call MPI_PUT(rightbuf,nbuf*(ny+2*nbuf)*(nz+2*nbuf),MPI_DOUBLE,right,offset,nbuf*(ny+2*nbuf)*(nz+2*nbuf),MPI_DOUBLE,winbuf1,ierr)
   endif
  call MPI_WIN_FENCE(0,winbuf1,ierr)
  if(left .ge. 0) then
    q(1-nbuf:0,1-nbuf:ny+nbuf,1-nbuf:nz+nbuf,i)=databuf1
  endif
enddo

do i=1,nvar
  leftbuf=q(1:nbuf,1-nbuf:ny+nbuf,1-nbuf:nz+nbuf,i)
  offset=0
  call MPI_WIN_FENCE(0,winbuf1,ierr)
    if(left .ge. 0) then
     call MPI_PUT(leftbuf,nbuf*(ny+2*nbuf)*(nz+2*nbuf),MPI_DOUBLE,left,offset,nbuf*(ny+2*nbuf)*(nz+2*nbuf),MPI_DOUBLE,winbuf1,ierr)
    endif
  call MPI_WIN_FENCE(0,winbuf1,ierr)
  if(right .ge. 0) then
    q(nx+1:nx+nbuf,1-nbuf:ny+nbuf,1-nbuf:nz+nbuf,i)=databuf1
  endif
enddo

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!! front and back boundaries !!!!!!

do i=1,nvar
  backbuf=q(1-nbuf:nx+nbuf,ny-nbuf+1:ny,1-nbuf:nz+nbuf,i)
  offset=0
  call MPI_WIN_FENCE(0,winbuf2,ierr)
    if (down .ge. 0) then
     call MPI_PUT(backbuf,nbuf*(nx+2*nbuf)*(nz+2*nbuf),MPI_DOUBLE,down,offset,nbuf*(nx+2*nbuf)*(nz+2*nbuf),MPI_DOUBLE,winbuf2,ierr)
    endif
  call MPI_WIN_FENCE(0,winbuf2,ierr)

  if(up .ge. 0) then
    q(1-nbuf:nx+nbuf,1-nbuf:0,1-nbuf:nz+nbuf,i)=databuf2
  endif

enddo

do i=1,nvar
  frontbuf=q(1-nbuf:nx+nbuf,1:nbuf,1-nbuf:nz+nbuf,i)
  offset=0
  call MPI_WIN_FENCE(0,winbuf2,ierr)
    if(up .ge. 0) then
     call MPI_PUT(frontbuf,nbuf*(nx+2*nbuf)*(nz+2*nbuf),MPI_DOUBLE,up,offset,nbuf*(nx+2*nbuf)*(nz+2*nbuf),MPI_DOUBLE,winbuf2,ierr)
    endif
  call MPI_WIN_FENCE(0,winbuf2,ierr)
  if(down .ge. 0) then
     q(1-nbuf:nx+nbuf,ny+1:ny+nbuf,1-nbuf:nz+nbuf,i)=databuf2
  endif

enddo

!print *,"exchgBdry2D.f03: myid=",myid,"den="
!do j=1-nbuf,ny+nbuf
!  print *,q(:,j,1)
!enddo

end subroutine exchgBdryMPI3D
