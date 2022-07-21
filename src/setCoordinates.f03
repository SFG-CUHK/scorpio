subroutine setCoordinates(this)
use gridModule
class(grid)::this
integer::ndim,i,j
double precision::dx,dlogr,logri,logro

ndim=this%ndim
!!!! coordType == 1 Cartesian coordinates
!!!! coordType == 2 cylindrical coordinate (logrithmic in r)
!!!! coordType == 3 cylindrical coordinate (uniform in r)


!!!! zwg try to understabd how to set mesh here!!!

if(this%coordType .eq. 1 .or. this%coordType .eq. 3) then
   do i=1, ndim
      dx=(this%rightBdry(i)-this%leftBdry(i))/dble(this%nMesh(i))
      do j=1-this%nbuf, this%nMesh(i)+this%nbuf
         this%dx(i)%coords(j)=dx
         this%xl(i)%coords(j)=this%leftBdry(i)+dble(j-1)*dx
         this%xr(i)%coords(j)=this%leftBdry(i)+dble(j  )*dx
         this%xc(i)%coords(j)=0.5d0*(this%xl(i)%coords(j)+this%xr(i)%coords(j))
      enddo
   enddo   
endif




if(this%coordType .eq. 2) then
!!!!! set r 
  logri=dlog(this%leftBdry(1))
  logro=dlog(this%rightBdry(1))
  dlogr=(logro-logri)/dble(this%nMesh(1))
  do j=1-this%nbuf,this%nMesh(1)+this%nbuf
    this%xl(1)%coords(j)=dexp(logri+dble(j-1)*dlogr)
    this%xr(1)%coords(j)=dexp(logri+dble(j  )*dlogr)
    this%dx(1)%coords(j)=this%xr(1)%coords(j)-this%xl(1)%coords(j)
    this%xc(1)%coords(j)=0.5d0*(this%xl(1)%coords(j)+this%xr(1)%coords(j))
  enddo

!!!!! set phi and z  
  do i=2, ndim  !! uniform in phi, z
      dx=(this%rightBdry(i)-this%leftBdry(i))/dble(this%nMesh(i))
      do j=1-this%nbuf, this%nMesh(i)+this%nbuf
         this%dx(i)%coords(j)=dx
         this%xl(i)%coords(j)=this%leftBdry(i)+dble(j-1)*dx
         this%xr(i)%coords(j)=this%leftBdry(i)+dble(j  )*dx
         this%xc(i)%coords(j)=0.5d0*(this%xl(i)%coords(j)+this%xr(i)%coords(j))
      enddo    
  enddo
endif

end subroutine setCoordinates
