!!! subroutines for sink particle
!!! zcao, 24-12-2019
subroutine AddToList(list, element)
use gridModule
IMPLICIT NONE
integer :: i, isize
type(sp_element), intent(in) :: element
type(sp_element), dimension(:), allocatable, intent(inout) :: list
type(sp_element), dimension(:), allocatable :: clist
if(allocated(list)) then
    isize = size(list)
    allocate(clist(isize+1))
    do i=1,isize          
    clist(i) = list(i)
    end do
    clist(isize+1) = element

    deallocate(list)
    call move_alloc(clist, list)

else
    allocate(list(1))
    list(1) = element
end if
end subroutine AddToList

subroutine createSP2D(this,q)
use gridModule
implicit none
class(grid)::this
double precision,dimension(1-this%nbuf:this%nMesh(1)+this%nbuf,1-this%nbuf:this%nMesh(2)+this%nbuf,this%nvar)::q
double precision,dimension(this%nMesh(1),this%nMesh(2))::den
double precision:: xc,yc,dx,dy
double precision:: rho_sp, mass_sp
type(sp_element)::sp_el
integer::nx_global,ny_global
integer::i,j
integer::nsp

rho_sp=this%rho_sp
xc=this%xc(1)%coords
yc=this%xc(2)%coords
dx=this%dx(1)%coords
dy=this%dx(2)%coords
nx_global=this%nMesh_global(1)
ny_global=this%nMesh_global(2)
den=q(1:this%nMesh(1),1:this%nMesh(2),1)
do j = 1,ny_global
  do i = 1,nx_global
    if (den(i,j) .gt. rho_sp) then
     this%num_sp=this%num_sp+1
     dden=den(i,j)-rho_sp
     dmass=dden*dx(i)*dy(j)
     sp_el%loc(:) = /xc(i),yc(j),0.d0/
     sp_el%mass = dmass
     call AddToList(this%sp_arr,sp_el)
     q(i,j,1)=rho_sp
    endif
  enddo
enddo
