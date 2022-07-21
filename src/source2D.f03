subroutine source2D(this,q,q1)
use gridModule
use testSuiteMPI
class(grid)::this
double precision,dimension(1-this%nbuf:this%nMesh(1)+this%nbuf,1-this%nbuf:this%nMesh(2)+this%nbuf,this%nvar)::q,q1
double precision,dimension(1-this%nbuf:this%nMesh(1)+this%nbuf)::xc
double precision,dimension(1-this%nbuf:this%nMesh(2)+this%nbuf)::yc
integer::nx,nvar,nbuf,gridID,eosType
integer::i

gridID=this%gridID
xc=this%xc(1)%coords
yc=this%xc(2)%coords

if(testOnOff) then
select case (gridID)
   case (44)
     call sourceRayleighTaylorInstabilityHD2D(this,q,q1)
   case (45)
     call sourceRayleighTaylorInstabilityMHD2D(this,q,q1)
end select
endif

end subroutine source2D

