subroutine setBdry1D(this,q)
use gridModule
use testSuiteMPI
implicit none
class(grid)::this
double precision,dimension(1-this%nbuf:this%nMesh(1)+this%nbuf,this%nvar)::q
integer::gridID

gridID=this%gridID

if(testOnOff) then
select case (gridID)
   case(1)
     call bdryIsoShockTube1D(this,q)
   case(2)
     call bdrySodShockTube1D(this,q)
   case(3)
     call  bdryWCShockTube1D(this,q)
   case(10)
     call bdryBrioWuShockTube1D(this,q)
   case(11)
     call bdryRJ2aShockTube1D(this,q)
   case(22)
     call bdryCShockTest1Dn(this,q)
   case(23)
     call bdryCShockTest1Di(this,q)
   case(24)
     call bdryWardleInstabilityInitn(this,q)
   case(25)
     call bdryWardleInstabilityIniti(this,q)
   case(46)
     call bdryIsoShockTubeMHD1D(this,q)
   case(50)
     call bdryPolyShockTubeHD1D(this,q)
end select
endif
end subroutine setBdry1D

