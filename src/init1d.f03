subroutine init1d(this,q)
use gridModule
use testSuiteMPI
class(grid)::this
double precision,dimension(1-this%nbuf:this%nMesh(1)+this%nbuf,this%nvar)::q
double precision,dimension(1-this%nbuf:this%nMesh(1)+this%nbuf)::xc
integer::nx,nvar,nbuf,gridID,eosType
integer::i

gridID=this%gridID
xc=this%xc(1)%coords

if(testOnOff) then
select case (gridID)
   case (1)
     call initIsoShockTube1D(this,q)
   case (2)
     call initSodShockTube1D(this,q)
   case (3)
     call  initWCShockTube1D(this,q)
   case (10)
     call initBrioWuShockTube1D(this,q) 
   case (11)
     call initRJ2aShockTube1D(this,q)
   case (22)
     call initCShockTest1Dn(this,q)
   case (23)
     call initCShockTest1Di(this,q)
   case (24)
     call initWardleInstabilityInitn(this,q)
   case (25)
     call initWardleInstabilityIniti(this,q)
   case (46)
     call initIsoShockTubeMHD1D(this,q)
   case (50)
     call initPolyShockTubeHD1D(this,q)
end select
endif

end subroutine init1d

