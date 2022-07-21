subroutine setBdry3D(this,q)
use gridModule
use testSuiteMPI
implicit none
class(grid)::this
double precision,dimension(1-this%nbuf:this%nMesh(1)+this%nbuf,1-this%nbuf:this%nMesh(2)+this%nbuf,&
& 1-this%nbuf:this%nMesh(3)+this%nbuf,this%nvar)::q
integer::gridID

gridID=this%gridID

if(testOnOff) then
select case (gridID)
   case(30)
     call bdrySodShockTube3D(this,q)
   case(31)
     call bdryIsoShockTube3D(this,q)
   case(32)
     call bdryBrioWuShockTube3D(this,q)
   case(33)
     call bdryFieldLoopAdvection(this,q)
   case(34)
     call bdryMHDBlastWave(this,q)
   case(35)
     call bdryselfgravityTest3D(this,q)
   case(48)
     call bdryIsoShockTubeMHD3D(this,q)
   case(49)
     call bdryIsoMHDsg3D(this,q)
   case(52)
     call bdryCShockTest3Dn(this,q)
   case(53)
     call bdryCShockTest3Di(this,q)
   case(199)
     call bdryAdiShockTubeMHD3D(this,q)
   case(202)
     call bdryHDRichtmyerMeshkovInstability3D(this,q)
   case(203)
     call bdryMHDRichtmyerMeshkovInstability3D(this,q)
   case(400)
     call bdryTurbulenceDriving3DHD(this,q)

end select
endif
end subroutine setBdry3D

