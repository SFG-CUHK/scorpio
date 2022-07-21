subroutine setBdry2D(this,q)
use gridModule
use testSuiteMPI
implicit none
class(grid)::this
double precision,dimension(1-this%nbuf:this%nMesh(1)+this%nbuf,1-this%nbuf:this%nMesh(2)+this%nbuf,this%nvar)::q
integer::gridID

gridID=this%gridID

if(testOnOff) then
select case (gridID)
   case(4)
     call bdryIsoShockTube2D(this,q)
   case(5)
     call bdrySodShockTube2D(this,q)
   case(6)
     call bdryWCShockTube2D(this,q)
   case(7)
     call bdryDoubleMachReflection(this,q)
   case(8)
     call bdryStrongRarefactionTest(this,q)
   case(9)
     call bdryShuOsherShocktubeTest(this,q)
   case(12)
     call bdryBrioWuShockTube2D(this,q)
   case(13)
     call bdryOrsagTangVortex(this,q)
   case(14)
     call bdrySphericalBlastWaveMHD2D(this,q)
   case(15)
     call bdryBFieldLoopTest2D(this,q)
   case(16)
     call bdryRotorMHD2D(this,q)
   case(17)
     call bdryCurrentSheetMHD2D(this,q)
   case(18)
     call bdryCPAlvenWave(this,q)
   case(19)
     call bdryLinearSlowWaveTestMHD2D(this,q)
   case(20)
     call bdryLinearFastWaveTestMHD2D(this,q)
   case(21)
     call bdrySelfgravityTest2D(this,q)
   case(26)
     call bdryWardleInstabilityn(this,q)
   case(27)
     call bdryWardleInstabilityi(this,q)
   case(122)
     call bdryCShockTest2Dn(this,q)
   case(123)
     call bdryCShockTest2Di(this,q)
   case(28)
     call bdryCoreCollapseAD2Dn(this,q)
   case(29)
     call bdryCoreCollapseAD2Di(this,q)
   case(36)
     call bdryHDBlastWavePolar2D(this,q)
   case(37)
     call bdryHDBlastWaveCart2D(this,q)
   case(38)
     call bdryHDBlastWavePolarIso2D(this,q)
   case(39)
     call bdryHDBlastWaveCartIso2D(this,q)
   case(40)
     call bdrySphericalBlastWaveMHDPolar2D(this,q)
   case(41)
     call bdryForceBalanceMHDPolar2D(this,q)
   case(42)
     call bdryKelvinHelmholtzInstabilityHD2D(this,q)
   case(43)
     call bdryKelvinHelmholtzInstabilityMHD2D(this,q)
   case(44)
     call bdryRayleighTaylorInstabilityHD2D(this,q)
   case(45)
     call bdryRayleighTaylorInstabilityMHD2D(this,q)
   case(47)
     call bdryOrsagTangVortexIso(this,q)
   case(51)
     call bdryPolyShockTubeHD2D(this,q)
   case(54)
     call bdrySphericalBlastWaveAD2Dn(this,q)
   case(55)
     call bdrySphericalBlastWaveAD2Di(this,q)
   case(200)
     call bdryHDRichtmyerMeshkovInstability2D(this,q)
   case(201)
     call bdryMHDRichtmyerMeshkovInstability2D(this,q)
   case(399)
     call bdryTurbulenceDriving2DHD(this,q)


end select
endif
end subroutine setBdry2D

