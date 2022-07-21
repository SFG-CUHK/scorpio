subroutine dt1D(this,q)
use mpi
use gridModule
implicit none
class(grid)::this
double precision,dimension(1-this%nbuf:this%nMesh(1)+this%nbuf,this%nvar)::q
double precision::rho,vx,snd,vtot,CFL,wavespd,dt_temp,pressure,gam,ene
double precision::vy,vz,bxl,byl,bzl,bxr,byr,bzr,bxc,byc,bzc,vsq,bsq,cfast,bmin
integer::i,nx,nbuf,eosType,solverType,ierr
double precision::global_dt

nx=this%nMesh(1)
nbuf=this%nbuf
snd=this%snd
CFL=this%CFL
eosType=this%eosType
solverType=this%solverType
gam=this%adiGamma

dt_temp=1.d10

 if(eosType .eq. 1) then
   select case (solverType)
    case (1)
      do i=1,nx
         rho=q(i,1)
          vx=q(i,2)/rho
        vtot=dsqrt(vx**2)
        wavespd=vtot+snd
        dt_temp=dmin1(dt_temp,this%dx(1)%coords(i)/wavespd*CFL)
      enddo
    case (4,5)
      do i=1,nx
          rho=q(i,1)
           vx=q(i,2)/rho
           vy=q(i,3)/rho
           vz=q(i,4)/rho
          bxl=q(i,5)
          byl=q(i,6)
          bzl=q(i,7)
          ene=0.d0
          bxr=q(i,9)
          byr=q(i,10)
          bzr=q(i,11)
          bxc=0.5d0*(bxl+bxr)
          byc=0.5d0*(byl+byr)
          bzc=0.5d0*(bzl+bzr)
          vsq=vx**2+vy**2+vz**2
          vtot=dsqrt(vsq)
          bsq=(bxc**2+byc**2+bzc**2)/rho
          bmin=dmin1(dmin1(dabs(bxc),dabs(byc)),dabs(bzc))
          cfast=dsqrt(0.5d0*(snd**2+bsq+dsqrt((snd**2+bsq)**2-4.d0*snd**2*bmin**2/rho)))
          wavespd=vtot+cfast
          dt_temp=dmin1(dt_temp,this%dx(1)%coords(i)/wavespd*CFL)
      enddo
    case default
         print *, 'dt1D.f03: dt undetermined....'
         stop
   end select
 elseif(eosType .eq. 2) then
    select case (solverType)
      case (2,3)
         do i=1,nx
            rho=q(i,1)
             vx=q(i,2)/rho
            ene=q(i,3)
           vtot=dsqrt(vx**2)
           pressure=(gam-1.d0)*(ene-0.5d0*rho*vx**2)
           wavespd=vtot+dsqrt(gam*pressure/rho)
           dt_temp=dmin1(dt_temp,this%dx(1)%coords(i)/wavespd*CFL)
         enddo
      case (4,5)
      do i=1,nx
          rho=q(i,1)
           vx=q(i,2)/rho
           vy=q(i,3)/rho
           vz=q(i,4)/rho
          bxl=q(i,5)
          byl=q(i,6)
          bzl=q(i,7)
          ene=q(i,8)
          bxr=q(i,9)
          byr=q(i,10)
          bzr=q(i,11)
          bxc=0.5d0*(bxl+bxr)
          byc=0.5d0*(byl+byr)
          bzc=0.5d0*(bzl+bzr)
          vsq=vx**2+vy**2+vz**2
          vtot=dsqrt(vsq)
          bsq=bxc**2+byc**2+bzc**2
          pressure=(gam-1.d0)*(ene-0.5d0*rho*vsq-0.5d0*bsq)
          bmin=dmin1(dmin1(dabs(bxc),dabs(byc)),dabs(bzc))
          cfast=dsqrt((gam*pressure+bsq+dsqrt((gam*pressure+bsq)**2.d0-4.d0*gam*pressure*bmin**2.d0))/(2.d0*rho))
          wavespd=vtot+cfast
          dt_temp=dmin1(dt_temp,this%dx(1)%coords(i)/wavespd*CFL)
      enddo
      case default
         print *, 'dt1D.f03: dt undetermined....'
         stop
    end select
 elseif(eosType .eq. 3) then
      do i=1,nx
         rho=q(i,1)
          vx=q(i,2)/rho
        vtot=dsqrt(vx**2)
        pressure=this%polyK*rho**this%polyGamma
        wavespd=vtot+dsqrt(this%polyGamma*pressure/rho)
        dt_temp=dmin1(dt_temp,this%dx(1)%coords(i)/wavespd*CFL)
      enddo    
 endif
   
call MPI_ALLREDUCE(dt_temp,global_dt,1,MPI_DOUBLE,MPI_MIN,MPI_COMM_WORLD,ierr)

if(global_dt .gt. this%toutput-this%t) then
  if(this%enable_ad .eqv. .false.) then
   global_dt = this%toutput-this%t
   this%writeFlag=.true.
   this%toutput=this%toutput+this%dtout
   this%fnum=this%fnum+1
  endif
elseif(global_dt .gt. this%tend-this%t) then
  if(this%enable_ad .eqv. .false.) then
   global_dt = this%tend-this%t
   this%writeFlag=.true.
   this%fnum=this%fnum+1
  endif
endif

!print *,"dt1D.f03, myid=",myid,"t=",this%t,"globaldt=",global_dt
   this%dt=global_dt
end subroutine
