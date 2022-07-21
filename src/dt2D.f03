subroutine dt2D(this,q)
use mpi
use gridModule
implicit none
class(grid)::this
double precision,dimension(1-this%nbuf:this%nMesh(1)+this%nbuf,1-this%nbuf:this%nMesh(2)+this%nbuf,this%nvar)::q
double precision::rho,vx,vy,vz,snd,vtot,CFL,wavespd,dt_temp,pressure,gam,ene
integer::i,j,nx,ny,nbuf,eosType,solverType,ierr,coordType
double precision::bxl,byl,bzl,bxr,byr,bzr,bxc,byc,bzc
double precision::vsq,bsq,bmin,cfast
double precision::global_dt
double precision::rc,dr,rdphi

nx=this%nMesh(1)
ny=this%nMesh(2)
nbuf=this%nbuf
snd=this%snd
CFL=this%CFL
eosType=this%eosType
solverType=this%solverType
gam=this%adiGamma
coordType=this%coordType

dt_temp=1.d10

 if(eosType .eq. 1) then
   if(coordType .eq. 1) then
     select case (solverType)
       case (1,2)
         do j=1,ny
           do i=1,nx
             rho=q(i,j,1)
              vx=q(i,j,2)/rho
              vy=q(i,j,3)/rho
            vtot=dsqrt(vx**2+vy**2)
            wavespd=vtot+snd
            dt_temp=dmin1(dt_temp,dmin1(this%dx(1)%coords(i),this%dx(2)%coords(j))/wavespd*CFL)
           enddo
         enddo
       case (4,5)
         do j=1,ny
           do i=1,nx
              rho=q(i,j,1)
               vx=q(i,j,2)/rho
               vy=q(i,j,3)/rho
               vz=q(i,j,4)/rho
              bxl=q(i,j,5)
              byl=q(i,j,6)
              bzl=q(i,j,7)
              ene=0.d0
              bxr=q(i,j,9)
              byr=q(i,j,10)
              bzr=q(i,j,11)
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
         enddo

       case default
         print *,"dt2D.f03: cannot determine dt"
         stop         
     end select
   endif  !! if coordType .eq. 1

   if(coordType .eq. 2 .or. coordType .eq. 3) then
     do j=1,ny
       do i=1,nx
            rc=this%xc(1)%coords(i)
            dr=this%dx(1)%coords(i)
         rdphi=rc*this%dx(2)%coords(j)
           rho=q(i,j,1)
            vx=q(i,j,2)/rho
            vy=q(i,j,3)/rho
          vtot=dsqrt(vx**2+vy**2)
          wavespd=vtot+snd
          dt_temp=dmin1(dt_temp,dmin1(dr,rdphi)/wavespd*CFL)
       enddo
     enddo
   endif

 elseif(eosType .eq. 2) then
   select case (solverType)
     case (2,3)

       if(coordType .eq. 1) then !! if Cartesian coordinate
         do j=1,ny
          do i=1,nx
           rho=q(i,j,1)
            vx=q(i,j,2)/rho
            vy=q(i,j,3)/rho
            ene=q(i,j,4)
           vtot=dsqrt(vx**2+vy**2)
           pressure=(gam-1.d0)*(ene-0.5d0*rho*(vx**2+vy**2))
           wavespd=vtot+dsqrt(gam*pressure/rho)
           dt_temp=dmin1(dt_temp,dmin1(this%dx(1)%coords(i),this%dx(2)%coords(j))/wavespd*CFL)
          enddo
         enddo
        endif !! if coordType .eq .1

      if(coordType .eq. 2 .or. coordType .eq. 3) then
        do j=1,ny
          do i=1,nx
            rho=q(i,j,1)
             vx=q(i,j,2)/rho  !! vr
             vy=q(i,j,3)/rho  !! vphi
            ene=q(i,j,4)
           vtot=dsqrt(vx**2+vy**2)
           pressure=(gam-1.d0)*(ene-0.5d0*rho*(vx**2+vy**2))
           wavespd=vtot+dsqrt(gam*pressure/rho)
           dt_temp=dmin1(dt_temp,dmin1(this%dx(1)%coords(i),this%xc(1)%coords(i)*this%dx(2)%coords(j))/wavespd*CFL)
          enddo
        enddo
      endif  !! if polar coordinates

     case (4,5)
         do j=1,ny
           do i=1,nx
              rho=q(i,j,1)
               vx=q(i,j,2)/rho
               vy=q(i,j,3)/rho
               vz=q(i,j,4)/rho
              bxl=q(i,j,5)
              byl=q(i,j,6)
              bzl=q(i,j,7)
              ene=q(i,j,8)
              bxr=q(i,j,9)
              byr=q(i,j,10)
              bzr=q(i,j,11)
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
         enddo
 
     case default
       print *,"dt2D.f03: cannot determine dt"
       stop
   end select
 elseif(eostype .eq. 3) then
         do j=1,ny
           do i=1,nx
             rho=q(i,j,1)
              vx=q(i,j,2)/rho
              vy=q(i,j,3)/rho
            vtot=dsqrt(vx**2+vy**2)
            pressure=this%polyK*rho**this%polyGamma
            wavespd=vtot+dsqrt(this%polyGamma*pressure/rho)
            dt_temp=dmin1(dt_temp,dmin1(this%dx(1)%coords(i),this%dx(2)%coords(j))/wavespd*CFL)
           enddo
         enddo   
 endif

 call MPI_ALLREDUCE(dt_temp,global_dt,1,MPI_DOUBLE,MPI_MIN,MPI_COMM_WORLD,ierr)
 !print *,"dt2D.f03,myid=",myid,"local_dt=",dt_temp,"global_dt=",global_dt !! the reason of printing this is weird...might be a bug of compiler....  
 
if(global_dt .gt. this%toutput-this%t) then
 if(this%enable_ad .eqv. .false. ) then
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

   this%dt=global_dt

end subroutine dt2D
