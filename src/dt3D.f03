subroutine dt3D(this,q)
use mpi
use gridModule
implicit none
class(grid)::this
double precision,dimension(1-this%nbuf:this%nMesh(1)+this%nbuf,1-this%nbuf:this%nMesh(2)+this%nbuf,& 
&                          1-this%nbuf:this%nMesh(3)+this%nbuf,this%nvar)::q
double precision::rho,vx,vy,vz,snd,vtot,CFL,wavespd,dt_temp,pressure,gam,ene
integer::i,j,k,nx,ny,nz,nbuf,eosType,solverType,ierr
double precision::bxl,byl,bzl,bxr,byr,bzr,bxc,byc,bzc
double precision::vsq,bsq,bmin,cfast
double precision::global_dt,sgftot,sgfx,sgfy,sgfz
double precision::dt_pressure,iene,sgpower
integer:: pos

nx=this%nMesh(1)
ny=this%nMesh(2)
nz=this%nMesh(3)
nbuf=this%nbuf
snd=this%snd
CFL=this%CFL
eosType=this%eosType
solverType=this%solverType
gam=this%adiGamma

dt_temp=1.d10
dt_pressure=1.d10

 if(eosType .eq. 1) then
  select case (solverType)
    case(1,2)
      do k=1,nz
       do j=1,ny
        do i=1,nx
           rho=q(i,j,k,1)
            vx=q(i,j,k,2)/rho
            vy=q(i,j,k,3)/rho
            vz=q(i,j,k,4)/rho
          vtot=dsqrt(vx**2+vy**2+vz**2)
          wavespd=vtot+snd
          dt_temp=dmin1(dt_temp,dmin1(dmin1(this%dx(1)%coords(i),this%dx(2)%coords(j)),&
                    & this%dx(3)%coords(k))/wavespd*CFL)
        enddo
       enddo
      enddo

    case(4,5)
       do k=1,nz
         do j=1,ny
           do i=1,nx
              rho=q(i,j,k,1)
               vx=q(i,j,k,2)/rho
               vy=q(i,j,k,3)/rho
               vz=q(i,j,k,4)/rho
              bxl=q(i,j,k,5)
              byl=q(i,j,k,6)
              bzl=q(i,j,k,7)
              ene=q(i,j,k,8)
              bxr=q(i,j,k,9)
              byr=q(i,j,k,10)
              bzr=q(i,j,k,11)
              bxc=0.5d0*(bxl+bxr)
              byc=0.5d0*(byl+byr)
              bzc=0.5d0*(bzl+bzr)
              vsq=vx**2+vy**2+vz**2
              vtot=dsqrt(vsq)
              bsq=(bxc**2+byc**2+bzc**2)/rho
              bmin=dmin1(dmin1(dabs(bxc),dabs(byc)),dabs(bzc))
              cfast=dsqrt(0.5d0*(snd**2+bsq+dsqrt((snd**2+bsq)**2-4.d0*snd**2*bmin**2/rho)))
              wavespd=vtot+cfast
      dt_temp=dmin1(dt_temp,dmin1(dmin1(this%dx(1)%coords(i),this%dx(2)%coords(j)),&
                    & this%dx(3)%coords(k))/wavespd*CFL)
           enddo
         enddo
       enddo
      

    case default
       print *,"dt3D.f03: cannot determine dt"
       stop    
  end select

 elseif(eosType .eq. 2) then
   select case (solverType)
     case (2,3)

        do k=1,nz
         do j=1,ny
          do i=1,nx
            rho=q(i,j,k,1)
             vx=q(i,j,k,2)/rho
             vy=q(i,j,k,3)/rho
             vz=q(i,j,k,4)/rho
             ene=q(i,j,k,5)
            vtot=dsqrt(vx**2+vy**2+vz**2)
            pressure=(gam-1.d0)*(ene-0.5d0*rho*(vx**2+vy**2+vz**2))
            wavespd=vtot+dsqrt(gam*pressure/rho)
        dt_temp=dmin1(dt_temp,dmin1(dmin1(this%dx(1)%coords(i),this%dx(2)%coords(j)),&
                    & this%dx(3)%coords(k))/wavespd*CFL)
  
          enddo
         enddo
        enddo

     case (4,5)
       do k=1,nz
         do j=1,ny
           do i=1,nx
              rho=q(i,j,k,1)
               vx=q(i,j,k,2)/rho
               vy=q(i,j,k,3)/rho
               vz=q(i,j,k,4)/rho
              bxl=q(i,j,k,5)
              byl=q(i,j,k,6)
              bzl=q(i,j,k,7)
              ene=q(i,j,k,8)
              bxr=q(i,j,k,9)
              byr=q(i,j,k,10)
              bzr=q(i,j,k,11)
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
      dt_temp=dmin1(dt_temp,dmin1(dmin1(this%dx(1)%coords(i),this%dx(2)%coords(j)),&
                    & this%dx(3)%coords(k))/wavespd*CFL)
           enddo
         enddo
       enddo
     case default
       print *,"dt3D.f03: cannot determine dt"
       stop
   end select
 endif !! if eosType 

!!!! added by hhwang to avoid large time-step for strong self-gravitating forces
 if(this%enable_sg) then
    do k=1,nz
      do j=1,ny
         do i=1,nx
          pos=(k-1)*nx*ny+(j-1)*nx+i
          sgfx=this%sgfx(pos)
          sgfy=this%sgfy(pos)
          sgfz=this%sgfz(pos)
          rho=q(i,j,k,1)
          vx=q(i,j,k,2)/rho
          vy=q(i,j,k,3)/rho
          vz=q(i,j,k,4)/rho
          vsq=vx**2+vy**2+vz**2
          vtot=dsqrt(vsq)

          sgftot=dsqrt(sgfx**2+sgfy**2+sgfz**2)
          dt_temp=dmin1(dt_temp,0.2d0*(-vtot/sgftot+dsqrt(vtot**2/sgftot**2+2.d0*dmin1(dmin1(this%dx(1)%coords(i)&
                       ,this%dx(2)%coords(j)),this%dx(3)%coords(k))/sgftot)))
        enddo
      enddo
    enddo
 endif

 call MPI_ALLREDUCE(dt_temp,global_dt,1,MPI_DOUBLE,MPI_MIN,MPI_COMM_WORLD,ierr)
 !print *,"dt3D.f03,myid=",myid,"local_dt=",dt_temp,"global_dt=",global_dt !! the reason of printing this is weird...might be a bug of compiler....  
 
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

end subroutine dt3D
