Module riemannSolverModule
use gridModule
use limiterModule
implicit none

abstract interface
   subroutine riemannSolver3D(this,q,q1,q2,dd)
     import grid
     implicit none
     class(grid)::this
     double precision,dimension(1-this%nbuf:this%nMesh(1)+this%nbuf,1-this%nbuf:this%nMesh(2)+this%nbuf,&
&    1-this%nbuf:this%nMesh(3)+this%nbuf,this%nvar)::q,q1,q2
     integer:: dd
   end subroutine riemannSolver3D

   subroutine riemannSolver2D(this,q,q1,q2,dd)
     import grid
     implicit none
     class(grid)::this
     double precision,dimension(1-this%nbuf:this%nMesh(1)+this%nbuf,1-this%nbuf:this%nMesh(2)+this%nbuf,this%nvar)::q,q1,q2
     integer::dd
   end subroutine riemannSolver2D

   subroutine riemannSolver1D(this,q,q1,q2,dd)
     import grid
     implicit none
     class(grid)::this
     double precision,dimension(1-this%nbuf:this%nMesh(1)+this%nbuf,this%nvar)::q,q1,q2
     integer::dd
   end subroutine riemannSolver1D

  subroutine fluxSolver(ql,qr,slopeL,slopeR,snd,flux,nvar)
     implicit none
     integer::nvar
     double precision::ql(nvar),qr(nvar),slopeL(nvar),slopeR(nvar),flux(nvar),snd
  end subroutine fluxSolver
end interface



contains


subroutine solverIsoMHD3D(this,q,q1,q2,dd)
class(grid)::this
double precision,dimension(1-this%nbuf:this%nMesh(1)+this%nbuf,1-this%nbuf:this%nMesh(2)+this%nbuf,&
                           1-this%nbuf:this%nMesh(3)+this%nbuf,this%nvar)::q,q1,q2
integer::nx,ny,nz,nbuf,nvar,dd
integer::i,j,k
integer::cx,cy,cz
integer::f1,f2,f3,f4,f5,f6,f7,f8,f9,f10,f11
double precision::g1,g2,g3,g4,g5,g6,g7,g8,g9,g10,g11
double precision::dt,dx,dtddx,gam
double precision::rhoL,rhoM,rhoR
double precision::vxL,vxM,vxR
double precision::vyL,vyM,vyR
double precision::vzL,vzM,vzR
double precision::bxL,bxM,bxR,bxtmp
double precision::byL,byM,byR
double precision::bzL,bzM,bzR
double precision::eneL,eneM,eneR
double precision::pL ,pM ,pR
double precision::pressure,bxc,byc,bzc,vxc,vyc,vzc,denc,enec
double precision,dimension(:,:,:,:),allocatable::SL
double precision,dimension(:,:,:,:),allocatable::Fx
double precision,dimension(:,:,:),allocatable::Ez,Ezc,Ex,Exc,Ey,Eyc
double precision::ql(this%nvar),qr(this%nvar),slopeL(this%nvar),slopeR(this%nvar),flux(this%nvar),coef(this%nvar)
procedure(limiter),pointer::slope=>null()
procedure(fluxSolver),pointer::fluxPtr=>null()
double precision::rho,vx,vy,vz,bx,by,bz
integer::chgFluxPtrCount

nx=this%nMesh(1)
ny=this%nMesh(2)
nz=this%nMesh(3)

nbuf=this%nbuf
nvar=this%nvar
dt=this%dt
gam=this%adiGamma
chgFluxPtrCount=0

select case(this%solverType)
   case (4)
     fluxPtr=>fluxHLLIsoMHD1D
   case (5)
     fluxPtr=>fluxHLLDIsoMHD1D
   case default
     print *,'riemannSolverModule.f03: no appropriate 2D MHD solver found!'
end select

select case (this%limiterType)
   case(0)
     slope=>zslop
   case(1)
     slope=>vslop
   case(2)
     slope=>fslop
   case(3)
     slope=>minmod
   case default
     print *, "riemannSolverModule.f03: no appropriate slope limiter found!"
end select

allocate(SL(1-nbuf:nx+nbuf,1-nbuf:ny+nbuf,1-nbuf:nz+nbuf,nvar))
allocate(Fx(1-nbuf:nx+nbuf,1-nbuf:ny+nbuf,1-nbuf:nz+nbuf,nvar))
allocate( Ez(1-nbuf:nx+nbuf,1-nbuf:ny+nbuf,1-nbuf:nz+nbuf))
allocate(Ezc(1-nbuf:nx+nbuf,1-nbuf:ny+nbuf,1-nbuf:nz+nbuf))
allocate( Ex(1-nbuf:nx+nbuf,1-nbuf:ny+nbuf,1-nbuf:nz+nbuf))
allocate(Exc(1-nbuf:nx+nbuf,1-nbuf:ny+nbuf,1-nbuf:nz+nbuf))
allocate( Ey(1-nbuf:nx+nbuf,1-nbuf:ny+nbuf,1-nbuf:nz+nbuf))
allocate(Eyc(1-nbuf:nx+nbuf,1-nbuf:ny+nbuf,1-nbuf:nz+nbuf))

!!!!! dd=1 ==>x, dd=2 ==>y, dd=3 ==>z !!!!!

if(dd .eq. 1) then
     cx=1
     cy=0
     cz=0
     f1=1
     f2=2
     f3=3
     f4=4
     f5=5
     f6=6
     f7=7
     f8=8
     f9=9
    f10=10
    f11=11
 coef(1)=1.d0
 coef(2)=1.d0
 coef(3)=1.d0
 coef(4)=1.d0
 coef(5)=1.d0
 coef(6)=1.d0
 coef(7)=1.d0
 coef(8)=1.d0
     g1=1.d0
     g2=1.d0
     g3=1.d0
     g4=1.d0
     g5=1.d0
     g6=1.d0
     g7=1.d0
     g8=1.d0
     g9=1.d0
    g10=1.d0
    g11=1.d0
elseif (dd .eq. 2) then
     cx=0
     cy=1
     cz=0
     f1=1
     f2=3
     f3=2
     f4=4
     f5=6
     f6=5
     f7=7
     f8=8
     f9=10
    f10=9
    f11=11
 coef(1)=1.d0
 coef(2)=-1.d0
 coef(3)=1.d0
 coef(4)=1.d0
 coef(5)=-1.d0
 coef(6)=1.d0
 coef(7)=1.d0
 coef(8)=1.d0
     g1=1.d0
     g2=1.d0
     g3=-1.d0
     g4=1.d0
     g5=1.d0
     g6=-1.d0
     g7=1.d0
     g8=1.d0
     g9=1.d0
    g10=-1.d0
    g11=1.d0
elseif (dd .eq. 3) then
     cx=0
     cy=0
     cz=1
     f1=1
     f2=4
     f3=3
     f4=2
     f5=7
     f6=6
     f7=5
     f8=8
     f9=11
    f10=10
    f11=9
 coef(1)=1.d0
 coef(2)=-1.d0
 coef(3)=1.d0
 coef(4)=1.d0
 coef(5)=-1.d0
 coef(6)=1.d0
 coef(7)=1.d0
 coef(8)=1.d0
     g1=1.d0
     g2=1.d0
     g3=1.d0
     g4=-1.d0
     g5=1.d0
     g6=1.d0
     g7=-1.d0
     g8=1.d0
     g9=1.d0
    g10=1.d0
    g11=-1.d0
endif

do k=0,nz+1
  do j=0,ny+1
    do i=0,nx+1
      rhoL=g1*q(i-cx,j-cy,k-cz,f1)
      rhoM=g1*q(i   ,j   ,k   ,f1)
      rhoR=g1*q(i+cx,j+cy,k+cz,f1)

       vxL=g2*q(i-cx,j-cy,k-cz,f2)/rhoL
       vxM=g2*q(i   ,j   ,k   ,f2)/rhoM
       vxR=g2*q(i+cx,j+cy,k+cz,f2)/rhoR

       vyL=g3*q(i-cx,j-cy,k-cz,f3)/rhoL
       vyM=g3*q(i   ,j   ,k   ,f3)/rhoM
       vyR=g3*q(i+cx,j+cy,k+cz,f3)/rhoR

       vzL=g4*q(i-cx,j-cy,k-cz,f4)/rhoL
       vzM=g4*q(i   ,j   ,k   ,f4)/rhoM
       vzR=g4*q(i+cx,j+cy,k+cz,f4)/rhoR

       bxL=0.5d0*(g5*q(i-cx,j-cy,k-cz,f5)+ g9*q(i-cx,j-cy,k-cz,f9 ))
       bxM=0.5d0*(g5*q(i   ,j   ,k   ,f5)+ g9*q(i   ,j   ,k   ,f9 ))
       bxR=0.5d0*(g5*q(i+cx,j+cy,k+cz,f5)+ g9*q(i+cx,j+cy,k+cz,f9 ))

       byL=0.5d0*(g6*q(i-cx,j-cy,k-cz,f6)+g10*q(i-cx,j-cy,k-cz,f10))
       byM=0.5d0*(g6*q(i   ,j   ,k   ,f6)+g10*q(i   ,j   ,k   ,f10))
       byR=0.5d0*(g6*q(i+cx,j+cy,k+cz,f6)+g10*q(i+cx,j+cy,k+cz,f10))

       bzL=0.5d0*(g7*q(i-cx,j-cy,k-cz,f7)+g11*q(i-cx,j-cy,k-cz,f11))
       bzM=0.5d0*(g7*q(i   ,j   ,k   ,f7)+g11*q(i   ,j   ,k   ,f11))
       bzR=0.5d0*(g7*q(i+cx,j+cy,k+cz,f7)+g11*q(i+cx,j+cy,k+cz,f11))

      SL(i,j,k,1)=slope(rhoL,rhoM,rhoR)
      SL(i,j,k,2)=slope( vxL, vxM, vxR)
      SL(i,j,k,3)=slope( vyL, vyM, vyR)
      SL(i,j,k,4)=slope( vzL, vzM, vzR)
      SL(i,j,k,5)=slope( bxL, bxM, bxR)
      SL(i,j,k,6)=slope( byL, byM, byR)
      SL(i,j,k,7)=slope( bzL, bzM, bzR)
      SL(i,j,k,8)=0.d0
    enddo !! end do i
  enddo !! end do j
enddo !! endo k

   SL(:,:,:,5)=0.d0

do k=0,nz+1
  do j=0,ny+1
    do i=0,nx+1
       ql(1)=g1*q(i,j,k,f1)
       ql(2)=g2*q(i,j,k,f2)/ql(1)
       ql(3)=g3*q(i,j,k,f3)/ql(1)
       ql(4)=g4*q(i,j,k,f4)/ql(1)
       ql(5)=0.5d0*(g5*q(i,j,k,f5)+ g9*q(i,j,k, f9))
       ql(6)=0.5d0*(g6*q(i,j,k,f6)+g10*q(i,j,k,f10))
       ql(7)=0.5d0*(g7*q(i,j,k,f7)+g11*q(i,j,k,f11))
       ql(8)=0.d0

       qr(1)=g1*q(i+cx,j+cy,k+cz,f1)
       qr(2)=g2*q(i+cx,j+cy,k+cz,f2)/qr(1)
       qr(3)=g3*q(i+cx,j+cy,k+cz,f3)/qr(1)
       qr(4)=g4*q(i+cx,j+cy,k+cz,f4)/qr(1)
       qr(5)=0.5d0*(g5*q(i+cx,j+cy,k+cz,f5)+ g9*q(i+cx,j+cy,k+cz, f9))
       qr(6)=0.5d0*(g6*q(i+cx,j+cy,k+cz,f6)+g10*q(i+cx,j+cy,k+cz,f10))
       qr(7)=0.5d0*(g7*q(i+cx,j+cy,k+cz,f7)+g11*q(i+cx,j+cy,k+cz,f11))
       qr(8)=0.d0

       ql(5)=g5*q(i,j,k,f9)
       qr(5)=g5*q(i+cx,j+cy,k+cz,f5)

       slopeL(:)=SL(i   ,j   ,k   ,:)
       slopeR(:)=SL(i+cx,j+cy,k+cz,:)

       call fluxPtr(ql,qr,slopeL,slopeR,this%snd,flux,this%nvar)
       Fx(i,j,k,:)=flux(:)  
 
    enddo !! end do i
  enddo !!end do j
enddo !! end do k   

Ez=0.d0
Ex=0.d0
Ey=0.d0
Ezc=0.d0
Exc=0.d0
Eyc=0.d0


!!!!! calculate EMF at the cell centers
!!!!! Ezc,Exc,Eyc

do k=1-nbuf,nz+nbuf
  do j=1-nbuf,ny+nbuf
    do i=1-nbuf,nx+nbuf
      rho=q(i,j,k,1)
       vx=q(i,j,k,2)/rho
       vy=q(i,j,k,3)/rho
       vz=q(i,j,k,4)/rho
       bx=0.5d0*(q(i,j,k,5)+q(i,j,k,9 ))
       by=0.5d0*(q(i,j,k,6)+q(i,j,k,10))
       bz=0.5d0*(q(i,j,k,7)+q(i,j,k,11))
       Ezc(i,j,k)=vy*bx-by*vx
       Exc(i,j,k)=vz*by-vy*bz
       Eyc(i,j,k)=vx*bz-vz*bx
    enddo !! enddo i
  enddo  !! enddo j
enddo  !! enddo k


!!!!!!!! Gardiner & Stone, JCP, 2005, 205, 509
if(dd .eq. 1) then

do k=0,nz
  do j=0,ny
    do i=0,nx
       Ez(i,j,k)= 0.5d0*(-Fx(i,j,k,f6)-Fx(i,j+1,k,f6))*coef(6)
       Ez(i,j,k)=Ez(i,j,k)-0.25d0*(Ezc(i,j,k)+Ezc(i,j+1,k))

       Ey(i,j,k)= 0.5d0*( Fx(i,j,k,f7)+Fx(i,j,k+1,f7))*coef(7)
       Ey(i,j,k)=Ey(i,j,k)-0.25d0*(Eyc(i,j,k)+Eyc(i,j,k+1))
    enddo
  enddo
enddo

elseif (dd .eq. 2) then
do k=0,nz
  do j=0,ny
    do i=0,nx
       Ez(i,j,k)=0.5d0*( Fx(i+1,j,k,f5)+Fx(i,j,k,f5))*coef(5)
       Ez(i,j,k)=Ez(i,j,k)-0.25d0*(Ezc(i+1,j,k)+Ezc(i+1,j+1,k))

       Ex(i,j,k)=0.5d0*(-Fx(i,j,k+1,f7)-Fx(i,j,k,f7))*coef(7)
       Ex(i,j,k)=Ex(i,j,k)-0.25d0*(Exc(i,j,k)+Exc(i,j+1,k))
    enddo
  enddo
enddo

elseif (dd .eq. 3) then
do k=0,nz
  do j=0,ny
    do i=0,nx
      Ey(i,j,k)=0.5d0*(-Fx(i+1,j,k,f5)-Fx(i,j,k,f5))*coef(5)
      Ey(i,j,k)=Ey(i,j,k)-0.25d0*(Eyc(i+1,j,k)+Eyc(i+1,j,k+1))

      Ex(i,j,k)=0.5d0*( Fx(i,j+1,k,f6)+Fx(i,j,k,f6))*coef(6)
      Ex(i,j,k)=Ex(i,j,k)-0.25d0*(Exc(i,j,k+1)+Exc(i,j+1,k+1))
    enddo
  enddo
enddo

endif !! if( dd .eq. 1 )


do k=1,nz
  do j=1,ny
    do i=1,nx
      dtddx=this%dt/(this%dx(1)%coords(i)*cx+this%dx(2)%coords(j)*cy+this%dx(3)%coords(k)*cz)
      q2(i,j,k,1)=q1(i,j,k,1)+coef(1)*(Fx(i-cx,j-cy,k-cz,f1)-Fx(i,j,k,f1))*dtddx
      q2(i,j,k,2)=q1(i,j,k,2)+coef(2)*(Fx(i-cx,j-cy,k-cz,f2)-Fx(i,j,k,f2))*dtddx
      q2(i,j,k,3)=q1(i,j,k,3)+coef(3)*(Fx(i-cx,j-cy,k-cz,f3)-Fx(i,j,k,f3))*dtddx
      q2(i,j,k,4)=q1(i,j,k,4)+coef(4)*(Fx(i-cx,j-cy,k-cz,f4)-Fx(i,j,k,f4))*dtddx
      q2(i,j,k,5)=q1(i,j,k,5)+(Ez(i-1,j-1,k)-Ez(i-1,j  ,k))*dtddx+(Ey(i-1,j,k)-Ey(i-1,j,k-1))*dtddx
      q2(i,j,k,6)=q1(i,j,k,6)+(Ez(i  ,j-1,k)-Ez(i-1,j-1,k))*dtddx+(Ex(i,j-1,k-1)-Ex(i,j-1,k))*dtddx
      q2(i,j,k,7)=q1(i,j,k,7)+(Ey(i-1,j,k-1)-Ey(i  ,j,k-1))*dtddx+(Ex(i,j,k-1)-Ex(i,j-1,k-1))*dtddx
      q2(i,j,k,8)=0.d0
      q2(i,j,k,9)=q1(i,j,k,9)+(Ez(i,j-1,k)-Ez(i,j,k))*dtddx+(Ey(i,j,k)-Ey(i,j,k-1))*dtddx
      q2(i,j,k,10)=q1(i,j,k,10)+(Ez(i,j,k)-Ez(i-1,j,k))*dtddx+(Ex(i,j,k-1)-Ex(i,j,k))*dtddx
      q2(i,j,k,11)=q1(i,j,k,11)+(Ey(i-1,j,k)-Ey(i,j,k))*dtddx+(Ex(i,j,k)-Ex(i,j-1,k))*dtddx
    enddo
  enddo
enddo


deallocate(SL,Fx,Ez,Ezc,Ex,Exc,Ey,Eyc)
end subroutine solverIsoMHD3D

subroutine solverIsoMHD2D(this,q,q1,q2,dd)
class(grid)::this
double precision,dimension(1-this%nbuf:this%nMesh(1)+this%nbuf,1-this%nbuf:this%nMesh(2)+this%nbuf,this%nvar)::q,q1,q2
integer::nx,ny,nbuf,nvar,dd
integer::i,j
integer::cx,cy
integer::f1,f2,f3,f4,f5,f6,f7,f8,f9,f10,f11
double precision::g1,g2,g3,g4,g5,g6,g7,g8,g9,g10,g11
double precision::dt,dx,dtddx,gam
double precision::rhoL,rhoM,rhoR
double precision::vxL,vxM,vxR
double precision::vyL,vyM,vyR
double precision::vzL,vzM,vzR
double precision::bxL,bxM,bxR,bxtmp
double precision::byL,byM,byR
double precision::bzL,bzM,bzR
double precision::eneL,eneM,eneR
double precision::pL ,pM ,pR
double precision::pressure,bxc,byc,bzc,vxc,vyc,vzc,denc,enec
double precision::ri,rm,ro,dr,dphi,fac1,fac2
double precision,dimension(:,:,:),allocatable::SL
double precision,dimension(:,:,:),allocatable::Fx
double precision,dimension(:,:),allocatable::Ez,Ezc
double precision::ql(this%nvar),qr(this%nvar),slopeL(this%nvar),slopeR(this%nvar),flux(this%nvar),coef(this%nvar)
procedure(limiter),pointer::slope=>null()
procedure(limiter3),pointer::slope3=>null()
procedure(fluxSolver),pointer::fluxPtr=>null()
double precision::rho,vx,vy,bx,by,vr,vphi,ene,ptot,vz,brc,bphic
integer::chgFluxPtrCount
integer::coordType

nx=this%nMesh(1)
ny=this%nMesh(2)
nbuf=this%nbuf
nvar=this%nvar
dt=this%dt
gam=this%adiGamma
chgFluxPtrCount=0
coordType=this%coordType

select case(this%solverType)
   case (4)
     fluxPtr=>fluxHLLIsoMHD1D
   case (5)
     fluxPtr=>fluxHLLDIsoMHD1D
   case default
     print *,'riemannSolverModule.f03: no appropriate 2D MHD solver found!'
end select

select case (this%limiterType)
   case(0)
     slope=>zslop
   case(1)
     slope=>vslop
   case(2)
     slope=>fslop
   case(3)
     slope=>minmod
   case default
     print *, "riemannSolverModule.f03: no appropriate slope limiter found!"
end select

if(coordType .eq. 2 .or. coordType .eq. 3) then
    select case (this%limiterType)
      case(0)
        slope3=>zslop3
      case(1)
        print *,"solverAdi2D: not yet implemented...."
        stop
      case(2)
        print *,"solverAdi2D: not yet implemented...."
        stop
      case(3)
        slope3=>minmod3
    end select  
endif

allocate(SL(1-nbuf:nx+nbuf,1-nbuf:ny+nbuf,nvar))
allocate(Fx(1-nbuf:nx+nbuf,1-nbuf:ny+nbuf,nvar))
allocate(Ez(1-nbuf:nx+nbuf,1-nbuf:ny+nbuf))
allocate(Ezc(1-nbuf:nx+nbuf,1-nbuf:ny+nbuf))

!!!!! dd=1 ==>x, dd=2 ==>y !!!!!

if(dd .eq. 1) then
     cx=1
     cy=0
     f1=1
     f2=2
     f3=3
     f4=4
     f5=5
     f6=6
     f7=7
     f8=8
     f9=9
    f10=10
    f11=11
 coef(1)=1.d0
 coef(2)=1.d0
 coef(3)=1.d0
 coef(4)=1.d0
 coef(5)=1.d0
 coef(6)=1.d0
 coef(7)=1.d0
 coef(8)=1.d0
     g1=1.d0
     g2=1.d0
     g3=1.d0
     g4=1.d0
     g5=1.d0
     g6=1.d0
     g7=1.d0
     g8=1.d0
     g9=1.d0
    g10=1.d0
    g11=1.d0
elseif (dd .eq. 2) then
     cx=0
     cy=1
     f1=1
     f2=3
     f3=2
     f4=4
     f5=6
     f6=5
     f7=7
     f8=8
     f9=10
    f10=9
    f11=11
 coef(1)=1.d0
 coef(2)=-1.d0
 coef(3)=1.d0
 coef(4)=1.d0
 coef(5)=-1.d0
 coef(6)=1.d0
 coef(7)=1.d0
 coef(8)=1.d0
     g1=1.d0
     g2=1.d0
     g3=-1.d0
     g4=1.d0
     g5=1.d0
     g6=-1.d0
     g7=1.d0
     g8=1.d0
     g9=1.d0
    g10=-1.d0
    g11=1.d0
endif

if(coordType .eq. 1) then  !! Cartesian coordinates
do j=0,ny+1
  do i=0,nx+1
    rhoL=g1*q(i-cx,j-cy,f1)
    rhoM=g1*q(i   ,j   ,f1)
    rhoR=g1*q(i+cx,j+cy,f1)

     vxL=g2*q(i-cx,j-cy,f2)/rhoL
     vxM=g2*q(i   ,j   ,f2)/rhoM
     vxR=g2*q(i+cx,j+cy,f2)/rhoR

     vyL=g3*q(i-cx,j-cy,f3)/rhoL
     vyM=g3*q(i   ,j   ,f3)/rhoM
     vyR=g3*q(i+cx,j+cy,f3)/rhoR

     vzL=g4*q(i-cx,j-cy,f4)/rhoL
     vzM=g4*q(i   ,j   ,f4)/rhoM
     vzR=g4*q(i+cx,j+cy,f4)/rhoR

     bxL=0.5d0*(g5*q(i-cx,j-cy,f5)+ g9*q(i-cx,j-cy,f9 ))
     bxM=0.5d0*(g5*q(i   ,j   ,f5)+ g9*q(i   ,j   ,f9 ))
     bxR=0.5d0*(g5*q(i+cx,j+cy,f5)+ g9*q(i+cx,j+cy,f9 ))

     byL=0.5d0*(g6*q(i-cx,j-cy,f6)+g10*q(i-cx,j-cy,f10))
     byM=0.5d0*(g6*q(i   ,j   ,f6)+g10*q(i   ,j   ,f10))
     byR=0.5d0*(g6*q(i+cx,j+cy,f6)+g10*q(i+cx,j+cy,f10))

     bzL=0.5d0*(g7*q(i-cx,j-cy,f7)+g11*q(i-cx,j-cy,f11))
     bzM=0.5d0*(g7*q(i   ,j   ,f7)+g11*q(i   ,j   ,f11))
     bzR=0.5d0*(g7*q(i+cx,j+cy,f7)+g11*q(i+cx,j+cy,f11))

    SL(i,j,1)=slope(rhoL,rhoM,rhoR)
    SL(i,j,2)=slope( vxL, vxM, vxR)
    SL(i,j,3)=slope( vyL, vyM, vyR)
    SL(i,j,4)=slope( vzL, vzM, vzR)
    SL(i,j,5)=slope( bxL, bxM, bxR)
    SL(i,j,6)=slope( byL, byM, byR)
    SL(i,j,7)=slope( bzL, bzM, bzR)
    SL(i,j,8)=0.d0

  enddo !! end do i
enddo !! end do j
endif !! if coordType .eq. 1

if(coordType .eq. 3 .or. coordType .eq. 2) then !! polar coordinates, uniform or log in r
  if(dd .eq. 1) then
    do j=0,ny+1
      do i=0,nx+1
        rhoL=g1*q(i-cx,j-cy,f1)
        rhoM=g1*q(i   ,j   ,f1)
        rhoR=g1*q(i+cx,j+cy,f1)

         vxL=g2*q(i-cx,j-cy,f2)/rhoL
         vxM=g2*q(i   ,j   ,f2)/rhoM
         vxR=g2*q(i+cx,j+cy,f2)/rhoR

         vyL=g3*q(i-cx,j-cy,f3)/rhoL
         vyM=g3*q(i   ,j   ,f3)/rhoM
         vyR=g3*q(i+cx,j+cy,f3)/rhoR

         vzL=g4*q(i-cx,j-cy,f4)/rhoL
         vzM=g4*q(i   ,j   ,f4)/rhoM
         vzR=g4*q(i+cx,j+cy,f4)/rhoR

          ri=this%xl(1)%coords(i-cx)
          rm=this%xc(1)%coords(i-cx)
          ro=this%xr(1)%coords(i-cx)
         bxL=0.5d0*(ri*g5*q(i-cx,j-cy,f5)+ ro*g9*q(i-cx,j-cy,f9 ))/rm

          ri=this%xl(1)%coords(i)
          rm=this%xc(1)%coords(i)
          ro=this%xr(1)%coords(i)
         bxM=0.5d0*(ri*g5*q(i   ,j   ,f5)+ ro*g9*q(i   ,j   ,f9 ))/rm

          ri=this%xl(1)%coords(i+cx)
          rm=this%xc(1)%coords(i+cx)
          ro=this%xr(1)%coords(i+cx)
         bxR=0.5d0*(ri*g5*q(i+cx,j+cy,f5)+ ro*g9*q(i+cx,j+cy,f9 ))/rm

         byL=0.5d0*(g6*q(i-cx,j-cy,f6)+g10*q(i-cx,j-cy,f10))
         byM=0.5d0*(g6*q(i   ,j   ,f6)+g10*q(i   ,j   ,f10))
         byR=0.5d0*(g6*q(i+cx,j+cy,f6)+g10*q(i+cx,j+cy,f10))

         bzL=0.5d0*(g7*q(i-cx,j-cy,f7)+g11*q(i-cx,j-cy,f11))
         bzM=0.5d0*(g7*q(i   ,j   ,f7)+g11*q(i   ,j   ,f11))
         bzR=0.5d0*(g7*q(i+cx,j+cy,f7)+g11*q(i+cx,j+cy,f11))

        SL(i,j,1)=slope(rhoL,rhoM,rhoR)
        SL(i,j,2)=slope( vxL, vxM, vxR)
        SL(i,j,3)=slope( vyL, vyM, vyR)
        SL(i,j,4)=slope( vzL, vzM, vzR)
        SL(i,j,5)=slope( bxL, bxM, bxR)
        SL(i,j,6)=slope( byL, byM, byR)
        SL(i,j,7)=slope( bzL, bzM, bzR)
        SL(i,j,8)=0.d0

        ri=this%xl(1)%coords(i-cx)
        rm=this%xc(1)%coords(i)
        ro=this%xr(1)%coords(i+cx)
        dr=this%dx(1)%coords(i)

        !SL(i,j,1)=slope3(rhoL,rhoM,rhoR,ri,rm,ro,dr)
        !SL(i,j,2)=slope3( vxL, vxM, vxR,ri,rm,ro,dr)
        !SL(i,j,3)=slope3( vyL, vyM, vyR,ri,rm,ro,dr)
        !SL(i,j,4)=slope3( vzL, vzM, vzR,ri,rm,ro,dr)
        !SL(i,j,5)=slope3( bxL, bxM, bxR,ri,rm,ro,dr)
        !SL(i,j,6)=slope3( byL, byM, byR,ri,rm,ro,dr)
        !SL(i,j,7)=slope3( bzL, bzM, bzR,ri,rm,ro,dr)
        !SL(i,j,8)=slope3(  pL,  pM,  pR,ri,rm,ro,dr)
       enddo !! end do i
    enddo !! end do j        
  endif !! dd .eq. 1 , in r direction

  if(dd .eq. 2) then
    do j=0,ny+1
      do i=0,nx+1
        rhoL=g1*q(i-cx,j-cy,f1)
        rhoM=g1*q(i   ,j   ,f1)
        rhoR=g1*q(i+cx,j+cy,f1)

         vxL=g2*q(i-cx,j-cy,f2)/rhoL
         vxM=g2*q(i   ,j   ,f2)/rhoM
         vxR=g2*q(i+cx,j+cy,f2)/rhoR

         vyL=g3*q(i-cx,j-cy,f3)/rhoL
         vyM=g3*q(i   ,j   ,f3)/rhoM
         vyR=g3*q(i+cx,j+cy,f3)/rhoR

         vzL=g4*q(i-cx,j-cy,f4)/rhoL
         vzM=g4*q(i   ,j   ,f4)/rhoM
         vzR=g4*q(i+cx,j+cy,f4)/rhoR

         bxL=0.5d0*(g5*q(i-cx,j-cy,f5)+ g9*q(i-cx,j-cy,f9 ))
         bxM=0.5d0*(g5*q(i   ,j   ,f5)+ g9*q(i   ,j   ,f9 ))
         bxR=0.5d0*(g5*q(i+cx,j+cy,f5)+ g9*q(i+cx,j+cy,f9 ))

          ri=this%xl(1)%coords(i) 
          rm=this%xc(1)%coords(i)
          ro=this%xr(1)%coords(i)
         byL=0.5d0*(ri*g6*q(i-cx,j-cy,f6)+ro*g10*q(i-cx,j-cy,f10))/rm
         byM=0.5d0*(ri*g6*q(i   ,j   ,f6)+ro*g10*q(i   ,j   ,f10))/rm
         byR=0.5d0*(ri*g6*q(i+cx,j+cy,f6)+ro*g10*q(i+cx,j+cy,f10))/rm

         bzL=0.5d0*(g7*q(i-cx,j-cy,f7)+g11*q(i-cx,j-cy,f11))
         bzM=0.5d0*(g7*q(i   ,j   ,f7)+g11*q(i   ,j   ,f11))
         bzR=0.5d0*(g7*q(i+cx,j+cy,f7)+g11*q(i+cx,j+cy,f11))

        SL(i,j,1)=slope(rhoL,rhoM,rhoR)
        SL(i,j,2)=slope( vxL, vxM, vxR)
        SL(i,j,3)=slope( vyL, vyM, vyR)
        SL(i,j,4)=slope( vzL, vzM, vzR)
        SL(i,j,5)=slope( bxL, bxM, bxR)
        SL(i,j,6)=slope( byL, byM, byR)
        SL(i,j,7)=slope( bzL, bzM, bzR)
        SL(i,j,8)=0.d0

      enddo !! end do i
    enddo !! end do j    
  endif  !! if dd .eq. 2, phi direction

endif !! coordType .eq. 2 or 3, polar coordinates

   SL(:,:,5)=0.d0
if(coordType .eq. 1) then
  do j=0,ny+1
    do i=0,nx+1
       ql(1)=g1*q(i,j,f1)
       ql(2)=g2*q(i,j,f2)/ql(1)
       ql(3)=g3*q(i,j,f3)/ql(1)
       ql(4)=g4*q(i,j,f4)/ql(1)
       ql(5)=0.5d0*(g5*q(i,j,f5)+ g9*q(i,j, f9))
       ql(6)=0.5d0*(g6*q(i,j,f6)+g10*q(i,j,f10))
       ql(7)=0.5d0*(g7*q(i,j,f7)+g11*q(i,j,f11))
       ql(8)=0.d0

       qr(1)=g1*q(i+cx,j+cy,f1)
       qr(2)=g2*q(i+cx,j+cy,f2)/qr(1)
       qr(3)=g3*q(i+cx,j+cy,f3)/qr(1)
       qr(4)=g4*q(i+cx,j+cy,f4)/qr(1)
       qr(5)=0.5d0*(g5*q(i+cx,j+cy,f5)+ g9*q(i+cx,j+cy, f9))
       qr(6)=0.5d0*(g6*q(i+cx,j+cy,f6)+g10*q(i+cx,j+cy,f10))
       qr(7)=0.5d0*(g7*q(i+cx,j+cy,f7)+g11*q(i+cx,j+cy,f11))
       qr(8)=0.d0

       ql(5)=g5*q(i,j,f9)
       qr(5)=g5*q(i+cx,j+cy,f5)

       slopeL(:)=SL(i   ,j   ,:)
       slopeR(:)=SL(i+cx,j+cy,:)

       call fluxPtr(ql,qr,slopeL,slopeR,this%snd,flux,this%nvar)
       Fx(i,j,:)=flux(:)  
 
    enddo !! end do i
  end do !!end do j   
endif !! coordType .eq. 1

if(coordType .eq. 2 .or. coordType .eq. 3) then
  if(dd .eq. 1) then
    do j=0,ny+1
      do i=0,nx+1
         ql(1)=g1*q(i,j,f1)
         ql(2)=g2*q(i,j,f2)/ql(1)
         ql(3)=g3*q(i,j,f3)/ql(1)
         ql(4)=g4*q(i,j,f4)/ql(1)
         ql(5)=0.5d0*(g5*q(i,j,f5)+ g9*q(i,j, f9))
         ql(6)=0.5d0*(g6*q(i,j,f6)+g10*q(i,j,f10))
         ql(7)=0.5d0*(g7*q(i,j,f7)+g11*q(i,j,f11))
         ql(8)=0.d0

         qr(1)=g1*q(i+cx,j+cy,f1)
         qr(2)=g2*q(i+cx,j+cy,f2)/qr(1)
         qr(3)=g3*q(i+cx,j+cy,f3)/qr(1)
         qr(4)=g4*q(i+cx,j+cy,f4)/qr(1)
         qr(5)=0.5d0*(g5*q(i+cx,j+cy,f5)+ g9*q(i+cx,j+cy, f9))
         qr(6)=0.5d0*(g6*q(i+cx,j+cy,f6)+g10*q(i+cx,j+cy,f10))
         qr(7)=0.5d0*(g7*q(i+cx,j+cy,f7)+g11*q(i+cx,j+cy,f11))
         qr(8)=0.d0

         ql(5)=g5*q(i,j,f9)
         qr(5)=g5*q(i+cx,j+cy,f5)

         slopeL(:)=SL(i   ,j   ,:)
         slopeR(:)=SL(i+cx,j+cy,:)

         call fluxPtr(ql,qr,slopeL,slopeR,this%snd,flux,this%nvar)
         Fx(i,j,:)=flux(:)

      enddo !! end do i
    end do !!end do j    
  endif !! dd .eq. 1
  
  if(dd .eq. 2) then
    do j=0,ny+1
      do i=0,nx+1
         ql(1)=g1*q(i,j,f1)
         ql(2)=g2*q(i,j,f2)/ql(1)
         ql(3)=g3*q(i,j,f3)/ql(1)
         ql(4)=g4*q(i,j,f4)/ql(1)
         ql(5)=0.5d0*(g5*q(i,j,f5)+ g9*q(i,j, f9))

         ri=this%xl(1)%coords(i)
         rm=this%xc(1)%coords(i)
         ro=this%xr(1)%coords(i)

         ql(6)=0.5d0*(ri*g6*q(i,j,f6)+ro*g10*q(i,j,f10))/rm
         ql(7)=0.5d0*(g7*q(i,j,f7)+g11*q(i,j,f11))
         ql(8)=0.d0

         qr(1)=g1*q(i+cx,j+cy,f1)
         qr(2)=g2*q(i+cx,j+cy,f2)/qr(1)
         qr(3)=g3*q(i+cx,j+cy,f3)/qr(1)
         qr(4)=g4*q(i+cx,j+cy,f4)/qr(1)
         qr(5)=0.5d0*(g5*q(i+cx,j+cy,f5)+ g9*q(i+cx,j+cy, f9))
         qr(6)=0.5d0*(ri*g6*q(i+cx,j+cy,f6)+ro*g10*q(i+cx,j+cy,f10))/rm
         qr(7)=0.5d0*(g7*q(i+cx,j+cy,f7)+g11*q(i+cx,j+cy,f11))
         qr(8)=0.d0

         ql(5)=g5*q(i,j,f9)
         qr(5)=g5*q(i+cx,j+cy,f5)

         slopeL(:)=SL(i   ,j   ,:)
         slopeR(:)=SL(i+cx,j+cy,:)

         call fluxPtr(ql,qr,slopeL,slopeR,this%snd,flux,this%nvar)
         Fx(i,j,:)=flux(:)

      enddo !! end do i
    end do !!end do j
  endif !! dd .eq. 2
endif  !! if coordType .eq. 2 or 3


Ez=0.d0
Ezc=0.d0

!!!!! calculate EMF at the cell centers
if(coordType .eq. 1) then
  do j=1-nbuf,ny+nbuf
    do i=1-nbuf,nx+nbuf
      rho=q(i,j,1)
       vx=q(i,j,2)/rho
       vy=q(i,j,3)/rho
       bx=0.5d0*(q(i,j,5)+q(i,j,9 ))
       by=0.5d0*(q(i,j,6)+q(i,j,10))
       Ezc(i,j)=vy*bx-by*vx
    enddo
  enddo
endif

if(coordType .eq. 2 .or. coordType .eq. 3) then
  do j=1-nbuf,ny+nbuf
    do i=1-nbuf,nx+nbuf
      rho=q(i,j,1)
       vx=q(i,j,2)/rho
       vy=q(i,j,3)/rho
       ri=this%xl(1)%coords(i)
       rm=this%xc(1)%coords(i)
       ro=this%xr(1)%coords(i)
       bx=0.5d0*(ri*q(i,j,5)+ro*q(i,j,9 ))/rm
       by=0.5d0*(q(i,j,6)+q(i,j,10))
       Ezc(i,j)=vy*bx-by*vx
    enddo
  enddo
endif

!!!!!!!! Gardiner & Stone, JCP, 2005, 205, 509
!!!!!!!! Skinner & Ostriker, ApJS, 2010, 188, 290
if(dd .eq. 1) then
  if(coordType .eq. 1) then
    do j=0,ny
      do i=0,nx
       Ez(i,j)= 0.5d0*(-Fx(i,j,f6)-Fx(i,j+1,f6))*coef(6)
       Ez(i,j)=Ez(i,j)-0.25d0*(Ezc(i,j)+Ezc(i,j+1))
      enddo
    enddo
  endif !! coordType .eq. 1
  if(coordType .eq. 3) then
    do j=0,ny
      do i=0,nx
        rm=this%xc(1)%coords(i)
        ro=this%xr(1)%coords(i)
        Ez(i,j)=0.25d0*(-2.d0*Fx(i,j,f6)-2.d0*Fx(i,j+1,f6))*coef(6)
        Ez(i,j)=Ez(i,j)-0.125d0*(1.d0+ro/rm)*(Ezc(i,j)+Ezc(i,j+1))
      enddo
    enddo
  endif !! 
  if(coordType .eq. 2) then
    print *,"riemannSolverModule.f03: not yet implemented..."
  endif
elseif (dd .eq. 2) then
  if(coordType .eq. 1) then
    do j=0,ny
      do i=0,nx
         Ez(i,j)=0.5d0*(Fx(i+1,j,f5)+Fx(i,j,f5))*coef(5)
         Ez(i,j)=Ez(i,j)-0.25d0*(Ezc(i+1,j)+Ezc(i+1,j+1))
      enddo
    enddo
  endif !! coordType .eq. 1
  if(coordType .eq. 3) then
    do j=0,ny
      do i=0,nx
        ri=this%xc(1)%coords(i)
        ro=this%xc(1)%coords(i+1)
        rm=this%xl(1)%coords(i+1)
        Ez(i,j)=0.25d0*((1.d0+rm/ri)*Fx(i,j,f5)+(1.d0+rm/ro)*Fx(i+1,j,f5))*coef(5)
        Ez(i,j)=Ez(i,j)-0.125d0*(1.d0+rm/ro)*(Ezc(i+1,j)+Ezc(i+1,j+1))
      enddo
    enddo
  endif 
  if(coordType .eq. 2) then
    print *,"riemannSolverModule.f03: not yet implemented...."
  endif
endif

if(coordType .eq. 1) then
  do j=1,ny
    do i=1,nx
      dtddx=this%dt/(this%dx(1)%coords(i)*cx+this%dx(2)%coords(j)*cy)
      q2(i,j,1)=q1(i,j,1)+coef(1)*(Fx(i-cx,j-cy,f1)-Fx(i,j,f1))*dtddx
      q2(i,j,2)=q1(i,j,2)+coef(2)*(Fx(i-cx,j-cy,f2)-Fx(i,j,f2))*dtddx
      q2(i,j,3)=q1(i,j,3)+coef(3)*(Fx(i-cx,j-cy,f3)-Fx(i,j,f3))*dtddx
      q2(i,j,4)=q1(i,j,4)+coef(4)*(Fx(i-cx,j-cy,f4)-Fx(i,j,f4))*dtddx
      q2(i,j,5)=q1(i,j,5)+(Ez(i-1,j-1)-Ez(i-1,j  ))*dtddx
      q2(i,j,6)=q1(i,j,6)+(Ez(i  ,j-1)-Ez(i-1,j-1))*dtddx
      q2(i,j,7)=q1(i,j,7)+coef(7)*(Fx(i-cx,j-cy,f7)-Fx(i,j,f7))*dtddx
      q2(i,j,8)=0.d0
      q2(i,j,9)=q1(i,j,  9)+(Ez(i,j-1)-Ez(i  ,j))*dtddx
      q2(i,j,10)=q1(i,j,10)+(Ez(i,j  )-Ez(i-1,j))*dtddx
      q2(i,j,11)=q2(i,j,7)
    enddo
  enddo
endif
if(coordType .eq. 3) then
  if(dd .eq. 1) then
    do j=1,ny
      do i=1,nx
        rm=this%xc(1)%coords(i)
        dr=this%dx(1)%coords(i)
       dphi=this%dx(2)%coords(j)
        fac1=rm-0.5d0*dr
        fac2=rm+0.5d0*dr   
        dtddx=this%dt/(rm*dr)
        rho=q(i,j,1)
         vr=q(i,j,2)/rho
       vphi=q(i,j,3)/rho
         vz=q(i,j,4)/rho
        brc=0.5d0*(fac1*q(i,j,5)+fac2*q(i,j,9))/rm
       bphic=0.5d0*(q(i,j,6)+q(i,j,10))
         bzc=0.5d0*(q(i,j,7)+q(i,j,11))
        pressure=this%snd**2*rho
        ptot=pressure+0.5d0*(brc**2+bphic**2+bzc**2)

        q2(i,j,1)=q1(i,j,1)+coef(1)*(Fx(i-cx,j-cy,f1)*fac1-Fx(i,j,f1)*fac2)*dtddx
        q2(i,j,2)=q1(i,j,2)+coef(2)*(Fx(i-cx,j-cy,f2)*fac1-Fx(i,j,f2)*fac2)*dtddx &
                                   +(rho*vphi**2+ptot-bphic**2)/rm*this%dt
        q2(i,j,3)=q1(i,j,3)+coef(3)*(Fx(i-cx,j-cy,f3)*fac1**2-Fx(i,j,f3)*fac2**2)*this%dt/(dr*rm**2)
        q2(i,j,4)=q1(i,j,4)+coef(4)*(Fx(i-cx,j-cy,f4)*fac1-Fx(i,j,f4)*fac2)*dtddx
        q2(i,j,5)=q1(i,j,5)+(Ez(i-1,j-1)-Ez(i-1,j  ))*this%dt/(fac1*dphi)
        q2(i,j,6)=q1(i,j,6)+(Ez(i  ,j-1)-Ez(i-1,j-1))*this%dt/dr
        q2(i,j,7)=q1(i,j,7)+coef(7)*(Fx(i-cx,j-cy,f7)*fac1-Fx(i,j,f7)*fac2)*dtddx
        q2(i,j,8)=0.d0
        q2(i,j,9)=q1(i,j,  9)+(Ez(i,j-1)-Ez(i  ,j))*this%dt/(fac2*dphi)
        q2(i,j,10)=q1(i,j,10)+(Ez(i,j  )-Ez(i-1,j))*this%dt/dr
        q2(i,j,11)=q2(i,j,7)
      enddo
    enddo  
  endif !! dd .eq. 1

  if(dd .eq. 2) then
    do j=1,ny
      do i=1,nx
        rm=this%xc(1)%coords(i)
        dr=this%dx(1)%coords(i)
       dphi=this%dx(2)%coords(j)
        fac1=rm-0.5d0*dr
        fac2=rm+0.5d0*dr
        dtddx=this%dt/(rm*dphi)

        q2(i,j,1)=q1(i,j,1)+coef(1)*(Fx(i-cx,j-cy,f1)-Fx(i,j,f1))*dtddx
        q2(i,j,2)=q1(i,j,2)+coef(2)*(Fx(i-cx,j-cy,f2)-Fx(i,j,f2))*dtddx
        q2(i,j,3)=q1(i,j,3)+coef(3)*(Fx(i-cx,j-cy,f3)-Fx(i,j,f3))*dtddx
        q2(i,j,4)=q1(i,j,4)+coef(4)*(Fx(i-cx,j-cy,f4)-Fx(i,j,f4))*dtddx
        q2(i,j,5)=q1(i,j,5)+(Ez(i-1,j-1)-Ez(i-1,j  ))*this%dt/(fac1*dphi)
        q2(i,j,6)=q1(i,j,6)+(Ez(i  ,j-1)-Ez(i-1,j-1))*this%dt/dr
        q2(i,j,7)=q1(i,j,7)+coef(7)*(Fx(i-cx,j-cy,f7)-Fx(i,j,f7))*dtddx
        q2(i,j,8)=0.d0
        q2(i,j,9)=q1(i,j,  9)+(Ez(i,j-1)-Ez(i  ,j))*this%dt/(fac2*dphi)
        q2(i,j,10)=q1(i,j,10)+(Ez(i,j  )-Ez(i-1,j))*this%dt/dr
        q2(i,j,11)=q2(i,j,7)
      enddo
    enddo    
  endif

endif
if(coordType .eq. 2) then
  print *,"riemannSolverModule.f03: not yet implemented....."
  stop
endif

deallocate(SL,Fx,Ez,Ezc)
end subroutine solverIsoMHD2D

subroutine solverIsoMHD1D(this,q,q1,q2,dd)
!!!!!!! when using isothermal MHD solver, energy equation is treated as a redundant variable
!!!!!!! it needs to be there, just not evolves
!!!!!!! the total number of variable should be eight

class(grid)::this
double precision,dimension(1-this%nbuf:this%nMesh(1)+this%nbuf,this%nvar)::q,q1,q2
integer::nx,nbuf,nvar,dd
integer::i
double precision::dt,dx,dtddx,gam
double precision::rhoL,rhoM,rhoR
double precision::vxL,vxM,vxR
double precision::vyL,vyM,vyR
double precision::vzL,vzM,vzR
double precision::bxL,bxM,bxR,bxtmp
double precision::byL,byM,byR
double precision::bzL,bzM,bzR
double precision,dimension(:,:),allocatable::SL
double precision,dimension(:,:),allocatable::Fx
double precision::ql(this%nvar),qr(this%nvar),slopeL(this%nvar),slopeR(this%nvar),flux(this%nvar)
procedure(limiter),pointer::slope=>null()
procedure(fluxSolver),pointer::fluxPtr=>null()

nx=this%nMesh(1)
nbuf=this%nbuf
nvar=this%nvar
dt=this%dt
gam=this%adiGamma
select case(this%solverType)
   case (4)
     fluxPtr=>fluxHLLIsoMHD1D
   case (5)
     fluxPtr=>fluxHLLDIsoMHD1D
   case default
     print *,'no appropriate solver found!'
     stop
end select

select case (this%limiterType)
   case(0)
     slope=>zslop
   case(1)
     slope=>vslop
   case(2)
     slope=>fslop
   case(3)
     slope=>minmod
end select

allocate(SL(1-nbuf:nx+nbuf,nvar))
allocate(Fx(1-nbuf:nx+nbuf,nvar))

do i=0,nx+1
   rhoL=q(i-1,1)
   rhoM=q(i  ,1)
   rhoR=q(i+1,1)

    vxL=q(i-1,2)/rhoL
    vxM=q(i  ,2)/rhoM
    vxR=q(i+1,2)/rhoR

    vyL=q(i-1,3)/rhoL
    vyM=q(i  ,3)/rhoM
    vyR=q(i+1,3)/rhoR

    vzL=q(i-1,4)/rhoL
    vzM=q(i  ,4)/rhoM
    vzR=q(i+1,4)/rhoR

    bxL=0.5d0*(q(i-1,5)+q(i-1, 9))
    bxM=0.5d0*(q(i  ,5)+q(i  , 9))    
    bxR=0.5d0*(q(i+1,5)+q(i+1, 9))
 
    byL=0.5d0*(q(i-1,6)+q(i-1,10))
    byM=0.5d0*(q(i  ,6)+q(i  ,10))
    byR=0.5d0*(q(i+1,6)+q(i+1,10))

    bzL=0.5d0*(q(i-1,7)+q(i-1,11))
    bzM=0.5d0*(q(i  ,7)+q(i  ,11))
    bzR=0.5d0*(q(i+1,7)+q(i+1,11))

   SL(i,1)=slope(rhoL,rhoM,rhoR)
   SL(i,2)=slope( vxL, vxM, vxR)
   SL(i,3)=slope( vyL, vyM, vyR)
   SL(i,4)=slope( vzL, vzM, vzR)
   SL(i,5)=0.d0
   SL(i,6)=slope( byL, byM, byR)
   SL(i,7)=slope( bzL, bzM, bzR)
   SL(i,8)=0.d0
enddo ! end i

do i=0,nx+1
   ql(1)=q(i,1)
   ql(2)=q(i,2)/ql(1)
   ql(3)=q(i,3)/ql(1)
   ql(4)=q(i,4)/ql(1)
   ql(5)=q(i,9)
   ql(6)=0.5d0*(q(i,6)+q(i,10)) 
   ql(7)=0.5d0*(q(i,7)+q(i,11))
   ql(8)=0.d0

   qr(1)=q(i+1,1)
   qr(2)=q(i+1,2)/qr(1)
   qr(3)=q(i+1,3)/qr(1)
   qr(4)=q(i+1,4)/qr(1)
   qr(5)=q(i+1,5)
   qr(6)=0.5d0*(q(i+1,6)+q(i+1,10))
   qr(7)=0.5d0*(q(i+1,7)+q(i+1,11))
   qr(8)=0.d0

   slopeL(:)=SL(i,:)
   slopeR(:)=SL(i+1,:)
   call fluxPtr(ql,qr,slopeL,slopeR,this%snd,flux,this%nvar)
   Fx(i,:)=flux(:)
enddo

do i=1,nx
   dtddx=dt/this%dx(1)%coords(i)
   q2(i,1)=q1(i,1)+(Fx(i-1,1)-Fx(i,1))*dtddx
   q2(i,2)=q1(i,2)+(Fx(i-1,2)-Fx(i,2))*dtddx
   q2(i,3)=q1(i,3)+(Fx(i-1,3)-Fx(i,3))*dtddx
   q2(i,4)=q1(i,4)+(Fx(i-1,4)-Fx(i,4))*dtddx
   q2(i,5)=q1(i,5)
   q2(i,6)=q1(i,6)+(Fx(i-1,6)-Fx(i,6))*dtddx
   q2(i,7)=q1(i,7)+(Fx(i-1,7)-Fx(i,7))*dtddx
   q2(i,8)=0.d0
   q2(i,9)=q1(i,9)
   q2(i,10)=q1(i,10)+(Fx(i-1,6)-Fx(i,6))*dtddx
   q2(i,11)=q1(i,11)+(Fx(i-1,7)-Fx(i,7))*dtddx
enddo


deallocate(SL,Fx)

end subroutine solverIsoMHD1D

subroutine fluxHLLIsoMHD1D(ql,qr,slopeL,slopeR,snd,flux,nvar)
!!! ref. A. Mignone, JCoPh, 2007

implicit none
integer::nvar
double precision::ql(NVAR),qr(NVAR),slopeL(NVAR),slopeR(NVAR),flux(NVAR)
double precision::rhoL,vxL,vyL,vzL,bxL,byL,bzL,fastL
double precision::rhoR,vxR,vyR,vzR,bxR,byR,bzR,fastR
double precision::SL,SR,SLs,SRs,snd
double precision::cfL,cfR,BsqL,BsqR,rhohll,mxhll,Frhohll,Fmxhll,pTL,pTR,bx
double precision::UL(nvar),UR(nvar)
double precision::FL(nvar),FR(nvar),FLs(nvar),FRs(nvar),FC(nvar)
double precision::rhos,vxs,vysL,vysR,vzsL,vzsR,bysL,bysR,bzsL,bzsR,vyc,vzc,byc,bzc,X

UL=ql+0.5d0*slopeL
UR=qr-0.5d0*slopeR

rhoL=UL(1)
vxL =UL(2)
vyL =UL(3)
vzL =UL(4)
bxL =0.5d0*(UL(5)+UR(5))
byL =UL(6)
bzL =UL(7)
BsqL=(bxL**2.d0+byL**2.d0+bzL**2.d0)/rhoL
cfL=dsqrt(0.5d0*(snd**2+BsqL+dsqrt((snd**2+BsqL)**2-4.d0*snd**2*bxL**2/rhoL)))
pTL=rhoL*snd**2+0.5d0*BsqL*rhoL

rhoR=UR(1)
vxR =UR(2)
vyR =UR(3)
vzR =UR(4)
bxR =0.5d0*(UL(5)+UR(5))
byR =UR(6)
bzR =UR(7)
BsqR=(bxR**2.d0+byR**2.d0+bzR**2.d0)/rhoR
cfR=dsqrt(0.5d0*(snd**2+BsqR+dsqrt((snd**2+BsqR)**2-4.d0*snd**2*BxR**2/rhoR)))
pTR=rhoR*snd**2+0.5d0*BsqR*rhoR

SL=dmin1(vxL,vxR)-dmax1(cfL,cfR)
SR=dmax1(vxL,vxR)+dmax1(cfL,cfR)

!SL=dmin1(vxL-cfL,vxR-cfR)
!SR=dmax1(vxL+cfR,vxR+cfR)

FL(1)=rhoL*vxL
FL(2)=rhoL*vxL**2.d0+pTL-bxL**2.d0
FL(3)=rhoL*vyL*vxL-bxL*byL
FL(4)=rhoL*vzL*vxL-bxL*bzL
FL(5)=0.d0
FL(6)=byL*vxL-bxL*vyL
FL(7)=bzL*vxL-bxL*vzL
FL(8)=0.d0

FR(1)=rhoR*vxR
FR(2)=rhoR*vxR**2.d0+pTR-bxR**2.d0
FR(3)=rhoR*vyR*vxR-bxR*byR
FR(4)=rhoR*vzR*vxR-bxR*bzR
FR(5)=0.d0
FR(6)=byR*vxR-bxR*vyR
FR(7)=bzR*vxR-bxR*vzR
FR(8)=0.d0

FC(1)=(SR*FL(1)-SL*FR(1)+SL*SR*(rhoR-rhoL))/(SR-SL)
FC(2)=(SR*FL(2)-SL*FR(2)+SL*SR*(rhoR*vxR-rhoL*vxL))/(SR-SL)
FC(3)=(SR*FL(3)-SL*FR(3)+SL*SR*(rhoR*vyR-rhoL*vyL))/(SR-SL)
FC(4)=(SR*FL(4)-SL*FR(4)+SL*SR*(rhoR*vzR-rhoL*vzL))/(SR-SL)
FC(5)=0.d0
FC(6)=(SR*FL(6)-SL*FR(6)+SL*SR*(byR-byL))/(SR-SL)
FC(7)=(SR*FL(7)-SL*FR(7)+SL*SR*(bzR-bzL))/(SR-SL)
FC(8)=0.d0

if(SL .gt. 0.d0) then
  flux=FL
elseif (SL .le. 0.d0 .and. SR .ge. 0.d0) then
  flux=FC
elseif (SR .lt. 0.d0) then
  flux=FR
endif

end subroutine fluxHLLIsoMHD1D

subroutine fluxHLLDIsoMHD1D(ql,qr,slopeL,slopeR,snd,flux,nvar)
!!! ref. A. Mignone, JCoPh, 2007

implicit none
integer::nvar
double precision::ql(NVAR),qr(NVAR),slopeL(NVAR),slopeR(NVAR),flux(NVAR)
double precision::rhoL,vxL,vyL,vzL,bxL,byL,bzL,fastL
double precision::rhoR,vxR,vyR,vzR,bxR,byR,bzR,fastR
double precision::SL,SR,SLs,SRs,snd
double precision::cfL,cfR,BsqL,BsqR,rhohll,mxhll,Frhohll,Fmxhll,pTL,pTR,bx
double precision::UL(nvar),UR(nvar)
double precision::FL(nvar),FR(nvar),FLs(nvar),FRs(nvar),FC(nvar)
double precision::rhos,vxs,vysL,vysR,vzsL,vzsR,bysL,bysR,bzsL,bzsR,vyc,vzc,byc,bzc,X

UL=ql+0.5d0*slopeL
UR=qr-0.5d0*slopeR

rhoL=UL(1)
vxL =UL(2)
vyL =UL(3)
vzL =UL(4)
bxL =0.5d0*(UL(5)+UR(5))
byL =UL(6)
bzL =UL(7)
BsqL=(bxL**2.d0+byL**2.d0+bzL**2.d0)/rhoL
cfL=dsqrt(0.5d0*(snd**2+BsqL+dsqrt((snd**2+BsqL)**2-4.d0*snd**2*bxL**2/rhoL)))
pTL=rhoL*snd**2+0.5d0*BsqL*rhoL

rhoR=UR(1)
vxR =UR(2)
vyR =UR(3)
vzR =UR(4)
bxR =0.5d0*(UL(5)+UR(5))
byR =UR(6)
bzR =UR(7)
BsqR=(bxR**2.d0+byR**2.d0+bzR**2.d0)/rhoR
cfR=dsqrt(0.5d0*(snd**2+BsqR+dsqrt((snd**2+BsqR)**2-4.d0*snd**2*BxR**2/rhoR)))
pTR=rhoR*snd**2+0.5d0*BsqR*rhoR

SL=dmin1(vxL,vxR)-dmax1(cfL,cfR)
SR=dmax1(vxL,vxR)+dmax1(cfL,cfR)

!SL=dmin1(vxL-cfL,vxR-cfR)
!SR=dmax1(vxL+cfR,vxR+cfR)
bx=0.5d0*(bxL+bxR)

rhohll=(SR*rhoR-SL*rhoL-rhoR*vxR+rhoL*vxL)/(SR-SL)
mxhll=(SR*rhoR*vxR-SL*rhoL*vxL-(rhoR*vxR**2+pTR-bxR**2)+(rhoL*vxL**2+pTL-bxL**2))/(SR-SL)
Frhohll=(SR*rhoL*vxL-SL*rhoR*vxR+SR*SL*(rhoR-rhoL))/(SR-SL)
vxs=Frhohll/rhohll
Fmxhll=(SR*(rhoL*vxL**2+pTL-BxL**2)-SL*(rhoR*vxR**2+pTR-BxR**2)+SR*SL*(rhoR*vxR-rhoL*vxL))/(SR-SL)

SLs=vxs-dabs(bx)/dsqrt(rhohll)
SRs=vxs+dabs(bx)/dsqrt(rhohll)

rhos=rhohll

vysL=vyL-bx*byL/rhos*(vxs-vxL)/((SL-SLs)*(SL-SRs))
vzsL=vzL-bx*bzL/rhos*(vxs-vxL)/((SL-SLs)*(SL-SRs))
vysR=vyR-bx*byR/rhos*(vxs-vxR)/((SR-SLs)*(SR-SRs))
vzsR=vzR-bx*bzR/rhos*(vxs-vxR)/((SR-SLs)*(SR-SRs))

bysL=byL/rhos*(rhoL*(SL-vxL)**2-bx**2)/((SL-SLs)*(SL-SRs))
bzsL=bzL/rhos*(rhoL*(SL-vxL)**2-bx**2)/((SL-SLs)*(SL-SRs))
bysR=byR/rhos*(rhoR*(SR-vxR)**2-bx**2)/((SR-SLs)*(SR-SRs))
bzsR=bzR/rhos*(rhoR*(SR-vxR)**2-bx**2)/((SR-SLs)*(SR-SRs))

if((SL-SLs)*(SL-SRs) .eq. 0) then
  vysL=vyL
  vzsL=vzL
  bysL=byL
  bzsL=bzL
endif

if((SR-SLs)*(SR-SRs) .eq. 0) then
  vysR=vyR
  vzsR=vzR
  bysR=byR
  bzsR=bzR
endif

FL(1)=rhoL*vxL
FL(2)=rhoL*vxL**2.d0+pTL-bxL**2.d0
FL(3)=rhoL*vyL*vxL-bxL*byL
FL(4)=rhoL*vzL*vxL-bxL*bzL
FL(5)=0.d0
FL(6)=byL*vxL-bxL*vyL
FL(7)=bzL*vxL-bxL*vzL
FL(8)=0.d0

FR(1)=rhoR*vxR
FR(2)=rhoR*vxR**2.d0+pTR-bxR**2.d0
FR(3)=rhoR*vyR*vxR-bxR*byR
FR(4)=rhoR*vzR*vxR-bxR*bzR
FR(5)=0.d0
FR(6)=byR*vxR-bxR*vyR
FR(7)=bzR*vxR-bxR*vzR
FR(8)=0.d0

! ============================
FLs(1)=rhos*vxs!FL(1)+SL*(rhos-rhoL)
FLs(2)=Fmxhll!FL(2)+SL*(mxhll-rhoL*vxL)
FLs(3)=FL(3)+SL*(rhos*vysL-rhoL*vyL)
FLs(4)=FL(4)+SL*(rhos*vzsL-rhoL*vzL)
FLs(5)=FL(5)+SL*(bx-bxL)
FLs(6)=FL(6)+SL*(bysL-byL)
FLs(7)=FL(7)+SL*(bzsL-bzL)
FLs(8)=0.d0

FRs(1)=rhos*vxs!FR(1)+SR*(rhos-rhoR)
FRs(2)=Fmxhll!FR(2)+SR*(mxhll-rhoR*vxR)
FRs(3)=FR(3)+SR*(rhos*vysR-rhoR*vyR)
FRs(4)=FR(4)+SR*(rhos*vzsR-rhoR*vzR)
FRs(5)=FR(5)+SR*(bx-bxR)
FRs(6)=FR(6)+SR*(bysR-byR)
FRs(7)=FR(7)+SR*(bzsR-bzR)
FRs(8)=0.d0
! ============================
X=dsqrt(rhos)*dsign(1.d0,bx)
vyc=0.5d0*(vysL+vysR)+0.5d0/rhos*X*(bysR-bysL)
vzc=0.5d0*(vzsL+vzsR)+0.5d0/rhos*X*(bzsR-bzsL)
byc=0.5d0*(bysL+bysR)+0.5d0*rhos/X*(vysR-vysL)
bzc=0.5d0*(bzsL+bzsR)+0.5d0*rhos/X*(vzsR-vzsL)

FC(1)=rhos*vxs
FC(2)=Fmxhll
FC(3)=rhos*vyc*vxs-bx*byc
FC(4)=rhos*vzc*vxs-bx*bzc
FC(5)=0.d0
FC(6)=byc*vxs-bx*vyc
FC(7)=bzc*vxs-bx*vzc
FC(8)=0.d0

if(SL .gt. 0.d0) then
  flux=FL
elseif (SL .le. 0.d0 .and. SLs .ge. 0.d0) then
  flux=FLs
elseif (SLs .le. 0.d0 .and. SRs .ge. 0.d0) then
  flux=FC
elseif (SRs .le. 0.d0 .and. SR .ge. 0.d0) then
  flux=FRs
elseif (SR .lt. 0.d0) then
  flux=FR
endif

end subroutine fluxHLLDIsoMHD1D
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!! Riemann solvers for 3D problems
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine solverAdiMHD3D(this,q,q1,q2,dd)
class(grid)::this
double precision,dimension(1-this%nbuf:this%nMesh(1)+this%nbuf,1-this%nbuf:this%nMesh(2)+this%nbuf,&
                           1-this%nbuf:this%nMesh(3)+this%nbuf,this%nvar)::q,q1,q2
integer::nx,ny,nz,nbuf,nvar,dd
integer::i,j,k
integer::cx,cy,cz
integer::f1,f2,f3,f4,f5,f6,f7,f8,f9,f10,f11
double precision::g1,g2,g3,g4,g5,g6,g7,g8,g9,g10,g11
double precision::dt,dx,dtddx,gam
double precision::rhoL,rhoM,rhoR
double precision::vxL,vxM,vxR
double precision::vyL,vyM,vyR
double precision::vzL,vzM,vzR
double precision::bxL,bxM,bxR,bxtmp
double precision::byL,byM,byR
double precision::bzL,bzM,bzR
double precision::eneL,eneM,eneR
double precision::pL ,pM ,pR
double precision::pressure,bxc,byc,bzc,vxc,vyc,vzc,denc,enec
double precision,dimension(:,:,:,:),allocatable::SL
double precision,dimension(:,:,:,:),allocatable::Fx
double precision,dimension(:,:,:),allocatable::Ez,Ezc,Ex,Exc,Ey,Eyc
double precision::ql(this%nvar),qr(this%nvar),slopeL(this%nvar),slopeR(this%nvar),flux(this%nvar),coef(this%nvar)
procedure(limiter),pointer::slope=>null()
procedure(fluxSolver),pointer::fluxPtr=>null()
double precision::rho,vx,vy,vz,bx,by,bz
integer::chgFluxPtrCount

nx=this%nMesh(1)
ny=this%nMesh(2)
nz=this%nMesh(3)

nbuf=this%nbuf
nvar=this%nvar
dt=this%dt
gam=this%adiGamma
chgFluxPtrCount=0

select case(this%solverType)
   case (4)
     fluxPtr=>fluxHLLAdiMHD1D
   case (5)
     fluxPtr=>fluxHLLDAdiMHD1D
   case default
     print *,'riemannSolverModule.f03: no appropriate 2D MHD solver found!'
end select

select case (this%limiterType)
   case(0)
     slope=>zslop
   case(1)
     slope=>vslop
   case(2)
     slope=>fslop
   case(3)
     slope=>minmod
   case default
     print *, "riemannSolverModule.f03: no appropriate slope limiter found!"
end select

allocate(SL(1-nbuf:nx+nbuf,1-nbuf:ny+nbuf,1-nbuf:nz+nbuf,nvar))
allocate(Fx(1-nbuf:nx+nbuf,1-nbuf:ny+nbuf,1-nbuf:nz+nbuf,nvar))
allocate( Ez(1-nbuf:nx+nbuf,1-nbuf:ny+nbuf,1-nbuf:nz+nbuf))
allocate(Ezc(1-nbuf:nx+nbuf,1-nbuf:ny+nbuf,1-nbuf:nz+nbuf))
allocate( Ex(1-nbuf:nx+nbuf,1-nbuf:ny+nbuf,1-nbuf:nz+nbuf))
allocate(Exc(1-nbuf:nx+nbuf,1-nbuf:ny+nbuf,1-nbuf:nz+nbuf))
allocate( Ey(1-nbuf:nx+nbuf,1-nbuf:ny+nbuf,1-nbuf:nz+nbuf))
allocate(Eyc(1-nbuf:nx+nbuf,1-nbuf:ny+nbuf,1-nbuf:nz+nbuf))

!!!!! dd=1 ==>x, dd=2 ==>y, dd=3 ==>z !!!!!

if(dd .eq. 1) then
     cx=1
     cy=0
     cz=0
     f1=1
     f2=2
     f3=3
     f4=4
     f5=5
     f6=6
     f7=7
     f8=8
     f9=9
    f10=10
    f11=11
 coef(1)=1.d0
 coef(2)=1.d0
 coef(3)=1.d0
 coef(4)=1.d0
 coef(5)=1.d0
 coef(6)=1.d0
 coef(7)=1.d0
 coef(8)=1.d0
     g1=1.d0
     g2=1.d0
     g3=1.d0
     g4=1.d0
     g5=1.d0
     g6=1.d0
     g7=1.d0
     g8=1.d0
     g9=1.d0
    g10=1.d0
    g11=1.d0
elseif (dd .eq. 2) then
     cx=0
     cy=1
     cz=0
     f1=1
     f2=3
     f3=2
     f4=4
     f5=6
     f6=5
     f7=7
     f8=8
     f9=10
    f10=9
    f11=11
 coef(1)=1.d0
 coef(2)=-1.d0
 coef(3)=1.d0
 coef(4)=1.d0
 coef(5)=-1.d0
 coef(6)=1.d0
 coef(7)=1.d0
 coef(8)=1.d0
     g1=1.d0
     g2=1.d0
     g3=-1.d0
     g4=1.d0
     g5=1.d0
     g6=-1.d0
     g7=1.d0
     g8=1.d0
     g9=1.d0
    g10=-1.d0
    g11=1.d0
elseif (dd .eq. 3) then
     cx=0
     cy=0
     cz=1
     f1=1
     f2=4
     f3=3
     f4=2
     f5=7
     f6=6
     f7=5
     f8=8
     f9=11
    f10=10
    f11=9
 coef(1)=1.d0
 coef(2)=-1.d0
 coef(3)=1.d0
 coef(4)=1.d0
 coef(5)=-1.d0
 coef(6)=1.d0
 coef(7)=1.d0
 coef(8)=1.d0
     g1=1.d0
     g2=1.d0
     g3=1.d0
     g4=-1.d0
     g5=1.d0
     g6=1.d0
     g7=-1.d0
     g8=1.d0
     g9=1.d0
    g10=1.d0
    g11=-1.d0
endif

do k=0,nz+1
  do j=0,ny+1
    do i=0,nx+1
      rhoL=g1*q(i-cx,j-cy,k-cz,f1)
      rhoM=g1*q(i   ,j   ,k   ,f1)
      rhoR=g1*q(i+cx,j+cy,k+cz,f1)

       vxL=g2*q(i-cx,j-cy,k-cz,f2)/rhoL
       vxM=g2*q(i   ,j   ,k   ,f2)/rhoM
       vxR=g2*q(i+cx,j+cy,k+cz,f2)/rhoR

       vyL=g3*q(i-cx,j-cy,k-cz,f3)/rhoL
       vyM=g3*q(i   ,j   ,k   ,f3)/rhoM
       vyR=g3*q(i+cx,j+cy,k+cz,f3)/rhoR

       vzL=g4*q(i-cx,j-cy,k-cz,f4)/rhoL
       vzM=g4*q(i   ,j   ,k   ,f4)/rhoM
       vzR=g4*q(i+cx,j+cy,k+cz,f4)/rhoR

       bxL=0.5d0*(g5*q(i-cx,j-cy,k-cz,f5)+ g9*q(i-cx,j-cy,k-cz,f9 ))
       bxM=0.5d0*(g5*q(i   ,j   ,k   ,f5)+ g9*q(i   ,j   ,k   ,f9 ))
       bxR=0.5d0*(g5*q(i+cx,j+cy,k+cz,f5)+ g9*q(i+cx,j+cy,k+cz,f9 ))

       byL=0.5d0*(g6*q(i-cx,j-cy,k-cz,f6)+g10*q(i-cx,j-cy,k-cz,f10))
       byM=0.5d0*(g6*q(i   ,j   ,k   ,f6)+g10*q(i   ,j   ,k   ,f10))
       byR=0.5d0*(g6*q(i+cx,j+cy,k+cz,f6)+g10*q(i+cx,j+cy,k+cz,f10))

       bzL=0.5d0*(g7*q(i-cx,j-cy,k-cz,f7)+g11*q(i-cx,j-cy,k-cz,f11))
       bzM=0.5d0*(g7*q(i   ,j   ,k   ,f7)+g11*q(i   ,j   ,k   ,f11))
       bzR=0.5d0*(g7*q(i+cx,j+cy,k+cz,f7)+g11*q(i+cx,j+cy,k+cz,f11))

      eneL=g8*q(i-cx,j-cy,k-cz,f8)/rhoL
      eneM=g8*q(i   ,j   ,k   ,f8)/rhoM
      eneR=g8*q(i+cx,j+cy,k+cz,f8)/rhoR

        pL=(gam-1.d0)*(rhoL*eneL-0.5d0*rhoL*(vxL**2+vyL**2+vzL**2)-0.5d0*(bxL**2+byL**2+bzL**2))
        pM=(gam-1.d0)*(rhoM*eneM-0.5d0*rhoM*(vxM**2+vyM**2+vzM**2)-0.5d0*(bxM**2+byM**2+bzM**2))
        pR=(gam-1.d0)*(rhoR*eneR-0.5d0*rhoR*(vxR**2+vyR**2+vzR**2)-0.5d0*(bxR**2+byR**2+bzR**2))

      SL(i,j,k,1)=slope(rhoL,rhoM,rhoR)
      SL(i,j,k,2)=slope( vxL, vxM, vxR)
      SL(i,j,k,3)=slope( vyL, vyM, vyR)
      SL(i,j,k,4)=slope( vzL, vzM, vzR)
      SL(i,j,k,5)=slope( bxL, bxM, bxR)
      SL(i,j,k,6)=slope( byL, byM, byR)
      SL(i,j,k,7)=slope( bzL, bzM, bzR)
      SL(i,j,k,8)=slope(  pL,  pM,  pR)
    enddo !! end do i
  enddo !! end do j
enddo !! endo k

   SL(:,:,:,5)=0.d0

do k=0,nz+1
  do j=0,ny+1
    do i=0,nx+1
       ql(1)=g1*q(i,j,k,f1)
       ql(2)=g2*q(i,j,k,f2)/ql(1)
       ql(3)=g3*q(i,j,k,f3)/ql(1)
       ql(4)=g4*q(i,j,k,f4)/ql(1)
       ql(5)=0.5d0*(g5*q(i,j,k,f5)+ g9*q(i,j,k, f9))
       ql(6)=0.5d0*(g6*q(i,j,k,f6)+g10*q(i,j,k,f10))
       ql(7)=0.5d0*(g7*q(i,j,k,f7)+g11*q(i,j,k,f11))
       ql(8)=g8*(gam-1.d0)*(q(i,j,k,f8)-0.5d0*ql(1)*(ql(2)**2+ql(3)**2+ql(4)**2)-0.5d0*(ql(5)**2+ql(6)**2+ql(7)**2))

       qr(1)=g1*q(i+cx,j+cy,k+cz,f1)
       qr(2)=g2*q(i+cx,j+cy,k+cz,f2)/qr(1)
       qr(3)=g3*q(i+cx,j+cy,k+cz,f3)/qr(1)
       qr(4)=g4*q(i+cx,j+cy,k+cz,f4)/qr(1)
       qr(5)=0.5d0*(g5*q(i+cx,j+cy,k+cz,f5)+ g9*q(i+cx,j+cy,k+cz, f9))
       qr(6)=0.5d0*(g6*q(i+cx,j+cy,k+cz,f6)+g10*q(i+cx,j+cy,k+cz,f10))
       qr(7)=0.5d0*(g7*q(i+cx,j+cy,k+cz,f7)+g11*q(i+cx,j+cy,k+cz,f11))
       qr(8)=g8*(gam-1.d0)*(q(i+cx,j+cy,k+cz,f8)-0.5d0*qr(1)*(qr(2)**2+qr(3)**2+qr(4)**2)-0.5d0*(qr(5)**2+qr(6)**2+qr(7)**2))

       ql(5)=g5*q(i,j,k,f9)
       qr(5)=g5*q(i+cx,j+cy,k+cz,f5)

       slopeL(:)=SL(i   ,j   ,k   ,:)
       slopeR(:)=SL(i+cx,j+cy,k+cz,:)

       call fluxPtr(ql,qr,slopeL,slopeR,this%adiGamma,flux,this%nvar)
       Fx(i,j,k,:)=flux(:)  
 
    enddo !! end do i
  enddo !!end do j
enddo !! end do k   

Ez=0.d0
Ex=0.d0
Ey=0.d0
Ezc=0.d0
Exc=0.d0
Eyc=0.d0


!!!!! calculate EMF at the cell centers
!!!!! Ezc,Exc,Eyc

do k=1-nbuf,nz+nbuf
  do j=1-nbuf,ny+nbuf
    do i=1-nbuf,nx+nbuf
      rho=q(i,j,k,1)
       vx=q(i,j,k,2)/rho
       vy=q(i,j,k,3)/rho
       vz=q(i,j,k,4)/rho
       bx=0.5d0*(q(i,j,k,5)+q(i,j,k,9 ))
       by=0.5d0*(q(i,j,k,6)+q(i,j,k,10))
       bz=0.5d0*(q(i,j,k,7)+q(i,j,k,11))
       Ezc(i,j,k)=vy*bx-by*vx
       Exc(i,j,k)=vz*by-vy*bz
       Eyc(i,j,k)=vx*bz-vz*bx
    enddo !! enddo i
  enddo  !! enddo j
enddo  !! enddo k


!!!!!!!! Gardiner & Stone, JCP, 2005, 205, 509
if(dd .eq. 1) then

do k=0,nz
  do j=0,ny
    do i=0,nx
       Ez(i,j,k)= 0.5d0*(-Fx(i,j,k,f6)-Fx(i,j+1,k,f6))*coef(6)
       Ez(i,j,k)=Ez(i,j,k)-0.25d0*(Ezc(i,j,k)+Ezc(i,j+1,k))

       Ey(i,j,k)= 0.5d0*( Fx(i,j,k,f7)+Fx(i,j,k+1,f7))*coef(7)
       Ey(i,j,k)=Ey(i,j,k)-0.25d0*(Eyc(i,j,k)+Eyc(i,j,k+1))
    enddo
  enddo
enddo

elseif (dd .eq. 2) then
do k=0,nz
  do j=0,ny
    do i=0,nx
       Ez(i,j,k)=0.5d0*( Fx(i+1,j,k,f5)+Fx(i,j,k,f5))*coef(5)
       Ez(i,j,k)=Ez(i,j,k)-0.25d0*(Ezc(i+1,j,k)+Ezc(i+1,j+1,k))

       Ex(i,j,k)=0.5d0*(-Fx(i,j,k+1,f7)-Fx(i,j,k,f7))*coef(7)
       Ex(i,j,k)=Ex(i,j,k)-0.25d0*(Exc(i,j,k)+Exc(i,j+1,k))
    enddo
  enddo
enddo

elseif (dd .eq. 3) then
do k=0,nz
  do j=0,ny
    do i=0,nx
      Ey(i,j,k)=0.5d0*(-Fx(i+1,j,k,f5)-Fx(i,j,k,f5))*coef(5)
      Ey(i,j,k)=Ey(i,j,k)-0.25d0*(Eyc(i+1,j,k)+Eyc(i+1,j,k+1))

      Ex(i,j,k)=0.5d0*( Fx(i,j+1,k,f6)+Fx(i,j,k,f6))*coef(6)
      Ex(i,j,k)=Ex(i,j,k)-0.25d0*(Exc(i,j,k+1)+Exc(i,j+1,k+1))
    enddo
  enddo
enddo

endif !! if( dd .eq. 1 )


do k=1,nz
  do j=1,ny
    do i=1,nx
      dtddx=this%dt/(this%dx(1)%coords(i)*cx+this%dx(2)%coords(j)*cy+this%dx(3)%coords(k)*cz)
      q2(i,j,k,1)=q1(i,j,k,1)+coef(1)*(Fx(i-cx,j-cy,k-cz,f1)-Fx(i,j,k,f1))*dtddx
      q2(i,j,k,2)=q1(i,j,k,2)+coef(2)*(Fx(i-cx,j-cy,k-cz,f2)-Fx(i,j,k,f2))*dtddx
      q2(i,j,k,3)=q1(i,j,k,3)+coef(3)*(Fx(i-cx,j-cy,k-cz,f3)-Fx(i,j,k,f3))*dtddx
      q2(i,j,k,4)=q1(i,j,k,4)+coef(4)*(Fx(i-cx,j-cy,k-cz,f4)-Fx(i,j,k,f4))*dtddx
      q2(i,j,k,5)=q1(i,j,k,5)+(Ez(i-1,j-1,k)-Ez(i-1,j  ,k))*dtddx+(Ey(i-1,j,k)-Ey(i-1,j,k-1))*dtddx
      q2(i,j,k,6)=q1(i,j,k,6)+(Ez(i  ,j-1,k)-Ez(i-1,j-1,k))*dtddx+(Ex(i,j-1,k-1)-Ex(i,j-1,k))*dtddx
      q2(i,j,k,7)=q1(i,j,k,7)+(Ey(i-1,j,k-1)-Ey(i  ,j,k-1))*dtddx+(Ex(i,j,k-1)-Ex(i,j-1,k-1))*dtddx
      q2(i,j,k,8)=q1(i,j,k,8)+coef(8)*(Fx(i-cx,j-cy,k-cz,f8)-Fx(i,j,k,f8))*dtddx
      q2(i,j,k,9)=q1(i,j,k,9)+(Ez(i,j-1,k)-Ez(i,j,k))*dtddx+(Ey(i,j,k)-Ey(i,j,k-1))*dtddx
      q2(i,j,k,10)=q1(i,j,k,10)+(Ez(i,j,k)-Ez(i-1,j,k))*dtddx+(Ex(i,j,k-1)-Ex(i,j,k))*dtddx
      q2(i,j,k,11)=q1(i,j,k,11)+(Ey(i-1,j,k)-Ey(i,j,k))*dtddx+(Ex(i,j,k)-Ex(i,j-1,k))*dtddx
    enddo
  enddo
enddo

select case(this%solverType)
  case(4)
  case(5)
   if(dd .eq. 3) then
    loopk:do k=1,nz
      do j=1,ny
        do i=1,nx
          denc=q2(i,j,k,1)
           vxc=q2(i,j,k,2)/denc
           vyc=q2(i,j,k,3)/denc
           vzc=q2(i,j,k,4)/denc
           bxc=0.5d0*(q2(i,j,k,5)+q2(i,j,k,9))
           byc=0.5d0*(q2(i,j,k,6)+q2(i,j,k,10))
           bzc=0.5d0*(q2(i,j,k,7)+q2(i,j,k,11))
           enec=q2(i,j,k,8)
          pressure=(gam-1.d0)*(enec-0.5d0*denc*(vxc**2+vyc**2+vzc**2)-0.5d0*(bxc**2+byc**2+bzc**2))

          if(pressure .lt. 0.d0) then
            print *,"riemannSolverModule.f03: negative pressure..."
            print *,"(myid,dd,i,j,k,pressure)=",myid,dd,i,j,k,pressure
            this%changeSolver=.true.
            exit loopk
          endif

        enddo
      enddo
    enddo loopk
     endif
  case default
   print *,"riemannSolverModule.f03: unknown MHD solver type..."
   stop
end select


deallocate(SL,Fx,Ez,Ezc,Ex,Exc,Ey,Eyc)
end subroutine solverAdiMHD3D


subroutine solverIso3D(this,q,q1,q2,dd)
class(grid)::this
double precision,dimension(1-this%nbuf:this%nMesh(1)+this%nbuf,1-this%nbuf:this%nMesh(2)+this%nbuf,&
                           1-this%nbuf:this%nMesh(3)+this%nbuf,this%nvar)::q,q1,q2
integer::dd
procedure(limiter),pointer::slope=>null()
procedure(fluxSolver),pointer::fluxPtr=>null()
double precision::dtddx
double precision::rhoL,rhoM,rhoR
double precision::vxL,vxM,vxR
double precision::vyL,vyM,vyR
double precision::vzL,vzM,vzR
double precision,dimension(:,:,:,:),allocatable::SL,Fx
double precision::flux(this%nvar),ql(this%nvar),qr(this%nvar),slopeL(this%nvar),slopeR(this%nvar)
integer::i,j,k,nbuf,nx,ny,nz,nvar
integer::cx,cy,cz,f1,f2,f3,f4
double precision:: coef(this%nvar),g1,g2,g3,g4

nx=this%nMesh(1)
ny=this%nMesh(2)
nz=this%nMesh(3)
nvar=this%nvar
nbuf=this%nbuf

if(this%solverType .eq. 1) then
   fluxPtr=>fluxExactIsoHD3D
endif
if(this%solverType .eq. 2) then
   fluxPtr=>fluxHLLIsoHD3D
endif

select case (this%limiterType)
   case(0)
     slope=>zslop
   case(1)
     slope=>vslop
   case(2)
     slope=>fslop
   case(3)
     slope=>minmod
end select

allocate(SL(1-nbuf:nx+nbuf,1-nbuf:ny+nbuf,1-nbuf:nz+nbuf,nvar))
allocate(Fx(1-nbuf:nx+nbuf,1-nbuf:ny+nbuf,1-nbuf:nz+nbuf,nvar))

if(dd .eq. 1) then
   cx=1
   cy=0
   cz=0
   f1=1
   f2=2
   f3=3
   f4=4
   coef(1)=1.d0
   coef(2)=1.d0
   coef(3)=1.d0
   coef(4)=1.d0
   g1=1.d0
   g2=1.d0
   g3=1.d0
   g4=1.d0
elseif (dd .eq. 2) then
   cx=0
   cy=1
   cz=0
   f1=1
   f2=3
   f3=2
   f4=4
   coef(1)=1.d0
   coef(2)=-1.d0
   coef(3)=1.d0
   coef(4)=1.d0
   g1=1.d0
   g2=1.d0
   g3=-1.d0
   g4=1.d0
elseif (dd .eq. 3) then
   cx=0
   cy=0
   cz=1
   f1=1
   f2=4
   f3=3
   f4=2
   coef(1)=1.d0
   coef(2)=-1.d0
   coef(3)=1.d0
   coef(4)=1.d0
   g1=1.d0
   g2=1.d0
   g3=1.d0
   g4=-1.d0
endif

do k=0,nz+1
  do j=0,ny+1
    do i=0,nx+1
      rhoL=g1*q(i-cx,j-cy,k-cz,f1)
      rhoM=g1*q(i   ,j   ,k   ,f1)
      rhoR=g1*q(i+cx,j+cy,k+cz,f1)

      vxL=g2*q(i-cx,j-cy,k-cz,f2)/rhoL
      vxM=g2*q(i   ,j   ,k   ,f2)/rhoM
      vxR=g2*q(i+cx,j+cy,k+cz,f2)/rhoR

      vyL=g3*q(i-cx,j-cy,k-cz,f3)/rhoL
      vyM=g3*q(i   ,j   ,k   ,f3)/rhoM
      vyR=g3*q(i+cx,j+cy,k+cz,f3)/rhoR
  
      vzL=g4*q(i-cx,j-cy,k-cz,f4)/rhoL
      vzM=g4*q(i   ,j   ,k   ,f4)/rhoM
      vzR=g4*q(i+cx,j+cy,k+cz,f4)/rhoR

      SL(i,j,k,1)=slope(rhoL,rhoM,rhoR)
      SL(i,j,k,2)=slope( vxL, vxM, vxR)
      SL(i,j,k,3)=slope( vyL, vyM, vyR)
      SL(i,j,k,4)=slope( vzL, vzM, vzR)
    enddo !! end i
  enddo ! end j
enddo ! end k

do k=0,nz+1
  do j=0,ny+1
    do i=0,nx+1
      ql(1)=g1*q(i,j,k,f1)
      ql(2)=g2*q(i,j,k,f2)/ql(1)
      ql(3)=g3*q(i,j,k,f3)/ql(1)
      ql(4)=g4*q(i,j,k,f4)/ql(1)

      qr(1)=g1*q(i+cx,j+cy,k+cz,f1)
      qr(2)=g2*q(i+cx,j+cy,k+cz,f2)/qr(1)
      qr(3)=g3*q(i+cx,j+cy,k+cz,f3)/qr(1)
      qr(4)=g4*q(i+cx,j+cy,k+cz,f4)/qr(1)

      slopeL(:)=SL(i,j,k,:)
      slopeR(:)=SL(i+cx,j+cy,k+cz,:)

      call fluxPtr(ql,qr,slopeL,slopeR,this%snd,flux,this%nvar)
      Fx(i,j,k,:)=flux(:)

    enddo ! end i
  enddo ! end j
enddo ! end k

do k=1,nz
  do j=1,ny
    do i=1,nx
      dtddx=this%dt/(this%dx(1)%coords(i)*cx+this%dx(2)%coords(j)*cy+this%dx(3)%coords(k)*cz)
      q2(i,j,k,1)=q1(i,j,k,1)+coef(1)*(Fx(i-cx,j-cy,k-cz,f1)-Fx(i,j,k,f1))*dtddx
      q2(i,j,k,2)=q1(i,j,k,2)+coef(2)*(Fx(i-cx,j-cy,k-cz,f2)-Fx(i,j,k,f2))*dtddx
      q2(i,j,k,3)=q1(i,j,k,3)+coef(3)*(Fx(i-cx,j-cy,k-cz,f3)-Fx(i,j,k,f3))*dtddx
      q2(i,j,k,4)=q1(i,j,k,4)+coef(4)*(Fx(i-cx,j-cy,k-cz,f4)-Fx(i,j,k,f4))*dtddx
    enddo
  enddo
enddo

deallocate(SL,Fx)
end subroutine solverIso3D


   !call rieSolvern(nthis,qn,qn,qn1,dd=1)
   !call rieSolvern(nthis,qn,qn1,qn1,dd=2)
   !call rieSolvern(nthis,qn,qn1,qn1,dd=3)



subroutine solverAdi3D(this,q,q1,q2,dd)
implicit none
class(grid)::this
double precision,dimension(1-this%nbuf:this%nMesh(1)+this%nbuf,1-this%nbuf:this%nMesh(2)+this%nbuf,&
&                          1-this%nbuf:this%nMesh(3)+this%nbuf,this%nvar)::q,q1,q2
integer::dd
procedure(limiter),pointer::slope=>null()
procedure(fluxSolver),pointer::fluxPtr=>null()
double precision::dtddx
double precision::rhoL,rhoM,rhoR
double precision::vxL,vxM,vxR
double precision::vyL,vyM,vyR
double precision::vzL,vzM,vzR
double precision::eneL,eneM,eneR
double precision::pL,pM,pR,gam
double precision,dimension(:,:,:,:),allocatable::SL,Fx
double precision::flux(this%nvar),ql(this%nvar),qr(this%nvar),slopeL(this%nvar),slopeR(this%nvar)
integer::i,j,k,nbuf,nx,ny,nz,nvar
integer::cx,cy,cz,f1,f2,f3,f4,f5
double precision:: coef(this%nvar),g1,g2,g3,g4,g5

nx=this%nMesh(1)
ny=this%nMesh(2)
nz=this%nMesh(3)

nvar=this%nvar
nbuf=this%nbuf
gam=this%adiGamma

if(this%solverType .eq. 1) then
   !fluxPtr=>fluxExactAdiHD2D
endif
if(this%solverType .eq. 2) then
   fluxPtr=>fluxHLLAdiHD3D
endif
if(this%solverType .eq. 3) then
   fluxPtr=>fluxHLLCAdiHD3D
endif

select case (this%limiterType)
   case(0)
     slope=>zslop
   case(1)
     slope=>vslop
   case(2)
     slope=>fslop
   case(3)
     slope=>minmod
end select

allocate(SL(1-nbuf:nx+nbuf,1-nbuf:ny+nbuf,1-nbuf:nz+nbuf,nvar))
allocate(Fx(1-nbuf:nx+nbuf,1-nbuf:ny+nbuf,1-nbuf:nz+nbuf,nvar))
if(dd .eq. 1) then
  cx=1
  cy=0
  cz=0
  f1=1
  f2=2
  f3=3
  f4=4
  f5=5
  coef(1)=1.d0
  coef(2)=1.d0
  coef(3)=1.d0
  coef(4)=1.d0
  coef(5)=1.d0
  g1=1.d0
  g2=1.d0
  g3=1.d0
  g4=1.d0
  g5=1.d0
elseif(dd .eq. 2) then
  cx=0
  cy=1
  cz=0
  f1=1
  f2=3
  f3=2
  f4=4
  f5=5
  coef(1)=1.d0
  coef(2)=-1.d0
  coef(3)=1.d0
  coef(4)=1.d0
  coef(5)=1.d0
  g1=1.d0
  g2=1.d0
  g3=-1.d0
  g4=1.d0
  g5=1.d0
elseif(dd .eq. 3) then
  cx=0
  cy=0
  cz=1
  f1=1
  f2=4
  f3=3
  f4=2
  f5=5
  coef(1)=1.d0
  coef(2)=-1.d0
  coef(3)=1.d0
  coef(4)=1.d0
  coef(5)=1.d0
  g1=1.d0
  g2=1.d0
  g3=1.d0
  g4=-1.d0
  g5=1.d0
endif

do k=0,nz+1
  do j=0,ny+1
    do i=0,nx+1
      rhoL=g1*q(i-cx,j-cy,k-cz,f1)
      rhoM=g1*q(i   ,j   ,k   ,f1)
      rhoR=g1*q(i+cx,j+cy,k+cz,f1)

      vxL=g2*q(i-cx,j-cy,k-cz,f2)/rhoL
      vxM=g2*q(i   ,j   ,k   ,f2)/rhoM
      vxR=g2*q(i+cx,j+cy,k+cz,f2)/rhoR

      vyL=g3*q(i-cx,j-cy,k-cz,f3)/rhoL
      vyM=g3*q(i   ,j   ,k   ,f3)/rhoM
      vyR=g3*q(i+cx,j+cy,k+cz,f3)/rhoR
  
      vzL=g4*q(i-cx,j-cy,k-cz,f4)/rhoL
      vzM=g4*q(i   ,j   ,k   ,f4)/rhoM
      vzR=g4*q(i+cx,j+cy,k+cz,f4)/rhoR

     eneL=g5*q(i-cx,j-cy,k-cz,f5)/rhoL
     eneM=g5*q(i   ,j   ,k   ,f5)/rhoM
     eneR=g5*q(i+cx,j+cy,k+cz,f5)/rhoR

       pL=(gam-1.d0)*(rhoL*eneL-0.5d0*rhoL*(vxL**2+vyL**2+vzL**2))
       pM=(gam-1.d0)*(rhoM*eneM-0.5d0*rhoM*(vxM**2+vyM**2+vzM**2))
       pR=(gam-1.d0)*(rhoR*eneR-0.5d0*rhoR*(vxR**2+vyR**2+vzR**2))

      SL(i,j,k,1)=slope(rhoL,rhoM,rhoR)
      SL(i,j,k,2)=slope( vxL, vxM, vxR)
      SL(i,j,k,3)=slope( vyL, vyM, vyR)
      SL(i,j,k,4)=slope( vzL, vzM, vzR)
      SL(i,j,k,5)=slope(  pL,  pM,  pR)
    enddo !! end i
  enddo ! end j
enddo ! end k

do k=0,nz+1
  do j=0,ny+1
    do i=0,nx+1
      ql(1)=g1*q(i,j,k,f1)
      ql(2)=g2*q(i,j,k,f2)/ql(1)
      ql(3)=g3*q(i,j,k,f3)/ql(1)
      ql(4)=g4*q(i,j,k,f4)/ql(1)
      ql(5)=g5*(gam-1.d0)*(q(i,j,k,f5)-0.5d0*ql(1)*(ql(2)**2+ql(3)**2+ql(4)**2))

      qr(1)=g1*q(i+cx,j+cy,k+cz,f1)
      qr(2)=g2*q(i+cx,j+cy,k+cz,f2)/qr(1)
      qr(3)=g3*q(i+cx,j+cy,k+cz,f3)/qr(1)
      qr(4)=g4*q(i+cx,j+cy,k+cz,f4)/qr(1)
      qr(5)=g5*(gam-1.d0)*(q(i+cx,j+cy,k+cz,f5)-0.5d0*qr(1)*(qr(2)**2+qr(3)**2+qr(4)**2))

      slopeL(:)=SL(i,j,k,:)
      slopeR(:)=SL(i+cx,j+cy,k+cz,:)

      call fluxPtr(ql,qr,slopeL,slopeR,gam,flux,this%nvar)
      Fx(i,j,k,:)=flux(:)

    enddo ! end i
  enddo ! end j
enddo ! end k

do k=1,nz
  do j=1,ny
    do i=1,nx
      dtddx=this%dt/(this%dx(1)%coords(i)*cx+this%dx(2)%coords(j)*cy+this%dx(3)%coords(k)*cz)
      q2(i,j,k,1)=q1(i,j,k,1)+coef(1)*(Fx(i-cx,j-cy,k-cz,f1)-Fx(i,j,k,f1))*dtddx
      q2(i,j,k,2)=q1(i,j,k,2)+coef(2)*(Fx(i-cx,j-cy,k-cz,f2)-Fx(i,j,k,f2))*dtddx
      q2(i,j,k,3)=q1(i,j,k,3)+coef(3)*(Fx(i-cx,j-cy,k-cz,f3)-Fx(i,j,k,f3))*dtddx
      q2(i,j,k,4)=q1(i,j,k,4)+coef(4)*(Fx(i-cx,j-cy,k-cz,f4)-Fx(i,j,k,f4))*dtddx
      q2(i,j,k,5)=q1(i,j,k,5)+coef(5)*(Fx(i-cx,j-cy,k-cz,f5)-Fx(i,j,k,f5))*dtddx
    enddo
  enddo
enddo

deallocate(SL,Fx)
end subroutine solverAdi3D

subroutine fluxHLLIsoHD3D(ql,qr,slopeL,slopeR,snd,flux,nvar)
implicit none
integer::nvar
double precision::ql(NVAR),qr(NVAR),slopeL(NVAR),slopeR(NVAR),flux(NVAR)
double precision, dimension(:), allocatable:: FL,FR,Fs,UL,UR
double precision:: spd,snd ! sound speed
double precision:: rhoL,vxL,vyL,vzL,pTL
double precision:: rhoR,vxR,vyR,vzR,pTR
double precision:: SL,SR

allocate(FL(NVAR))
allocate(FR(NVAR))
allocate(Fs(NVAR))
allocate(UL(NVAR))
allocate(UR(NVAR))

UL=ql+0.5d0*slopeL
UR=qr-0.5d0*slopeR

spd=snd
rhoL=UL(1)
vxL =UL(2)
vyL =UL(3)
vzL =UL(4)
pTL =rhoL*spd**2.d0

rhoR=UR(1)
vxR =UR(2)
vyR =UR(3)
vzR =UR(4)
pTR =rhoR*spd**2.d0

SL=dmin1(vxL-spd,vxR-spd)
SR=dmax1(vxL+spd,vxR+spd)

FL(1)=rhoL*vxL
FL(2)=rhoL*vxL**2.d0+pTL
FL(3)=rhoL*vxL*vyL
FL(4)=rhoL*vxL*vzL

FR(1)=rhoR*vxR
FR(2)=rhoR*vxR**2.d0+pTR
FR(3)=rhoR*vxR*vyR
FR(4)=rhoR*vxR*vzR

Fs(1)=(SR*FL(1)-SL*FR(1)+SR*SL*(UR(1)-UL(1)))/(SR-SL)
Fs(2)=(SR*FL(2)-SL*FR(2)+SR*SL*(UR(1)*UR(2)-UL(1)*UL(2)))/(SR-SL)
Fs(3)=(SR*FL(3)-SL*FR(3)+SR*SL*(UR(1)*UR(3)-UL(1)*UL(3)))/(SR-SL)
Fs(4)=(SR*FL(4)-SL*FR(4)+SR*SL*(UR(1)*UR(4)-UL(1)*UL(4)))/(SR-SL)

if(SL .gt. 0.d0) then
   flux=FL
elseif(SL .le. 0.d0 .and. SR .ge. 0.d0) then
   flux=Fs
else
   flux=FR
endif


deallocate(FL)
deallocate(FR)
deallocate(Fs)
deallocate(UL)
deallocate(UR)
end subroutine fluxHLLIsoHD3D


subroutine fluxExactIsoHD3D(ql,qr,slopeL,slopeR,snd,flux,nvar)
implicit none
integer::nvar
double precision::ql(nvar),qr(nvar),qm(nvar),slopeL(nvar),slopeR(nvar),flux(nvar),snd
double precision:: rhoL,rhoR,vxL,vxR,vyL,vyR,vzL,vzR
double precision::c,shkR,rarR,csq,vstar,rhostar,rho0,rhosR,srhosR,aisrhosR,sigma2, &
                   rhosL,srhosL,aisrhosL,sigma1,aa,bb,uu,srhostar
integer:: Intr,n
double precision:: tolerr = 1.d-10

ql(1) = ql(1)+0.5d0*slopeL(1)
ql(2) = ql(2)+0.5d0*slopeL(2)
ql(3) = ql(3)+0.5d0*slopeL(3)
ql(4) = ql(4)+0.5d0*slopeL(4)

qr(1) = qr(1)-0.5d0*slopeR(1)
qr(2) = qr(2)-0.5d0*slopeR(2)
qr(3) = qr(3)-0.5d0*slopeR(3)
qr(4) = qr(4)-0.5d0*slopeR(4)

c = snd
Intr = 10

if (ql(1).eq.qr(1).and.ql(2).eq.qr(2)) then
     qm(1) = qr(1)
     qm(2) = qr(2)
if(ql(2).ge. 0.   ) then
     qm(3) = ql(3)
     qm(4) = ql(4)
else
     qm(3) = qr(3)
     qm(4) = qr(4)
endif
goto 1000
endif

      rhoL = ql(1)
      rhoR = qr(1)
       vxL = ql(2)
       vxR = qr(2)
       vyL = ql(3)
       vyR = qr(3)
       vzL = ql(4)
       vzR = qr(4)

      shkR = -c*(rhoR-rhoL)/dsqrt(rhoR*rhoL)
      rarR = -c*dlog(rhoR/rhoL)
       csq =  c*c

      if(rhoL .ge. rhoR) then
        if ((vxR-vxL) .gt. rarR) then
          goto 100
        elseif ((vxR-vxL) .gt. -shkR) then
          goto 200
        else
          goto 400
        endif
      else
        if ((vxR-vxL) .lt. shkR ) then
          goto 400
        elseif ((vxR-vxL) .lt. -rarR) then
          goto 300
        else
          goto 100
        endif
       endif
!cccccc 1-rarefaction wave, 2-rarefaction wave
!cccccc 1-rarefaction wave, 2-rarefaction wave
!cccccc 1-rarefaction wave, 2-rarefaction wave
100   vstar   = 0.5d0*(vxL+vxR)+0.5d0*c*dlog(rhoL/rhoR)
      rhostar = rhoL*dexp(-(vstar-vxL)/c)

      if ((vxL-c) .ge. 0.) then
         qm(2) = vxL
         qm(1) = rhoL
      elseif ((vstar-c) .ge. 0.) then
         qm(2) =  c
         qm(1) =  rhoL*dexp((vxL-c)/c)
      elseif ((vstar+c) .ge. 0.) then
         qm(2) =  vstar
         qm(1) =  rhostar
      elseif ((  vxR+c) .ge. 0.) then
         qm(2) =  -c
         qm(1) =  rhoR*dexp(-(vxR+c)/c)
      else
         qm(2) =   vxR
         qm(1) =  rhoR
      endif
      if ( qm(2) .ge. 0.) then
         qm(3) =  ql(3)
         qm(4) =  ql(4)
      else
         qm(3) =  qr(3)
         qm(4) =  qr(4)
      endif
      goto 1000
!cccccc1-rarefaction wave, 2-shock wave
!cccccc1-rarefaction wave, 2-shock wave
!cccccc1-rarefaction wave, 2-shock wave
200   rho0 = rhoR
      do n = 1, INtr
        rhosR    = rho0/rhoR
        srhosR   = dsqrt(rhosR)
        aisrhosR = 1.d0/srhosR
        rhostar  = rho0 - &
       (vxR-vxL+c*(srhosR-aisrhosR)+c*dlog(rho0/rhoL)) &
       /(0.5d0*c/rho0*(srhosR+aisrhosR)+c/rho0)
        if (dabs((rho0-rhostar)/rho0) .lt. tolerr) then
            goto 210
        endif
        rho0 = rhostar
      enddo
      write(6,*) 'Newton Interation is divergent at x-dir 1R2S'
      write(6,*) 'rhostar,rhoL=',rhostar,rhoL
      write(6,*) 'uxL    , uxR=',vxL,      vxR
      stop
210   vstar  = vxL - c*dlog(rhostar/rhoL)
      sigma2 = vxR + c*dsqrt(rhostar/rhoR)
      if ( sigma2 .le. 0.) then
         qm(2) =  vxR
         qm(1) = rhoR
      elseif ((vstar-c) .le. 0.) then
         qm(2) = vstar
         qm(1) = rhostar
      elseif ((vxL-c) .le. 0.) then
         qm(2) =  c
         qm(1) =  rhoL*dexp((vxL-c)/c)
      else
         qm(2) = vxL
         qm(1) = rhoL
      endif
      if ( qm(2) .ge. 0.) then
         qm(3) = ql(3)
         qm(4) = ql(4)
      else
         qm(3) = qr(3)
         qm(4) = qr(4)
      endif
      goto 1000
!cccccc 1-shock wave, 2-rarefaction wave
!cccccc 1-shock wave, 2-rarefaction wave
!cccccc 1-shock wave, 2-rarefaction wave
300   rho0 = rhoL
      do n = 1, INtr
        rhosL    = rho0/rhoL
        srhosL   = dsqrt(rhosL)
        aisrhosL = 1.d0/srhosL
        rhostar  = rho0- &
      (vxR-vxL+c*(srhosL-aisrhosL)+c*dlog(rho0/rhoR)) &
       /(0.5d0*c/rho0*(srhosL+aisrhosL)+c/rho0)
        if (dabs((rho0-rhostar)/rho0) .lt. tolerr) then
            goto 310
        endif
        rho0 = rhostar
      enddo
      write(6,*) 'Newton Interation is divergent at x-dir 1S2R'
      write(6,*) 'rhostar,rhoL=',rhostar,rhoL
      write(6,*) 'uxL    , uxR=',vxL,      vxR
      stop
310   vstar  = vxR + c*dlog(rhostar/rhoR)
      sigma1 = vxL - c*dsqrt(rhostar/rhoL)
      if (sigma1 .ge. 0.) then
         qm(2) = vxL
         qm(1) = rhoL
      elseif ((vstar+c) .ge. 0.) then
         qm(2) = vstar
         qm(1) = rhostar
      elseif ((vxR+c) .ge. 0.) then
         qm(2) =  -c
         qm(1) =  rhoR*dexp(-(vxR+c)/c)
      else
         qm(2) = vxR
         qm(1) = rhoR
      endif
      if ( qm(2) .ge. 0.) then
         qm(3) = ql(3)
         qm(4) = ql(4)
      else
         qm(3) = qr(3)
         qm(4) = qr(4)
      endif
      goto 1000
!       1-shock wave, 2-shock wave
!       1-shock wave, 2-shock wave
!       1-shock wave, 2-shock wave
400   aa       = 1.d0/dsqrt(rhoL)+1.d0/dsqrt(rhoR)
      bb       =    dsqrt(rhoL)+   dsqrt(rhoR)
      uu       = vxR-vxL
      srhostar = (-uu+dsqrt(uu**2+4.d0*csq*aa*bb))/(2.d0*c*aa)
      rhostar  = srhostar*srhostar
      vstar    = vxL - c*(dsqrt(rhostar/rhoL)-dsqrt(rhoL/rhostar))
      sigma1   = vxL - c*dsqrt(rhostar/rhoL)
      sigma2   = vxR + c*dsqrt(rhostar/rhoR)
      if     (sigma1 .ge. 0.) then
       qm(2) = vxL
       qm(1) = rhoL
      elseif (sigma2 .ge. 0.) then
       qm(2) = vstar
       qm(1) = rhostar
      else
       qm(2) = vxR
       qm(1) = rhoR
      endif
         if  (qm(2) .ge. 0.) then
       qm(3) = ql(3)
       qm(4) = ql(4)
      else
       qm(3) = qr(3)
       qm(4) = qr(4)
      endif
      goto 1000

1000  flux(1) = qm(1)*(   qm(2)          )
      flux(2) = qm(1)*( qm(2)*qm(2) + c*c)
      flux(3) = qm(1)*( qm(2)*qm(3)      )
      flux(4) = qm(1)*( qm(2)*qm(4)      )

end subroutine fluxExactIsoHD3D

subroutine fluxHLLAdiHD3D(ql,qr,slopeL,slopeR,gam,flux,nvar)
implicit none
integer::nvar
double precision::ql(NVAR),qr(NVAR),slopeL(NVAR),slopeR(NVAR),flux(NVAR)
double precision, dimension(:), allocatable:: FL,FR,Fs,UL,UR
double precision:: spd,gam ! sound speed
double precision:: rhoL,vxL,vyL,vzL,eneL,spdL,pTL,pL
double precision:: rhoR,vxR,vyR,vzR,eneR,spdR,pTR,pR
double precision:: SL,SR

allocate(FL(NVAR))
allocate(FR(NVAR))
allocate(Fs(NVAR))
allocate(UL(NVAR))
allocate(UR(NVAR))

UL=ql+0.5d0*slopeL
UR=qr-0.5d0*slopeR

rhoL=UL(1)
vxL =UL(2)
vyL =UL(3)
vzL =UL(4)
pL  =UL(5)
eneL = (pL/(gam-1.d0)+0.5d0*rhoL*(vxL**2+vyL**2+vzL**2))/rhoL
pTL = rhoL*(eneL-0.5d0*(vxL**2.d0+vyL**2.d0+vzL**2.d0))*(gam-1.d0)
spdL = dsqrt(gam*pTL/rhoL)

rhoR=UR(1)
vxR =UR(2)
vyR =UR(3)
vzR =UR(4)
pR  =UR(5)
eneR = (pR/(gam-1.d0)+0.5d0*rhoR*(vxR**2+vyR**2+vzR**2))/rhoR
pTR = rhoR*(eneR-0.5d0*(vxR**2.d0+vyR**2.d0+vzR**2.d0))*(gam-1.d0)
spdR = dsqrt(gam*pTR/rhoR)

SL=dmin1(vxL-spdL,vxR-spdR)
SR=dmax1(vxL+spdL,vxR+spdR)

FL(1)=rhoL*vxL
FL(2)=rhoL*vxL**2.d0+pTL
FL(3)=rhoL*vxL*vyL
FL(4)=rhoL*vxL*vzL
FL(5)=(rhoL*eneL+pTL)*vxL

FR(1)=rhoR*vxR
FR(2)=rhoR*vxR**2.d0+pTR
FR(3)=rhoR*vxR*vyR
FR(4)=rhoR*vxR*vzR
FR(5)=(rhoR*eneR+pTR)*vxR

Fs(1)=(SR*FL(1)-SL*FR(1)+SR*SL*(UR(1)-UL(1)))/(SR-SL)
Fs(2)=(SR*FL(2)-SL*FR(2)+SR*SL*(UR(1)*UR(2)-UL(1)*UL(2)))/(SR-SL)
Fs(3)=(SR*FL(3)-SL*FR(3)+SR*SL*(UR(1)*UR(3)-UL(1)*UL(3)))/(SR-SL)
Fs(4)=(SR*FL(4)-SL*FR(4)+SR*SL*(UR(1)*UR(4)-UL(1)*UL(4)))/(SR-SL)
!Fs(5)=(SR*FL(5)-SL*FR(5)+SR*SL*(UR(1)*eneR-UL(1)*eneL))/(SR-SL)
Fs(5)=(SR*FL(5)-SL*FR(5)+SR*SL*(UR(1)*UR(5)-UL(1)*UL(5)))/(SR-SL)

if(SL .gt. 0.d0) then
   flux=FL
elseif(SL .le. 0.d0 .and. SR .ge. 0.d0) then
   flux=Fs
else
   flux=FR
endif


deallocate(FL)
deallocate(FR)
deallocate(Fs)
deallocate(UL)
deallocate(UR)
end subroutine fluxHLLAdiHD3D

subroutine fluxHLLCAdiHD3D(ql,qr,slopeL,slopeR,gam,flux,nvar)
implicit none
integer::nvar
double precision::flux(nvar),ql(nvar),qr(nvar),slopeL(nvar),slopeR(nvar), &
                  FL(nvar),FR(nvar),FsL(nvar), FsR(nvar),Fs(nvar)
double precision::rhoL,rhoR,vxL,vxR,seneL,seneR,eneL,eneR,pTL,pTR,pxL,pxR, &
             spdL,spdR,SL,SR,SM,pstarL,pstarR,rhosL,rhosR,pxsL,pxsR,enesL,enesR, &
             vyL,vyR,pyL,pyR,pysL,pysR,pL,pR,gam, &
             vzL,vzR,pzL,pzR,pzsL,pzsR

rhoL =ql(1)+0.5d0*slopeL(1)
vxL  =ql(2)+0.5d0*slopeL(2)
vyL  =ql(3)+0.5d0*slopeL(3)
vzL  =ql(4)+0.5d0*slopeL(4)
pL   =ql(5)+0.5d0*slopeL(5)
eneL =pL/(gam-1.d0)+0.5d0*rhoL*(vxL**2+vyL**2+vzL**2)
!seneL=ql(4)+0.5d0*slopeL(4)
!eneL =seneL*rhoL

rhoR =qr(1)-0.5d0*slopeR(1)
vxR  =qr(2)-0.5d0*slopeR(2)
vyR  =qr(3)-0.5d0*slopeR(3)
vzR  =qr(4)-0.5d0*slopeR(4)
pR   =qr(5)-0.5d0*slopeR(5)
eneR =pR/(gam-1.d0)+0.5d0*rhoR*(vxR**2+vyR**2+vzR**2)
!seneR=qr(4)-0.5d0*slopeR(4)
!eneR =seneR*rhoR

pTL = (gam-1.d0)*(eneL-0.5d0*rhoL*(vxL**2+vyL**2+vzL**2))
pTR = (gam-1.d0)*(eneR-0.5d0*rhoR*(vxR**2+vyR**2+vzR**2))

!if(pTL .le. 0.d0 .or. pTR .le. 0.d0) then
if(pL .le. 0.d0 .or. pR .le. 0.d0) then
rhoL =ql(1)
vxL  =ql(2)
vyL  =ql(3)
vzL  =ql(4)
pL   =ql(5)
eneL =pL/(gam-1.d0)+0.5d0*rhoL*(vxL**2+vyL**2+vzL**2)
!seneL=ql(4)
!eneL =seneL*rhoL
rhoR =qr(1)
vxR  =qr(2)
vyR  =qr(3)
vzR  =qr(4)
pR   =qr(5)
eneR =pR/(gam-1.d0)+0.5d0*rhoR*(vxR**2+vyR**2+vzR**2)
!seneR=qr(4)
!eneR =seneR*rhoR

pTL = (gam-1.d0)*(eneL-0.5d0*rhoL*(vxL**2+vyL**2+vzL**2))
pTR = (gam-1.d0)*(eneR-0.5d0*rhoR*(vxR**2+vyR**2+vzR**2))
!write(*,*) "negative pressure!"
endif

spdL = dsqrt(gam*pTL/rhoL)
pxL = vxL*rhoL
pyL = vyL*rhoL
pzL = vzL*rhoL

spdR = dsqrt(gam*pTR/rhoR)
pxR = vxR*rhoR
pyR = vyR*rhoR
pzR = vzR*rhoR

SL = dmin1(vxL-spdL,vxR-spdR)
SR = dmax1(vxL+spdL,vxR+spdR)
SM = ((SR-vxR)*rhoR*vxR-(SL-vxL)*rhoL*vxL-pTR+pTL)  &
        /((SR-vxR)*rhoR-(SL-vxL)*rhoL)
pstarL = rhoL*(vxL-SL)*(vxL-SM)+pTL
pstarR = rhoR*(vxR-SR)*(vxR-SM)+pTR

rhosL = rhoL*(SL-vxL)/(SL-SM)
pxsL = ((SL-vxL)*pxL+(pstarL-pTL))/(SL-SM)
pysL = rhosL*vyL
pzsL = rhosL*vzL
enesL = ((SL-vxL)*eneL-pTL*vxL+pstarL*SM)/(SL-SM)

rhosR = rhoR*(SR-vxR)/(SR-SM)
pxsR = ((SR-vxR)*pxR+(pstarR-pTR))/(SR-SM)
pysR = rhosR*vyR
pzsR = rhosR*vzR
enesR = ((SR-vxR)*eneR-pTR*vxR+pstarR*SM)/(SR-SM)

FL(1) = rhoL*vxL
FL(2) = rhoL*vxL**2.d0 +pTL
FL(3) = rhoL*vxL*vyL
FL(4) = rhoL*vxL*vzL
FL(5) = (eneL+pTL)*vxL

FR(1) = rhoR*vxR
FR(2) = rhoR*vxR**2.d0 +pTR
FR(3) = rhoR*vxR*vyR
FR(4) = rhoR*vxR*vzR
FR(5) = (eneR+pTR)*vxR

FsL(1)= FL(1) + SL*(rhosL-rhoL)
FsL(2)= FL(2) + SL*(pxsL-pxL)
FsL(3)= FL(3) + SL*(pysL-pyL)
FsL(4)= FL(4) + SL*(pzsL-pzL)
FsL(5)= FL(5) + SL*(enesL-eneL)

FsR(1)= FR(1) + SR*(rhosR-rhoR)
FsR(2)= FR(2) + SR*(pxsR-pxR)
FsR(3)= FR(3) + SR*(pysR-pyR)
FsR(4)= FR(4) + SR*(pzsR-pzR)
FsR(5)= FR(5) + SR*(enesR-eneR)

IF (SL .gt. 0.d0)then
       flux=FL
ELSEIF (SL .le. 0.d0 .and. SM .ge. 0.d0)then
       flux=FsL
ELSEIF (SM .le. 0.d0 .and. SR .ge. 0.d0)then
       flux=FsR
ELSEIF (SR .lt. 0.d0) then
       flux=FR
ENDIF
end subroutine fluxHLLCAdiHD3D
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!! Riemann solvers for 2D problems
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine solverIso2D(this,q,q1,q2,dd)
class(grid)::this
double precision,dimension(1-this%nbuf:this%nMesh(1)+this%nbuf,1-this%nbuf:this%nMesh(2)+this%nbuf,this%nvar)::q,q1,q2
integer::dd
procedure(limiter),pointer::slope=>null()
procedure(limiter3),pointer::slope3=>null()
procedure(fluxSolver),pointer::fluxPtr=>null()
double precision::dtddx
double precision::rhoL,rhoM,rhoR
double precision::vxL,vxM,vxR
double precision::vyL,vyM,vyR
double precision,dimension(:,:,:),allocatable::SL,Fx
double precision::flux(this%nvar),ql(this%nvar),qr(this%nvar),slopeL(this%nvar),slopeR(this%nvar)
integer::i,j,nbuf,nx,ny,nvar
integer::cx,cy,f1,f2,f3
double precision:: coef(this%nvar),g1,g2,g3
double precision:: ri,rm,ro,dr,dphi,fac1,fac2,vphi,vr,snd,den
integer::coordType

nx=this%nMesh(1)
ny=this%nMesh(2)
nvar=this%nvar
nbuf=this%nbuf
coordType=this%coordType

if(this%solverType .eq. 1) then
   fluxPtr=>fluxExactIsoHD2D
endif
if(this%solverType .eq. 2) then
   fluxPtr=>fluxHLLIsoHD2D
endif

select case (this%limiterType)
   case(0)
     slope=>zslop
   case(1)
     slope=>vslop
   case(2)
     slope=>fslop
   case(3)
     slope=>minmod
end select

  if(coordType .eq. 2 .or. coordType .eq. 3) then
    select case (this%limiterType)
      case(0)
        slope3=>zslop3
      case(1)
        print *,"solverIso2D: not yet implemented...."
        stop
      case(2)
        print *,"solverIso2D: not yet implemented...."
        stop
      case(3)
        slope3=>minmod3
    end select
  endif

allocate(SL(1-nbuf:nx+nbuf,1-nbuf:ny+nbuf,nvar))
allocate(Fx(1-nbuf:nx+nbuf,1-nbuf:ny+nbuf,nvar))

if(dd .eq. 1) then
   cx=1
   cy=0
   f1=1
   f2=2
   f3=3
   coef(1)=1.d0
   coef(2)=1.d0
   coef(3)=1.d0
   g1=1.d0
   g2=1.d0
   g3=1.d0
elseif (dd .eq. 2) then
   cx=0
   cy=1
   f1=1
   f2=3
   f3=2
   coef(1)=1.d0
   coef(2)=-1.d0
   coef(3)=1.d0
   g1=1.d0
   g2=1.d0
   g3=-1.d0
endif


if(coordType .eq. 1) then !! Cartesian coordinates
do j=0,ny+1
  do i=0,nx+1
    rhoL=g1*q(i-cx,j-cy,f1)
    rhoM=g1*q(i   ,j   ,f1)
    rhoR=g1*q(i+cx,j+cy,f1)

    vxL=g2*q(i-cx,j-cy,f2)/rhoL
    vxM=g2*q(i   ,j   ,f2)/rhoM
    vxR=g2*q(i+cx,j+cy,f2)/rhoR

    vyL=g3*q(i-cx,j-cy,f3)/rhoL
    vyM=g3*q(i   ,j   ,f3)/rhoM
    vyR=g3*q(i+cx,j+cy,f3)/rhoR

    SL(i,j,1)=slope(rhoL,rhoM,rhoR)
    SL(i,j,2)=slope( vxL, vxM, vxR)
    SL(i,j,3)=slope( vyL, vyM, vyR)
  enddo !! end i
enddo ! end j
endif ! if coordType .eq. 1

if(coordType .eq. 2 .or. coordType .eq. 3) then !! polar coordinates, uniform and log in r
  if(dd .eq. 1) then
    do j=0,ny+1
      do i=0,nx+1
        rhoL=g1*q(i-cx,j-cy,f1)
        rhoM=g1*q(i   ,j   ,f1)
        rhoR=g1*q(i+cx,j+cy,f1)

        vxL=g2*q(i-cx,j-cy,f2)/rhoL
        vxM=g2*q(i   ,j   ,f2)/rhoM
        vxR=g2*q(i+cx,j+cy,f2)/rhoR

        vyL=g3*q(i-cx,j-cy,f3)/rhoL
        vyM=g3*q(i   ,j   ,f3)/rhoM
        vyR=g3*q(i+cx,j+cy,f3)/rhoR

         ri=this%xc(1)%coords(i-cx)
         rm=this%xc(1)%coords(i   )
         ro=this%xc(1)%coords(i+cx)
         dr=this%dx(1)%coords(i   )

        !SL(i,j,1)=slope(rhoL,rhoM,rhoR)
        !SL(i,j,2)=slope( vxL, vxM, vxR)
        !SL(i,j,3)=slope( vyL, vyM, vyR)

        SL(i,j,1)=slope3(rhoL,rhoM,rhoR,ri,rm,ro,dr)
        SL(i,j,2)=slope3( vxL, vxM, vxR,ri,rm,ro,dr)
        SL(i,j,3)=slope3( vyL, vyM, vyR,ri,rm,ro,dr)

      enddo
    enddo
  endif  !! if dd .eq. 1
  
  if(dd .eq. 2) then
    do j=0,ny+1
      do i=0,nx+1
        rhoL=g1*q(i-cx,j-cy,f1)
        rhoM=g1*q(i   ,j   ,f1)
        rhoR=g1*q(i+cx,j+cy,f1)

        vxL=g2*q(i-cx,j-cy,f2)/rhoL
        vxM=g2*q(i   ,j   ,f2)/rhoM
        vxR=g2*q(i+cx,j+cy,f2)/rhoR

        vyL=g3*q(i-cx,j-cy,f3)/rhoL
        vyM=g3*q(i   ,j   ,f3)/rhoM
        vyR=g3*q(i+cx,j+cy,f3)/rhoR


        SL(i,j,1)=slope(rhoL,rhoM,rhoR)
        SL(i,j,2)=slope( vxL, vxM, vxR)
        SL(i,j,3)=slope( vyL, vyM, vyR)
      enddo
    enddo
  endif !! if dd. eq. 2

endif

do j=0,ny+1
  do i=0,nx+1
    ql(1)=g1*q(i,j,f1)
    ql(2)=g2*q(i,j,f2)/ql(1)
    ql(3)=g3*q(i,j,f3)/ql(1)

    qr(1)=g1*q(i+cx,j+cy,f1)
    qr(2)=g2*q(i+cx,j+cy,f2)/qr(1)
    qr(3)=g3*q(i+cx,j+cy,f3)/qr(1)

    slopeL(:)=SL(i,j,:)
    slopeR(:)=SL(i+cx,j+cy,:)

    call fluxPtr(ql,qr,slopeL,slopeR,this%snd,flux,this%nvar)
    Fx(i,j,:)=flux(:)

  enddo ! end i
enddo ! end j

if(coordType .eq. 1) then
  do j=1,ny
    do i=1,nx
      dtddx=this%dt/(this%dx(1)%coords(i)*cx+this%dx(2)%coords(j)*cy)
      q2(i,j,1)=q1(i,j,1)+coef(1)*(Fx(i-cx,j-cy,f1)-Fx(i,j,f1))*dtddx
      q2(i,j,2)=q1(i,j,2)+coef(2)*(Fx(i-cx,j-cy,f2)-Fx(i,j,f2))*dtddx
      q2(i,j,3)=q1(i,j,3)+coef(3)*(Fx(i-cx,j-cy,f3)-Fx(i,j,f3))*dtddx
    enddo
  enddo
endif !! coordType .eq. 1

if(coordType .eq. 3 .or. coordType .eq. 2) then
  if(dd .eq. 1) then  !! radial direction
    do j=1,ny
      do i=1,nx
        rm=this%xc(1)%coords(i)
        dr=this%dx(1)%coords(i)
        fac1=rm-0.5d0*dr
        fac2=rm+0.5d0*dr
        dtddx=this%dt/(rm*dr)
        den=q(i,j,1)
         vr=q(i,j,2)/den
       vphi=q(i,j,3)/den
        snd=this%snd

        q2(i,j,1)=q1(i,j,1)+coef(1)*(Fx(i-cx,j-cy,f1)*fac1-Fx(i,j,f1)*fac2)*dtddx
        q2(i,j,2)=q1(i,j,2)+coef(2)*(Fx(i-cx,j-cy,f2)*fac1-Fx(i,j,f2)*fac2)*dtddx+(den*vphi**2+den*snd**2)/rm*this%dt
        q2(i,j,3)=q1(i,j,3)+coef(3)*(Fx(i-cx,j-cy,f3)*fac1**2-Fx(i,j,f3)*fac2**2)*this%dt/(dr*rm**2)
        !q2(i,j,3)=q1(i,j,3)+coef(3)*(Fx(i-cx,j-cy,f3)*fac1-Fx(i,j,f3)*fac2)*dtddx-(den*vr*vphi)/rm*this%dt
      enddo
    enddo
  endif  !!! dd .eq. 1

  if(dd .eq. 2) then  !! phi direction
    do j=1,ny
      do i=1,nx
        rm=this%xc(1)%coords(i)
        dphi=this%dx(2)%coords(j)
        dtddx=this%dt/(dphi*rm)
        q2(i,j,1)=q1(i,j,1)+coef(1)*(Fx(i-cx,j-cy,f1)-Fx(i,j,f1))*dtddx
        q2(i,j,2)=q1(i,j,2)+coef(2)*(Fx(i-cx,j-cy,f2)-Fx(i,j,f2))*dtddx
        q2(i,j,3)=q1(i,j,3)+coef(3)*(Fx(i-cx,j-cy,f3)-Fx(i,j,f3))*dtddx
      enddo
    enddo
  endif  !!! dd .eq. 2
endif !!! coordType .eq. 3

deallocate(SL,Fx)

end subroutine solverIso2D

subroutine solverPoly2D(this,q,q1,q2,dd)
class(grid)::this
double precision,dimension(1-this%nbuf:this%nMesh(1)+this%nbuf,1-this%nbuf:this%nMesh(2)+this%nbuf,this%nvar)::q,q1,q2
integer::dd
procedure(limiter),pointer::slope=>null()
procedure(limiter3),pointer::slope3=>null()
procedure(fluxSolver),pointer::fluxPtr=>null()
double precision::dtddx
double precision::rhoL,rhoM,rhoR
double precision::vxL,vxM,vxR
double precision::vyL,vyM,vyR
double precision,dimension(:,:,:),allocatable::SL,Fx
double precision::flux(this%nvar),ql(this%nvar),qr(this%nvar),slopeL(this%nvar),slopeR(this%nvar)
integer::i,j,nbuf,nx,ny,nvar
integer::cx,cy,f1,f2,f3
double precision:: coef(this%nvar),g1,g2,g3
double precision:: ri,rm,ro,dr,dphi,fac1,fac2,vphi,vr,snd,den
integer::coordType

nx=this%nMesh(1)
ny=this%nMesh(2)
nvar=this%nvar
nbuf=this%nbuf
coordType=this%coordType

if(this%solverType .eq. 6) then
   !fluxPtr=>fluxHLLIsoHD2D
endif

select case (this%limiterType)
   case(0)
     slope=>zslop
   case(1)
     slope=>vslop
   case(2)
     slope=>fslop
   case(3)
     slope=>minmod
end select

  if(coordType .eq. 2 .or. coordType .eq. 3) then
    select case (this%limiterType)
      case(0)
        slope3=>zslop3
      case(1)
        print *,"solverIso2D: not yet implemented...."
        stop
      case(2)
        print *,"solverIso2D: not yet implemented...."
        stop
      case(3)
        slope3=>minmod3
    end select
  endif

allocate(SL(1-nbuf:nx+nbuf,1-nbuf:ny+nbuf,nvar))
allocate(Fx(1-nbuf:nx+nbuf,1-nbuf:ny+nbuf,nvar))

if(dd .eq. 1) then
   cx=1
   cy=0
   f1=1
   f2=2
   f3=3
   coef(1)=1.d0
   coef(2)=1.d0
   coef(3)=1.d0
   g1=1.d0
   g2=1.d0
   g3=1.d0
elseif (dd .eq. 2) then
   cx=0
   cy=1
   f1=1
   f2=3
   f3=2
   coef(1)=1.d0
   coef(2)=-1.d0
   coef(3)=1.d0
   g1=1.d0
   g2=1.d0
   g3=-1.d0
endif


if(coordType .eq. 1) then !! Cartesian coordinates
do j=0,ny+1
  do i=0,nx+1
    rhoL=g1*q(i-cx,j-cy,f1)
    rhoM=g1*q(i   ,j   ,f1)
    rhoR=g1*q(i+cx,j+cy,f1)

    vxL=g2*q(i-cx,j-cy,f2)/rhoL
    vxM=g2*q(i   ,j   ,f2)/rhoM
    vxR=g2*q(i+cx,j+cy,f2)/rhoR

    vyL=g3*q(i-cx,j-cy,f3)/rhoL
    vyM=g3*q(i   ,j   ,f3)/rhoM
    vyR=g3*q(i+cx,j+cy,f3)/rhoR

    SL(i,j,1)=slope(rhoL,rhoM,rhoR)
    SL(i,j,2)=slope( vxL, vxM, vxR)
    SL(i,j,3)=slope( vyL, vyM, vyR)
  enddo !! end i
enddo ! end j
endif ! if coordType .eq. 1

if(coordType .eq. 2 .or. coordType .eq. 3) then !! polar coordinates, uniform and log in r
  if(dd .eq. 1) then
    do j=0,ny+1
      do i=0,nx+1
        rhoL=g1*q(i-cx,j-cy,f1)
        rhoM=g1*q(i   ,j   ,f1)
        rhoR=g1*q(i+cx,j+cy,f1)

        vxL=g2*q(i-cx,j-cy,f2)/rhoL
        vxM=g2*q(i   ,j   ,f2)/rhoM
        vxR=g2*q(i+cx,j+cy,f2)/rhoR

        vyL=g3*q(i-cx,j-cy,f3)/rhoL
        vyM=g3*q(i   ,j   ,f3)/rhoM
        vyR=g3*q(i+cx,j+cy,f3)/rhoR

         ri=this%xc(1)%coords(i-cx)
         rm=this%xc(1)%coords(i   )
         ro=this%xc(1)%coords(i+cx)
         dr=this%dx(1)%coords(i   )

        !SL(i,j,1)=slope(rhoL,rhoM,rhoR)
        !SL(i,j,2)=slope( vxL, vxM, vxR)
        !SL(i,j,3)=slope( vyL, vyM, vyR)

        SL(i,j,1)=slope3(rhoL,rhoM,rhoR,ri,rm,ro,dr)
        SL(i,j,2)=slope3( vxL, vxM, vxR,ri,rm,ro,dr)
        SL(i,j,3)=slope3( vyL, vyM, vyR,ri,rm,ro,dr)

      enddo
    enddo
  endif  !! if dd .eq. 1
  
  if(dd .eq. 2) then
    do j=0,ny+1
      do i=0,nx+1
        rhoL=g1*q(i-cx,j-cy,f1)
        rhoM=g1*q(i   ,j   ,f1)
        rhoR=g1*q(i+cx,j+cy,f1)

        vxL=g2*q(i-cx,j-cy,f2)/rhoL
        vxM=g2*q(i   ,j   ,f2)/rhoM
        vxR=g2*q(i+cx,j+cy,f2)/rhoR

        vyL=g3*q(i-cx,j-cy,f3)/rhoL
        vyM=g3*q(i   ,j   ,f3)/rhoM
        vyR=g3*q(i+cx,j+cy,f3)/rhoR


        SL(i,j,1)=slope(rhoL,rhoM,rhoR)
        SL(i,j,2)=slope( vxL, vxM, vxR)
        SL(i,j,3)=slope( vyL, vyM, vyR)
      enddo
    enddo
  endif !! if dd. eq. 2

endif

do j=0,ny+1
  do i=0,nx+1
    ql(1)=g1*q(i,j,f1)
    ql(2)=g2*q(i,j,f2)/ql(1)
    ql(3)=g3*q(i,j,f3)/ql(1)

    qr(1)=g1*q(i+cx,j+cy,f1)
    qr(2)=g2*q(i+cx,j+cy,f2)/qr(1)
    qr(3)=g3*q(i+cx,j+cy,f3)/qr(1)

    slopeL(:)=SL(i,j,:)
    slopeR(:)=SL(i+cx,j+cy,:)

    !call fluxPtr(ql,qr,slopeL,slopeR,this%snd,flux,this%nvar)
    call fluxHLLPolyHD2D(ql,qr,slopeL,slopeR,this%polyK,this%polyGamma,flux,this%nvar)
    Fx(i,j,:)=flux(:)

  enddo ! end i
enddo ! end j

if(coordType .eq. 1) then
  do j=1,ny
    do i=1,nx
      dtddx=this%dt/(this%dx(1)%coords(i)*cx+this%dx(2)%coords(j)*cy)
      q2(i,j,1)=q1(i,j,1)+coef(1)*(Fx(i-cx,j-cy,f1)-Fx(i,j,f1))*dtddx
      q2(i,j,2)=q1(i,j,2)+coef(2)*(Fx(i-cx,j-cy,f2)-Fx(i,j,f2))*dtddx
      q2(i,j,3)=q1(i,j,3)+coef(3)*(Fx(i-cx,j-cy,f3)-Fx(i,j,f3))*dtddx
    enddo
  enddo
endif !! coordType .eq. 1

if(coordType .eq. 3 .or. coordType .eq. 2) then
  if(dd .eq. 1) then  !! radial direction
    do j=1,ny
      do i=1,nx
        rm=this%xc(1)%coords(i)
        dr=this%dx(1)%coords(i)
        fac1=rm-0.5d0*dr
        fac2=rm+0.5d0*dr
        dtddx=this%dt/(rm*dr)
        den=q(i,j,1)
         vr=q(i,j,2)/den
       vphi=q(i,j,3)/den
        snd=this%snd

        q2(i,j,1)=q1(i,j,1)+coef(1)*(Fx(i-cx,j-cy,f1)*fac1-Fx(i,j,f1)*fac2)*dtddx
        q2(i,j,2)=q1(i,j,2)+coef(2)*(Fx(i-cx,j-cy,f2)*fac1-Fx(i,j,f2)*fac2)*dtddx+(den*vphi**2+den*snd**2)/rm*this%dt
        q2(i,j,3)=q1(i,j,3)+coef(3)*(Fx(i-cx,j-cy,f3)*fac1**2-Fx(i,j,f3)*fac2**2)*this%dt/(dr*rm**2)
        !q2(i,j,3)=q1(i,j,3)+coef(3)*(Fx(i-cx,j-cy,f3)*fac1-Fx(i,j,f3)*fac2)*dtddx-(den*vr*vphi)/rm*this%dt
      enddo
    enddo
  endif  !!! dd .eq. 1

  if(dd .eq. 2) then  !! phi direction
    do j=1,ny
      do i=1,nx
        rm=this%xc(1)%coords(i)
        dphi=this%dx(2)%coords(j)
        dtddx=this%dt/(dphi*rm)
        q2(i,j,1)=q1(i,j,1)+coef(1)*(Fx(i-cx,j-cy,f1)-Fx(i,j,f1))*dtddx
        q2(i,j,2)=q1(i,j,2)+coef(2)*(Fx(i-cx,j-cy,f2)-Fx(i,j,f2))*dtddx
        q2(i,j,3)=q1(i,j,3)+coef(3)*(Fx(i-cx,j-cy,f3)-Fx(i,j,f3))*dtddx
      enddo
    enddo
  endif  !!! dd .eq. 2
endif !!! coordType .eq. 3

deallocate(SL,Fx)

end subroutine solverPoly2D


subroutine solverAdiMHD2D(this,q,q1,q2,dd)
class(grid)::this
double precision,dimension(1-this%nbuf:this%nMesh(1)+this%nbuf,1-this%nbuf:this%nMesh(2)+this%nbuf,this%nvar)::q,q1,q2
integer::nx,ny,nbuf,nvar,dd
integer::i,j
integer::cx,cy
integer::f1,f2,f3,f4,f5,f6,f7,f8,f9,f10,f11
double precision::g1,g2,g3,g4,g5,g6,g7,g8,g9,g10,g11
double precision::dt,dx,dtddx,gam
double precision::rhoL,rhoM,rhoR
double precision::vxL,vxM,vxR
double precision::vyL,vyM,vyR
double precision::vzL,vzM,vzR
double precision::bxL,bxM,bxR,bxtmp
double precision::byL,byM,byR
double precision::bzL,bzM,bzR
double precision::eneL,eneM,eneR
double precision::pL ,pM ,pR
double precision::pressure,bxc,byc,bzc,vxc,vyc,vzc,denc,enec
double precision::ri,rm,ro,dr,dphi,fac1,fac2
double precision,dimension(:,:,:),allocatable::SL
double precision,dimension(:,:,:),allocatable::Fx
double precision,dimension(:,:),allocatable::Ez,Ezc
double precision::ql(this%nvar),qr(this%nvar),slopeL(this%nvar),slopeR(this%nvar),flux(this%nvar),coef(this%nvar)

procedure(limiter),pointer::slope=>null()
procedure(limiter3),pointer::slope3=>null()
procedure(fluxSolver),pointer::fluxPtr=>null()

double precision::rho,vx,vy,bx,by,vr,vphi,ene,ptot,vz,brc,bphic
integer::chgFluxPtrCount
integer::coordType

nx=this%nMesh(1)
ny=this%nMesh(2)

nbuf=this%nbuf
nvar=this%nvar

dt=this%dt
gam=this%adiGamma
chgFluxPtrCount=0
coordType=this%coordType

select case(this%solverType)
   case (4)
     fluxPtr=>fluxHLLAdiMHD1D
   case (5)
     fluxPtr=>fluxHLLDAdiMHD1D
   case default
     print *,'riemannSolverModule.f03: no appropriate 2D MHD solver found!'
end select

select case (this%limiterType)
   case(0)
     slope=>zslop
   case(1)
     slope=>vslop
   case(2)
     slope=>fslop
   case(3)
     slope=>minmod
   case default
     print *, "riemannSolverModule.f03: no appropriate slope limiter found!"
end select

if(coordType .eq. 2 .or. coordType .eq. 3) then
    select case (this%limiterType)
      case(0)
        slope3=>zslop3
      case(1)
        print *,"solverAdi2D: not yet implemented...."
        stop
      case(2)
        print *,"solverAdi2D: not yet implemented...."
        stop
      case(3)
        slope3=>minmod3
    end select  
endif

allocate(SL(1-nbuf:nx+nbuf,1-nbuf:ny+nbuf,nvar))
allocate(Fx(1-nbuf:nx+nbuf,1-nbuf:ny+nbuf,nvar))
allocate(Ez(1-nbuf:nx+nbuf,1-nbuf:ny+nbuf))
allocate(Ezc(1-nbuf:nx+nbuf,1-nbuf:ny+nbuf))

!!!!! dd=1 ==>x, dd=2 ==>y !!!!!

if(dd .eq. 1) then
     cx=1
     cy=0
     f1=1
     f2=2
     f3=3
     f4=4
     f5=5
     f6=6
     f7=7
     f8=8
     f9=9
    f10=10
    f11=11
 coef(1)=1.d0
 coef(2)=1.d0
 coef(3)=1.d0
 coef(4)=1.d0
 coef(5)=1.d0
 coef(6)=1.d0
 coef(7)=1.d0
 coef(8)=1.d0
     g1=1.d0
     g2=1.d0
     g3=1.d0
     g4=1.d0
     g5=1.d0
     g6=1.d0
     g7=1.d0
     g8=1.d0
     g9=1.d0
    g10=1.d0
    g11=1.d0
elseif (dd .eq. 2) then
     cx=0
     cy=1
     f1=1
     f2=3
     f3=2
     f4=4
     f5=6
     f6=5
     f7=7
     f8=8
     f9=10
    f10=9
    f11=11
 coef(1)=1.d0
 coef(2)=-1.d0
 coef(3)=1.d0
 coef(4)=1.d0
 coef(5)=-1.d0
 coef(6)=1.d0
 coef(7)=1.d0
 coef(8)=1.d0
     g1=1.d0
     g2=1.d0
     g3=-1.d0
     g4=1.d0
     g5=1.d0
     g6=-1.d0
     g7=1.d0
     g8=1.d0
     g9=1.d0
    g10=-1.d0
    g11=1.d0
endif

if(coordType .eq. 1) then  !! Cartesian coordinates
do j=0,ny+1
  do i=0,nx+1
    rhoL=g1*q(i-cx,j-cy,f1)
    rhoM=g1*q(i   ,j   ,f1)
    rhoR=g1*q(i+cx,j+cy,f1)

     vxL=g2*q(i-cx,j-cy,f2)/rhoL
     vxM=g2*q(i   ,j   ,f2)/rhoM
     vxR=g2*q(i+cx,j+cy,f2)/rhoR

     vyL=g3*q(i-cx,j-cy,f3)/rhoL
     vyM=g3*q(i   ,j   ,f3)/rhoM
     vyR=g3*q(i+cx,j+cy,f3)/rhoR

     vzL=g4*q(i-cx,j-cy,f4)/rhoL
     vzM=g4*q(i   ,j   ,f4)/rhoM
     vzR=g4*q(i+cx,j+cy,f4)/rhoR

     bxL=0.5d0*(g5*q(i-cx,j-cy,f5)+ g9*q(i-cx,j-cy,f9 ))
     bxM=0.5d0*(g5*q(i   ,j   ,f5)+ g9*q(i   ,j   ,f9 ))
     bxR=0.5d0*(g5*q(i+cx,j+cy,f5)+ g9*q(i+cx,j+cy,f9 ))

     byL=0.5d0*(g6*q(i-cx,j-cy,f6)+g10*q(i-cx,j-cy,f10))
     byM=0.5d0*(g6*q(i   ,j   ,f6)+g10*q(i   ,j   ,f10))
     byR=0.5d0*(g6*q(i+cx,j+cy,f6)+g10*q(i+cx,j+cy,f10))

     bzL=0.5d0*(g7*q(i-cx,j-cy,f7)+g11*q(i-cx,j-cy,f11))
     bzM=0.5d0*(g7*q(i   ,j   ,f7)+g11*q(i   ,j   ,f11))
     bzR=0.5d0*(g7*q(i+cx,j+cy,f7)+g11*q(i+cx,j+cy,f11))

    eneL=g8*q(i-cx,j-cy,f8)/rhoL
    eneM=g8*q(i   ,j   ,f8)/rhoM
    eneR=g8*q(i+cx,j+cy,f8)/rhoR

      pL=(gam-1.d0)*(rhoL*eneL-0.5d0*rhoL*(vxL**2+vyL**2+vzL**2)-0.5d0*(bxL**2+byL**2+bzL**2))
      pM=(gam-1.d0)*(rhoM*eneM-0.5d0*rhoM*(vxM**2+vyM**2+vzM**2)-0.5d0*(bxM**2+byM**2+bzM**2))
      pR=(gam-1.d0)*(rhoR*eneR-0.5d0*rhoR*(vxR**2+vyR**2+vzR**2)-0.5d0*(bxR**2+byR**2+bzR**2))

    SL(i,j,1)=slope(rhoL,rhoM,rhoR)
    SL(i,j,2)=slope( vxL, vxM, vxR)
    SL(i,j,3)=slope( vyL, vyM, vyR)
    SL(i,j,4)=slope( vzL, vzM, vzR)
    SL(i,j,5)=slope( bxL, bxM, bxR)
    SL(i,j,6)=slope( byL, byM, byR)
    SL(i,j,7)=slope( bzL, bzM, bzR)
    SL(i,j,8)=slope(  pL,  pM,  pR)

  enddo !! end do i
enddo !! end do j
endif !! if coordType .eq. 1

if(coordType .eq. 3 .or. coordType .eq. 2) then !! polar coordinates, uniform or log in r
  if(dd .eq. 1) then
    do j=0,ny+1
      do i=0,nx+1
        rhoL=g1*q(i-cx,j-cy,f1)
        rhoM=g1*q(i   ,j   ,f1)
        rhoR=g1*q(i+cx,j+cy,f1)

         vxL=g2*q(i-cx,j-cy,f2)/rhoL
         vxM=g2*q(i   ,j   ,f2)/rhoM
         vxR=g2*q(i+cx,j+cy,f2)/rhoR

         vyL=g3*q(i-cx,j-cy,f3)/rhoL
         vyM=g3*q(i   ,j   ,f3)/rhoM
         vyR=g3*q(i+cx,j+cy,f3)/rhoR

         vzL=g4*q(i-cx,j-cy,f4)/rhoL
         vzM=g4*q(i   ,j   ,f4)/rhoM
         vzR=g4*q(i+cx,j+cy,f4)/rhoR

          ri=this%xl(1)%coords(i-cx)
          rm=this%xc(1)%coords(i-cx)
          ro=this%xr(1)%coords(i-cx)
         bxL=0.5d0*(ri*g5*q(i-cx,j-cy,f5)+ ro*g9*q(i-cx,j-cy,f9 ))/rm

          ri=this%xl(1)%coords(i)
          rm=this%xc(1)%coords(i)
          ro=this%xr(1)%coords(i)
         bxM=0.5d0*(ri*g5*q(i   ,j   ,f5)+ ro*g9*q(i   ,j   ,f9 ))/rm

          ri=this%xl(1)%coords(i+cx)
          rm=this%xc(1)%coords(i+cx)
          ro=this%xr(1)%coords(i+cx)
         bxR=0.5d0*(ri*g5*q(i+cx,j+cy,f5)+ ro*g9*q(i+cx,j+cy,f9 ))/rm

         byL=0.5d0*(g6*q(i-cx,j-cy,f6)+g10*q(i-cx,j-cy,f10))
         byM=0.5d0*(g6*q(i   ,j   ,f6)+g10*q(i   ,j   ,f10))
         byR=0.5d0*(g6*q(i+cx,j+cy,f6)+g10*q(i+cx,j+cy,f10))

         bzL=0.5d0*(g7*q(i-cx,j-cy,f7)+g11*q(i-cx,j-cy,f11))
         bzM=0.5d0*(g7*q(i   ,j   ,f7)+g11*q(i   ,j   ,f11))
         bzR=0.5d0*(g7*q(i+cx,j+cy,f7)+g11*q(i+cx,j+cy,f11))

        eneL=g8*q(i-cx,j-cy,f8)/rhoL
        eneM=g8*q(i   ,j   ,f8)/rhoM
        eneR=g8*q(i+cx,j+cy,f8)/rhoR

          pL=(gam-1.d0)*(rhoL*eneL-0.5d0*rhoL*(vxL**2+vyL**2+vzL**2)-0.5d0*(bxL**2+byL**2+bzL**2))
          pM=(gam-1.d0)*(rhoM*eneM-0.5d0*rhoM*(vxM**2+vyM**2+vzM**2)-0.5d0*(bxM**2+byM**2+bzM**2))
          pR=(gam-1.d0)*(rhoR*eneR-0.5d0*rhoR*(vxR**2+vyR**2+vzR**2)-0.5d0*(bxR**2+byR**2+bzR**2))

        SL(i,j,1)=slope(rhoL,rhoM,rhoR)
        SL(i,j,2)=slope( vxL, vxM, vxR)
        SL(i,j,3)=slope( vyL, vyM, vyR)
        SL(i,j,4)=slope( vzL, vzM, vzR)
        SL(i,j,5)=slope( bxL, bxM, bxR)
        SL(i,j,6)=slope( byL, byM, byR)
        SL(i,j,7)=slope( bzL, bzM, bzR)
        SL(i,j,8)=slope(  pL,  pM,  pR)

        ri=this%xl(1)%coords(i-cx)
        rm=this%xc(1)%coords(i)
        ro=this%xr(1)%coords(i+cx)
        dr=this%dx(1)%coords(i)

        !SL(i,j,1)=slope3(rhoL,rhoM,rhoR,ri,rm,ro,dr)
        !SL(i,j,2)=slope3( vxL, vxM, vxR,ri,rm,ro,dr)
        !SL(i,j,3)=slope3( vyL, vyM, vyR,ri,rm,ro,dr)
        !SL(i,j,4)=slope3( vzL, vzM, vzR,ri,rm,ro,dr)
        !SL(i,j,5)=slope3( bxL, bxM, bxR,ri,rm,ro,dr)
        !SL(i,j,6)=slope3( byL, byM, byR,ri,rm,ro,dr)
        !SL(i,j,7)=slope3( bzL, bzM, bzR,ri,rm,ro,dr)
        !SL(i,j,8)=slope3(  pL,  pM,  pR,ri,rm,ro,dr)
       enddo !! end do i
    enddo !! end do j        
  endif !! dd .eq. 1 , in r direction

  if(dd .eq. 2) then
    do j=0,ny+1
      do i=0,nx+1
        rhoL=g1*q(i-cx,j-cy,f1)
        rhoM=g1*q(i   ,j   ,f1)
        rhoR=g1*q(i+cx,j+cy,f1)

         vxL=g2*q(i-cx,j-cy,f2)/rhoL
         vxM=g2*q(i   ,j   ,f2)/rhoM
         vxR=g2*q(i+cx,j+cy,f2)/rhoR

         vyL=g3*q(i-cx,j-cy,f3)/rhoL
         vyM=g3*q(i   ,j   ,f3)/rhoM
         vyR=g3*q(i+cx,j+cy,f3)/rhoR

         vzL=g4*q(i-cx,j-cy,f4)/rhoL
         vzM=g4*q(i   ,j   ,f4)/rhoM
         vzR=g4*q(i+cx,j+cy,f4)/rhoR

         bxL=0.5d0*(g5*q(i-cx,j-cy,f5)+ g9*q(i-cx,j-cy,f9 ))
         bxM=0.5d0*(g5*q(i   ,j   ,f5)+ g9*q(i   ,j   ,f9 ))
         bxR=0.5d0*(g5*q(i+cx,j+cy,f5)+ g9*q(i+cx,j+cy,f9 ))

          ri=this%xl(1)%coords(i) 
          rm=this%xc(1)%coords(i)
          ro=this%xr(1)%coords(i)
         byL=0.5d0*(ri*g6*q(i-cx,j-cy,f6)+ro*g10*q(i-cx,j-cy,f10))/rm
         byM=0.5d0*(ri*g6*q(i   ,j   ,f6)+ro*g10*q(i   ,j   ,f10))/rm
         byR=0.5d0*(ri*g6*q(i+cx,j+cy,f6)+ro*g10*q(i+cx,j+cy,f10))/rm

         bzL=0.5d0*(g7*q(i-cx,j-cy,f7)+g11*q(i-cx,j-cy,f11))
         bzM=0.5d0*(g7*q(i   ,j   ,f7)+g11*q(i   ,j   ,f11))
         bzR=0.5d0*(g7*q(i+cx,j+cy,f7)+g11*q(i+cx,j+cy,f11))

        eneL=g8*q(i-cx,j-cy,f8)/rhoL
        eneM=g8*q(i   ,j   ,f8)/rhoM
        eneR=g8*q(i+cx,j+cy,f8)/rhoR

          pL=(gam-1.d0)*(rhoL*eneL-0.5d0*rhoL*(vxL**2+vyL**2+vzL**2)-0.5d0*(bxL**2+byL**2+bzL**2))
          pM=(gam-1.d0)*(rhoM*eneM-0.5d0*rhoM*(vxM**2+vyM**2+vzM**2)-0.5d0*(bxM**2+byM**2+bzM**2))
          pR=(gam-1.d0)*(rhoR*eneR-0.5d0*rhoR*(vxR**2+vyR**2+vzR**2)-0.5d0*(bxR**2+byR**2+bzR**2))

        SL(i,j,1)=slope(rhoL,rhoM,rhoR)
        SL(i,j,2)=slope( vxL, vxM, vxR)
        SL(i,j,3)=slope( vyL, vyM, vyR)
        SL(i,j,4)=slope( vzL, vzM, vzR)
        SL(i,j,5)=slope( bxL, bxM, bxR)
        SL(i,j,6)=slope( byL, byM, byR)
        SL(i,j,7)=slope( bzL, bzM, bzR)
        SL(i,j,8)=slope(  pL,  pM,  pR)

      enddo !! end do i
    enddo !! end do j    
  endif  !! if dd .eq. 2, phi direction

endif !! coordType .eq. 2 or 3, polar coordinates

   SL(:,:,5)=0.d0
if(coordType .eq. 1) then
  do j=0,ny+1
    do i=0,nx+1
       ql(1)=g1*q(i,j,f1)
       ql(2)=g2*q(i,j,f2)/ql(1)
       ql(3)=g3*q(i,j,f3)/ql(1)
       ql(4)=g4*q(i,j,f4)/ql(1)
       ql(5)=0.5d0*(g5*q(i,j,f5)+ g9*q(i,j, f9))
       ql(6)=0.5d0*(g6*q(i,j,f6)+g10*q(i,j,f10))
       ql(7)=0.5d0*(g7*q(i,j,f7)+g11*q(i,j,f11))
       ql(8)=g8*(gam-1.d0)*(q(i,j,f8)-0.5d0*ql(1)*(ql(2)**2+ql(3)**2+ql(4)**2)-0.5d0*(ql(5)**2+ql(6)**2+ql(7)**2))

       qr(1)=g1*q(i+cx,j+cy,f1)
       qr(2)=g2*q(i+cx,j+cy,f2)/qr(1)
       qr(3)=g3*q(i+cx,j+cy,f3)/qr(1)
       qr(4)=g4*q(i+cx,j+cy,f4)/qr(1)
       qr(5)=0.5d0*(g5*q(i+cx,j+cy,f5)+ g9*q(i+cx,j+cy, f9))
       qr(6)=0.5d0*(g6*q(i+cx,j+cy,f6)+g10*q(i+cx,j+cy,f10))
       qr(7)=0.5d0*(g7*q(i+cx,j+cy,f7)+g11*q(i+cx,j+cy,f11))
       qr(8)=g8*(gam-1.d0)*(q(i+cx,j+cy,f8)-0.5d0*qr(1)*(qr(2)**2+qr(3)**2+qr(4)**2)-0.5d0*(qr(5)**2+qr(6)**2+qr(7)**2))

       ql(5)=g5*q(i,j,f9)
       qr(5)=g5*q(i+cx,j+cy,f5)

       slopeL(:)=SL(i   ,j   ,:)
       slopeR(:)=SL(i+cx,j+cy,:)

       call fluxPtr(ql,qr,slopeL,slopeR,this%adiGamma,flux,this%nvar)
       Fx(i,j,:)=flux(:)  
 
    enddo !! end do i
  end do !!end do j   
endif !! coordType .eq. 1

if(coordType .eq. 2 .or. coordType .eq. 3) then
  if(dd .eq. 1) then
    do j=0,ny+1
      do i=0,nx+1
         ql(1)=g1*q(i,j,f1)
         ql(2)=g2*q(i,j,f2)/ql(1)
         ql(3)=g3*q(i,j,f3)/ql(1)
         ql(4)=g4*q(i,j,f4)/ql(1)
         ql(5)=0.5d0*(g5*q(i,j,f5)+ g9*q(i,j, f9))
         ql(6)=0.5d0*(g6*q(i,j,f6)+g10*q(i,j,f10))
         ql(7)=0.5d0*(g7*q(i,j,f7)+g11*q(i,j,f11))
         ql(8)=g8*(gam-1.d0)*(q(i,j,f8)-0.5d0*ql(1)*(ql(2)**2+ql(3)**2+ql(4)**2)-0.5d0*(ql(5)**2+ql(6)**2+ql(7)**2))

         qr(1)=g1*q(i+cx,j+cy,f1)
         qr(2)=g2*q(i+cx,j+cy,f2)/qr(1)
         qr(3)=g3*q(i+cx,j+cy,f3)/qr(1)
         qr(4)=g4*q(i+cx,j+cy,f4)/qr(1)
         qr(5)=0.5d0*(g5*q(i+cx,j+cy,f5)+ g9*q(i+cx,j+cy, f9))
         qr(6)=0.5d0*(g6*q(i+cx,j+cy,f6)+g10*q(i+cx,j+cy,f10))
         qr(7)=0.5d0*(g7*q(i+cx,j+cy,f7)+g11*q(i+cx,j+cy,f11))
         qr(8)=g8*(gam-1.d0)*(q(i+cx,j+cy,f8)-0.5d0*qr(1)*(qr(2)**2+qr(3)**2+qr(4)**2)-0.5d0*(qr(5)**2+qr(6)**2+qr(7)**2))

         ql(5)=g5*q(i,j,f9)
         qr(5)=g5*q(i+cx,j+cy,f5)

         slopeL(:)=SL(i   ,j   ,:)
         slopeR(:)=SL(i+cx,j+cy,:)

         call fluxPtr(ql,qr,slopeL,slopeR,this%adiGamma,flux,this%nvar)
         Fx(i,j,:)=flux(:)

      enddo !! end do i
    end do !!end do j    
  endif !! dd .eq. 1
  
  if(dd .eq. 2) then
    do j=0,ny+1
      do i=0,nx+1
         ql(1)=g1*q(i,j,f1)
         ql(2)=g2*q(i,j,f2)/ql(1)
         ql(3)=g3*q(i,j,f3)/ql(1)
         ql(4)=g4*q(i,j,f4)/ql(1)
         ql(5)=0.5d0*(g5*q(i,j,f5)+ g9*q(i,j, f9))

         ri=this%xl(1)%coords(i)
         rm=this%xc(1)%coords(i)
         ro=this%xr(1)%coords(i)

         ql(6)=0.5d0*(ri*g6*q(i,j,f6)+ro*g10*q(i,j,f10))/rm
         ql(7)=0.5d0*(g7*q(i,j,f7)+g11*q(i,j,f11))
         ql(8)=g8*(gam-1.d0)*(q(i,j,f8)-0.5d0*ql(1)*(ql(2)**2+ql(3)**2+ql(4)**2)-0.5d0*(ql(5)**2+ql(6)**2+ql(7)**2))

         qr(1)=g1*q(i+cx,j+cy,f1)
         qr(2)=g2*q(i+cx,j+cy,f2)/qr(1)
         qr(3)=g3*q(i+cx,j+cy,f3)/qr(1)
         qr(4)=g4*q(i+cx,j+cy,f4)/qr(1)
         qr(5)=0.5d0*(g5*q(i+cx,j+cy,f5)+ g9*q(i+cx,j+cy, f9))
         qr(6)=0.5d0*(ri*g6*q(i+cx,j+cy,f6)+ro*g10*q(i+cx,j+cy,f10))/rm
         qr(7)=0.5d0*(g7*q(i+cx,j+cy,f7)+g11*q(i+cx,j+cy,f11))
         qr(8)=g8*(gam-1.d0)*(q(i+cx,j+cy,f8)-0.5d0*qr(1)*(qr(2)**2+qr(3)**2+qr(4)**2)-0.5d0*(qr(5)**2+qr(6)**2+qr(7)**2))

         ql(5)=g5*q(i,j,f9)
         qr(5)=g5*q(i+cx,j+cy,f5)

         slopeL(:)=SL(i   ,j   ,:)
         slopeR(:)=SL(i+cx,j+cy,:)

         call fluxPtr(ql,qr,slopeL,slopeR,this%adiGamma,flux,this%nvar)
         Fx(i,j,:)=flux(:)

      enddo !! end do i
    end do !!end do j
  endif !! dd .eq. 2
endif  !! if coordType .eq. 2 or 3


Ez=0.d0
Ezc=0.d0

!!!!! calculate EMF at the cell centers
if(coordType .eq. 1) then
  do j=1-nbuf,ny+nbuf
    do i=1-nbuf,nx+nbuf
      rho=q(i,j,1)
       vx=q(i,j,2)/rho
       vy=q(i,j,3)/rho
       bx=0.5d0*(q(i,j,5)+q(i,j,9 ))
       by=0.5d0*(q(i,j,6)+q(i,j,10))
       Ezc(i,j)=vy*bx-by*vx
    enddo
  enddo
endif

if(coordType .eq. 2 .or. coordType .eq. 3) then
  do j=1-nbuf,ny+nbuf
    do i=1-nbuf,nx+nbuf
      rho=q(i,j,1)
       vx=q(i,j,2)/rho
       vy=q(i,j,3)/rho
       ri=this%xl(1)%coords(i)
       rm=this%xc(1)%coords(i)
       ro=this%xr(1)%coords(i)
       bx=0.5d0*(ri*q(i,j,5)+ro*q(i,j,9 ))/rm
       by=0.5d0*(q(i,j,6)+q(i,j,10))
       Ezc(i,j)=vy*bx-by*vx
    enddo
  enddo
endif

!!!!!!!! Gardiner & Stone, JCP, 2005, 205, 509
!!!!!!!! Skinner & Ostriker, ApJS, 2010, 188, 290
if(dd .eq. 1) then
  if(coordType .eq. 1) then
    do j=0,ny
      do i=0,nx
       Ez(i,j)= 0.5d0*(-Fx(i,j,f6)-Fx(i,j+1,f6))*coef(6)
       Ez(i,j)=Ez(i,j)-0.25d0*(Ezc(i,j)+Ezc(i,j+1))
       !Ez(i,j)= 0.25d0*(-Fx(i,j,f6)-Fx(i,j+1,f6))*coef(6)

       !if(q(i,j,3) .gt. 0.d0 .and. q(i,j+1,3) .gt. 0.d0) then
       !  Ez(i,j)=Ez(i,j)+0.125d0*(-Fx(i,j,f6)*coef(6)-Ezc(i,j))
       !elseif(q(i,j,3) .lt. 0.d0 .and. q(i,j+1,3) .lt. 0.d0) then
       !  Ez(i,j)=Ez(i,j)+0.125d0*(-Fx(i,j+1,f6)*coef(6)-Ezc(i,j+1))
       !else
       !  Ez(i,j)=Ez(i,j)+0.0625d0*(-Fx(i,j,f6)*coef(6)-Ezc(i,j)-Fx(i,j+1,f6)*coef(6)-Ezc(i,j+1))
       !endif
       !
       !if(q(i+1,j,3) .gt. 0.d0 .and. q(i+1,j+1,3) .gt. 0.d0) then
       !   Ez(i,j)=Ez(i,j)-0.125d0*(Ezc(i+1,j)+Fx(i,j,f6)*coef(6))
       !elseif(q(i+1,j,3) .lt. 0.d0 .and. q(i+1,j+1,3) .lt. 0.d0) then
       !   Ez(i,j)=Ez(i,j)-0.125d0*(Ezc(i+1,j+1)+Fx(i,j+1,f6)*coef(6))
       !else
       !   Ez(i,j)=Ez(i,j)-0.0625d0*(Ezc(i+1,j)+Fx(i,j,f6)*coef(6)+Ezc(i+1,j+1)+Fx(i,j+1,f6)*coef(6))
       !endif
      enddo
    enddo
  endif !! coordType .eq. 1
  if(coordType .eq. 3) then
    do j=0,ny
      do i=0,nx
        rm=this%xc(1)%coords(i)
        ro=this%xr(1)%coords(i)
        Ez(i,j)=0.25d0*(-2.d0*Fx(i,j,f6)-2.d0*Fx(i,j+1,f6))*coef(6)
        Ez(i,j)=Ez(i,j)-0.125d0*(1.d0+ro/rm)*(Ezc(i,j)+Ezc(i,j+1))
      enddo
    enddo
  endif !! 
  if(coordType .eq. 2) then
    print *,"riemannSolverModule.f03: not yet implemented..."
  endif
elseif (dd .eq. 2) then
  if(coordType .eq. 1) then
    do j=0,ny
      do i=0,nx
         Ez(i,j)=0.5d0*(Fx(i+1,j,f5)+Fx(i,j,f5))*coef(5)
         Ez(i,j)=Ez(i,j)-0.25d0*(Ezc(i+1,j)+Ezc(i+1,j+1))
       !Ez(i,j)=0.25d0*(Fx(i+1,j,f5)+Fx(i,j,f5))*coef(5)

       !if(q(i,j,2) .gt. 0.d0 .and. q(i+1,j,2) .gt. 0.d0) then
       !  Ez(i,j)=Ez(i,j)+0.125d0*(Fx(i,j,f5)*coef(5)-Ezc(i,j))
       !elseif(q(i,j,2) .lt. 0.d0 .and. q(i+1,j,2) .lt. 0.d0) then
       !  Ez(i,j)=Ez(i,j)+0.125d0*(Fx(i+1,j,f5)*coef(5)-Ezc(i+1,j))
       !else
       !  Ez(i,j)=Ez(i,j)+0.0625d0*(Fx(i,j,f5)*coef(5)-Ezc(i,j)+Fx(i+1,j,f5)*coef(5)-Ezc(i+1,j))
       !endif       
       !
       !if(q(i,j+1,2) .gt. 0.d0 .and. q(i+1,j+1,2) .gt. 0.d0) then
       !  Ez(i,j)=Ez(i,j)-0.125d0*(Ezc(i,j+1)-Fx(i,j,f5)*coef(5))
       !elseif(q(i,j+1,2) .lt. 0.d0 .and. q(i+1,j+1,2) .lt. 0.d0) then
       !  Ez(i,j)=Ez(i,j)-0.125d0*(Ezc(i+1,j+1)-Fx(i+1,j,f5)*coef(5))
       !else
       !  Ez(i,j)=Ez(i,j)-0.0625d0*(Ezc(i,j+1)-Fx(i,j,f5)*coef(5)+Ezc(i+1,j+1)-Fx(i+1,j,f5)*coef(5))
       !endif
      enddo
    enddo
  endif !! coordType .eq. 1
  if(coordType .eq. 3) then
    do j=0,ny
      do i=0,nx
        ri=this%xc(1)%coords(i)
        ro=this%xc(1)%coords(i+1)
        rm=this%xl(1)%coords(i+1)
        Ez(i,j)=0.25d0*((1.d0+rm/ri)*Fx(i,j,f5)+(1.d0+rm/ro)*Fx(i+1,j,f5))*coef(5)
        Ez(i,j)=Ez(i,j)-0.125d0*(1.d0+rm/ro)*(Ezc(i+1,j)+Ezc(i+1,j+1))
      enddo
    enddo
  endif 
  if(coordType .eq. 2) then
    print *,"riemannSolverModule.f03: not yet implemented...."
  endif
endif

if(coordType .eq. 1) then
  do j=1,ny
    do i=1,nx
      dtddx=this%dt/(this%dx(1)%coords(i)*cx+this%dx(2)%coords(j)*cy)
      q2(i,j,1)=q1(i,j,1)+coef(1)*(Fx(i-cx,j-cy,f1)-Fx(i,j,f1))*dtddx
      q2(i,j,2)=q1(i,j,2)+coef(2)*(Fx(i-cx,j-cy,f2)-Fx(i,j,f2))*dtddx
      q2(i,j,3)=q1(i,j,3)+coef(3)*(Fx(i-cx,j-cy,f3)-Fx(i,j,f3))*dtddx
      q2(i,j,4)=q1(i,j,4)+coef(4)*(Fx(i-cx,j-cy,f4)-Fx(i,j,f4))*dtddx
      q2(i,j,5)=q1(i,j,5)+(Ez(i-1,j-1)-Ez(i-1,j  ))*dtddx
      q2(i,j,6)=q1(i,j,6)+(Ez(i  ,j-1)-Ez(i-1,j-1))*dtddx
      q2(i,j,7)=q1(i,j,7)+coef(7)*(Fx(i-cx,j-cy,f7)-Fx(i,j,f7))*dtddx
      q2(i,j,8)=q1(i,j,8)+coef(8)*(Fx(i-cx,j-cy,f8)-Fx(i,j,f8))*dtddx
      q2(i,j,9)=q1(i,j,  9)+(Ez(i,j-1)-Ez(i  ,j))*dtddx
      q2(i,j,10)=q1(i,j,10)+(Ez(i,j  )-Ez(i-1,j))*dtddx
      q2(i,j,11)=q2(i,j,7)
    enddo
  enddo
endif
if(coordType .eq. 3) then
  if(dd .eq. 1) then
    do j=1,ny
      do i=1,nx
        rm=this%xc(1)%coords(i)
        dr=this%dx(1)%coords(i)
       dphi=this%dx(2)%coords(j)
        fac1=rm-0.5d0*dr
        fac2=rm+0.5d0*dr   
        dtddx=this%dt/(rm*dr)
        rho=q(i,j,1)
         vr=q(i,j,2)/rho
       vphi=q(i,j,3)/rho
         vz=q(i,j,4)/rho
        brc=0.5d0*(fac1*q(i,j,5)+fac2*q(i,j,9))/rm
       bphic=0.5d0*(q(i,j,6)+q(i,j,10))
         bzc=0.5d0*(q(i,j,7)+q(i,j,11))
        ene=q(i,j,8)
        pressure=(ene-0.5d0*rho*(vr**2+vphi**2+vz**2)-0.5d0*(brc**2+bphic**2+bzc**2))*(gam-1.d0)
        ptot=pressure+0.5d0*(brc**2+bphic**2+bzc**2)

        q2(i,j,1)=q1(i,j,1)+coef(1)*(Fx(i-cx,j-cy,f1)*fac1-Fx(i,j,f1)*fac2)*dtddx
        q2(i,j,2)=q1(i,j,2)+coef(2)*(Fx(i-cx,j-cy,f2)*fac1-Fx(i,j,f2)*fac2)*dtddx &
                                   +(rho*vphi**2+ptot-bphic**2)/rm*this%dt
        q2(i,j,3)=q1(i,j,3)+coef(3)*(Fx(i-cx,j-cy,f3)*fac1**2-Fx(i,j,f3)*fac2**2)*this%dt/(dr*rm**2)
        q2(i,j,4)=q1(i,j,4)+coef(4)*(Fx(i-cx,j-cy,f4)*fac1-Fx(i,j,f4)*fac2)*dtddx
        q2(i,j,5)=q1(i,j,5)+(Ez(i-1,j-1)-Ez(i-1,j  ))*this%dt/(fac1*dphi)
        q2(i,j,6)=q1(i,j,6)+(Ez(i  ,j-1)-Ez(i-1,j-1))*this%dt/dr
        q2(i,j,7)=q1(i,j,7)+coef(7)*(Fx(i-cx,j-cy,f7)*fac1-Fx(i,j,f7)*fac2)*dtddx
        q2(i,j,8)=q1(i,j,8)+coef(8)*(Fx(i-cx,j-cy,f8)*fac1-Fx(i,j,f8)*fac2)*dtddx
        q2(i,j,9)=q1(i,j,  9)+(Ez(i,j-1)-Ez(i  ,j))*this%dt/(fac2*dphi)
        q2(i,j,10)=q1(i,j,10)+(Ez(i,j  )-Ez(i-1,j))*this%dt/dr
        q2(i,j,11)=q2(i,j,7)
      enddo
    enddo  
  endif !! dd .eq. 1

  if(dd .eq. 2) then
    do j=1,ny
      do i=1,nx
        rm=this%xc(1)%coords(i)
        dr=this%dx(1)%coords(i)
       dphi=this%dx(2)%coords(j)
        fac1=rm-0.5d0*dr
        fac2=rm+0.5d0*dr
        dtddx=this%dt/(rm*dphi)

        q2(i,j,1)=q1(i,j,1)+coef(1)*(Fx(i-cx,j-cy,f1)-Fx(i,j,f1))*dtddx
        q2(i,j,2)=q1(i,j,2)+coef(2)*(Fx(i-cx,j-cy,f2)-Fx(i,j,f2))*dtddx
        q2(i,j,3)=q1(i,j,3)+coef(3)*(Fx(i-cx,j-cy,f3)-Fx(i,j,f3))*dtddx
        q2(i,j,4)=q1(i,j,4)+coef(4)*(Fx(i-cx,j-cy,f4)-Fx(i,j,f4))*dtddx
        q2(i,j,5)=q1(i,j,5)+(Ez(i-1,j-1)-Ez(i-1,j  ))*this%dt/(fac1*dphi)
        q2(i,j,6)=q1(i,j,6)+(Ez(i  ,j-1)-Ez(i-1,j-1))*this%dt/dr
        q2(i,j,7)=q1(i,j,7)+coef(7)*(Fx(i-cx,j-cy,f7)-Fx(i,j,f7))*dtddx
        q2(i,j,8)=q1(i,j,8)+coef(8)*(Fx(i-cx,j-cy,f8)-Fx(i,j,f8))*dtddx
        q2(i,j,9)=q1(i,j,  9)+(Ez(i,j-1)-Ez(i  ,j))*this%dt/(fac2*dphi)
        q2(i,j,10)=q1(i,j,10)+(Ez(i,j  )-Ez(i-1,j))*this%dt/dr
        q2(i,j,11)=q2(i,j,7)
      enddo
    enddo    
  endif

endif
if(coordType .eq. 2) then
  print *,"riemannSolverModule.f03: not yet implemented....."
  stop
endif

select case(this%solverType)
  case(4)
  case(5)
   if(dd .eq. 2) then
    loopj:do j=1,ny
      do i=1,nx
        denc=q2(i,j,1)
         vxc=q2(i,j,2)/denc
         vyc=q2(i,j,3)/denc
         vzc=q2(i,j,4)/denc
         bxc=0.5d0*(q2(i,j,5)+q2(i,j,9))
         byc=0.5d0*(q2(i,j,6)+q2(i,j,10))
         bzc=0.5d0*(q2(i,j,7)+q2(i,j,11))
         enec=q2(i,j,8)
        pressure=(gam-1.d0)*(enec-0.5d0*denc*(vxc**2+vyc**2+vzc**2)-0.5d0*(bxc**2+byc**2+bzc**2))
        if(pressure .lt. 0.d0) then
          print *,"riemannSolverModule.f03: negative pressure..."
          print *,"(myid,dd,i,j,pressure)=",myid,dd,i,j,pressure
          this%changeSolver=.true.
          exit loopj
        endif
      enddo
    enddo loopj
   endif
  case default
   print *,"riemannSolverModule.f03: unknown MHD solver type..."
   stop
end select


deallocate(SL,Fx,Ez,Ezc)
end subroutine solverAdiMHD2D

subroutine solverAdi2D(this,q,q1,q2,dd)
class(grid)::this
double precision,dimension(1-this%nbuf:this%nMesh(1)+this%nbuf,1-this%nbuf:this%nMesh(2)+this%nbuf,this%nvar)::q,q1,q2
integer::dd
procedure(limiter),pointer::slope=>null()
procedure(limiter3),pointer::slope3=>null()
procedure(fluxSolver),pointer::fluxPtr=>null()
double precision::dtddx
double precision::rhoL,rhoM,rhoR
double precision::vxL,vxM,vxR
double precision::vyL,vyM,vyR
double precision::eneL,eneM,eneR
double precision::pL,pM,pR,gam
double precision,dimension(:,:,:),allocatable::SL,Fx
double precision::flux(this%nvar),ql(this%nvar),qr(this%nvar),slopeL(this%nvar),slopeR(this%nvar)
integer::i,j,nbuf,nx,ny,nvar
integer::cx,cy,f1,f2,f3,f4
double precision:: coef(this%nvar),g1,g2,g3,g4
double precision:: ri,rm,ro,dr,dphi,fac1,fac2,vphi,vr,ene,ptot,den
integer:: coordType

nx=this%nMesh(1)
ny=this%nMesh(2)
nvar=this%nvar
nbuf=this%nbuf
gam=this%adiGamma
coordType=this%coordType

if(this%solverType .eq. 1) then
   !fluxPtr=>fluxExactAdiHD2D  !! required to be implemented
endif
if(this%solverType .eq. 2) then
   fluxPtr=>fluxHLLAdiHD2D
endif
if(this%solverType .eq. 3) then
   fluxPtr=>fluxHLLCAdiHD2D
endif

  select case (this%limiterType)
    case(0)
      slope=>zslop
    case(1)
      slope=>vslop
    case(2)
      slope=>fslop
    case(3)
      slope=>minmod
  end select

  if(coordType .eq. 2 .or. coordType .eq. 3) then
    select case (this%limiterType)
      case(0)
        slope3=>zslop3
      case(1)
        print *,"solverAdi2D: not yet implemented...."
        stop
      case(2)
        print *,"solverAdi2D: not yet implemented...."
        stop
      case(3)
        slope3=>minmod3
    end select
  endif

allocate(SL(1-nbuf:nx+nbuf,1-nbuf:ny+nbuf,nvar))
allocate(Fx(1-nbuf:nx+nbuf,1-nbuf:ny+nbuf,nvar))
if(dd .eq. 1) then
  cx=1
  cy=0
  f1=1
  f2=2
  f3=3
  f4=4
  coef(1)=1.d0
  coef(2)=1.d0
  coef(3)=1.d0
  coef(4)=1.d0
  g1=1.d0
  g2=1.d0
  g3=1.d0
  g4=1.d0
elseif(dd .eq. 2) then
  cx=0
  cy=1
  f1=1
  f2=3
  f3=2
  f4=4
  coef(1)=1.d0
  coef(2)=-1.d0
  coef(3)=1.d0
  coef(4)=1.d0
  g1=1.d0
  g2=1.d0
  g3=-1.d0
  g4=1.d0
endif


if(coordType .eq. 1) then !! Cartesian coordinates
do j=0,ny+1
  do i=0,nx+1
    rhoL=g1*q(i-cx,j-cy,f1)
    rhoM=g1*q(i   ,j   ,f1)
    rhoR=g1*q(i+cx,j+cy,f1)

    vxL=g2*q(i-cx,j-cy,f2)/rhoL
    vxM=g2*q(i   ,j   ,f2)/rhoM
    vxR=g2*q(i+cx,j+cy,f2)/rhoR

    vyL=g3*q(i-cx,j-cy,f3)/rhoL
    vyM=g3*q(i   ,j   ,f3)/rhoM
    vyR=g3*q(i+cx,j+cy,f3)/rhoR

   eneL=g4*q(i-cx,j-cy,f4)/rhoL
   eneM=g4*q(i   ,j   ,f4)/rhoM
   eneR=g4*q(i+cx,j+cy,f4)/rhoR

     pL=(gam-1.d0)*(rhoL*eneL-0.5d0*rhoL*(vxL**2+vyL**2))
     pM=(gam-1.d0)*(rhoM*eneM-0.5d0*rhoM*(vxM**2+vyM**2))
     pR=(gam-1.d0)*(rhoR*eneR-0.5d0*rhoR*(vxR**2+vyR**2))

    SL(i,j,1)=slope(rhoL,rhoM,rhoR)
    SL(i,j,2)=slope( vxL, vxM, vxR)
    SL(i,j,3)=slope( vyL, vyM, vyR)
    SL(i,j,4)=slope(  pL,  pM,  pR)
  enddo !! end i
enddo ! end j
endif !! if coordType = 1

if(coordType .eq. 3 .or. coordType .eq. 2) then !! polar coodinates, uniform in r or log in r
  if(dd .eq. 1) then
    do j=0,ny+1
      do i=0,nx+1
        rhoL=g1*q(i-cx,j-cy,f1)
        rhoM=g1*q(i   ,j   ,f1)
        rhoR=g1*q(i+cx,j+cy,f1)

        vxL=g2*q(i-cx,j-cy,f2)/rhoL  !! vr
        vxM=g2*q(i   ,j   ,f2)/rhoM
        vxR=g2*q(i+cx,j+cy,f2)/rhoR

        vyL=g3*q(i-cx,j-cy,f3)/rhoL  !! vphi
        vyM=g3*q(i   ,j   ,f3)/rhoM
        vyR=g3*q(i+cx,j+cy,f3)/rhoR

       eneL=g4*q(i-cx,j-cy,f4)/rhoL  
       eneM=g4*q(i   ,j   ,f4)/rhoM
       eneR=g4*q(i+cx,j+cy,f4)/rhoR

       pL=(gam-1.d0)*(rhoL*eneL-0.5d0*rhoL*(vxL**2+vyL**2))
       pM=(gam-1.d0)*(rhoM*eneM-0.5d0*rhoM*(vxM**2+vyM**2))
       pR=(gam-1.d0)*(rhoR*eneR-0.5d0*rhoR*(vxR**2+vyR**2))

       ri=this%xc(1)%coords(i-cx)
       rm=this%xc(1)%coords(i   )
       ro=this%xc(1)%coords(i+cx)
       dr=this%dx(1)%coords(i   )

       SL(i,j,1)=slope(rhoL,rhoM,rhoR)
       SL(i,j,2)=slope( vxL, vxM, vxR)
       SL(i,j,3)=slope( vyL, vyM, vyR)
       SL(i,j,4)=slope(  pL,  pM,  pR)

       !SL(i,j,1)=slope3(rhoL,rhoM,rhoR,ri,rm,ro,dr)
       !SL(i,j,2)=slope3( vxL, vxM, vxR,ri,rm,ro,dr)
       !SL(i,j,3)=slope3( vyL, vyM, vyR,ri,rm,ro,dr)
       !SL(i,j,4)=slope3(  pL,  pM,  pR,ri,rm,ro,dr)
      enddo
    enddo
  endif  !! if dd .eq. 1

  if(dd .eq. 2) then
    do j=0,ny+1
      do i=0,nx+1
        rhoL=g1*q(i-cx,j-cy,f1)
        rhoM=g1*q(i   ,j   ,f1)
        rhoR=g1*q(i+cx,j+cy,f1)

        vxL=g2*q(i-cx,j-cy,f2)/rhoL
        vxM=g2*q(i   ,j   ,f2)/rhoM
        vxR=g2*q(i+cx,j+cy,f2)/rhoR

        vyL=g3*q(i-cx,j-cy,f3)/rhoL
        vyM=g3*q(i   ,j   ,f3)/rhoM
        vyR=g3*q(i+cx,j+cy,f3)/rhoR

       eneL=g4*q(i-cx,j-cy,f4)/rhoL
       eneM=g4*q(i   ,j   ,f4)/rhoM
       eneR=g4*q(i+cx,j+cy,f4)/rhoR

         pL=(gam-1.d0)*(rhoL*eneL-0.5d0*rhoL*(vxL**2+vyL**2))
         pM=(gam-1.d0)*(rhoM*eneM-0.5d0*rhoM*(vxM**2+vyM**2))
         pR=(gam-1.d0)*(rhoR*eneR-0.5d0*rhoR*(vxR**2+vyR**2))

        SL(i,j,1)=slope(rhoL,rhoM,rhoR)
        SL(i,j,2)=slope( vxL, vxM, vxR)
        SL(i,j,3)=slope( vyL, vyM, vyR)
        SL(i,j,4)=slope(  pL,  pM,  pR)
      enddo !! end i
    enddo ! end j    
  endif  !! if dd .eq. 2

endif  !! if coordType .eq. 3


do j=0,ny+1
  do i=0,nx+1
    ql(1)=g1*q(i,j,f1)
    ql(2)=g2*q(i,j,f2)/ql(1)
    ql(3)=g3*q(i,j,f3)/ql(1)
    ql(4)=g4*(gam-1.d0)*(q(i,j,f4)-0.5d0*ql(1)*(ql(2)**2+ql(3)**2))

    qr(1)=g1*q(i+cx,j+cy,f1)
    qr(2)=g2*q(i+cx,j+cy,f2)/qr(1)
    qr(3)=g3*q(i+cx,j+cy,f3)/qr(1)
    qr(4) =g4*(gam-1.d0)*(q(i+cx,j+cy,f4)-0.5d0*qr(1)*(qr(2)**2+qr(3)**2))

    slopeL(:)=SL(i,j,:)
    slopeR(:)=SL(i+cx,j+cy,:)

    call fluxPtr(ql,qr,slopeL,slopeR,gam,flux,this%nvar)
    Fx(i,j,:)=flux(:)

  enddo ! end i
enddo ! end j

if(coordType .eq. 1) then  
  do j=1,ny
    do i=1,nx
      dtddx=this%dt/(this%dx(1)%coords(i)*cx+this%dx(2)%coords(j)*cy)
      q2(i,j,1)=q1(i,j,1)+coef(1)*(Fx(i-cx,j-cy,f1)-Fx(i,j,f1))*dtddx
      q2(i,j,2)=q1(i,j,2)+coef(2)*(Fx(i-cx,j-cy,f2)-Fx(i,j,f2))*dtddx
      q2(i,j,3)=q1(i,j,3)+coef(3)*(Fx(i-cx,j-cy,f3)-Fx(i,j,f3))*dtddx
      q2(i,j,4)=q1(i,j,4)+coef(4)*(Fx(i-cx,j-cy,f4)-Fx(i,j,f4))*dtddx
    enddo
  enddo
endif

if(coordType .eq. 3 .or. coordType .eq. 2) then
  if(dd .eq. 1) then  !! radial direction
    do j=1,ny
      do i=1,nx
        rm=this%xc(1)%coords(i)
        dr=this%dx(1)%coords(i)
        fac1=rm-0.5d0*dr
        fac2=rm+0.5d0*dr
        dtddx=this%dt/(rm*dr)
        den=q(i,j,1)
         vr=q(i,j,2)/den
       vphi=q(i,j,3)/den
        ene=q(i,j,4)
        ptot= (ene-0.5d0*den*(vr**2+vphi**2))*(gam-1.d0)

        q2(i,j,1)=q1(i,j,1)+coef(1)*(Fx(i-cx,j-cy,f1)*fac1-Fx(i,j,f1)*fac2)*dtddx
        q2(i,j,2)=q1(i,j,2)+coef(2)*(Fx(i-cx,j-cy,f2)*fac1-Fx(i,j,f2)*fac2)*dtddx+(den*vphi**2+ptot)/rm*this%dt
        q2(i,j,3)=q1(i,j,3)+coef(3)*(Fx(i-cx,j-cy,f3)*fac1**2-Fx(i,j,f3)*fac2**2)*this%dt/(dr*rm**2)
        !q2(i,j,3)=q1(i,j,3)+coef(3)*(Fx(i-cx,j-cy,f3)*fac1-Fx(i,j,f3)*fac2)*dtddx-(den*vr*vphi)/rm*this%dt
        q2(i,j,4)=q1(i,j,4)+coef(4)*(Fx(i-cx,j-cy,f4)*fac1-Fx(i,j,f4)*fac2)*dtddx
      enddo
    enddo 
  endif  !!! dd .eq. 1

  if(dd .eq. 2) then  !! phi direction
    do j=1,ny
      do i=1,nx
        rm=this%xc(1)%coords(i)
        dphi=this%dx(2)%coords(j)
        dtddx=this%dt/(dphi*rm)
        q2(i,j,1)=q1(i,j,1)+coef(1)*(Fx(i-cx,j-cy,f1)-Fx(i,j,f1))*dtddx
        q2(i,j,2)=q1(i,j,2)+coef(2)*(Fx(i-cx,j-cy,f2)-Fx(i,j,f2))*dtddx 
        q2(i,j,3)=q1(i,j,3)+coef(3)*(Fx(i-cx,j-cy,f3)-Fx(i,j,f3))*dtddx
        q2(i,j,4)=q1(i,j,4)+coef(4)*(Fx(i-cx,j-cy,f4)-Fx(i,j,f4))*dtddx
      enddo
    enddo    
  endif  !!! dd .eq. 2
endif !!! coordType .eq. 3

deallocate(SL,Fx)
end subroutine solverAdi2D

subroutine fluxHLLCAdiHD2D(ql,qr,slopeL,slopeR,gam,flux,nvar)
implicit none
integer::nvar
double precision::flux(nvar),ql(nvar),qr(nvar),slopeL(nvar),slopeR(nvar), &
                  FL(nvar),FR(nvar),FsL(nvar), FsR(nvar),Fs(nvar)
double precision::rhoL,rhoR,vxL,vxR,seneL,seneR,eneL,eneR,pTL,pTR,pxL,pxR, &
             spdL,spdR,SL,SR,SM,pstarL,pstarR,rhosL,rhosR,pxsL,pxsR,enesL,enesR, &
             vyL,vyR,pyL,pyR,pysL,pysR,pL,pR,gam

rhoL =ql(1)+0.5d0*slopeL(1)
vxL  =ql(2)+0.5d0*slopeL(2)
vyL  =ql(3)+0.5d0*slopeL(3)
pL   =ql(4)+0.5d0*slopeL(4)
eneL =pL/(gam-1.d0)+0.5d0*rhoL*(vxL**2+vyL**2)
!seneL=ql(4)+0.5d0*slopeL(4)
!eneL =seneL*rhoL

rhoR =qr(1)-0.5d0*slopeR(1)
vxR  =qr(2)-0.5d0*slopeR(2)
vyR  =qr(3)-0.5d0*slopeR(3)
pR   =qr(4)-0.5d0*slopeR(4)
eneR =pR/(gam-1.d0)+0.5d0*rhoR*(vxR**2+vyR**2)
!seneR=qr(4)-0.5d0*slopeR(4)
!eneR =seneR*rhoR

pTL = (gam-1.d0)*(eneL-0.5d0*rhoL*(vxL**2+vyL**2))
pTR = (gam-1.d0)*(eneR-0.5d0*rhoR*(vxR**2+vyR**2))

!if(pTL .le. 0.d0 .or. pTR .le. 0.d0) then
if(pL .le. 0.d0 .or. pR .le. 0.d0) then
rhoL =ql(1)
vxL  =ql(2)
vyL  =ql(3)
pL   =ql(4)
eneL =pL/(gam-1.d0)+0.5d0*rhoL*(vxL**2+vyL**2)
!seneL=ql(4)
!eneL =seneL*rhoL
rhoR =qr(1)
vxR  =qr(2)
vyR  =qr(3)
pR   =qr(4)
eneR =pR/(gam-1.d0)+0.5d0*rhoR*(vxR**2+vyR**2)
!seneR=qr(4)
!eneR =seneR*rhoR

pTL = (gam-1.d0)*(eneL-0.5d0*rhoL*(vxL**2+vyL**2))
pTR = (gam-1.d0)*(eneR-0.5d0*rhoR*(vxR**2+vyR**2))
!write(*,*) "negative pressure!"
endif

spdL = dsqrt(gam*pTL/rhoL)
pxL = vxL*rhoL
pyL = vyL*rhoL

spdR = dsqrt(gam*pTR/rhoR)
pxR = vxR*rhoR
pyR = vyR*rhoR

SL = dmin1(vxL-spdL,vxR-spdR)
SR = dmax1(vxL+spdL,vxR+spdR)
SM = ((SR-vxR)*rhoR*vxR-(SL-vxL)*rhoL*vxL-pTR+pTL)  &
        /((SR-vxR)*rhoR-(SL-vxL)*rhoL)
pstarL = rhoL*(vxL-SL)*(vxL-SM)+pTL
pstarR = rhoR*(vxR-SR)*(vxR-SM)+pTR

rhosL = rhoL*(SL-vxL)/(SL-SM)
pxsL = ((SL-vxL)*pxL+(pstarL-pTL))/(SL-SM)
pysL = rhosL*vyL
enesL = ((SL-vxL)*eneL-pTL*vxL+pstarL*SM)/(SL-SM)

rhosR = rhoR*(SR-vxR)/(SR-SM)
pxsR = ((SR-vxR)*pxR+(pstarR-pTR))/(SR-SM)
pysR = rhosR*vyR
enesR = ((SR-vxR)*eneR-pTR*vxR+pstarR*SM)/(SR-SM)

FL(1) = rhoL*vxL
FL(2) = rhoL*vxL**2.d0 +pTL
FL(3) = rhoL*vxL*vyL
FL(4) = (eneL+pTL)*vxL

FR(1) = rhoR*vxR
FR(2) = rhoR*vxR**2.d0 +pTR
FR(3) = rhoR*vxR*vyR
FR(4) = (eneR+pTR)*vxR

FsL(1)= FL(1) + SL*(rhosL-rhoL)
FsL(2)= FL(2) + SL*(pxsL-pxL)
FsL(3)= FL(3) + SL*(pysL-pyL)
FsL(4)= FL(4) + SL*(enesL-eneL)

FsR(1)= FR(1) + SR*(rhosR-rhoR)
FsR(2)= FR(2) + SR*(pxsR-pxR)
FsR(3)= FR(3) + SR*(pysR-pyR)
FsR(4)= FR(4) + SR*(enesR-eneR)

IF (SL .gt. 0.d0)then
       flux=FL
ELSEIF (SL .le. 0.d0 .and. SM .ge. 0.d0)then
       flux=FsL
ELSEIF (SM .le. 0.d0 .and. SR .ge. 0.d0)then
       flux=FsR
ELSEIF (SR .lt. 0.d0) then
       flux=FR
ENDIF
end subroutine fluxHLLCAdiHD2D

subroutine fluxHLLAdiHD2D(ql,qr,slopeL,slopeR,gam,flux,nvar)
implicit none
integer::nvar
double precision::ql(NVAR),qr(NVAR),slopeL(NVAR),slopeR(NVAR),flux(NVAR)
double precision, dimension(:), allocatable:: FL,FR,Fs,UL,UR
double precision:: spd,gam ! sound speed
double precision:: rhoL,vxL,vyL,eneL,spdL,pTL,pL
double precision:: rhoR,vxR,vyR,eneR,spdR,pTR,pR
double precision:: SL,SR

allocate(FL(NVAR))
allocate(FR(NVAR))
allocate(Fs(NVAR))
allocate(UL(NVAR))
allocate(UR(NVAR))

UL=ql+0.5d0*slopeL
UR=qr-0.5d0*slopeR

rhoL=UL(1)
vxL =UL(2)
vyL =UL(3)
pL  =UL(4)
eneL = (pL/(gam-1.d0)+0.5d0*rhoL*(vxL**2+vyL**2))/rhoL
pTL = rhoL*(eneL-0.5d0*(vxL**2.d0+vyL**2.d0))*(gam-1.d0)
spdL = dsqrt(gam*pTL/rhoL)

rhoR=UR(1)
vxR =UR(2)
vyR =UR(3)
pR  =UR(4)
eneR = (pR/(gam-1.d0)+0.5d0*rhoR*(vxR**2+vyR**2))/rhoR
pTR = rhoR*(eneR-0.5d0*(vxR**2.d0+vyR**2.d0))*(gam-1.d0)
spdR = dsqrt(gam*pTR/rhoR)

SL=dmin1(vxL-spdL,vxR-spdR)
SR=dmax1(vxL+spdL,vxR+spdR)

FL(1)=rhoL*vxL
FL(2)=rhoL*vxL**2.d0+pTL
FL(3)=rhoL*vxL*vyL
FL(4)=(rhoL*eneL+pTL)*vxL

FR(1)=rhoR*vxR
FR(2)=rhoR*vxR**2.d0+pTR
FR(3)=rhoR*vxR*vyR
FR(4)=(rhoR*eneR+pTR)*vxR

Fs(1)=(SR*FL(1)-SL*FR(1)+SR*SL*(UR(1)-UL(1)))/(SR-SL)
Fs(2)=(SR*FL(2)-SL*FR(2)+SR*SL*(UR(1)*UR(2)-UL(1)*UL(2)))/(SR-SL)
Fs(3)=(SR*FL(3)-SL*FR(3)+SR*SL*(UR(1)*UR(3)-UL(1)*UL(3)))/(SR-SL)
Fs(4)=(SR*FL(4)-SL*FR(4)+SR*SL*(UR(1)*UR(4)-UL(1)*UL(4)))/(SR-SL)

if(SL .gt. 0.d0) then
   flux=FL
elseif(SL .le. 0.d0 .and. SR .ge. 0.d0) then
   flux=Fs
else
   flux=FR
endif


deallocate(FL)
deallocate(FR)
deallocate(Fs)
deallocate(UL)
deallocate(UR)
end subroutine fluxHLLAdiHD2D

subroutine fluxHLLIsoHD2D(ql,qr,slopeL,slopeR,snd,flux,nvar)
implicit none
integer::nvar
double precision::ql(NVAR),qr(NVAR),slopeL(NVAR),slopeR(NVAR),flux(NVAR)
double precision, dimension(:), allocatable:: FL,FR,Fs,UL,UR
double precision:: spd,snd ! sound speed
double precision:: rhoL,vxL,vyL,pTL
double precision:: rhoR,vxR,vyR,pTR
double precision:: SL,SR

allocate(FL(NVAR))
allocate(FR(NVAR))
allocate(Fs(NVAR))
allocate(UL(NVAR))
allocate(UR(NVAR))

UL=ql+0.5d0*slopeL
UR=qr-0.5d0*slopeR

spd=snd
rhoL=UL(1)
vxL =UL(2)
vyL =UL(3)
pTL =rhoL*spd**2.d0

rhoR=UR(1)
vxR =UR(2)
vyR =UR(3)
pTR =rhoR*spd**2.d0

SL=dmin1(vxL-spd,vxR-spd)
SR=dmax1(vxL+spd,vxR+spd)

FL(1)=rhoL*vxL
FL(2)=rhoL*vxL**2.d0+pTL
FL(3)=rhoL*vxL*vyL

FR(1)=rhoR*vxR
FR(2)=rhoR*vxR**2.d0+pTR
FR(3)=rhoR*vxR*vyR

Fs(1)=(SR*FL(1)-SL*FR(1)+SR*SL*(UR(1)-UL(1)))/(SR-SL)
Fs(2)=(SR*FL(2)-SL*FR(2)+SR*SL*(UR(1)*UR(2)-UL(1)*UL(2)))/(SR-SL)
Fs(3)=(SR*FL(3)-SL*FR(3)+SR*SL*(UR(1)*UR(3)-UL(1)*UL(3)))/(SR-SL)

if(SL .gt. 0.d0) then
   flux=FL
elseif(SL .le. 0.d0 .and. SR .ge. 0.d0) then
   flux=Fs
else
   flux=FR
endif


deallocate(FL)
deallocate(FR)
deallocate(Fs)
deallocate(UL)
deallocate(UR)
end subroutine fluxHLLIsoHD2D

subroutine fluxHLLPolyHD2D(ql,qr,slopeL,slopeR,polyK,polyGamma,flux,nvar)
implicit none
integer::nvar
double precision::polyK, polyGamma
double precision::ql(NVAR),qr(NVAR),slopeL(NVAR),slopeR(NVAR),flux(NVAR)
double precision, dimension(:), allocatable:: FL,FR,Fs,UL,UR
double precision:: spdL,spdR ! sound speed
double precision:: rhoL,vxL,vyL,pTL
double precision:: rhoR,vxR,vyR,pTR
double precision:: SL,SR

allocate(FL(NVAR))
allocate(FR(NVAR))
allocate(Fs(NVAR))
allocate(UL(NVAR))
allocate(UR(NVAR))

UL=ql+0.5d0*slopeL
UR=qr-0.5d0*slopeR

rhoL=UL(1)
vxL =UL(2)
vyL =UL(3)
pTL =polyK*rhoL**polyGamma
spdL=dsqrt(polyGamma*pTL/rhoL)

rhoR=UR(1)
vxR =UR(2)
vyR =UR(3)
pTR =polyK*rhoR**polyGamma
spdR=dsqrt(polyGamma*pTR/rhoR)

SL=dmin1(vxL-spdL,vxR-spdR)
SR=dmax1(vxL+spdL,vxR+spdR)

FL(1)=rhoL*vxL
FL(2)=rhoL*vxL**2.d0+pTL
FL(3)=rhoL*vxL*vyL

FR(1)=rhoR*vxR
FR(2)=rhoR*vxR**2.d0+pTR
FR(3)=rhoR*vxR*vyR

Fs(1)=(SR*FL(1)-SL*FR(1)+SR*SL*(UR(1)-UL(1)))/(SR-SL)
Fs(2)=(SR*FL(2)-SL*FR(2)+SR*SL*(UR(1)*UR(2)-UL(1)*UL(2)))/(SR-SL)
Fs(3)=(SR*FL(3)-SL*FR(3)+SR*SL*(UR(1)*UR(3)-UL(1)*UL(3)))/(SR-SL)

if(SL .gt. 0.d0) then
   flux=FL
elseif(SL .le. 0.d0 .and. SR .ge. 0.d0) then
   flux=Fs
else
   flux=FR
endif


deallocate(FL)
deallocate(FR)
deallocate(Fs)
deallocate(UL)
deallocate(UR)
end subroutine fluxHLLPolyHD2D


subroutine fluxExactIsoHD2D(ql,qr,slopeL,slopeR,snd,flux,nvar)
implicit none
integer::nvar
double precision::ql(nvar),qr(nvar),qm(nvar),slopeL(nvar),slopeR(nvar),flux(nvar),snd
double precision:: rhoL,rhoR,vxL,vxR,vyL,vyR
double precision::c,shkR,rarR,csq,vstar,rhostar,rho0,rhosR,srhosR,aisrhosR,sigma2, &
                   rhosL,srhosL,aisrhosL,sigma1,aa,bb,uu,srhostar
integer:: Intr,n
double precision:: tolerr = 1.d-10

ql(1) = ql(1)+0.5d0*slopeL(1)
ql(2) = ql(2)+0.5d0*slopeL(2)
ql(3) = ql(3)+0.5d0*slopeL(3)

qr(1) = qr(1)-0.5d0*slopeR(1)
qr(2) = qr(2)-0.5d0*slopeR(2)
qr(3) = qr(3)-0.5d0*slopeR(3)

c = snd
Intr = 10

if (ql(1).eq.qr(1).and.ql(2).eq.qr(2)) then
     qm(1) = qr(1)
     qm(2) = qr(2)
if(ql(2).ge. 0.   ) then
     qm(3) = ql(3)
else
     qm(3) = qr(3)
endif
goto 1000
endif

      rhoL = ql(1)
      rhoR = qr(1)
       vxL = ql(2)
       vxR = qr(2)
       vyL = ql(3)
       vyR = qr(3)


      shkR = -c*(rhoR-rhoL)/dsqrt(rhoR*rhoL)
      rarR = -c*dlog(rhoR/rhoL)
       csq =  c*c

      if(rhoL .ge. rhoR) then
        if ((vxR-vxL) .gt. rarR) then
          goto 100
        elseif ((vxR-vxL) .gt. -shkR) then
          goto 200
        else
          goto 400
        endif
      else
        if ((vxR-vxL) .lt. shkR ) then
          goto 400
        elseif ((vxR-vxL) .lt. -rarR) then
          goto 300
        else
          goto 100
        endif
       endif
!cccccc 1-rarefaction wave, 2-rarefaction wave
!cccccc 1-rarefaction wave, 2-rarefaction wave
!cccccc 1-rarefaction wave, 2-rarefaction wave
100   vstar   = 0.5d0*(vxL+vxR)+0.5d0*c*dlog(rhoL/rhoR)
      rhostar = rhoL*dexp(-(vstar-vxL)/c)

      if ((vxL-c) .ge. 0.) then
         qm(2) = vxL
         qm(1) = rhoL
      elseif ((vstar-c) .ge. 0.) then
         qm(2) =  c
         qm(1) =  rhoL*dexp((vxL-c)/c)
      elseif ((vstar+c) .ge. 0.) then
         qm(2) =  vstar
         qm(1) =  rhostar
      elseif ((  vxR+c) .ge. 0.) then
         qm(2) =  -c
         qm(1) =  rhoR*dexp(-(vxR+c)/c)
      else
         qm(2) =   vxR
         qm(1) =  rhoR
      endif
      if ( qm(2) .ge. 0.) then
         qm(3) =  ql(3)
      else
         qm(3) =  qr(3)
      endif
      goto 1000
!cccccc1-rarefaction wave, 2-shock wave
!cccccc1-rarefaction wave, 2-shock wave
!cccccc1-rarefaction wave, 2-shock wave
200   rho0 = rhoR
      do n = 1, INtr
        rhosR    = rho0/rhoR
        srhosR   = dsqrt(rhosR)
        aisrhosR = 1.d0/srhosR
        rhostar  = rho0 - &
       (vxR-vxL+c*(srhosR-aisrhosR)+c*dlog(rho0/rhoL)) &
       /(0.5d0*c/rho0*(srhosR+aisrhosR)+c/rho0)
        if (dabs((rho0-rhostar)/rho0) .lt. tolerr) then
            goto 210
        endif
        rho0 = rhostar
      enddo
      write(6,*) 'Newton Interation is divergent at x-dir 1R2S'
      write(6,*) 'rhostar,rhoL=',rhostar,rhoL
      write(6,*) 'uxL    , uxR=',vxL,      vxR
      stop
210   vstar  = vxL - c*dlog(rhostar/rhoL)
      sigma2 = vxR + c*dsqrt(rhostar/rhoR)
      if ( sigma2 .le. 0.) then
         qm(2) =  vxR
         qm(1) = rhoR
      elseif ((vstar-c) .le. 0.) then
         qm(2) = vstar
         qm(1) = rhostar
      elseif ((vxL-c) .le. 0.) then
         qm(2) =  c
         qm(1) =  rhoL*dexp((vxL-c)/c)
      else
         qm(2) = vxL
         qm(1) = rhoL
      endif
      if ( qm(2) .ge. 0.) then
         qm(3) = ql(3)
      else
         qm(3) = qr(3)
      endif
      goto 1000
!cccccc 1-shock wave, 2-rarefaction wave
!cccccc 1-shock wave, 2-rarefaction wave
!cccccc 1-shock wave, 2-rarefaction wave
300   rho0 = rhoL
      do n = 1, INtr
        rhosL    = rho0/rhoL
        srhosL   = dsqrt(rhosL)
        aisrhosL = 1.d0/srhosL
        rhostar  = rho0- &
      (vxR-vxL+c*(srhosL-aisrhosL)+c*dlog(rho0/rhoR)) &
       /(0.5d0*c/rho0*(srhosL+aisrhosL)+c/rho0)
        if (dabs((rho0-rhostar)/rho0) .lt. tolerr) then
            goto 310
        endif
        rho0 = rhostar
      enddo
      write(6,*) 'Newton Interation is divergent at x-dir 1S2R'
      write(6,*) 'rhostar,rhoL=',rhostar,rhoL
      write(6,*) 'uxL    , uxR=',vxL,      vxR
      stop
310   vstar  = vxR + c*dlog(rhostar/rhoR)
      sigma1 = vxL - c*dsqrt(rhostar/rhoL)
      if (sigma1 .ge. 0.) then
         qm(2) = vxL
         qm(1) = rhoL
      elseif ((vstar+c) .ge. 0.) then
         qm(2) = vstar
         qm(1) = rhostar
      elseif ((vxR+c) .ge. 0.) then
         qm(2) =  -c
         qm(1) =  rhoR*dexp(-(vxR+c)/c)
      else
         qm(2) = vxR
         qm(1) = rhoR
      endif
      if ( qm(2) .ge. 0.) then
         qm(3) = ql(3)
      else
         qm(3) = qr(3)
      endif
      goto 1000
!       1-shock wave, 2-shock wave
!       1-shock wave, 2-shock wave
!       1-shock wave, 2-shock wave
400   aa       = 1.d0/dsqrt(rhoL)+1.d0/dsqrt(rhoR)
      bb       =    dsqrt(rhoL)+   dsqrt(rhoR)
      uu       = vxR-vxL
      srhostar = (-uu+dsqrt(uu**2+4.d0*csq*aa*bb))/(2.d0*c*aa)
      rhostar  = srhostar*srhostar
      vstar    = vxL - c*(dsqrt(rhostar/rhoL)-dsqrt(rhoL/rhostar))
      sigma1   = vxL - c*dsqrt(rhostar/rhoL)
      sigma2   = vxR + c*dsqrt(rhostar/rhoR)
      if     (sigma1 .ge. 0.) then
       qm(2) = vxL
       qm(1) = rhoL
      elseif (sigma2 .ge. 0.) then
       qm(2) = vstar
       qm(1) = rhostar
      else
       qm(2) = vxR
       qm(1) = rhoR
      endif
         if  (qm(2) .ge. 0.) then
       qm(3) = ql(3)
      else
       qm(3) = qr(3)
      endif
      goto 1000

1000  flux(1) = qm(1)*(   qm(2)          )
      flux(2) = qm(1)*( qm(2)*qm(2) + c*c)
      flux(3) = qm(1)*( qm(2)*qm(3)      )

end subroutine fluxExactIsoHD2D


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!! Riemann solvers for 1D problems
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine fluxHLLCAdiHD1D(ql,qr,slopeL,slopeR,gam,flux,nvar)
integer::nvar
double precision::flux(nvar),ql(nvar),qr(nvar),slopeL(nvar),slopeR(nvar), &
                  FL(nvar),FR(nvar),FsL(nvar), FsR(nvar),Fs(nvar)
double precision::rhoL,rhoR,vxL,vxR,seneL,seneR,eneL,eneR,pTL,pTR,pxL,pxR, &
             spdL,spdR,SL,SR,SM,pstarL,pstarR,rhosL,rhosR,pxsL,pxsR,enesL,enesR, &
             vyL,vyR,pyL,pyR,pysL,pysR,pL,pR,gam


rhoL =ql(1)+0.5d0*slopeL(1)
vxL  =ql(2)+0.5d0*slopeL(2)
pL   =ql(3)+0.5d0*slopeL(3)
eneL =pL/(gam-1.d0)+0.5d0*rhoL*(vxL**2+vyL**2)
!seneL=ql(4)+0.5d0*slopeL(4)
!eneL =seneL*rhoL

rhoR =qr(1)-0.5d0*slopeR(1)
vxR  =qr(2)-0.5d0*slopeR(2)
pR   =qr(3)-0.5d0*slopeR(3)
eneR =pR/(gam-1.d0)+0.5d0*rhoR*(vxR**2+vyR**2)
!seneR=qr(4)-0.5d0*slopeR(4)
!eneR =seneR*rhoR

pTL = (gam-1.d0)*(eneL-0.5d0*rhoL*(vxL**2+vyL**2))
pTR = (gam-1.d0)*(eneR-0.5d0*rhoR*(vxR**2+vyR**2))

!if(pTL .le. 0.d0 .or. pTR .le. 0.d0) then
if(pL .le. 0.d0 .or. pR .le. 0.d0) then
rhoL =ql(1)
vxL  =ql(2)
pL   =ql(3)
eneL =pL/(gam-1.d0)+0.5d0*rhoL*(vxL**2+vyL**2)
!seneL=ql(4)
!eneL =seneL*rhoL

rhoR =qr(1)
vxR  =qr(2)
pR   =qr(3)
eneR =pR/(gam-1.d0)+0.5d0*rhoR*(vxR**2+vyR**2)
!seneR=qr(4)
!eneR =seneR*rhoR

pTL = (gam-1.d0)*(eneL-0.5d0*rhoL*(vxL**2+vyL**2))
pTR = (gam-1.d0)*(eneR-0.5d0*rhoR*(vxR**2+vyR**2))
!write(*,*) "negative pressure!"
endif

spdL = dsqrt(gam*pTL/rhoL)
pxL = vxL*rhoL
pyL = vyL*rhoL

spdR = dsqrt(gam*pTR/rhoR)
pxR = vxR*rhoR
pyR = vyR*rhoR

SL = dmin1(vxL-spdL,vxR-spdR)
SR = dmax1(vxL+spdL,vxR+spdR)
SM = ((SR-vxR)*rhoR*vxR-(SL-vxL)*rhoL*vxL-pTR+pTL)  &
        /((SR-vxR)*rhoR-(SL-vxL)*rhoL)
pstarL = rhoL*(vxL-SL)*(vxL-SM)+pTL
pstarR = rhoR*(vxR-SR)*(vxR-SM)+pTR

rhosL = rhoL*(SL-vxL)/(SL-SM)
pxsL = ((SL-vxL)*pxL+(pstarL-pTL))/(SL-SM)
pysL = rhosL*vyL
enesL = ((SL-vxL)*eneL-pTL*vxL+pstarL*SM)/(SL-SM)

rhosR = rhoR*(SR-vxR)/(SR-SM)
pxsR = ((SR-vxR)*pxR+(pstarR-pTR))/(SR-SM)
pysR = rhosR*vyR
enesR = ((SR-vxR)*eneR-pTR*vxR+pstarR*SM)/(SR-SM)

FL(1) = rhoL*vxL
FL(2) = rhoL*vxL**2.d0 +pTL
FL(3) = (eneL+pTL)*vxL

FR(1) = rhoR*vxR
FR(2) = rhoR*vxR**2.d0 +pTR
FR(3) = (eneR+pTR)*vxR

FsL(1)= FL(1) + SL*(rhosL-rhoL)
FsL(2)= FL(2) + SL*(pxsL-pxL)
FsL(3)= FL(3) + SL*(enesL-eneL)

FsR(1)= FR(1) + SR*(rhosR-rhoR)
FsR(2)= FR(2) + SR*(pxsR-pxR)
FsR(3)= FR(3) + SR*(enesR-eneR)

IF (SL .gt. 0.d0)then
       flux=FL
ELSEIF (SL .le. 0.d0 .and. SM .ge. 0.d0)then
       flux=FsL
ELSEIF (SM .le. 0.d0 .and. SR .ge. 0.d0)then
       flux=FsR
ELSEIF (SR .lt. 0.d0) then
       flux=FR
ENDIF
end subroutine fluxHLLCAdiHD1D



subroutine fluxHLLAdiHD1D(ql,qr,slopeL,slopeR,gam,flux,nvar)
integer::nvar
double precision::ql(NVAR),qr(NVAR),slopeL(NVAR),slopeR(NVAR),flux(NVAR)
double precision, dimension(:), allocatable:: FL,FR,Fs,UL,UR
double precision:: spd,gam ! sound speed
double precision:: rhoL,vxL,vyL,eneL,spdL,pTL,pL
double precision:: rhoR,vxR,vyR,eneR,spdR,pTR,pR
double precision:: SL,SR

allocate(FL(NVAR))
allocate(FR(NVAR))
allocate(Fs(NVAR))
allocate(UL(NVAR))
allocate(UR(NVAR))

UL=ql+0.5d0*slopeL
UR=qr-0.5d0*slopeR

rhoL=UL(1)
vxL =UL(2)
pL  =UL(3)
eneL = (pL/(gam-1.d0)+0.5d0*rhoL*(vxL**2))/rhoL
pTL = rhoL*(eneL-0.5d0*(vxL**2.d0))*(gam-1.d0)
spdL = dsqrt(gam*pTL/rhoL)

rhoR=UR(1)
vxR =UR(2)
pR  =UR(3)
eneR = (pR/(gam-1.d0)+0.5d0*rhoR*(vxR**2))/rhoR
pTR = rhoR*(eneR-0.5d0*(vxR**2.d0))*(gam-1.d0)
spdR = dsqrt(gam*pTR/rhoR)

SL=dmin1(vxL-spdL,vxR-spdR)
SR=dmax1(vxL+spdL,vxR+spdR)

FL(1)=rhoL*vxL
FL(2)=rhoL*vxL**2.d0+pTL
FL(3)=(rhoL*eneL+pTL)*vxL

FR(1)=rhoR*vxR
FR(2)=rhoR*vxR**2.d0+pTR
FR(3)=(rhoR*eneR+pTR)*vxR

Fs(1)=(SR*FL(1)-SL*FR(1)+SR*SL*(UR(1)-UL(1)))/(SR-SL)
Fs(2)=(SR*FL(2)-SL*FR(2)+SR*SL*(UR(1)*UR(2)-UL(1)*UL(2)))/(SR-SL)
Fs(3)=(SR*FL(3)-SL*FR(3)+SR*SL*(UR(1)*UR(3)-UL(1)*UL(3)))/(SR-SL)

if(SL .gt. 0.d0) then
   flux=FL
elseif(SL .le. 0.d0 .and. SR .ge. 0.d0) then
   flux=Fs
else
   flux=FR
endif


deallocate(FL)
deallocate(FR)
deallocate(Fs)
deallocate(UL)
deallocate(UR)
end subroutine fluxHLLAdiHD1D

subroutine solverAdiMHD1D(this,q,q1,q2,dd)
class(grid)::this
double precision,dimension(1-this%nbuf:this%nMesh(1)+this%nbuf,this%nvar)::q,q1,q2
integer::nx,nbuf,nvar,dd
integer::i
double precision::dt,dx,dtddx,gam
double precision::rhoL,rhoM,rhoR
double precision::vxL,vxM,vxR
double precision::vyL,vyM,vyR
double precision::vzL,vzM,vzR
double precision::bxL,bxM,bxR,bxtmp
double precision::byL,byM,byR
double precision::bzL,bzM,bzR
double precision::eneL,eneM,eneR
double precision::pL ,pM ,pR
double precision,dimension(:,:),allocatable::SL
double precision,dimension(:,:),allocatable::Fx
double precision::ql(this%nvar),qr(this%nvar),slopeL(this%nvar),slopeR(this%nvar),flux(this%nvar)
procedure(limiter),pointer::slope=>null()
procedure(fluxSolver),pointer::fluxPtr=>null()

nx=this%nMesh(1)
nbuf=this%nbuf
nvar=this%nvar
dt=this%dt
gam=this%adiGamma
select case(this%solverType)
   case (4)
     fluxPtr=>fluxHLLAdiMHD1D
   case (5)
     fluxPtr=>fluxHLLDAdiMHD1D
   case default
     print *,'no appropriate solver found!'
end select

select case (this%limiterType)
   case(0)
     slope=>zslop
   case(1)
     slope=>vslop
   case(2)
     slope=>fslop
   case(3)
     slope=>minmod
end select

allocate(SL(1-nbuf:nx+nbuf,nvar))
allocate(Fx(1-nbuf:nx+nbuf,nvar))

do i=0,nx+1
   rhoL=q(i-1,1)
   rhoM=q(i  ,1)
   rhoR=q(i+1,1)

    vxL=q(i-1,2)/rhoL
    vxM=q(i  ,2)/rhoM
    vxR=q(i+1,2)/rhoR

    vyL=q(i-1,3)/rhoL
    vyM=q(i  ,3)/rhoM
    vyR=q(i+1,3)/rhoR

    vzL=q(i-1,4)/rhoL
    vzM=q(i  ,4)/rhoM
    vzR=q(i+1,4)/rhoR

    bxL=0.5d0*(q(i-1,5)+q(i-1, 9))
    bxM=0.5d0*(q(i  ,5)+q(i  , 9))    
    bxR=0.5d0*(q(i+1,5)+q(i+1, 9))
 
    byL=0.5d0*(q(i-1,6)+q(i-1,10))
    byM=0.5d0*(q(i  ,6)+q(i  ,10))
    byR=0.5d0*(q(i+1,6)+q(i+1,10))

    bzL=0.5d0*(q(i-1,7)+q(i-1,11))
    bzM=0.5d0*(q(i  ,7)+q(i  ,11))
    bzR=0.5d0*(q(i+1,7)+q(i+1,11))

   eneL=q(i-1,8)/rhoL
   eneM=q(i  ,8)/rhoM
   eneR=q(i+1,8)/rhoR

     pL=(gam-1.d0)*(rhoL*eneL-0.5d0*rhoL*(vxL**2+vyL**2+vzL**2)-0.5d0*(bxL**2+byL**2+bzL**2))
     pM=(gam-1.d0)*(rhoM*eneM-0.5d0*rhoM*(vxM**2+vyM**2+vzM**2)-0.5d0*(bxM**2+byM**2+bzM**2))
     pR=(gam-1.d0)*(rhoR*eneR-0.5d0*rhoR*(vxR**2+vyR**2+vzR**2)-0.5d0*(bxR**2+byR**2+bzR**2))

   SL(i,1)=slope(rhoL,rhoM,rhoR)
   SL(i,2)=slope( vxL, vxM, vxR)
   SL(i,3)=slope( vyL, vyM, vyR)
   SL(i,4)=slope( vzL, vzM, vzR)
   SL(i,5)=0.d0
   SL(i,6)=slope( byL, byM, byR)
   SL(i,7)=slope( bzL, bzM, bzR)
   SL(i,8)=slope(  pL,  pM,  pR)
enddo ! end i

do i=0,nx+1
   ql(1)=q(i,1)
   ql(2)=q(i,2)/ql(1)
   ql(3)=q(i,3)/ql(1)
   ql(4)=q(i,4)/ql(1)
   ql(5)=q(i,9)
   bxtmp=0.5d0*(q(i,5)+q(i, 9))
   ql(6)=0.5d0*(q(i,6)+q(i,10)) 
   ql(7)=0.5d0*(q(i,7)+q(i,11))
   ql(8)=(gam-1.d0)*(q(i,8)-0.5d0*ql(1)*(ql(2)**2+ql(3)**2+ql(4)**2)-0.5d0*(bxtmp**2+ql(6)**2+ql(7)**2))

   qr(1)=q(i+1,1)
   qr(2)=q(i+1,2)/qr(1)
   qr(3)=q(i+1,3)/qr(1)
   qr(4)=q(i+1,4)/qr(1)
   qr(5)=q(i+1,5)
   bxtmp=0.5d0*(q(i+1,5)+q(i+1, 9))
   qr(6)=0.5d0*(q(i+1,6)+q(i+1,10))
   qr(7)=0.5d0*(q(i+1,7)+q(i+1,11))
   qr(8)=(gam-1.d0)*(q(i+1,8)-0.5d0*qr(1)*(qr(2)**2+qr(3)**2+qr(4)**2)-0.5d0*(bxtmp**2+qr(6)**2+qr(7)**2))

   slopeL(:)=SL(i,:)
   slopeR(:)=SL(i+1,:)
   call fluxPtr(ql,qr,slopeL,slopeR,this%adiGamma,flux,this%nvar)
   Fx(i,:)=flux(:)
enddo

do i=1,nx
   dtddx=dt/this%dx(1)%coords(i)
   q2(i,1)=q1(i,1)+(Fx(i-1,1)-Fx(i,1))*dtddx
   q2(i,2)=q1(i,2)+(Fx(i-1,2)-Fx(i,2))*dtddx
   q2(i,3)=q1(i,3)+(Fx(i-1,3)-Fx(i,3))*dtddx
   q2(i,4)=q1(i,4)+(Fx(i-1,4)-Fx(i,4))*dtddx
   q2(i,5)=q1(i,5)
   q2(i,6)=q1(i,6)+(Fx(i-1,6)-Fx(i,6))*dtddx
   q2(i,7)=q1(i,7)+(Fx(i-1,7)-Fx(i,7))*dtddx
   q2(i,8)=q1(i,8)+(Fx(i-1,8)-Fx(i,8))*dtddx
   q2(i,9)=q1(i,9)
   q2(i,10)=q1(i,10)+(Fx(i-1,6)-Fx(i,6))*dtddx
   q2(i,11)=q1(i,11)+(Fx(i-1,7)-Fx(i,7))*dtddx
enddo


deallocate(SL,Fx)

end subroutine solverAdiMHD1D

subroutine solverAdi1D(this,q,q1,q2,dd)
class(grid)::this
double precision,dimension(1-this%nbuf:this%nMesh(1)+this%nbuf,this%nvar)::q,q1,q2
integer::nx,nbuf,nvar,dd
integer::i
double precision::dt,dx,dtddx,gam
double precision::rhoL,rhoM,rhoR
double precision::vxL,vxM,vxR
double precision::eneL,eneM,eneR
double precision::pL,pM,pR
double precision,dimension(:,:),allocatable::SL
double precision,dimension(:,:),allocatable::Fx
double precision::ql(this%nvar),qr(this%nvar),slopeL(this%nvar),slopeR(this%nvar),flux(this%nvar)
procedure(limiter),pointer::slope=>null()
procedure(fluxSolver),pointer::fluxPtr=>null()

nx=this%nMesh(1)
nbuf=this%nbuf
nvar=this%nvar
dt=this%dt
gam=this%adiGamma
if(this%solverType .eq. 1) then
   print *,"riemannSolverModule.f03: not yet implemented!"
   stop
   !fluxPtr=>fluxExactAdiHD1D
endif
if(this%solverType .eq. 2) then
   fluxPtr=>fluxHLLAdiHD1D
endif
if(this%solverType .eq. 3) then
   fluxPtr=>fluxHLLCAdiHD1D
endif

select case (this%limiterType)
   case(0)
     slope=>zslop
   case(1)
     slope=>vslop
   case(2)
     slope=>fslop
   case(3)
     slope=>minmod
end select

allocate(SL(1-nbuf:nx+nbuf,nvar))
allocate(Fx(1-nbuf:nx+nbuf,nvar))
do i=0,nx+1
   rhoL=q(i-1,1)
   rhoM=q(i  ,1)
   rhoR=q(i+1,1)

    vxL=q(i-1,2)/rhoL
    vxM=q(i  ,2)/rhoM
    vxR=q(i+1,2)/rhoR

   eneL=q(i-1,3)/rhoL
   eneM=q(i  ,3)/rhoM
   eneR=q(i+1,3)/rhoR

     pL=(gam-1.d0)*(rhoL*eneL-0.5d0*rhoL*vxL**2)
     pM=(gam-1.d0)*(rhoM*eneM-0.5d0*rhoM*vxM**2)
     pR=(gam-1.d0)*(rhoR*eneR-0.5d0*rhoR*vxR**2)

    SL(i,1)=slope(rhoL,rhoM,rhoR)
    SL(i,2)=slope( vxL, vxM, vxR)
    SL(i,3)=slope(  pL,  pM,  pR)
enddo

do i=0,nx+1
    ql(1)=q(i,1)
    ql(2)=q(i,2)/ql(1)
    ql(3)=(gam-1.d0)*(q(i,3)-0.5d0*ql(1)*ql(2)**2)

    qr(1)=q(i+1,1)
    qr(2)=q(i+1,2)/qr(1)
    qr(3)=(gam-1.d0)*(q(i+1,3)-0.5d0*qr(1)*qr(2)**2)

    slopeL(:)=SL(i,:)
    slopeR(:)=SL(i+1,:)
       call fluxPtr(ql,qr,slopeL,slopeR,this%adiGamma,flux,this%nvar)
    Fx(i,:)=flux(:)
enddo

do i=1,nx
   dtddx=dt/this%dx(1)%coords(i)
   q2(i,1)=q1(i,1)+(Fx(i-1,1)-Fx(i,1))*dtddx
   q2(i,2)=q1(i,2)+(Fx(i-1,2)-Fx(i,2))*dtddx
   q2(i,3)=q1(i,3)+(Fx(i-1,3)-Fx(i,3))*dtddx
enddo

deallocate(SL,Fx)

end subroutine solverAdi1D

subroutine solverPoly1D(this,q,q1,q2,dd)  !!!! add by hhwang, 2018, Feb. 1
class(grid)::this
double precision,dimension(1-this%nbuf:this%nMesh(1)+this%nbuf,this%nvar)::q,q1,q2
integer::nx,nbuf,nvar,dd
integer::i
double precision::dt,dx,dtddx
double precision::rhoL,rhoM,rhoR
double precision::vxL,vxM,vxR
double precision,dimension(:,:),allocatable::SL
double precision,dimension(:,:),allocatable::Fx
double precision::ql(this%nvar),qr(this%nvar),slopeL(this%nvar),slopeR(this%nvar),flux(this%nvar)
procedure(limiter),pointer::slope=>null()
procedure(fluxSolver),pointer::fluxPtr=>null()

nx=this%nMesh(1)
nbuf=this%nbuf
nvar=this%nvar
dt=this%dt
if(this%solverType .eq. 6) then 
   !print *, "warning...this solver does not conform the fluxSolver"
   !fluxPtr=>fluxHLLPolyHD1D
endif

select case (this%limiterType)
   case(1)
     slope=>vslop
end select

allocate(SL(1-nbuf:nx+nbuf,nvar))
allocate(Fx(1-nbuf:nx+nbuf,nvar))

do i=0,nx+1
   rhoL=q(i-1,1)
   rhoM=q(i  ,1)
   rhoR=q(i+1,1)

    vxL=q(i-1,2)/rhoL
    vxM=q(i  ,2)/rhoM
    vxR=q(i+1,2)/rhoR

    SL(i,1)=slope(rhoL,rhoM,rhoR)
    SL(i,2)=slope( vxL, vxM, vxR)
enddo

do i=0,nx+1
    ql(1)=q(i,1)
    ql(2)=q(i,2)/ql(1)

    qr(1)=q(i+1,1)
    qr(2)=q(i+1,2)/qr(1)

    slopeL(:)=SL(i,:)
    slopeR(:)=SL(i+1,:)
       call fluxHLLPolyHD1D(ql,qr,slopeL,slopeR,this%polyK,this%polyGamma,flux,this%nvar)
    Fx(i,:)=flux(:)
enddo

do i=1,nx
   dtddx=dt/this%dx(1)%coords(i)
   q2(i,1)=q1(i,1)+(Fx(i-1,1)-Fx(i,1))*dtddx
   q2(i,2)=q1(i,2)+(Fx(i-1,2)-Fx(i,2))*dtddx
enddo

deallocate(SL,Fx)
end subroutine solverPoly1D

subroutine solverIso1D(this,q,q1,q2,dd)
class(grid)::this
double precision,dimension(1-this%nbuf:this%nMesh(1)+this%nbuf,this%nvar)::q,q1,q2
integer::nx,nbuf,nvar,dd
integer::i
double precision::dt,dx,dtddx
double precision::rhoL,rhoM,rhoR
double precision::vxL,vxM,vxR
double precision,dimension(:,:),allocatable::SL
double precision,dimension(:,:),allocatable::Fx
double precision::ql(this%nvar),qr(this%nvar),slopeL(this%nvar),slopeR(this%nvar),flux(this%nvar)
procedure(limiter),pointer::slope=>null()
procedure(fluxSolver),pointer::fluxPtr=>null()

nx=this%nMesh(1)
nbuf=this%nbuf
nvar=this%nvar
dt=this%dt
if(this%solverType .eq. 1) then 
   fluxPtr=>fluxExactIsoHD1D
endif
if(this%solverType .eq. 2) then
   fluxPtr=>fluxHLLIsoHD1D
endif

select case (this%limiterType)
   case(1)
     slope=>vslop
end select

allocate(SL(1-nbuf:nx+nbuf,nvar))
allocate(Fx(1-nbuf:nx+nbuf,nvar))

do i=0,nx+1
   rhoL=q(i-1,1)
   rhoM=q(i  ,1)
   rhoR=q(i+1,1)

    vxL=q(i-1,2)/rhoL
    vxM=q(i  ,2)/rhoM
    vxR=q(i+1,2)/rhoR

    SL(i,1)=slope(rhoL,rhoM,rhoR)
    SL(i,2)=slope( vxL, vxM, vxR)
enddo

do i=0,nx+1
    ql(1)=q(i,1)
    ql(2)=q(i,2)/ql(1)

    qr(1)=q(i+1,1)
    qr(2)=q(i+1,2)/qr(1)

    slopeL(:)=SL(i,:)
    slopeR(:)=SL(i+1,:)
       call fluxPtr(ql,qr,slopeL,slopeR,this%snd,flux,this%nvar)
    Fx(i,:)=flux(:)
enddo

do i=1,nx
   dtddx=dt/this%dx(1)%coords(i)
   q2(i,1)=q1(i,1)+(Fx(i-1,1)-Fx(i,1))*dtddx
   q2(i,2)=q1(i,2)+(Fx(i-1,2)-Fx(i,2))*dtddx
enddo

deallocate(SL,Fx)
end subroutine solverIso1D

subroutine fluxHLLDAdiMHD1D(ql,qr,slopeL,slopeR,gam,flux,nvar)
implicit none
integer::nvar
double precision::ql(NVAR),qr(NVAR),slopeL(NVAR),slopeR(NVAR),flux(NVAR)
double precision::rhoL,vxL,vyL,vzL,bxL,byL,bzL,eneL,vsqL,BsqL,pL,pTL,fastL
double precision::rhoR,vxR,vyR,vzR,bxR,byR,bzR,eneR,vsqR,BsqR,pR,pTR,fastR
double precision::SM,SL,SR,SLs,SRs,pTs
double precision::vxLs,vxLss,vxRs,vxRss,vyss,vzss
double precision::pTLs,pTLss,pTRs,pTRss,byss,bzss
double precision::rhoLs,rhoRs
double precision::gam

double precision,dimension(:),allocatable::UL,UR,ULs,URs,ULss,URss
double precision,dimension(:),allocatable::FL,FR,FLs,FRs,FLss,FRss
double precision:: tempsL,tempsR,tempss

allocate(UL(NVAR))
allocate(UR(NVAR))
allocate(ULs(NVAR))
allocate(URs(NVAR))
allocate(ULss(NVAR))
allocate(URss(NVAR))
allocate(FL(NVAR))
allocate(FR(NVAR))
allocate(FLs(NVAR))
allocate(FRs(NVAR))
allocate(FLss(NVAR))
allocate(FRss(NVAR))

UL=ql+0.5d0*slopeL
UR=qr-0.5d0*slopeR

rhoL=UL(1)
vxL =UL(2)
vyL =UL(3)
vzL =UL(4)
bxL =0.5d0*(UL(5)+UR(5))
byL =UL(6)
bzL =UL(7)
pL = UL(8)
!eneL=UL(8)
vsqL=vxL**2.d0+vyL**2.d0+vzL**2.d0
BsqL=bxL**2.d0+byL**2.d0+bzL**2.d0
eneL = (pL/(gam-1.d0)+0.5d0*rhoL*vsqL+0.5d0*BsqL)/rhoL
!pL  =(gam-1.d0)*(rhoL*eneL-0.5d0*rhoL*vsqL-0.5d0*BsqL)
pTL =pL+0.5d0*BsqL
fastL=dsqrt((gam*pL+BsqL+dsqrt((gam*pL+BsqL)**2.d0-4.d0*gam*pL*bxL**2))/(2.d0*rhoL))

rhoR=UR(1)
vxR =UR(2)
vyR =UR(3)
vzR =UR(4)
bxR =0.5d0*(UL(5)+UR(5))
byR =UR(6)
bzR =UR(7)
!eneR=UR(8)
pR = UR(8)
vsqR=vxR**2.d0+vyR**2.d0+vzR**2.d0
BsqR=bxR**2.d0+byR**2.d0+bzR**2.d0
eneR = (pR/(gam-1.d0)+0.5d0*rhoR*vsqR+0.5d0*BsqR)/rhoR
!pR  =(gam-1.d0)*(rhoR*eneR-0.5d0*rhoR*vsqR-0.5d0*BsqR)
pTR =pR+0.5d0*BsqR
fastR=dsqrt((gam*pR+BsqR+dsqrt((gam*pR+BsqR)**2.d0-4.d0*gam*pR*bxR**2))/(2.d0*rhoR))

SL=dmin1(vxL,vxR)-dmax1(fastL,fastR)
SR=dmax1(vxL,vxR)+dmax1(fastL,fastR)

SM=((SR-vxR)*rhoR*vxR-(SL-vxL)*rhoL*vxL-pTR+pTL)/((SR-vxR)*rhoR-(SL-vxL)*rhoL)
pTs=((SR-vxR)*rhoR*pTL-(SL-vxL)*rhoL*pTR+rhoL*rhoR*(SR-vxR)*(SL-vxL)*(vxR-vxL))/((SR-vxR)*rhoR-(SL-vxL)*rhoL)

vxLs =SM
vxLss=SM
vxRss=SM
vxRs =SM

pTLs =pTs
pTLss=pTs
pTRs =pTs
pTRss=pTs

FL(1)=rhoL*vxL
FL(2)=rhoL*vxL**2.d0+pTL-bxL**2.d0
FL(3)=rhoL*vyL*vxL-bxL*byL
FL(4)=rhoL*vzL*vxL-bxL*bzL
FL(5)=0.d0
FL(6)=byL*vxL-bxL*vyL
FL(7)=bzL*vxL-bxL*vzL
FL(8)=(rhoL*eneL+pTL)*vxL-bxL*(vxL*bxL+vyL*byL+vzL*bzL)

FR(1)=rhoR*vxR
FR(2)=rhoR*vxR**2.d0+pTR-bxR**2.d0
FR(3)=rhoR*vyR*vxR-bxR*byR
FR(4)=rhoR*vzR*vxR-bxR*bzR
FR(5)=0.d0
FR(6)=byR*vxR-bxR*vyR
FR(7)=bzR*vxR-bxR*vzR
FR(8)=(rhoR*eneR+pTR)*vxR-bxR*(vxR*bxR+vyR*byR+vzR*bzR)

! ============================
if(SL .eq. SM) then
ULs(1)=rhoL
else
ULs(1)=rhoL*(SL-vxL)/(SL-SM)
endif
ULs(2)=vxLs
tempsL=rhoL*(SL-vxL)*(SL-SM)-bxL**2.d0
if(tempsL .eq. 0.d0) then
ULs(3)=vyL
ULs(4)=vzL
ULs(5)=bxL
ULs(6)=byL
ULs(7)=bzL
else
ULs(3)=vyL-bxL*byL*(SM-vxL)/tempsL
ULs(4)=vzL-bxL*bzL*(SM-vxL)/tempsL
ULs(5)=bxL
ULs(6)=byL*(rhoL*(SL-vxL)**2.d0-bxL**2)/tempsL
ULs(7)=bzL*(rhoL*(SL-vxL)**2.d0-bxL**2)/tempsL
endif
ULs(8)=((SL-vxL)*rhoL*eneL-pTL*vxL+pTs*SM+bxL*(vxL*bxL+vyL*byL+vzL*bzL-ULs(2)*ULs(5)-ULs(3)*ULs(6)-ULs(4)*ULs(7)))/(SL-SM)/ULs(1)

if(SR .eq. SM) then
URs(1)=rhoR
else
URs(1)=rhoR*(SR-vxR)/(SR-SM)
endif
URs(2)=vxRs
tempsR=rhoR*(SR-vxR)*(SR-SM)-bxR**2.d0
if(tempsR .eq. 0.d0) then
URs(3)=vyR
URs(4)=vzR
URs(5)=bxR
URs(6)=byR
URs(7)=bzR
else
URs(3)=vyR-bxR*byR*(SM-vxR)/tempsR
URs(4)=vzR-bxR*bzR*(SM-vxR)/tempsR
URs(5)=bxR
URs(6)=byR*(rhoR*(SR-vxR)**2.d0-bxR**2)/tempsR
URs(7)=bzR*(rhoR*(SR-vxR)**2.d0-bxR**2)/tempsR
endif
URs(8)=((SR-vxR)*rhoR*eneR-pTR*vxR+pTs*SM+bxR*(vxR*bxR+vyR*byR+vzR*bzR-URs(2)*URs(5)-URs(3)*URs(6)-URs(4)*URs(7)))/(SR-SM)/URs(1)
!==========================
tempss=dsqrt(ULs(1))+dsqrt(URs(1))
vyss=(dsqrt(ULs(1))*ULs(3)+dsqrt(URs(1))*URs(3)+(URs(6)-ULs(6))*dsign(1.d0,bxL))/tempss
vzss=(dsqrt(ULs(1))*ULs(4)+dsqrt(URs(1))*URs(4)+(URs(7)-ULs(7))*dsign(1.d0,bxL))/tempss
byss=(dsqrt(ULs(1))*URs(6)+dsqrt(URs(1))*ULs(6)+dsqrt(URs(1)*ULs(1))*(URs(3)-ULs(3))*dsign(1.d0,bxL))/tempss
bzss=(dsqrt(ULs(1))*URs(7)+dsqrt(URs(1))*ULs(7)+dsqrt(URs(1)*ULs(1))*(URs(4)-ULs(4))*dsign(1.d0,bxL))/tempss

ULss(1)=ULs(1)
ULss(2)=vxLss
if(BxL .eq. 0.d0) then
ULss(3)=ULs(3)
ULss(4)=ULs(4)
ULss(5)=0.d0
ULss(6)=ULs(6)
ULss(7)=ULs(7)
ULss(8)=ULs(8)
else
ULss(3)=vyss
ULss(4)=vzss
ULss(5)=bxL
ULss(6)=byss
ULss(7)=bzss
ULss(8)=(ULs(8)*ULs(1)-dsqrt(ULs(1))*(ULs(2)*ULs(5)+ULs(3)*ULs(6)+ULs(4)*ULs(7)&
 & -ULss(2)*ULss(5)-ULss(3)*ULss(6)-ULss(4)*ULss(7))*dsign(1.d0,bxL))/ULss(1)
endif

URss(1)=URs(1)
URss(2)=vxRss
if(BxR .eq. 0.d0) then
URss(3)=URs(3)
URss(4)=URs(4)
URss(5)=0.d0
URss(6)=URs(6)
URss(7)=URs(7)
URss(8)=URs(8)
else
URss(3)=vyss
URss(4)=vzss
URss(5)=bxR
URss(6)=byss
URss(7)=bzss
URss(8)=(URs(8)*URs(1)+dsqrt(URs(1))*(URs(2)*URs(5)+URs(3)*URs(6)+URs(4)*URs(7)&
&-URss(2)*URss(5)-URss(3)*URss(6)-URss(4)*URss(7))*dsign(1.d0,bxL))/URss(1)
endif

FLs(1)=FL(1)+SL*(ULs(1)-UL(1))
FLs(2)=FL(2)+SL*(ULs(1)*ULs(2)-UL(1)*UL(2))
FLs(3)=FL(3)+SL*(ULs(1)*ULs(3)-UL(1)*UL(3))
FLs(4)=FL(4)+SL*(ULs(1)*ULs(4)-UL(1)*UL(4))
FLs(5)=FL(5)+SL*(ULs(5)-UL(5))
FLs(6)=FL(6)+SL*(ULs(6)-UL(6))
FLs(7)=FL(7)+SL*(ULs(7)-UL(7))
FLs(8)=FL(8)+SL*(ULs(1)*ULs(8)-UL(1)*eneL)

FRs(1)=FR(1)+SR*(URs(1)-UR(1))
FRs(2)=FR(2)+SR*(URs(1)*URs(2)-UR(1)*UR(2))
FRs(3)=FR(3)+SR*(URs(1)*URs(3)-UR(1)*UR(3))
FRs(4)=FR(4)+SR*(URs(1)*URs(4)-UR(1)*UR(4))
FRs(5)=FR(5)+SR*(URs(5)-UR(5))
FRs(6)=FR(6)+SR*(URs(6)-UR(6))
FRs(7)=FR(7)+SR*(URs(7)-UR(7))
FRs(8)=FR(8)+SR*(URs(1)*URs(8)-UR(1)*eneR)

SLs=SM-dabs(bxL)/dsqrt(ULs(1))
SRs=SM+dabs(bxR)/dsqrt(URs(1))

FLss(1)=FL(1)+SLs*ULss(1)-(SLs-SL)*ULs(1)-SL*UL(1)
FLss(2)=FL(2)+SLs*(ULss(1)*ULss(2))-(SLs-SL)*(ULs(1)*ULs(2))-SL*(UL(1)*UL(2))
FLss(3)=FL(3)+SLs*(ULss(1)*ULss(3))-(SLs-SL)*(ULs(1)*ULs(3))-SL*(UL(1)*UL(3))
FLss(4)=FL(4)+SLs*(ULss(1)*ULss(4))-(SLs-SL)*(ULs(1)*ULs(4))-SL*(UL(1)*UL(4))
FLss(5)=FL(5)+SLs*ULss(5)-(SLs-SL)*ULs(5)-SL*UL(5)
FLss(6)=FL(6)+SLs*ULss(6)-(SLs-SL)*ULs(6)-SL*UL(6)
FLss(7)=FL(7)+SLs*ULss(7)-(SLs-SL)*ULs(7)-SL*UL(7)
FLss(8)=FL(8)+SLs*(ULss(1)*ULss(8))-(SLs-SL)*(ULs(1)*ULs(8))-SL*(UL(1)*eneL)

FRss(1)=FR(1)+SRs*URss(1)-(SRs-SR)*URs(1)-SR*UR(1)
FRss(2)=FR(2)+SRs*(URss(1)*URss(2))-(SRs-SR)*(URs(1)*URs(2))-SR*(UR(1)*UR(2))
FRss(3)=FR(3)+SRs*(URss(1)*URss(3))-(SRs-SR)*(URs(1)*URs(3))-SR*(UR(1)*UR(3))
FRss(4)=FR(4)+SRs*(URss(1)*URss(4))-(SRs-SR)*(URs(1)*URs(4))-SR*(UR(1)*UR(4))
FRss(5)=FR(5)+SRs*URss(5)-(SRs-SR)*URs(5)-SR*UR(5)
FRss(6)=FR(6)+SRs*URss(6)-(SRs-SR)*URs(6)-SR*UR(6)
FRss(7)=FR(7)+SRs*URss(7)-(SRs-SR)*URs(7)-SR*UR(7)
FRss(8)=FR(8)+SRs*(URss(1)*URss(8))-(SRs-SR)*(URs(1)*URs(8))-SR*(UR(1)*eneR)

if(SL .gt. 0.d0) then
  flux=FL
elseif (SL .le. 0.d0 .and. SLs .ge. 0.d0) then
  flux=FLs
elseif (SLs .le. 0.d0 .and. SM .ge. 0.d0) then
  flux=FLss
elseif (SM .le. 0.d0 .and. SRs .ge. 0.d0) then
  flux=FRss
elseif (SRs .le. 0.d0 .and. SR .ge. 0.d0) then
  flux=FRs
elseif (SR .lt. 0.d0) then
  flux=FR
endif

deallocate(UL)
deallocate(UR)
deallocate(ULs)
deallocate(URs)
deallocate(ULss)
deallocate(URss)
deallocate(FL)
deallocate(FR)
deallocate(FLs)
deallocate(FRs)
deallocate(FLss)
deallocate(FRss)
end subroutine fluxHLLDAdiMHD1D

subroutine fluxHLLAdiMHD1D(ql,qr,slopeL,slopeR,gam,flux,nvar)
implicit none
integer::nvar
double precision::ql(NVAR),qr(NVAR),slopeL(NVAR),slopeR(NVAR),flux(NVAR)
double precision, dimension(:), allocatable:: FL,FR,Fs,UL,UR
double precision:: gam ! sound speed
double precision:: rhoL,vxL,vyL,vzL,bxL,byL,bzL,eneL,sndfL,pL,BsqL,vsqL,pTL
double precision:: rhoR,vxR,vyR,vzR,bxR,byR,bzR,eneR,sndfR,pR,BsqR,vsqR,pTR
double precision:: SL,SR

allocate(FL(NVAR))
allocate(FR(NVAR))
allocate(Fs(NVAR))
allocate(UL(NVAR))
allocate(UR(NVAR))

UL=ql+0.5d0*slopeL
UR=qr-0.5d0*slopeR

rhoL=UL(1)
vxL =UL(2)
vyL =UL(3)
vzL =UL(4)
bxL =0.5d0*(UL(5)+UR(5))
byL =UL(6)
bzL =UL(7)
pL = UL(8)
!eneL = UL(8)
vsqL = vxL**2.d0+vyL**2.d0+vzL**2.d0
BsqL = bxL**2.d0+byL**2.d0+bzL**2.d0
eneL = (pL/(gam-1.d0)+0.5d0*rhoL*vsqL+0.5d0*BsqL)/rhoL
!pL = (gam-1.d0)*(rhoL*eneL-0.5d0*rhoL*vsqL-0.5d0*BsqL) ! gaseous pressure
sndfL = dsqrt((gam*pL+BsqL+dsqrt((gam*pL+BsqL)**2.d0-4.d0*gam*pL*bxL**2.d0))/(2.d0*rhoL)) ! fast wave
pTL = pL+0.5d0*BsqL

rhoR=UR(1)
vxR =UR(2)
vyR =UR(3)
vzR =UR(4)
bxR =0.5d0*(UL(5)+UR(5))
byR =UR(6)
bzR =UR(7)
!eneR=UR(8)
pR=UR(8)
vsqR = vxR**2.d0+vyR**2.d0+vzR**2.d0
BsqR = bxR**2.d0+byR**2.d0+bzR**2.d0
eneR = (pR/(gam-1.d0)+0.5d0*rhoR*vsqR+0.5d0*BsqR)/rhoR
!pR = (gam-1.d0)*(rhoR*eneR-0.5d0*rhoR*vsqR-0.5d0*BsqR) ! gaseous pressure
sndfR = dsqrt((gam*pR+BsqR+dsqrt((gam*pR+BsqR)**2.d0-4.d0*gam*pR*bxR**2.d0))/(2.d0*rhoR)) ! fast wave
pTR = pR+0.5d0*BsqR

SL=dmin1(vxL,vxR)-dmax1(sndfL,sndfR)
SR=dmax1(vxL,vxR)+dmax1(sndfL,sndfR)

FL(1)=rhoL*vxL
FL(2)=rhoL*vxL**2.d0+pTL-bxL**2.d0
FL(3)=rhoL*vxL*vyL-bxL*byL
FL(4)=rhoL*vxL*vzL-bxL*bzL
FL(5)=0.d0
FL(6)=byL*vxL-bxL*vyL
FL(7)=bzL*vxL-bxL*vzL
FL(8)=(rhoL*eneL+pTL)*vxL-bxL*(vxL*bxL+vyL*byL+vzL*bzL)

FR(1)=rhoR*vxR
FR(2)=rhoR*vxR**2.d0+pTR-bxR**2.d0
FR(3)=rhoR*vxR*vyR-bxR*byR
FR(4)=rhoR*vxR*vzR-bxR*bzR
FR(5)=0.d0
FR(6)=byR*vxR-bxR*vyR
FR(7)=bzR*vxR-bxR*vzR
FR(8)=(rhoR*eneR+pTR)*vxR-bxR*(vxR*bxR+vyR*byR+vzR*bzR)

Fs(1)=(SR*FL(1)-SL*FR(1)+SR*SL*(UR(1)-UL(1)))/(SR-SL)
Fs(2)=(SR*FL(2)-SL*FR(2)+SR*SL*(UR(1)*UR(2)-UL(1)*UL(2)))/(SR-SL)
Fs(3)=(SR*FL(3)-SL*FR(3)+SR*SL*(UR(1)*UR(3)-UL(1)*UL(3)))/(SR-SL)
Fs(4)=(SR*FL(4)-SL*FR(4)+SR*SL*(UR(1)*UR(4)-UL(1)*UL(4)))/(SR-SL)
Fs(5)=0.d0
Fs(6)=(SR*FL(6)-SL*FR(6)+SR*SL*(UR(6)-UL(6)))/(SR-SL)
Fs(7)=(SR*FL(7)-SL*FR(7)+SR*SL*(UR(7)-UL(7)))/(SR-SL)
Fs(8)=(SR*FL(8)-SL*FR(8)+SR*SL*(UR(1)*eneR-UL(1)*eneL))/(SR-SL)

if(SL .gt. 0.d0) then
   flux=FL
elseif(SL .le. 0.d0 .and. SR .ge. 0.d0) then
   flux=Fs
else
   flux=FR
endif


deallocate(FL)
deallocate(FR)
deallocate(Fs)
deallocate(UL)
deallocate(UR)

end subroutine fluxHLLAdiMHD1D

subroutine fluxHLLPolyHD1D(ql,qr,slopeL,slopeR,polyK,polyGamma,flux,nvar)
implicit none
integer::nvar
double precision::polyK,polyGamma
double precision::ql(NVAR),qr(NVAR),slopeL(NVAR),slopeR(NVAR),flux(NVAR)
double precision, dimension(:), allocatable:: FL,FR,Fs,UL,UR
double precision:: spdL,spdR ! sound speed
double precision:: rhoL,vxL,vyL,pTL
double precision:: rhoR,vxR,vyR,pTR
double precision:: SL,SR

allocate(FL(NVAR))
allocate(FR(NVAR))
allocate(Fs(NVAR))
allocate(UL(NVAR))
allocate(UR(NVAR))

UL=ql+0.5d0*slopeL
UR=qr-0.5d0*slopeR

rhoL=UL(1)
vxL =UL(2)
pTL =polyK*rhoL**polyGamma
spdL=dsqrt(polyGamma*pTL/rhoL)

rhoR=UR(1)
vxR =UR(2)
pTR =polyK*rhoR**polyGamma
spdR=dsqrt(polyGamma*pTR/rhoR)

SL=dmin1(vxL-spdL,vxR-spdR)
SR=dmax1(vxL+spdL,vxR+spdR)

FL(1)=rhoL*vxL
FL(2)=rhoL*vxL**2.d0+pTL

FR(1)=rhoR*vxR
FR(2)=rhoR*vxR**2.d0+pTR

Fs(1)=(SR*FL(1)-SL*FR(1)+SR*SL*(UR(1)-UL(1)))/(SR-SL)
Fs(2)=(SR*FL(2)-SL*FR(2)+SR*SL*(UR(1)*UR(2)-UL(1)*UL(2)))/(SR-SL)

if(SL .gt. 0.d0) then
   flux=FL
elseif(SL .le. 0.d0 .and. SR .ge. 0.d0) then
   flux=Fs
else
   flux=FR
endif


deallocate(FL)
deallocate(FR)
deallocate(Fs)
deallocate(UL)
deallocate(UR)
end subroutine fluxHLLPolyHD1D



subroutine fluxHLLIsoHD1D(ql,qr,slopeL,slopeR,snd,flux,nvar)
implicit none
integer::nvar
double precision::ql(NVAR),qr(NVAR),slopeL(NVAR),slopeR(NVAR),flux(NVAR)
double precision, dimension(:), allocatable:: FL,FR,Fs,UL,UR
double precision:: snd,spd ! sound speed
double precision:: rhoL,vxL,vyL,pTL
double precision:: rhoR,vxR,vyR,pTR
double precision:: SL,SR

allocate(FL(NVAR))
allocate(FR(NVAR))
allocate(Fs(NVAR))
allocate(UL(NVAR))
allocate(UR(NVAR))

UL=ql+0.5d0*slopeL
UR=qr-0.5d0*slopeR

spd=snd
rhoL=UL(1)
vxL =UL(2)
pTL =rhoL*spd**2.d0

rhoR=UR(1)
vxR =UR(2)
pTR =rhoR*spd**2.d0

SL=dmin1(vxL-spd,vxR-spd)
SR=dmax1(vxL+spd,vxR+spd)

FL(1)=rhoL*vxL
FL(2)=rhoL*vxL**2.d0+pTL

FR(1)=rhoR*vxR
FR(2)=rhoR*vxR**2.d0+pTR

Fs(1)=(SR*FL(1)-SL*FR(1)+SR*SL*(UR(1)-UL(1)))/(SR-SL)
Fs(2)=(SR*FL(2)-SL*FR(2)+SR*SL*(UR(1)*UR(2)-UL(1)*UL(2)))/(SR-SL)

if(SL .gt. 0.d0) then
   flux=FL
elseif(SL .le. 0.d0 .and. SR .ge. 0.d0) then
   flux=Fs
else
   flux=FR
endif


deallocate(FL)
deallocate(FR)
deallocate(Fs)
deallocate(UL)
deallocate(UR)
end subroutine fluxHLLIsoHD1D

subroutine fluxExactIsoHD1D(ql,qr,slopeL,slopeR,snd,flux,nvar)
implicit none
integer::nvar
double precision::ql(nvar),qr(nvar),qm(nvar),slopeL(nvar),slopeR(nvar),flux(nvar),snd
double precision:: rhoL,rhoR,vxL,vxR,vyL,vyR
double precision::c,shkR,rarR,csq,vstar,rhostar,rho0,rhosR,srhosR,aisrhosR,sigma2, &
                   rhosL,srhosL,aisrhosL,sigma1,aa,bb,uu,srhostar
integer:: Intr,n
double precision:: tolerr = 1.d-10

ql(1) = ql(1)+0.5d0*slopeL(1)
ql(2) = ql(2)+0.5d0*slopeL(2)

qr(1) = qr(1)-0.5d0*slopeR(1)
qr(2) = qr(2)-0.5d0*slopeR(2)

c = snd
Intr = 10


if (ql(1).eq.qr(1).and.ql(2).eq.qr(2)) then
     qm(1) = qr(1)
     qm(2) = qr(2)
goto 1000
endif

      rhoL = ql(1)
      rhoR = qr(1)
       vxL = ql(2)
       vxR = qr(2)


      shkR = -c*(rhoR-rhoL)/dsqrt(rhoR*rhoL)
      rarR = -c*dlog(rhoR/rhoL)
       csq =  c*c

      if(rhoL .ge. rhoR) then
        if ((vxR-vxL) .gt. rarR) then
          goto 100
        elseif ((vxR-vxL) .gt. -shkR) then
          goto 200
        else
          goto 400
        endif
      else
        if ((vxR-vxL) .lt. shkR ) then
          goto 400
        elseif ((vxR-vxL) .lt. -rarR) then
          goto 300
        else
          goto 100
        endif
      endif
!cccccc 1-rarefaction wave, 2-rarefaction wave
!cccccc 1-rarefaction wave, 2-rarefaction wave
!cccccc 1-rarefaction wave, 2-rarefaction wave
100   vstar   = 0.5d0*(vxL+vxR)+0.5d0*c*dlog(rhoL/rhoR)
      rhostar = rhoL*dexp(-(vstar-vxL)/c)

      if ((vxL-c) .ge. 0.) then
         qm(2) = vxL
         qm(1) = rhoL
      elseif ((vstar-c) .ge. 0.) then
         qm(2) =  c
         qm(1) =  rhoL*dexp((vxL-c)/c)
      elseif ((vstar+c) .ge. 0.) then
         qm(2) =  vstar
         qm(1) =  rhostar
      elseif ((  vxR+c) .ge. 0.) then
         qm(2) =  -c
         qm(1) =  rhoR*dexp(-(vxR+c)/c)
      else
         qm(2) =   vxR
         qm(1) =  rhoR
      endif
      goto 1000
!cccccc1-rarefaction wave, 2-shock wave
!cccccc1-rarefaction wave, 2-shock wave
!cccccc1-rarefaction wave, 2-shock wave
200   rho0 = rhoR
      do n = 1, INtr
        rhosR    = rho0/rhoR
        srhosR   = dsqrt(rhosR)
        aisrhosR = 1.d0/srhosR
        rhostar  = rho0 - &
       (vxR-vxL+c*(srhosR-aisrhosR)+c*dlog(rho0/rhoL)) &
       /(0.5d0*c/rho0*(srhosR+aisrhosR)+c/rho0)
        if (dabs((rho0-rhostar)/rho0) .lt. tolerr) then
            goto 210
        endif
        rho0 = rhostar
      enddo
      write(6,*) 'Newton Interation is divergent at x-dir 1R2S'
      write(6,*) 'rhostar,rhoL=',rhostar,rhoL
      write(6,*) 'uxL    , uxR=',vxL,      vxR
      stop
210   vstar  = vxL - c*dlog(rhostar/rhoL)
      sigma2 = vxR + c*dsqrt(rhostar/rhoR)
      if ( sigma2 .le. 0.) then
         qm(2) =  vxR
         qm(1) = rhoR
      elseif ((vstar-c) .le. 0.) then
         qm(2) = vstar
         qm(1) = rhostar
      elseif ((vxL-c) .le. 0.) then
         qm(2) =  c
         qm(1) =  rhoL*dexp((vxL-c)/c)
      else
         qm(2) = vxL
         qm(1) = rhoL
      endif
      goto 1000
!cccccc 1-shock wave, 2-rarefaction wave
!cccccc 1-shock wave, 2-rarefaction wave
!cccccc 1-shock wave, 2-rarefaction wave
300   rho0 = rhoL
      do n = 1, INtr
        rhosL    = rho0/rhoL
        srhosL   = dsqrt(rhosL)
        aisrhosL = 1.d0/srhosL
        rhostar  = rho0- &
      (vxR-vxL+c*(srhosL-aisrhosL)+c*dlog(rho0/rhoR)) &
       /(0.5d0*c/rho0*(srhosL+aisrhosL)+c/rho0)
        if (dabs((rho0-rhostar)/rho0) .lt. tolerr) then
            goto 310
        endif
        rho0 = rhostar
      enddo
      write(6,*) 'Newton Interation is divergent at x-dir 1S2R'
      write(6,*) 'rhostar,rhoL=',rhostar,rhoL
      write(6,*) 'uxL    , uxR=',vxL,      vxR
      stop
310   vstar  = vxR + c*dlog(rhostar/rhoR)
      sigma1 = vxL - c*dsqrt(rhostar/rhoL)
      if (sigma1 .ge. 0.) then
         qm(2) = vxL
         qm(1) = rhoL
      elseif ((vstar+c) .ge. 0.) then
         qm(2) = vstar
         qm(1) = rhostar
      elseif ((vxR+c) .ge. 0.) then
         qm(2) =  -c
         qm(1) =  rhoR*dexp(-(vxR+c)/c)
      else
         qm(2) = vxR
         qm(1) = rhoR
      endif
      goto 1000
!       1-shock wave, 2-shock wave
!       1-shock wave, 2-shock wave
!       1-shock wave, 2-shock wave
400   aa       = 1.d0/dsqrt(rhoL)+1.d0/dsqrt(rhoR)
      bb       =    dsqrt(rhoL)+   dsqrt(rhoR)
      uu       = vxR-vxL
      srhostar = (-uu+dsqrt(uu**2+4.d0*csq*aa*bb))/(2.d0*c*aa)
      rhostar  = srhostar*srhostar
      vstar    = vxL - c*(dsqrt(rhostar/rhoL)-dsqrt(rhoL/rhostar))
      sigma1   = vxL - c*dsqrt(rhostar/rhoL)
      sigma2   = vxR + c*dsqrt(rhostar/rhoR)
      if     (sigma1 .ge. 0.) then
       qm(2) = vxL
       qm(1) = rhoL
      elseif (sigma2 .ge. 0.) then
       qm(2) = vstar
       qm(1) = rhostar
      else
       qm(2) = vxR
       qm(1) = rhoR
      endif
      goto 1000

1000  flux(1) = qm(1)*(   qm(2)          )
      flux(2) = qm(1)*( qm(2)*qm(2) + c*c)
end subroutine fluxExactIsoHD1D
end module riemannSolverModule
