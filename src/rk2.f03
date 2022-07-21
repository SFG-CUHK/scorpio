subroutine rk2ADsg_3D(nthis,qn,qn1,qn2,ithis,qi,qi1,qi2)
use gridModule
use riemannSolverModule
use mpi
implicit none
type(grid)::nthis,ithis
double precision,dimension(1-nthis%nbuf:nthis%nMesh(1)+nthis%nbuf,1-nthis%nbuf:nthis%nMesh(2)+nthis%nbuf,&
                           1-nthis%nbuf:nthis%nMesh(3)+nthis%nbuf,nthis%nvar)::qn,qn1,qn2
double precision,dimension(:),allocatable::totalden
double precision,dimension(1-ithis%nbuf:ithis%nMesh(1)+ithis%nbuf,1-ithis%nbuf:ithis%nMesh(2)+ithis%nbuf,&
                           1-ithis%nbuf:ithis%nMesh(3)+ithis%nbuf,ithis%nvar)::qi,qi1,qi2
procedure(riemannSolver3D),pointer::rieSolvern,rieSolveri
integer::nvar,nx,ny,nz,nbuf,eosTypen,eosTypei,SolverTypen,SolverTypei,solverTypeOriginaln,solverTypeOriginali
double precision::sgfx,sgfy,sgfz,den,momx,momy,momz
double precision::t,dt
integer::ierr,datasize,i,j,k,pos
logical::global_changeSolvern,global_changeSolveri

nbuf=nthis%nbuf
nx=nthis%nMesh(1)
ny=nthis%nMesh(2)
nz=nthis%nMesh(3)
nvar=nthis%nvar
allocate(totalden((nx+2*nbuf)*(ny+2*nbuf)*(nz+2*nbuf)*nvar))
eosTypen=nthis%eosType
eosTypei=ithis%eosType
solverTypen=nthis%solverType
solverTypei=ithis%solverType
solverTypeOriginaln=0
solverTypeOriginali=0

if(eosTypen .eq. 1) then
  select case (solverTypen)
    case (1,2)
      rieSolvern=>solverIso3D
    case default
      print *, "rk2.f03: no suitable 3D solver found!"
      stop
  end select
elseif(eosTypen .eq. 2) then
  select case (solverTypen)
    case (2,3)
      rieSolvern=>solverAdi3D
    case (4,5)
      rieSolvern=>solverAdiMHD3D
    case default
      print *, "rk2.f03: no suitable 3D solver found!"
      stop
  end select
endif

if(eosTypei .eq. 1) then
  select case (solverTypei)
    case (1,2)
      rieSolveri=>solverIso3D
    case default
      print *, "rk2.f03: no suitable 3D solver found!"
      stop
  end select
elseif(eosTypei .eq. 2) then
  select case (solverTypei)
    case (2,3)
      rieSolveri=>solverAdi3D
    case (4,5)
      rieSolveri=>solverAdiMHD3D
    case default
      print *, "rk2.f03: no suitable 3D solver found!"
      stop
  end select
endif

t=nthis%t
datasize=(nthis%nMesh(1)+2*nthis%nbuf)*(nthis%nMesh(2)+2*nthis%nbuf)*(nthis%nMesh(3)+2*nthis%nbuf)

do

!!!!!!!!!!! zwg: step1 of RK2 

   call rieSolvern(nthis,qn,qn,qn1,dd=1)
   call rieSolvern(nthis,qn,qn1,qn1,dd=2)
   call rieSolvern(nthis,qn,qn1,qn1,dd=3)

   call rieSolveri(ithis,qi,qi,qi1,dd=1)
   call rieSolveri(ithis,qi,qi1,qi1,dd=2)
   call rieSolveri(ithis,qi,qi1,qi1,dd=3)

   totalden(1:datasize)=nthis%q(1:datasize)+ithis%q(1:datasize)
   call nthis%calcSelfgravity(totalden)
   
   dt=nthis%dt
   !!!!!! apply gravity !!!!!!!
   do k=1,nz
     do j=1,ny
       do i=1,nx
         pos=(k-1)*nx*ny+(j-1)*nx+i
         den = qn(i,j,k,1)
         momx = qn(i,j,k,2)
         momy = qn(i,j,k,3)
         momz = qn(i,j,k,4)
         sgfx = nthis%sgfx(pos)
         sgfy = nthis%sgfy(pos)
         sgfz = nthis%sgfz(pos)

         qn1(i,j,k,2)=qn1(i,j,k,2)+den*sgfx*dt
         qn1(i,j,k,3)=qn1(i,j,k,3)+den*sgfy*dt
         qn1(i,j,k,4)=qn1(i,j,k,4)+den*sgfz*dt
         qn1(i,j,k,5)=qn1(i,j,k,5)+(momx*sgfx+momy*sgfy+momz*sgfz)*dt

         den = qi(i,j,k,1)
         momx = qi(i,j,k,2)
         momy = qi(i,j,k,3)
         momz = qi(i,j,k,4)

         qi1(i,j,k,2)=qi1(i,j,k,2)+den*sgfx*dt
         qi1(i,j,k,3)=qi1(i,j,k,3)+den*sgfy*dt
         qi1(i,j,k,4)=qi1(i,j,k,4)+den*sgfz*dt
         qi1(i,j,k,8)=qi1(i,j,k,8)+(momx*sgfx+momy*sgfy+momz*sgfz)*dt

       enddo
     enddo
   enddo 
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!
   call evolveAD3D(nthis,ithis,nthis%q1,ithis%q1)

   call nthis%exchangeBdryMPI(nthis%q1,nthis%winq1)
   call nthis%setBoundary(nthis%q1)
   call ithis%exchangeBdryMPI(ithis%q1,ithis%winq1)
   call ithis%setBoundary(ithis%q1)


   !!!!!!!!!!! zwg: step2  of RK2

   call rieSolvern(nthis,qn1,qn1,qn2,dd=1)
   call rieSolvern(nthis,qn1,qn2,qn2,dd=2)
   call rieSolvern(nthis,qn1,qn2,qn2,dd=3)

   call rieSolveri(ithis,qi1,qi1,qi2,dd=1)
   call rieSolveri(ithis,qi1,qi2,qi2,dd=2)
   call rieSolveri(ithis,qi1,qi2,qi2,dd=3)

   totalden(1:datasize)=nthis%q1(1:datasize)+ithis%q1(1:datasize)
   call nthis%calcSelfgravity(totalden)
   !!!!!!! apply gravity !!!!!!!
   do k=1,nz
     do j=1,ny
       do i=1,nx
         pos=(k-1)*nx*ny+(j-1)*nx+i
         den = qn1(i,j,k,1)
         momx = qn1(i,j,k,2)
         momy = qn1(i,j,k,3)
         momz = qn1(i,j,k,4)
         sgfx = nthis%sgfx(pos)
         sgfy = nthis%sgfy(pos)
         sgfz = nthis%sgfz(pos)

         qn2(i,j,k,2)=qn2(i,j,k,2)+den*sgfx*dt
         qn2(i,j,k,3)=qn2(i,j,k,3)+den*sgfy*dt
         qn2(i,j,k,4)=qn2(i,j,k,4)+den*sgfz*dt
         qn2(i,j,k,5)=qn2(i,j,k,5)+(momx*sgfx+momy*sgfy+momz*sgfz)*dt

         den = qi1(i,j,k,1)
         momx = qi1(i,j,k,2)
         momy = qi1(i,j,k,3)
         momz = qi1(i,j,k,4)

         qi2(i,j,k,2)=qi2(i,j,k,2)+den*sgfx*dt
         qi2(i,j,k,3)=qi2(i,j,k,3)+den*sgfy*dt
         qi2(i,j,k,4)=qi2(i,j,k,4)+den*sgfz*dt
         qi2(i,j,k,8)=qi2(i,j,k,8)+(momx*sgfx+momy*sgfy+momz*sgfz)*dt

       enddo
     enddo
   enddo
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!

   call evolveAD3D(nthis,ithis,nthis%q2,ithis%q2)

   call nthis%exchangeBdryMPI(nthis%q2,nthis%winq2)
   call nthis%setBoundary(nthis%q2)
   call ithis%exchangeBdryMPI(ithis%q2,ithis%winq2)
   call ithis%setBoundary(ithis%q2)

   call MPI_ALLREDUCE(nthis%changeSolver,global_changeSolvern,1,MPI_LOGICAL,MPI_LOR,MPI_COMM_WORLD,ierr)
   call MPI_ALLREDUCE(ithis%changeSolver,global_changeSolveri,1,MPI_LOGICAL,MPI_LOR,MPI_COMM_WORLD,ierr)

   exit
   if(global_changeSolveri .eqv. .false.) then
      exit
   else
      select case(ithis%solverType)
        case(5)
          print *,"rk2.f03: myid",myid,"change to use HLL MHD solver"
          solverTypeOriginali=5
          ithis%solverType=4
          ithis%changeSolver=.false.
        case(4)
          print *,"rk2.f03: myid=",myid,"HLL MHD solver encounters negative pressure!!"
          stop
      end select
   endif
enddo

select case(solverTypeOriginali)
    case(5)
      print *,"rk2.f03: myid=",myid,"change back to HLLD MHD solver"
      ithis%solverType=solverTypeOriginali
      solverTypeOriginali=0
    case default
end select

qn(1:nx,1:ny,1:nz,:)=0.5d0*(qn(1:nx,1:ny,1:nz,:)+qn2(1:nx,1:ny,1:nz,:))
call nthis%exchangeBdryMPI(nthis%q,nthis%winq)
call nthis%setBoundary(nthis%q)
nthis%t=nthis%t+nthis%dt

qi(1:nx,1:ny,1:nz,:)=0.5d0*(qi(1:nx,1:ny,1:nz,:)+qi2(1:nx,1:ny,1:nz,:))
call ithis%exchangeBdryMPI(ithis%q,ithis%winq)
call ithis%setBoundary(ithis%q)
ithis%t=ithis%t+ithis%dt

deallocate(totalden)
end subroutine rk2ADsg_3D



subroutine rk2ADsg_3D_zwgV1(nthis,qn,qn1,qn2,ithis,qi,qi1,qi2)
use gridModule
use riemannSolverModule
use mpi
implicit none
type(grid)::nthis,ithis
double precision,dimension(1-nthis%nbuf:nthis%nMesh(1)+nthis%nbuf,1-nthis%nbuf:nthis%nMesh(2)+nthis%nbuf,&
                           1-nthis%nbuf:nthis%nMesh(3)+nthis%nbuf,nthis%nvar)::qn,qn1,qn2,qn0
double precision,dimension(:),allocatable::totalden
double precision,dimension(1-ithis%nbuf:ithis%nMesh(1)+ithis%nbuf,1-ithis%nbuf:ithis%nMesh(2)+ithis%nbuf,&
                           1-ithis%nbuf:ithis%nMesh(3)+ithis%nbuf,ithis%nvar)::qi,qi1,qi2,qi0
procedure(riemannSolver3D),pointer::rieSolvern,rieSolveri
integer::nvar,nx,ny,nz,nbuf,eosTypen,eosTypei,SolverTypen,SolverTypei,solverTypeOriginaln,solverTypeOriginali
double precision::sgfx,sgfy,sgfz,den,momx,momy,momz
double precision::t,dt
integer::ierr,datasize,i,j,k,pos
logical::global_changeSolvern,global_changeSolveri

nbuf=nthis%nbuf
nx=nthis%nMesh(1)
ny=nthis%nMesh(2)
nz=nthis%nMesh(3)
nvar=nthis%nvar
allocate(totalden((nx+2*nbuf)*(ny+2*nbuf)*(nz+2*nbuf)*nvar))
eosTypen=nthis%eosType
eosTypei=ithis%eosType
solverTypen=nthis%solverType
solverTypei=ithis%solverType
solverTypeOriginaln=0
solverTypeOriginali=0

if(eosTypen .eq. 1) then
  select case (solverTypen)
    case (1,2)
      rieSolvern=>solverIso3D
    case default
      print *, "rk2.f03: no suitable 3D solver found!"
      stop
  end select
elseif(eosTypen .eq. 2) then
  select case (solverTypen)
    case (2,3)
      rieSolvern=>solverAdi3D
    case (4,5)
      rieSolvern=>solverAdiMHD3D
    case default
      print *, "rk2.f03: no suitable 3D solver found!"
      stop
  end select
endif

if(eosTypei .eq. 1) then
  select case (solverTypei)
    case (1,2)
      rieSolveri=>solverIso3D
    case default
      print *, "rk2.f03: no suitable 3D solver found!"
      stop
  end select
elseif(eosTypei .eq. 2) then
  select case (solverTypei)
    case (2,3)
      rieSolveri=>solverAdi3D
    case (4,5)
      rieSolveri=>solverAdiMHD3D
    case default
      print *, "rk2.f03: no suitable 3D solver found!"
      stop
  end select
endif

t=nthis%t
datasize=(nthis%nMesh(1)+2*nthis%nbuf)*(nthis%nMesh(2)+2*nthis%nbuf)*(nthis%nMesh(3)+2*nthis%nbuf)



do


!!!!!!!!!!! zwg: step1 of RK2 



    !qn0(1-nthis%nbuf:nthis%nMesh(1)+nthis%nbuf,1-nthis%nbuf:nthis%nMesh(2)+nthis%nbuf,:) &
 !& =qn(1-nthis%nbuf:nthis%nMesh(1)+nthis%nbuf,1-nthis%nbuf:nthis%nMesh(2)+nthis%nbuf,:)

   ! qi0(1-ithis%nbuf:ithis%nMesh(1)+ithis%nbuf,1-ithis%nbuf:ithis%nMesh(2)+ithis%nbuf,:) &
 !& = qi(1-ithis%nbuf:ithis%nMesh(1)+ithis%nbuf,1-ithis%nbuf:ithis%nMesh(2)+ithis%nbuf,:)

   !qn0(1:nx,1:ny,1:nz,:)=qn(1:nx,1:ny,1:nz,:)
   !qi0(1:nx,1:ny,1:nz,:)=qi(1:nx,1:ny,1:nz,:)
   qn0=qn
   qi0=qi


   qn1=qn
   qi1=qi

    !qn1(1-nthis%nbuf:nthis%nMesh(1)+nthis%nbuf,1-nthis%nbuf:nthis%nMesh(2)+nthis%nbuf,:) &
 !& =qn(1-nthis%nbuf:nthis%nMesh(1)+nthis%nbuf,1-nthis%nbuf:nthis%nMesh(2)+nthis%nbuf,:)

    !qi1(1-ithis%nbuf:ithis%nMesh(1)+ithis%nbuf,1-ithis%nbuf:ithis%nMesh(2)+ithis%nbuf,:) &
 !& =qi(1-ithis%nbuf:ithis%nMesh(1)+ithis%nbuf,1-ithis%nbuf:ithis%nMesh(2)+ithis%nbuf,:)

   call evolveAD3D(nthis,ithis,nthis%q1,ithis%q1)

   call nthis%exchangeBdryMPI(nthis%q1,nthis%winq1)
   call nthis%setBoundary(nthis%q1)
   call ithis%exchangeBdryMPI(ithis%q1,ithis%winq1)
   call ithis%setBoundary(ithis%q1)

   qn=qn1
   qi=qi1

    !qn(1-nthis%nbuf:nthis%nMesh(1)+nthis%nbuf,1-nthis%nbuf:nthis%nMesh(2)+nthis%nbuf,:) &
 !& =qn1(1-nthis%nbuf:nthis%nMesh(1)+nthis%nbuf,1-nthis%nbuf:nthis%nMesh(2)+nthis%nbuf,:)

     !qi(1-ithis%nbuf:ithis%nMesh(1)+ithis%nbuf,1-ithis%nbuf:ithis%nMesh(2)+ithis%nbuf,:) &
 !& =qi1(1-ithis%nbuf:ithis%nMesh(1)+ithis%nbuf,1-ithis%nbuf:ithis%nMesh(2)+ithis%nbuf,:)






   call rieSolvern(nthis,qn,qn,qn1,dd=1)
   call rieSolvern(nthis,qn,qn1,qn1,dd=2)
   call rieSolvern(nthis,qn,qn1,qn1,dd=3)

   call rieSolveri(ithis,qi,qi,qi1,dd=1)
   call rieSolveri(ithis,qi,qi1,qi1,dd=2)
   call rieSolveri(ithis,qi,qi1,qi1,dd=3)

   totalden(1:datasize)=nthis%q(1:datasize)+ithis%q(1:datasize)
   call nthis%calcSelfgravity(totalden)
   
   dt=nthis%dt
   !!!!!! apply gravity !!!!!!!
   do k=1,nz
     do j=1,ny
       do i=1,nx
         pos=(k-1)*nx*ny+(j-1)*nx+i
         den = qn(i,j,k,1)
         momx = qn(i,j,k,2)
         momy = qn(i,j,k,3)
         momz = qn(i,j,k,4)
         sgfx = nthis%sgfx(pos)
         sgfy = nthis%sgfy(pos)
         sgfz = nthis%sgfz(pos)

         qn1(i,j,k,2)=qn1(i,j,k,2)+den*sgfx*dt
         qn1(i,j,k,3)=qn1(i,j,k,3)+den*sgfy*dt
         qn1(i,j,k,4)=qn1(i,j,k,4)+den*sgfz*dt
         qn1(i,j,k,5)=qn1(i,j,k,5)+(momx*sgfx+momy*sgfy+momz*sgfz)*dt

         den = qi(i,j,k,1)
         momx = qi(i,j,k,2)
         momy = qi(i,j,k,3)
         momz = qi(i,j,k,4)

         qi1(i,j,k,2)=qi1(i,j,k,2)+den*sgfx*dt
         qi1(i,j,k,3)=qi1(i,j,k,3)+den*sgfy*dt
         qi1(i,j,k,4)=qi1(i,j,k,4)+den*sgfz*dt
         qi1(i,j,k,8)=qi1(i,j,k,8)+(momx*sgfx+momy*sgfy+momz*sgfz)*dt

       enddo
     enddo
   enddo 
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!



   !!!!!!!!!!! zwg: step2  of RK2

   call evolveAD3D(nthis,ithis,nthis%q1,ithis%q1)

   call nthis%exchangeBdryMPI(nthis%q1,nthis%winq1)
   call nthis%setBoundary(nthis%q1)
   call ithis%exchangeBdryMPI(ithis%q1,ithis%winq1)
   call ithis%setBoundary(ithis%q1)


   call rieSolvern(nthis,qn1,qn1,qn2,dd=1)
   call rieSolvern(nthis,qn1,qn2,qn2,dd=2)
   call rieSolvern(nthis,qn1,qn2,qn2,dd=3)

   call rieSolveri(ithis,qi1,qi1,qi2,dd=1)
   call rieSolveri(ithis,qi1,qi2,qi2,dd=2)
   call rieSolveri(ithis,qi1,qi2,qi2,dd=3)

   totalden(1:datasize)=nthis%q1(1:datasize)+ithis%q1(1:datasize)
   call nthis%calcSelfgravity(totalden)
   !!!!!!! apply gravity !!!!!!!
   do k=1,nz
     do j=1,ny
       do i=1,nx
         pos=(k-1)*nx*ny+(j-1)*nx+i
         den = qn1(i,j,k,1)
         momx = qn1(i,j,k,2)
         momy = qn1(i,j,k,3)
         momz = qn1(i,j,k,4)
         sgfx = nthis%sgfx(pos)
         sgfy = nthis%sgfy(pos)
         sgfz = nthis%sgfz(pos)

         qn2(i,j,k,2)=qn2(i,j,k,2)+den*sgfx*dt
         qn2(i,j,k,3)=qn2(i,j,k,3)+den*sgfy*dt
         qn2(i,j,k,4)=qn2(i,j,k,4)+den*sgfz*dt
         qn2(i,j,k,5)=qn2(i,j,k,5)+(momx*sgfx+momy*sgfy+momz*sgfz)*dt

         den = qi1(i,j,k,1)
         momx = qi1(i,j,k,2)
         momy = qi1(i,j,k,3)
         momz = qi1(i,j,k,4)

         qi2(i,j,k,2)=qi2(i,j,k,2)+den*sgfx*dt
         qi2(i,j,k,3)=qi2(i,j,k,3)+den*sgfy*dt
         qi2(i,j,k,4)=qi2(i,j,k,4)+den*sgfz*dt
         qi2(i,j,k,8)=qi2(i,j,k,8)+(momx*sgfx+momy*sgfy+momz*sgfz)*dt

       enddo
     enddo
   enddo
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!

   !!!!  call evolveAD3D(nthis,ithis,nthis%q2,ithis%q2)

   call nthis%exchangeBdryMPI(nthis%q2,nthis%winq2)
   call nthis%setBoundary(nthis%q2)
   call ithis%exchangeBdryMPI(ithis%q2,ithis%winq2)
   call ithis%setBoundary(ithis%q2)

   call MPI_ALLREDUCE(nthis%changeSolver,global_changeSolvern,1,MPI_LOGICAL,MPI_LOR,MPI_COMM_WORLD,ierr)
   call MPI_ALLREDUCE(ithis%changeSolver,global_changeSolveri,1,MPI_LOGICAL,MPI_LOR,MPI_COMM_WORLD,ierr)

   exit
   if(global_changeSolveri .eqv. .false.) then
      exit
   else
      select case(ithis%solverType)
        case(5)
          print *,"rk2.f03: myid",myid,"change to use HLL MHD solver"
          solverTypeOriginali=5
          ithis%solverType=4
          ithis%changeSolver=.false.
        case(4)
          print *,"rk2.f03: myid=",myid,"HLL MHD solver encounters negative pressure!!"
          stop
      end select
   endif
enddo

select case(solverTypeOriginali)
    case(5)
      print *,"rk2.f03: myid=",myid,"change back to HLLD MHD solver"
      ithis%solverType=solverTypeOriginali
      solverTypeOriginali=0
    case default
end select

qn(1:nx,1:ny,1:nz,:)=0.5d0*(qn0(1:nx,1:ny,1:nz,:)+qn2(1:nx,1:ny,1:nz,:))
call nthis%exchangeBdryMPI(nthis%q,nthis%winq)
call nthis%setBoundary(nthis%q)
nthis%t=nthis%t+nthis%dt

qi(1:nx,1:ny,1:nz,:)=0.5d0*(qi0(1:nx,1:ny,1:nz,:)+qi2(1:nx,1:ny,1:nz,:))
call ithis%exchangeBdryMPI(ithis%q,ithis%winq)
call ithis%setBoundary(ithis%q)
ithis%t=ithis%t+ithis%dt

deallocate(totalden)
end subroutine rk2ADsg_3D_zwgV1




!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!!!!!!!!!!! rk2 procedure for 3D MHD with self-gravity!!!!!
!!!subroutine rk2MHDsg_3D(nthis,qn,qn1,qn2,ithis,qi,qi1,qi2)
subroutine rk2MHDsg_3D(ithis,qi,qi1,qi2)
use gridModule
use riemannSolverModule
use mpi
implicit none
!type(grid)::nthis,ithis
type(grid):: ithis

!double precision,dimension(1-nthis%nbuf:nthis%nMesh(1)+nthis%nbuf,1-nthis%nbuf:nthis%nMesh(2)+nthis%nbuf,&
!                           1-nthis%nbuf:nthis%nMesh(3)+nthis%nbuf,nthis%nvar)::qn,qn1,qn2

double precision,dimension(:),allocatable::totalden

double precision,dimension(1-ithis%nbuf:ithis%nMesh(1)+ithis%nbuf,1-ithis%nbuf:ithis%nMesh(2)+ithis%nbuf,&
                           1-ithis%nbuf:ithis%nMesh(3)+ithis%nbuf,ithis%nvar)::qi,qi1,qi2
!procedure(riemannSolver3D),pointer::rieSolvern,rieSolveri
procedure(riemannSolver3D),pointer::rieSolveri

!integer::nvar,nx,ny,nz,nbuf,eosTypen,eosTypei,SolverTypen,SolverTypei,solverTypeOriginaln,solverTypeOriginali
 integer::nvar,nx,ny,nz,nbuf,eosTypei,SolverTypei,solverTypeOriginali

double precision::sgfx,sgfy,sgfz,den,momx,momy,momz

double precision::t,dt
integer::ierr,datasize,i,j,k,pos

!logical::global_changeSolvern,global_changeSolveri
logical::global_changeSolveri

!nbuf=nthis%nbuf
!nx=nthis%nMesh(1)
!ny=nthis%nMesh(2)
!nz=nthis%nMesh(3)
!nvar=nthis%nvar

nbuf=ithis%nbuf
nx=ithis%nMesh(1)
ny=ithis%nMesh(2)
nz=ithis%nMesh(3)
nvar=ithis%nvar

allocate(totalden((nx+2*nbuf)*(ny+2*nbuf)*(nz+2*nbuf)*nvar))

!eosTypen=nthis%eosType
eosTypei=ithis%eosType

!solverTypen=nthis%solverType
solverTypei=ithis%solverType

!solverTypeOriginaln=0
solverTypeOriginali=0

!if(eosTypen .eq. 1) then
  !select case (solverTypen)
   ! case (1,2)
     ! rieSolvern=>solverIso3D
   ! case default
      !print *, "rk2.f03: no suitable 3D solver found!"
      !stop
  !end select
!elseif(eosTypen .eq. 2) then
  !select case (solverTypen)
    !case (2,3)
      !rieSolvern=>solverAdi3D
    !case (4,5)
     ! rieSolvern=>solverAdiMHD3D
    !case default
      !print *, "rk2.f03: no suitable 3D solver found!"
      !stop
  !end select
!endif

if(eosTypei .eq. 1) then
  select case (solverTypei)
    case (4,5)
      rieSolveri=>solverIsoMHD3D
    case default
      print *, "rk2.f03: no suitable 3D solver found!"
      stop
  end select
elseif(eosTypei .eq. 2) then
  select case (solverTypei)
    case (4,5)
      rieSolveri=>solverAdiMHD3D
    case default
      print *, "rk2.f03: no suitable 3D solver found!"
      stop
  end select
endif

!t=nthis%t
 t=ithis%t
!datasize=(nthis%nMesh(1)+2*nthis%nbuf)*(nthis%nMesh(2)+2*nthis%nbuf)*(nthis%nMesh(3)+2*nthis%nbuf)
 datasize=(ithis%nMesh(1)+2*ithis%nbuf)*(ithis%nMesh(2)+2*ithis%nbuf)*(ithis%nMesh(3)+2*ithis%nbuf)

do

   !!! RK2 step1
   !call rieSolvern(nthis,qn,qn,qn1,dd=1)
   !call rieSolvern(nthis,qn,qn1,qn1,dd=2)
   !call rieSolvern(nthis,qn,qn1,qn1,dd=3)

   call rieSolveri(ithis,qi,qi,qi1,dd=1)
   call rieSolveri(ithis,qi,qi1,qi1,dd=2)
   call rieSolveri(ithis,qi,qi1,qi1,dd=3)

   !totalden(1:datasize)=nthis%q(1:datasize)+ithis%q(1:datasize)
    totalden(1:datasize)=ithis%q(1:datasize)
   !call nthis%calcSelfgravity(totalden)
    call ithis%calcSelfgravity(totalden)
   
   !dt=nthis%dt

   dt=ithis%dt

   !!!!!! apply gravity !!!!!!!
   do k=1,nz
     do j=1,ny
       do i=1,nx
         pos=(k-1)*nx*ny+(j-1)*nx+i

         !den = qn(i,j,k,1)
         !momx = qn(i,j,k,2)
         !momy = qn(i,j,k,3)
         !momz = qn(i,j,k,4)

         !sgfx = nthis%sgfx(pos)
         !sgfy = nthis%sgfy(pos)
         !sgfz = nthis%sgfz(pos)

         !qn1(i,j,k,2)=qn1(i,j,k,2)+den*sgfx*dt
         !qn1(i,j,k,3)=qn1(i,j,k,3)+den*sgfy*dt
         !qn1(i,j,k,4)=qn1(i,j,k,4)+den*sgfz*dt
         !qn1(i,j,k,5)=qn1(i,j,k,5)+(momx*sgfx+momy*sgfy+momz*sgfz)*dt



         sgfx = ithis%sgfx(pos)
         sgfy = ithis%sgfy(pos)
         sgfz = ithis%sgfz(pos)

         den = qi(i,j,k,1)
         momx = qi(i,j,k,2)
         momy = qi(i,j,k,3)
         momz = qi(i,j,k,4)

         qi1(i,j,k,2)=qi1(i,j,k,2)+den*sgfx*dt
         qi1(i,j,k,3)=qi1(i,j,k,3)+den*sgfy*dt
         qi1(i,j,k,4)=qi1(i,j,k,4)+den*sgfz*dt
         qi1(i,j,k,8)=qi1(i,j,k,8)+(momx*sgfx+momy*sgfy+momz*sgfz)*dt

       enddo
     enddo
   enddo 
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !call evolveAD3D(nthis,ithis,nthis%q1,ithis%q1)

   !call nthis%exchangeBdryMPI(nthis%q1,nthis%winq1)
   !call nthis%setBoundary(nthis%q1)

   call ithis%exchangeBdryMPI(ithis%q1,ithis%winq1)
   call ithis%setBoundary(ithis%q1)


   !!! RK2 step2

   !call rieSolvern(nthis,qn1,qn1,qn2,dd=1)
   !call rieSolvern(nthis,qn1,qn2,qn2,dd=2)
   !call rieSolvern(nthis,qn1,qn2,qn2,dd=3)

   call rieSolveri(ithis,qi1,qi1,qi2,dd=1)
   call rieSolveri(ithis,qi1,qi2,qi2,dd=2)
   call rieSolveri(ithis,qi1,qi2,qi2,dd=3)

   !totalden(1:datasize)=nthis%q1(1:datasize)+ithis%q1(1:datasize)
    totalden(1:datasize)=ithis%q1(1:datasize)

   !call nthis%calcSelfgravity(totalden)
    call ithis%calcSelfgravity(totalden)

   !!!!!!! apply gravity !!!!!!!
   do k=1,nz
     do j=1,ny
       do i=1,nx
         !pos=(k-1)*nx*ny+(j-1)*nx+i
         !den = qn1(i,j,k,1)
         !momx = qn1(i,j,k,2)
         !momy = qn1(i,j,k,3)
         !momz = qn1(i,j,k,4)
         !sgfx = nthis%sgfx(pos)
         !sgfy = nthis%sgfy(pos)
         !sgfz = nthis%sgfz(pos)

         !qn2(i,j,k,2)=qn2(i,j,k,2)+den*sgfx*dt
         !qn2(i,j,k,3)=qn2(i,j,k,3)+den*sgfy*dt
         !qn2(i,j,k,4)=qn2(i,j,k,4)+den*sgfz*dt
         !qn2(i,j,k,5)=qn2(i,j,k,5)+(momx*sgfx+momy*sgfy+momz*sgfz)*dt
 
         pos=(k-1)*nx*ny+(j-1)*nx+i
         sgfx = ithis%sgfx(pos)
         sgfy = ithis%sgfy(pos)
         sgfz = ithis%sgfz(pos)

         den = qi1(i,j,k,1)
         momx = qi1(i,j,k,2)
         momy = qi1(i,j,k,3)
         momz = qi1(i,j,k,4)

         qi2(i,j,k,2)=qi2(i,j,k,2)+den*sgfx*dt
         qi2(i,j,k,3)=qi2(i,j,k,3)+den*sgfy*dt
         qi2(i,j,k,4)=qi2(i,j,k,4)+den*sgfz*dt
         qi2(i,j,k,8)=qi2(i,j,k,8)+(momx*sgfx+momy*sgfy+momz*sgfz)*dt

       enddo
     enddo
   enddo
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!

   !call evolveAD3D(nthis,ithis,nthis%q2,ithis%q2)

   !call nthis%exchangeBdryMPI(nthis%q2,nthis%winq2)
   !call nthis%setBoundary(nthis%q2)

   call ithis%exchangeBdryMPI(ithis%q2,ithis%winq2)
   call ithis%setBoundary(ithis%q2)

   !call MPI_ALLREDUCE(nthis%changeSolver,global_changeSolvern,1,MPI_LOGICAL,MPI_LOR,MPI_COMM_WORLD,ierr)
   call MPI_ALLREDUCE(ithis%changeSolver,global_changeSolveri,1,MPI_LOGICAL,MPI_LOR,MPI_COMM_WORLD,ierr)

   exit
   if(global_changeSolveri .eqv. .false.) then
      exit
   else
      select case(ithis%solverType)
        case(5)
          print *,"rk2.f03: myid",myid,"change to use HLL MHD solver"
          solverTypeOriginali=5
          ithis%solverType=4
          ithis%changeSolver=.false.
        case(4)
          print *,"rk2.f03: myid=",myid,"HLL MHD solver encounters negative pressure!!"
          stop
      end select
   endif
enddo

select case(solverTypeOriginali)
    case(5)
      print *,"rk2.f03: myid=",myid,"change back to HLLD MHD solver"
      ithis%solverType=solverTypeOriginali
      solverTypeOriginali=0
    case default
end select

!qn(1:nx,1:ny,1:nz,:)=0.5d0*(qn(1:nx,1:ny,1:nz,:)+qn2(1:nx,1:ny,1:nz,:))
!call nthis%exchangeBdryMPI(nthis%q,nthis%winq)
!call nthis%setBoundary(nthis%q)
!nthis%t=nthis%t+nthis%dt

qi(1:nx,1:ny,1:nz,:)=0.5d0*(qi(1:nx,1:ny,1:nz,:)+qi2(1:nx,1:ny,1:nz,:))
call ithis%exchangeBdryMPI(ithis%q,ithis%winq)
call ithis%setBoundary(ithis%q)
ithis%t=ithis%t+ithis%dt

deallocate(totalden)
end subroutine rk2MHDsg_3D



!!!!!!!!!! rk2 procedure for 3D ambipolar diffusion !!!!!
subroutine rk2AD_3D(nthis,qn,qn1,qn2,ithis,qi,qi1,qi2)
use gridModule
use riemannSolverModule
use mpi
implicit none
type(grid)::nthis,ithis
double precision,dimension(1-nthis%nbuf:nthis%nMesh(1)+nthis%nbuf,1-nthis%nbuf:nthis%nMesh(2)+nthis%nbuf,&
                           1-nthis%nbuf:nthis%nMesh(3)+nthis%nbuf,nthis%nvar)::qn,qn1,qn2

double precision,dimension(1-ithis%nbuf:ithis%nMesh(1)+ithis%nbuf,1-ithis%nbuf:ithis%nMesh(2)+ithis%nbuf,&
                           1-ithis%nbuf:ithis%nMesh(3)+ithis%nbuf,ithis%nvar)::qi,qi1,qi2
procedure(riemannSolver3D),pointer::rieSolvern,rieSolveri
integer::nvar,nx,ny,nz,nbuf,eosTypen,eosTypei,SolverTypen,SolverTypei,solverTypeOriginaln,solverTypeOriginali
double precision::t,dt
integer::ierr,datasize,i,j,k,pos
logical::global_changeSolvern,global_changeSolveri

nbuf=nthis%nbuf
nx=nthis%nMesh(1)
ny=nthis%nMesh(2)
nz=nthis%nMesh(3)
nvar=nthis%nvar
eosTypen=nthis%eosType
eosTypei=ithis%eosType
solverTypen=nthis%solverType
solverTypei=ithis%solverType
solverTypeOriginaln=0
solverTypeOriginali=0

if(eosTypen .eq. 1) then
  select case (solverTypen)
    case (1,2)
      rieSolvern=>solverIso3D
    case default
      print *, "rk2.f03: no suitable 3D solver found!"
      stop
  end select
elseif(eosTypen .eq. 2) then
  select case (solverTypen)
    case (2,3)
      rieSolvern=>solverAdi3D
    case (4,5)
      rieSolvern=>solverAdiMHD3D
    case default
      print *, "rk2.f03: no suitable 3D solver found!"
      stop
  end select
endif

if(eosTypei .eq. 1) then
  select case (solverTypei)
    case (1,2)
      rieSolveri=>solverIso3D
    case default
      print *, "rk2.f03: no suitable 3D solver found!"
      stop
  end select
elseif(eosTypei .eq. 2) then
  select case (solverTypei)
    case (2,3)
      rieSolveri=>solverAdi3D
    case (4,5)
      rieSolveri=>solverAdiMHD3D
    case default
      print *, "rk2.f03: no suitable 3D solver found!"
      stop
  end select
endif

t=nthis%t

do
   call rieSolvern(nthis,qn,qn,qn1,dd=1)
   call rieSolvern(nthis,qn,qn1,qn1,dd=2)
   call rieSolvern(nthis,qn,qn1,qn1,dd=3)

   call rieSolveri(ithis,qi,qi,qi1,dd=1)
   call rieSolveri(ithis,qi,qi1,qi1,dd=2)
   call rieSolveri(ithis,qi,qi1,qi1,dd=3)

   call evolveAD3D(nthis,ithis,nthis%q1,ithis%q1)

   call nthis%exchangeBdryMPI(nthis%q1,nthis%winq1)
   call nthis%setBoundary(nthis%q1)
   call ithis%exchangeBdryMPI(ithis%q1,ithis%winq1)
   call ithis%setBoundary(ithis%q1)

   call rieSolvern(nthis,qn1,qn1,qn2,dd=1)
   call rieSolvern(nthis,qn1,qn2,qn2,dd=2)
   call rieSolvern(nthis,qn1,qn2,qn2,dd=3)

   call rieSolveri(ithis,qi1,qi1,qi2,dd=1)
   call rieSolveri(ithis,qi1,qi2,qi2,dd=2)
   call rieSolveri(ithis,qi1,qi2,qi2,dd=3)   

   call evolveAD3D(nthis,ithis,nthis%q2,ithis%q2)

   call nthis%exchangeBdryMPI(nthis%q2,nthis%winq2)
   call nthis%setBoundary(nthis%q2)
   call ithis%exchangeBdryMPI(ithis%q2,ithis%winq2)
   call ithis%setBoundary(ithis%q2)

   call MPI_ALLREDUCE(nthis%changeSolver,global_changeSolvern,1,MPI_LOGICAL,MPI_LOR,MPI_COMM_WORLD,ierr)
   call MPI_ALLREDUCE(ithis%changeSolver,global_changeSolveri,1,MPI_LOGICAL,MPI_LOR,MPI_COMM_WORLD,ierr)

   exit
   if(global_changeSolveri .eqv. .false.) then
      exit
   else
      select case(ithis%solverType)
        case(5)
          print *,"rk2.f03: myid",myid,"change to use HLL MHD solver"
          solverTypeOriginali=5
          ithis%solverType=4
          ithis%changeSolver=.false.
        case(4)
          print *,"rk2.f03: myid=",myid,"HLL MHD solver encounters negative pressure!!"
          stop
      end select
   endif
enddo

select case(solverTypeOriginali)
    case(5)
      print *,"rk2.f03: myid=",myid,"change back to HLLD MHD solver"
      ithis%solverType=solverTypeOriginali
      solverTypeOriginali=0
    case default
end select

qn(1:nx,1:ny,1:nz,:)=0.5d0*(qn(1:nx,1:ny,1:nz,:)+qn2(1:nx,1:ny,1:nz,:))
call nthis%exchangeBdryMPI(nthis%q,nthis%winq)
call nthis%setBoundary(nthis%q)
nthis%t=nthis%t+nthis%dt

qi(1:nx,1:ny,1:nz,:)=0.5d0*(qi(1:nx,1:ny,1:nz,:)+qi2(1:nx,1:ny,1:nz,:))
call ithis%exchangeBdryMPI(ithis%q,ithis%winq)
call ithis%setBoundary(ithis%q)
ithis%t=ithis%t+ithis%dt  

end subroutine rk2AD_3D


subroutine rk2AD_3D_HSHSMD(nthis,qn,qn1,qn2,ithis,qi,qi1,qi2)
use gridModule
use riemannSolverModule
use mpi
implicit none
type(grid)::nthis,ithis
double precision,dimension(1-nthis%nbuf:nthis%nMesh(1)+nthis%nbuf,1-nthis%nbuf:nthis%nMesh(2)+nthis%nbuf,&
                           1-nthis%nbuf:nthis%nMesh(3)+nthis%nbuf,nthis%nvar)::qn,qn1,qn2

double precision,dimension(1-ithis%nbuf:ithis%nMesh(1)+ithis%nbuf,1-ithis%nbuf:ithis%nMesh(2)+ithis%nbuf,&
                           1-ithis%nbuf:ithis%nMesh(3)+ithis%nbuf,ithis%nvar)::qi,qi1,qi2
procedure(riemannSolver3D),pointer::rieSolvern,rieSolveri
integer::nvar,nx,ny,nz,nbuf,eosTypen,eosTypei,SolverTypen,SolverTypei,solverTypeOriginaln,solverTypeOriginali
double precision::t,dt
integer::ierr,datasize,i,j,k,pos
logical::global_changeSolvern,global_changeSolveri

nbuf=nthis%nbuf
nx=nthis%nMesh(1)
ny=nthis%nMesh(2)
nz=nthis%nMesh(3)
nvar=nthis%nvar
eosTypen=nthis%eosType
eosTypei=ithis%eosType
solverTypen=nthis%solverType
solverTypei=ithis%solverType
solverTypeOriginaln=0
solverTypeOriginali=0

if(eosTypen .eq. 1) then
  select case (solverTypen)
    case (1,2)
      rieSolvern=>solverIso3D
    case default
      print *, "rk2.f03: no suitable 3D solver found!"
      stop
  end select
elseif(eosTypen .eq. 2) then
  select case (solverTypen)
    case (2,3)
      rieSolvern=>solverAdi3D
    case (4,5)
      rieSolvern=>solverAdiMHD3D
    case default
      print *, "rk2.f03: no suitable 3D solver found!"
      stop
  end select
endif

if(eosTypei .eq. 1) then
  select case (solverTypei)
    case (1,2)
      rieSolveri=>solverIso3D
    case default
      print *, "rk2.f03: no suitable 3D solver found!"
      stop
  end select
elseif(eosTypei .eq. 2) then
  select case (solverTypei)
    case (2,3)
      rieSolveri=>solverAdi3D
    case (4,5)
      rieSolveri=>solverAdiMHD3D
    case default
      print *, "rk2.f03: no suitable 3D solver found!"
      stop
  end select
endif

t=nthis%t

do
   call rieSolvern(nthis,qn,qn,qn1,dd=1)
   call rieSolvern(nthis,qn,qn1,qn1,dd=2)
   call rieSolvern(nthis,qn,qn1,qn1,dd=3)

   call rieSolveri(ithis,qi,qi,qi1,dd=1)
   call rieSolveri(ithis,qi,qi1,qi1,dd=2)
   call rieSolveri(ithis,qi,qi1,qi1,dd=3)
    
   !!!!!!!!!!!zwg modified!!!!!!!!!!!!!!!!!!!!
   !!!call evolveAD3D(nthis,ithis,nthis%q1,ithis%q1)
      call evolveAD3D_MD(nthis,ithis,nthis%q,ithis%q,nthis%q1,ithis%q1)

   call nthis%exchangeBdryMPI(nthis%q1,nthis%winq1)
   call nthis%setBoundary(nthis%q1)
   call ithis%exchangeBdryMPI(ithis%q1,ithis%winq1)
   call ithis%setBoundary(ithis%q1)

   call rieSolvern(nthis,qn1,qn1,qn2,dd=1)
   call rieSolvern(nthis,qn1,qn2,qn2,dd=2)
   call rieSolvern(nthis,qn1,qn2,qn2,dd=3)

   call rieSolveri(ithis,qi1,qi1,qi2,dd=1)
   call rieSolveri(ithis,qi1,qi2,qi2,dd=2)
   call rieSolveri(ithis,qi1,qi2,qi2,dd=3)   

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !!!call evolveAD3D(nthis,ithis,nthis%q2,ithis%q2)
      call evolveAD3D_MD(nthis,ithis,nthis%q1,ithis%q1,nthis%q2,ithis%q2)

   call nthis%exchangeBdryMPI(nthis%q2,nthis%winq2)
   call nthis%setBoundary(nthis%q2)
   call ithis%exchangeBdryMPI(ithis%q2,ithis%winq2)
   call ithis%setBoundary(ithis%q2)

   call MPI_ALLREDUCE(nthis%changeSolver,global_changeSolvern,1,MPI_LOGICAL,MPI_LOR,MPI_COMM_WORLD,ierr)
   call MPI_ALLREDUCE(ithis%changeSolver,global_changeSolveri,1,MPI_LOGICAL,MPI_LOR,MPI_COMM_WORLD,ierr)

   exit
   if(global_changeSolveri .eqv. .false.) then
      exit
   else
      select case(ithis%solverType)
        case(5)
          print *,"rk2.f03: myid",myid,"change to use HLL MHD solver"
          solverTypeOriginali=5
          ithis%solverType=4
          ithis%changeSolver=.false.
        case(4)
          print *,"rk2.f03: myid=",myid,"HLL MHD solver encounters negative pressure!!"
          stop
      end select
   endif
enddo

select case(solverTypeOriginali)
    case(5)
      print *,"rk2.f03: myid=",myid,"change back to HLLD MHD solver"
      ithis%solverType=solverTypeOriginali
      solverTypeOriginali=0
    case default
end select

qn(1:nx,1:ny,1:nz,:)=0.5d0*(qn(1:nx,1:ny,1:nz,:)+qn2(1:nx,1:ny,1:nz,:))
call nthis%exchangeBdryMPI(nthis%q,nthis%winq)
call nthis%setBoundary(nthis%q)
nthis%t=nthis%t+nthis%dt

qi(1:nx,1:ny,1:nz,:)=0.5d0*(qi(1:nx,1:ny,1:nz,:)+qi2(1:nx,1:ny,1:nz,:))
call ithis%exchangeBdryMPI(ithis%q,ithis%winq)
call ithis%setBoundary(ithis%q)
ithis%t=ithis%t+ithis%dt  

end subroutine rk2AD_3D_HSHSMD





!!!!!!!!!!! rk2 procedure for 3D MHD (without sg)!!!!!
!!!subroutine rk2MHD_3D(nthis,qn,qn1,qn2,ithis,qi,qi1,qi2)
subroutine rk2MHD_3D(ithis,qi,qi1,qi2)
use gridModule
use riemannSolverModule
use mpi
implicit none
!type(grid)::nthis,ithis
type(grid):: ithis

!double precision,dimension(1-nthis%nbuf:nthis%nMesh(1)+nthis%nbuf,1-nthis%nbuf:nthis%nMesh(2)+nthis%nbuf,&
!                           1-nthis%nbuf:nthis%nMesh(3)+nthis%nbuf,nthis%nvar)::qn,qn1,qn2

double precision,dimension(1-ithis%nbuf:ithis%nMesh(1)+ithis%nbuf,1-ithis%nbuf:ithis%nMesh(2)+ithis%nbuf,&
                           1-ithis%nbuf:ithis%nMesh(3)+ithis%nbuf,ithis%nvar)::qi,qi1,qi2
!procedure(riemannSolver3D),pointer::rieSolvern,rieSolveri
procedure(riemannSolver3D),pointer::rieSolveri

!integer::nvar,nx,ny,nz,nbuf,eosTypen,eosTypei,SolverTypen,SolverTypei,solverTypeOriginaln,solverTypeOriginali
 integer::nvar,nx,ny,nz,nbuf,eosTypei,SolverTypei,solverTypeOriginali
double precision::t,dt
integer::ierr,datasize,i,j,k,pos
!logical::global_changeSolvern,global_changeSolveri
 logical::global_changeSolveri

nbuf=ithis%nbuf
nx=ithis%nMesh(1)
ny=ithis%nMesh(2)
nz=ithis%nMesh(3)
nvar=ithis%nvar

!eosTypen=nthis%eosType
eosTypei=ithis%eosType
!solverTypen=nthis%solverType
solverTypei=ithis%solverType

!solverTypeOriginaln=0
solverTypeOriginali=0

!if(eosTypen .eq. 1) then
  !select case (solverTypen)
    !case (1,2)
      !rieSolvern=>solverIso3D
    !case default
      !print *, "rk2.f03: no suitable 3D solver found!"
      !stop
  !end select
!elseif(eosTypen .eq. 2) then
  !select case (solverTypen)
    !case (2,3)
      !rieSolvern=>solverAdi3D
    !case (4,5)
      !rieSolvern=>solverAdiMHD3D
    !case default
      !print *, "rk2.f03: no suitable 3D solver found!"
      !stop
  !end select
!endif

if(eosTypei .eq. 1) then
  select case (solverTypei)
    case (4,5)
      rieSolveri=>solverIsoMHD3D
    case default
      print *, "rk2.f03: no suitable 3D solver found!"
      stop
  end select
elseif(eosTypei .eq. 2) then
  select case (solverTypei)
    case (4,5)
      rieSolveri=>solverAdiMHD3D
    case default
      print *, "rk2.f03: no suitable 3D solver found!"
      stop
  end select
endif




t=ithis%t

do

   !!! RK2 step1
   !call rieSolvern(nthis,qn,qn,qn1,dd=1)
   !call rieSolvern(nthis,qn,qn1,qn1,dd=2)
   !call rieSolvern(nthis,qn,qn1,qn1,dd=3)

   call rieSolveri(ithis,qi,qi,qi1,dd=1)
   call rieSolveri(ithis,qi,qi1,qi1,dd=2)
   call rieSolveri(ithis,qi,qi1,qi1,dd=3)

   !call evolveAD3D(nthis,ithis,nthis%q1,ithis%q1)

   !call nthis%exchangeBdryMPI(nthis%q1,nthis%winq1)
   !call nthis%setBoundary(nthis%q1)

   call ithis%exchangeBdryMPI(ithis%q1,ithis%winq1)
   call ithis%setBoundary(ithis%q1)


  !!! RK2 step2
   !call rieSolvern(nthis,qn1,qn1,qn2,dd=1)
   !call rieSolvern(nthis,qn1,qn2,qn2,dd=2)
   !call rieSolvern(nthis,qn1,qn2,qn2,dd=3)

   call rieSolveri(ithis,qi1,qi1,qi2,dd=1)
   call rieSolveri(ithis,qi1,qi2,qi2,dd=2)
   call rieSolveri(ithis,qi1,qi2,qi2,dd=3)   

   !call evolveAD3D(nthis,ithis,nthis%q2,ithis%q2)

   !call nthis%exchangeBdryMPI(nthis%q2,nthis%winq2)
   !call nthis%setBoundary(nthis%q2)

   call ithis%exchangeBdryMPI(ithis%q2,ithis%winq2)
   call ithis%setBoundary(ithis%q2)

   !call MPI_ALLREDUCE(nthis%changeSolver,global_changeSolvern,1,MPI_LOGICAL,MPI_LOR,MPI_COMM_WORLD,ierr)
   call MPI_ALLREDUCE(ithis%changeSolver,global_changeSolveri,1,MPI_LOGICAL,MPI_LOR,MPI_COMM_WORLD,ierr)

   exit
   if(global_changeSolveri .eqv. .false.) then
      exit
   else
      select case(ithis%solverType)
        case(5)
          print *,"rk2.f03: myid",myid,"change to use HLL MHD solver"
          solverTypeOriginali=5
          ithis%solverType=4
          ithis%changeSolver=.false.
        case(4)
          print *,"rk2.f03: myid=",myid,"HLL MHD solver encounters negative pressure!!"
          stop
      end select
   endif
enddo

select case(solverTypeOriginali)
    case(5)
      print *,"rk2.f03: myid=",myid,"change back to HLLD MHD solver"
      ithis%solverType=solverTypeOriginali
      solverTypeOriginali=0
    case default
end select

!qn(1:nx,1:ny,1:nz,:)=0.5d0*(qn(1:nx,1:ny,1:nz,:)+qn2(1:nx,1:ny,1:nz,:))
!call nthis%exchangeBdryMPI(nthis%q,nthis%winq)
!call nthis%setBoundary(nthis%q)
!nthis%t=nthis%t+nthis%dt

qi(1:nx,1:ny,1:nz,:)=0.5d0*(qi(1:nx,1:ny,1:nz,:)+qi2(1:nx,1:ny,1:nz,:))
call ithis%exchangeBdryMPI(ithis%q,ithis%winq)
call ithis%setBoundary(ithis%q)
ithis%t=ithis%t+ithis%dt  


end subroutine rk2MHD_3D







subroutine rk2_3D(this,q,q1,q2)
use gridModule
use riemannSolverModule
use mpi
implicit none
class(grid)::this
double precision,dimension(1-this%nbuf:this%nMesh(1)+this%nbuf,1-this%nbuf:this%nMesh(2)+this%nbuf,&
&                          1-this%nbuf:this%nMesh(3)+this%nbuf,this%nvar)::q,q1,q2
procedure(riemannSolver3D),pointer::rieSolver
integer::nvar,nx,ny,nz,nbuf,ndim,eosType,solverType,solverTypeOriginal
double precision::t,dt,pressure,den,vx,vy,vz,bx,by,bz,gam
integer::ierr,i,j,k,pos
logical::global_changeSolver,global_neg_pressure
!! eostype: 1=isothermal, 2=adiabatic
!! solverType: 1=exactHD, 2=HLLHD, 3=HLLCHD, 4=HLLMHD, 5=HLLDMHD)
eostype=this%eosType
solverType=this%solverType
solverTypeOriginal=0
if(eosType .eq. 1) then
  select case (solverType)
     case (1,2)
        rieSolver=>solverIso3D
     case (4,5)
        rieSolver=>solverIsoMHD3D
     case default
        print *, "rk2.f03: no suitable 3D solver found!"
        stop
  end select
elseif(eosType .eq. 2) then
  select case (solverType)
     case (2,3)
        rieSolver=>solverAdi3D
     case (4,5)
        rieSolver=>solverAdiMHD3D
     case default
        print *, "rk2.f03: no suitable 3D solver found!"
        stop
  end select 
endif

nvar=this%nvar
nx=this%nMesh(1)
ny=this%nMesh(2)
nz=this%nMesh(3)
nbuf=this%nbuf
t=this%t
!!!!!!! evolve one time step using RK2
!!! RK1
do 
   call rieSolver(this,q,q,q1,dd=1)
   call rieSolver(this,q,q1,q1,dd=2)
   call rieSolver(this,q,q1,q1,dd=3)

   !!!!! added on Feb. 22, 2017 by hhwang
   if(this%enable_sg) then
     call this%calcSelfgravity(this%q)
     dt=this%dt
     if(eostype .eq. 1) then !! if isothermal HD gas
      select case (solverType)
        case(1,2,4,5)
           do k=1,nz
             do j=1,ny
               do i=1,nx
                 pos=(k-1)*nx*ny+(j-1)*nx+i
                 q1(i,j,k,2)=q1(i,j,k,2)+q(i,j,k,1)*this%sgfx(pos)*dt !! momx
                 q1(i,j,k,3)=q1(i,j,k,3)+q(i,j,k,1)*this%sgfy(pos)*dt !! momy
                 q1(i,j,k,4)=q1(i,j,k,4)+q(i,j,k,1)*this%sgfz(pos)*dt !! momz
               enddo
             enddo
           enddo
        case default
          print *,'rk2_3D: cannot recognize solver'
          stop          
      end select
     elseif(eostype .eq. 2) then !! adiabatic
       if(solverType .eq. 2 .or. solverType .eq. 3) then !! HD solvers
         do k=1,nz
           do j=1,ny
             do i=1,nx
               pos=(k-1)*nx*ny+(j-1)*nx+i
               q1(i,j,k,2)=q1(i,j,k,2)+q(i,j,k,1)*this%sgfx(pos)*dt !! momx
               q1(i,j,k,3)=q1(i,j,k,3)+q(i,j,k,1)*this%sgfy(pos)*dt !! momy
               q1(i,j,k,4)=q1(i,j,k,4)+q(i,j,k,1)*this%sgfz(pos)*dt !! momz
               q1(i,j,k,5)=q1(i,j,k,5)+(q(i,j,k,2)*this%sgfx(pos)&
             +q(i,j,k,3)*this%sgfy(pos)+q(i,j,k,4)*this%sgfz(pos))*dt !! energy
             enddo
           enddo
         enddo
       elseif(solverType .eq. 4 .or. solverType .eq. 5) then !! MHD solvers 
         do k=1,nz
           do j=1,ny
             do i=1,nx
               pos=(k-1)*nx*ny+(j-1)*nx+i
               q1(i,j,k,2)=q1(i,j,k,2)+q(i,j,k,1)*this%sgfx(pos)*dt !! momx
               q1(i,j,k,3)=q1(i,j,k,3)+q(i,j,k,1)*this%sgfy(pos)*dt !! momy
               q1(i,j,k,4)=q1(i,j,k,4)+q(i,j,k,1)*this%sgfz(pos)*dt !! momz
               q1(i,j,k,8)=q1(i,j,k,8)+(q(i,j,k,2)*this%sgfx(pos)&
             +q(i,j,k,3)*this%sgfy(pos)+q(i,j,k,4)*this%sgfz(pos))*dt !! energy
             enddo
           enddo
         enddo
       endif
     endif !! end eostype
     
   endif !! end enable_sg: added by hhwang on Feb. 22, 2017
   
   call this%exchangeBdryMPI(this%q1,this%winq1)
   call this%setBoundary(this%q1)
!!! RK2
   call rieSolver(this,q1,q1,q2,dd=1)
   call rieSolver(this,q1,q2,q2,dd=2)
   call rieSolver(this,q1,q2,q2,dd=3)

   !!!!! added on Feb. 22, 2017 by hhwang
   if(this%enable_sg) then

     call this%calcSelfgravity(this%q1)
     dt=this%dt
     if(eostype .eq. 1) then !! if isothermal HD gas
       select case (solverType)
         case(1,2,4,5)
           do k=1,nz
             do j=1,ny
               do i=1,nx
                 pos=(k-1)*nx*ny+(j-1)*nx+i
                 q2(i,j,k,2)=q2(i,j,k,2)+q1(i,j,k,1)*this%sgfx(pos)*dt !! momx
                 q2(i,j,k,3)=q2(i,j,k,3)+q1(i,j,k,1)*this%sgfy(pos)*dt !! momy
                 q2(i,j,k,4)=q2(i,j,k,4)+q1(i,j,k,1)*this%sgfz(pos)*dt !! momz
               enddo
             enddo
           enddo
         case default
          print *,'rk2_3D: cannot recognize solver'
          stop
       end select
     elseif(eostype .eq. 2) then !! adiabatic
       if(solverType .eq. 2 .or. solverType .eq. 3) then !! HD solvers
         do k=1,nz
           do j=1,ny
             do i=1,nx
               pos=(k-1)*nx*ny+(j-1)*nx+i
               q2(i,j,k,2)=q2(i,j,k,2)+q1(i,j,k,1)*this%sgfx(pos)*dt !! momx
               q2(i,j,k,3)=q2(i,j,k,3)+q1(i,j,k,1)*this%sgfy(pos)*dt !! momy
               q2(i,j,k,4)=q2(i,j,k,4)+q1(i,j,k,1)*this%sgfz(pos)*dt !! momz
               q2(i,j,k,5)=q2(i,j,k,5)+(q1(i,j,k,2)*this%sgfx(pos)&
            +q1(i,j,k,3)*this%sgfy(pos)+q1(i,j,k,4)*this%sgfz(pos))*dt !! energy
             enddo
           enddo
         enddo
       elseif(solverType .eq. 4 .or. solverType .eq. 5) then !! MHD solvers 
         do k=1,nz
           do j=1,ny
             do i=1,nx
               pos=(k-1)*nx*ny+(j-1)*nx+i
               q2(i,j,k,2)=q2(i,j,k,2)+q1(i,j,k,1)*this%sgfx(pos)*dt !! momx
               q2(i,j,k,3)=q2(i,j,k,3)+q1(i,j,k,1)*this%sgfy(pos)*dt !! momy
               q2(i,j,k,4)=q2(i,j,k,4)+q1(i,j,k,1)*this%sgfz(pos)*dt !! momz
               q2(i,j,k,8)=q2(i,j,k,8)+(q1(i,j,k,2)*this%sgfx(pos)&
            +q1(i,j,k,3)*this%sgfy(pos)+q1(i,j,k,4)*this%sgfz(pos))*dt !! energy
             enddo
           enddo
         enddo
       endif
     endif !! end eostype

   endif !! end enable_sg: added by hhwang on Feb. 22, 2017

   !!!!! check if negative pressure, added by hhwang Feb. 27, 2017 !!!!!
   pressure=1.d0
   gam=this%adiGamma
   if(eostype .eq. 2) then  !! adiabatic
     if(solverType .eq. 2 .or. solverType .eq. 3) then !! HD solvers
       do k=1,nz
         do j=1,ny
           do i=1,nx
             den=q2(i,j,k,1)
              vx=q2(i,j,k,2)/den
              vy=q2(i,j,k,3)/den
              vz=q2(i,j,k,4)/den
             pressure=dmin1(pressure,(gam-1.d0)*(q2(i,j,k,5)-0.5d0*den*(vx**2+vy**2+vz**2)))
           enddo
         enddo
       enddo 
     endif  !! solverType 2 or 3
     if(solverType .eq. 4 .or. solverType .eq. 5) then !! MHD solvers
       do k=1,nz
         do j=1,ny
           do i=1,nx
              den=q2(i,j,k,1)
              vx=q2(i,j,k,2)/den
              vy=q2(i,j,k,3)/den
              vz=q2(i,j,k,4)/den
              bx=0.5d0*(q2(i,j,k,5)+q2(i,j,k, 9))
              by=0.5d0*(q2(i,j,k,6)+q2(i,j,k,10))
              bz=0.5d0*(q2(i,j,k,7)+q2(i,j,k,11))
              pressure=dmin1(pressure,(gam-1.d0)*(q2(i,j,k,8)-0.5d0*den*(vx**2+vy**2+vz**2)&
                                                             -0.5d0*(bx**2+by**2+bz**2)))
           enddo
         enddo
       enddo       
     endif !! solverType 4 or 5    
   endif !! end eosType .eq. 2
   if(pressure .lt. 0.d0) then
     this%neg_pressure=.true.
   endif
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   call this%exchangeBdryMPI(this%q2,this%winq2)
   call this%setBoundary(this%q2)

   call MPI_ALLREDUCE(this%changeSolver,global_changeSolver,1,MPI_LOGICAL,MPI_LOR,MPI_COMM_WORLD,ierr)
   call MPI_ALLREDUCE(this%neg_pressure,global_neg_pressure,1,MPI_LOGICAL,MPI_LOR,MPI_COMM_WORLD,ierr)
   !!! check the global status of changeSolver
  
   if(global_changeSolver .eqv. .false.) then
      if(global_neg_pressure .eqv. .false.) then
        exit
      else
        print *,"negative pressure..reduce dt from ", this%dt, "to", this%dt*0.5d0
        print *,"negative pressure:", pressure
        this%dt=0.5d0*this%dt
        this%neg_pressure=.false.
      endif      
   else
      select case(this%solverType)
        case(5)
          print *,"rk2.f03: myid",myid,"change to use HLL MHD solver"
          solverTypeOriginal=5
          this%solverType=4
          this%changeSolver=.false.
        case(4)
          print *,"rk2.f03: myid=",myid,"HLL MHD solver encounters negative pressure!!"
          stop
      end select
   endif   

enddo

   select case(solverTypeOriginal)
      case(5)
        print *,"rk2.f03: myid=",myid,"change back to HLLD MHD solver"
        this%solverType=solverTypeOriginal
        solverTypeOriginal=0
      case default
   end select   

   q(1:nx,1:ny,1:nz,:)=0.5d0*(q(1:nx,1:ny,1:nz,:)+q2(1:nx,1:ny,1:nz,:))
   call this%exchangeBdryMPI(this%q,this%winq)
   call this%setBoundary(this%q)
   this%t=this%t+this%dt
end subroutine rk2_3D


subroutine rk2ADsg_2D(nthis,qn,qn1,qn2,ithis,qi,qi1,qi2)
use gridModule
use riemannSolverModule
use mpi
implicit none
type(grid)::nthis,ithis
double precision,dimension(1-nthis%nbuf:nthis%nMesh(1)+nthis%nbuf,1-nthis%nbuf:nthis%nMesh(2)+nthis%nbuf,nthis%nvar)::qn,qn1,qn2
double precision,dimension(:),allocatable::totalden
double precision,dimension(1-ithis%nbuf:ithis%nMesh(1)+ithis%nbuf,1-ithis%nbuf:ithis%nMesh(2)+ithis%nbuf,ithis%nvar)::qi,qi1,qi2
procedure(riemannSolver2D),pointer::rieSolvern,rieSolveri
integer::nvar,nx,ny,nbuf,eosTypen,eosTypei,SolverTypen,SolverTypei,solverTypeOriginaln,solverTypeOriginali
double precision::t,dt
integer::ierr,datasize,i,j,pos
logical::global_changeSolvern,global_changeSolveri

nbuf=nthis%nbuf
nx=nthis%nMesh(1)
ny=nthis%nMesh(2)
nvar=nthis%nvar
allocate(totalden((nx+2*nbuf)*(ny+2*nbuf)*nvar))
eosTypen=nthis%eosType
eosTypei=ithis%eosType
solverTypen=nthis%solverType
solverTypei=ithis%solverType
solverTypeOriginaln=0
solverTypeOriginali=0

if(eosTypen .eq. 1) then
  select case (solverTypen)
    case (1,2)
      rieSolvern=>solverIso2D
    case default
      print *, "rk2.f03: no suitable 2D solver found!"
      stop
  end select
elseif(eosTypen .eq. 2) then
  select case (solverTypen)
    case (2,3)
      rieSolvern=>solverAdi2D
    case (4,5)
      rieSolvern=>solverAdiMHD2D
    case default
      print *, "rk2.f03: no suitable 2D solver found!"
      stop
  end select
endif

if(eosTypei .eq. 1) then
  select case (solverTypei)
    case (1,2)
      rieSolveri=>solverIso2D
    case default
      print *, "rk2.f03: no suitable 2D solver found!"
      stop
  end select
elseif(eosTypei .eq. 2) then
  select case (solverTypei)
    case (2,3)
      rieSolveri=>solverAdi2D
    case (4,5)
      rieSolveri=>solverAdiMHD2D
    case default
      print *, "rk2.f03: no suitable 2D solver found!"
      stop
  end select
endif

t=nthis%t
datasize=(nthis%nMesh(1)+2*nthis%nbuf)*(nthis%nMesh(2)+2*nthis%nbuf)

do 
   call rieSolvern(nthis,qn,qn,qn1,dd=1)
   call rieSolvern(nthis,qn,qn1,qn1,dd=2)
   call rieSolveri(ithis,qi,qi,qi1,dd=1)
   call rieSolveri(ithis,qi,qi1,qi1,dd=2)
   totalden(1:datasize)=nthis%q(1:datasize)+ithis%q(1:datasize)
   call nthis%calcSelfgravity(totalden)
    
   dt=nthis%dt 
   do j=1,ny
     do i=1,nx
       pos=(j-1)*nx+i
       qn1(i,j,2)=qn1(i,j,2)+qn(i,j,1)*nthis%sgfx(pos)*dt
       qn1(i,j,3)=qn1(i,j,3)+qn(i,j,1)*nthis%sgfy(pos)*dt
       qn1(i,j,4)=qn1(i,j,4)+(qn(i,j,2)*nthis%sgfx(pos)+qn(i,j,3)*nthis%sgfy(pos))*dt
       qi1(i,j,2)=qi1(i,j,2)+qi(i,j,1)*nthis%sgfx(pos)*dt
       qi1(i,j,3)=qi1(i,j,3)+qi(i,j,1)*nthis%sgfy(pos)*dt
       qi1(i,j,8)=qi1(i,j,8)+(qi(i,j,2)*nthis%sgfx(pos)+qi(i,j,3)*nthis%sgfy(pos))*dt
     enddo
   enddo
   call evolveAD2D(nthis,ithis,nthis%q1,ithis%q1)
 
   call nthis%exchangeBdryMPI(nthis%q1,nthis%winq1)
   call nthis%setBoundary(nthis%q1)
   call ithis%exchangeBdryMPI(ithis%q1,ithis%winq1)
   call ithis%setBoundary(ithis%q1)

   call rieSolvern(nthis,qn1,qn1,qn2,dd=1)
   call rieSolvern(nthis,qn1,qn2,qn2,dd=2)
   call rieSolveri(ithis,qi1,qi1,qi2,dd=1)
   call rieSolveri(ithis,qi1,qi2,qi2,dd=2)
   totalden(1:datasize)=nthis%q1(1:datasize)+ithis%q1(1:datasize)
   call nthis%calcSelfgravity(totalden)
   do j=1,ny
     do i=1,nx
       pos=(j-1)*nx+i
       qn2(i,j,2)=qn2(i,j,2)+qn1(i,j,1)*nthis%sgfx(pos)*dt
       qn2(i,j,3)=qn2(i,j,3)+qn1(i,j,1)*nthis%sgfy(pos)*dt
       qn2(i,j,4)=qn2(i,j,4)+(qn1(i,j,2)*nthis%sgfx(pos)+qn1(i,j,3)*nthis%sgfy(pos))*dt
       qi2(i,j,2)=qi2(i,j,2)+qi1(i,j,1)*nthis%sgfx(pos)*dt
       qi2(i,j,3)=qi2(i,j,3)+qi1(i,j,1)*nthis%sgfy(pos)*dt
       qi2(i,j,8)=qi2(i,j,8)+(qi1(i,j,2)*nthis%sgfx(pos)+qi1(i,j,3)*nthis%sgfy(pos))*dt
     enddo
   enddo
   call evolveAD2D(nthis,ithis,nthis%q2,ithis%q2)

   call nthis%exchangeBdryMPI(nthis%q2,nthis%winq2)
   call nthis%setBoundary(nthis%q2)
   call ithis%exchangeBdryMPI(ithis%q2,ithis%winq2)
   call ithis%setBoundary(ithis%q2)

   call MPI_ALLREDUCE(nthis%changeSolver,global_changeSolvern,1,MPI_LOGICAL,MPI_LOR,MPI_COMM_WORLD,ierr)
   call MPI_ALLREDUCE(ithis%changeSolver,global_changeSolveri,1,MPI_LOGICAL,MPI_LOR,MPI_COMM_WORLD,ierr)

   if(global_changeSolveri .eqv. .false.) then
      exit
   else
      select case(ithis%solverType)
        case(5)
          print *,"rk2.f03: myid",myid,"change to use HLL MHD solver"
          solverTypeOriginali=5
          ithis%solverType=4
          ithis%changeSolver=.false.
        case(4)
          print *,"rk2.f03: myid=",myid,"HLL MHD solver encounters negative pressure!!"
          stop
      end select
   endif
enddo

   select case(solverTypeOriginali)
      case(5)
        print *,"rk2.f03: myid=",myid,"change back to HLLD MHD solver"
        ithis%solverType=solverTypeOriginali
        solverTypeOriginali=0
      case default
   end select

   qn(1:nx,1:ny,:)=0.5d0*(qn(1:nx,1:ny,:)+qn2(1:nx,1:ny,:))
   call nthis%exchangeBdryMPI(nthis%q,nthis%winq)
   call nthis%setBoundary(nthis%q)
   nthis%t=nthis%t+nthis%dt
   
   qi(1:nx,1:ny,:)=0.5d0*(qi(1:nx,1:ny,:)+qi2(1:nx,1:ny,:))
   call ithis%exchangeBdryMPI(ithis%q,ithis%winq)
   call ithis%setBoundary(ithis%q)
   ithis%t=ithis%t+ithis%dt

deallocate(totalden)
end subroutine rk2ADsg_2D

subroutine rk2AD_2D(nthis,qn,qn1,qn2,ithis,qi,qi1,qi2)
use gridModule
use riemannSolverModule
use mpi
implicit none
type(grid)::nthis,ithis
double precision,dimension(1-nthis%nbuf:nthis%nMesh(1)+nthis%nbuf,1-nthis%nbuf:nthis%nMesh(2)+nthis%nbuf,nthis%nvar)::qn,qn1,qn2
double precision,dimension(:),allocatable::totalden
double precision,dimension(1-ithis%nbuf:ithis%nMesh(1)+ithis%nbuf,1-ithis%nbuf:ithis%nMesh(2)+ithis%nbuf,ithis%nvar)::qi,qi1,qi2
procedure(riemannSolver2D),pointer::rieSolvern,rieSolveri
integer::nvar,nx,ny,nbuf,eosTypen,eosTypei,SolverTypen,SolverTypei,solverTypeOriginaln,solverTypeOriginali
double precision::t,dt
integer::ierr,datasize,i,j,pos
logical::global_changeSolvern,global_changeSolveri

nbuf=nthis%nbuf
nx=nthis%nMesh(1)
ny=nthis%nMesh(2)
nvar=nthis%nvar
eosTypen=nthis%eosType
eosTypei=ithis%eosType
solverTypen=nthis%solverType
solverTypei=ithis%solverType
solverTypeOriginaln=0
solverTypeOriginali=0

if(eosTypen .eq. 1) then
  select case (solverTypen)
    case (1,2)
      rieSolvern=>solverIso2D
    case default
      print *, "rk2.f03: no suitable 2D solver found!"
      stop
  end select
elseif(eosTypen .eq. 2) then
  select case (solverTypen)
    case (2,3)
      rieSolvern=>solverAdi2D
    case (4,5)
      rieSolvern=>solverAdiMHD2D
    case default
      print *, "rk2.f03: no suitable 2D solver found!"
      stop
  end select
endif

if(eosTypei .eq. 1) then
  select case (solverTypei)
    case (1,2)
      rieSolveri=>solverIso2D
    case default
      print *, "rk2.f03: no suitable 2D solver found!"
      stop
  end select
elseif(eosTypei .eq. 2) then
  select case (solverTypei)
    case (2,3)
      rieSolveri=>solverAdi2D
    case (4,5)
      rieSolveri=>solverAdiMHD2D
    case default
      print *, "rk2.f03: no suitable 2D solver found!"
      stop
  end select
endif

t=nthis%t

do 
   call rieSolvern(nthis,qn,qn,qn1,dd=1)
   call rieSolvern(nthis,qn,qn1,qn1,dd=2)
   call rieSolveri(ithis,qi,qi,qi1,dd=1)
   call rieSolveri(ithis,qi,qi1,qi1,dd=2)
    
   call evolveAD2D(nthis,ithis,nthis%q1,ithis%q1)
 
   call nthis%exchangeBdryMPI(nthis%q1,nthis%winq1)
   call nthis%setBoundary(nthis%q1)
   call ithis%exchangeBdryMPI(ithis%q1,ithis%winq1)
   call ithis%setBoundary(ithis%q1)

   call rieSolvern(nthis,qn1,qn1,qn2,dd=1)
   call rieSolvern(nthis,qn1,qn2,qn2,dd=2)
   call rieSolveri(ithis,qi1,qi1,qi2,dd=1)
   call rieSolveri(ithis,qi1,qi2,qi2,dd=2)
   call evolveAD2D(nthis,ithis,nthis%q2,ithis%q2)

   call nthis%exchangeBdryMPI(nthis%q2,nthis%winq2)
   call nthis%setBoundary(nthis%q2)
   call ithis%exchangeBdryMPI(ithis%q2,ithis%winq2)
   call ithis%setBoundary(ithis%q2)

   call MPI_ALLREDUCE(nthis%changeSolver,global_changeSolvern,1,MPI_LOGICAL,MPI_LOR,MPI_COMM_WORLD,ierr)
   call MPI_ALLREDUCE(ithis%changeSolver,global_changeSolveri,1,MPI_LOGICAL,MPI_LOR,MPI_COMM_WORLD,ierr)

   exit  !!! this will ignore the check for negative pressure

   if(global_changeSolveri .eqv. .false.) then
      exit
   else
      select case(ithis%solverType)
        case(5)
          print *,"rk2.f03: myid",myid,"change to use HLL MHD solver"
          solverTypeOriginali=5
          ithis%solverType=4
          ithis%changeSolver=.false.
        case(4)
          print *,"rk2.f03: myid=",myid,"HLL MHD solver encounters negative pressure!!"
          stop
      end select
   endif
enddo

   select case(solverTypeOriginali)
      case(5)
        print *,"rk2.f03: myid=",myid,"change back to HLLD MHD solver"
        ithis%solverType=solverTypeOriginali
        solverTypeOriginali=0
      case default
   end select

   qn(1:nx,1:ny,:)=0.5d0*(qn(1:nx,1:ny,:)+qn2(1:nx,1:ny,:))
   call nthis%exchangeBdryMPI(nthis%q,nthis%winq)
   call nthis%setBoundary(nthis%q)
   nthis%t=nthis%t+nthis%dt
   
   qi(1:nx,1:ny,:)=0.5d0*(qi(1:nx,1:ny,:)+qi2(1:nx,1:ny,:))
   call ithis%exchangeBdryMPI(ithis%q,ithis%winq)
   call ithis%setBoundary(ithis%q)
   ithis%t=ithis%t+ithis%dt

end subroutine rk2AD_2D






subroutine rk2AD_2D_HSHSMD(nthis,qn,qn1,qn2,ithis,qi,qi1,qi2)
use gridModule
use riemannSolverModule
use mpi
implicit none
type(grid)::nthis,ithis
double precision,dimension(1-nthis%nbuf:nthis%nMesh(1)+nthis%nbuf,1-nthis%nbuf:nthis%nMesh(2)+nthis%nbuf,nthis%nvar)::qn,qn1,qn2
double precision,dimension(:),allocatable::totalden
double precision,dimension(1-ithis%nbuf:ithis%nMesh(1)+ithis%nbuf,1-ithis%nbuf:ithis%nMesh(2)+ithis%nbuf,ithis%nvar)::qi,qi1,qi2
procedure(riemannSolver2D),pointer::rieSolvern,rieSolveri
integer::nvar,nx,ny,nbuf,eosTypen,eosTypei,SolverTypen,SolverTypei,solverTypeOriginaln,solverTypeOriginali
double precision::t,dt
integer::ierr,datasize,i,j,pos
logical::global_changeSolvern,global_changeSolveri

nbuf=nthis%nbuf
nx=nthis%nMesh(1)
ny=nthis%nMesh(2)
nvar=nthis%nvar
eosTypen=nthis%eosType
eosTypei=ithis%eosType
solverTypen=nthis%solverType
solverTypei=ithis%solverType
solverTypeOriginaln=0
solverTypeOriginali=0

if(eosTypen .eq. 1) then
  select case (solverTypen)
    case (1,2)
      rieSolvern=>solverIso2D
    case default
      print *, "rk2.f03: no suitable 2D solver found!"
      stop
  end select
elseif(eosTypen .eq. 2) then
  select case (solverTypen)
    case (2,3)
      rieSolvern=>solverAdi2D
    case (4,5)
      rieSolvern=>solverAdiMHD2D
    case default
      print *, "rk2.f03: no suitable 2D solver found!"
      stop
  end select
endif

if(eosTypei .eq. 1) then
  select case (solverTypei)
    case (1,2)
      rieSolveri=>solverIso2D
    case default
      print *, "rk2.f03: no suitable 2D solver found!"
      stop
  end select
elseif(eosTypei .eq. 2) then
  select case (solverTypei)
    case (2,3)
      rieSolveri=>solverAdi2D
    case (4,5)
      rieSolveri=>solverAdiMHD2D
    case default
      print *, "rk2.f03: no suitable 2D solver found!"
      stop
  end select
endif

t=nthis%t

do 
   call rieSolvern(nthis,qn,qn,qn1,dd=1)
   call rieSolvern(nthis,qn,qn1,qn1,dd=2)
   call rieSolveri(ithis,qi,qi,qi1,dd=1)
   call rieSolveri(ithis,qi,qi1,qi1,dd=2)
    
   !!!!!!!!!!!!!!!!zwg modified
   !!!call evolveAD2D(nthis,ithis,nthis%q1,ithis%q1)
      call evolveAD2D_MD(nthis,ithis,nthis%q,ithis%q,nthis%q1,ithis%q1)
 
   call nthis%exchangeBdryMPI(nthis%q1,nthis%winq1)
   call nthis%setBoundary(nthis%q1)
   call ithis%exchangeBdryMPI(ithis%q1,ithis%winq1)
   call ithis%setBoundary(ithis%q1)

   call rieSolvern(nthis,qn1,qn1,qn2,dd=1)
   call rieSolvern(nthis,qn1,qn2,qn2,dd=2)
   call rieSolveri(ithis,qi1,qi1,qi2,dd=1)
   call rieSolveri(ithis,qi1,qi2,qi2,dd=2)

   !!!!!!!!!!!!!!!!!zwg modified
   !!!call evolveAD2D(nthis,ithis,nthis%q2,ithis%q2)
      call evolveAD2D_MD(nthis,ithis,nthis%q1,ithis%q1,nthis%q2,ithis%q2)

   call nthis%exchangeBdryMPI(nthis%q2,nthis%winq2)
   call nthis%setBoundary(nthis%q2)
   call ithis%exchangeBdryMPI(ithis%q2,ithis%winq2)
   call ithis%setBoundary(ithis%q2)

   call MPI_ALLREDUCE(nthis%changeSolver,global_changeSolvern,1,MPI_LOGICAL,MPI_LOR,MPI_COMM_WORLD,ierr)
   call MPI_ALLREDUCE(ithis%changeSolver,global_changeSolveri,1,MPI_LOGICAL,MPI_LOR,MPI_COMM_WORLD,ierr)

   exit  !!! this will ignore the check for negative pressure

   if(global_changeSolveri .eqv. .false.) then
      exit
   else
      select case(ithis%solverType)
        case(5)
          print *,"rk2.f03: myid",myid,"change to use HLL MHD solver"
          solverTypeOriginali=5
          ithis%solverType=4
          ithis%changeSolver=.false.
        case(4)
          print *,"rk2.f03: myid=",myid,"HLL MHD solver encounters negative pressure!!"
          stop
      end select
   endif
enddo

   select case(solverTypeOriginali)
      case(5)
        print *,"rk2.f03: myid=",myid,"change back to HLLD MHD solver"
        ithis%solverType=solverTypeOriginali
        solverTypeOriginali=0
      case default
   end select

   qn(1:nx,1:ny,:)=0.5d0*(qn(1:nx,1:ny,:)+qn2(1:nx,1:ny,:))
   call nthis%exchangeBdryMPI(nthis%q,nthis%winq)
   call nthis%setBoundary(nthis%q)
   nthis%t=nthis%t+nthis%dt
   
   qi(1:nx,1:ny,:)=0.5d0*(qi(1:nx,1:ny,:)+qi2(1:nx,1:ny,:))
   call ithis%exchangeBdryMPI(ithis%q,ithis%winq)
   call ithis%setBoundary(ithis%q)
   ithis%t=ithis%t+ithis%dt

end subroutine rk2AD_2D_HSHSMD



subroutine rk2_2D(this,q,q1,q2)
use gridModule
use riemannSolverModule
use mpi
implicit none
class(grid)::this
double precision,dimension(1-this%nbuf:this%nMesh(1)+this%nbuf,1-this%nbuf:this%nMesh(2)+this%nbuf,this%nvar)::q,q1,q2
procedure(riemannSolver2D),pointer::rieSolver
integer::nvar,nx,ny,nbuf,ndim,eosType,solverType,solverTypeOriginal
double precision::t,dt
integer::ierr,i,j,pos
logical::global_changeSolver
!! eostype: 1=isothermal, 2=adiabatic
!! solverType: 1=exactHD, 2=HLLHD, 3=HLLCHD, 4=HLLMHD, 5=HLLDMHD)

eostype=this%eosType
solverType=this%solverType
solverTypeOriginal=0
if(eosType .eq. 1) then
  select case (solverType)
     case (1,2)
        rieSolver=>solverIso2D
     case (4,5)
        rieSolver=>solverIsoMHD2D
     case default
        print *, "rk2.f03: no suitable 2D solver found!"
        stop
  end select
elseif(eosType .eq. 2) then
  select case (solverType)
     case (2,3)
        rieSolver=>solverAdi2D
     case (4,5)
        rieSolver=>solverAdiMHD2D
     case default
        print *, "rk2.f03: no suitable 2D solver found!"
        stop
  end select
elseif(eosType .eq. 3) then
  select case (solverType)
     case(6)
       rieSolver=>solverPoly2D
     case default
       print *, "rk2.f03: no suitable 2D solver found!"
       stop
  end select
endif

nvar=this%nvar
nx=this%nMesh(1)
ny=this%nMesh(2)
nbuf=this%nbuf
t=this%t
!!!!!!! evolve one time step using RK2
!!! RK1
do 

   !write(*,*)"RK2RK2RK2RK2:::AAAAAAAAAAAAAAAAAAAAAAAAAAA!!!!!!!"

   call rieSolver(this,q,q,q1,dd=1)
   call rieSolver(this,q,q1,q1,dd=2)
   call source2D(this,q,q1)  !!! adding the externl source term

   !!!! added on Sep. 18, 2017 by hhwang
   if(this%enable_sg) then
     call this%calcSelfgravity(this%q)
     dt=this%dt
     if(eostype .eq. 1 .or. eostype .eq. 3) then !! if isothermal HD gas
       select case(solverType)
         case(1,2,4,5,6)
           do j=1,ny
             do i=1,nx
               pos=(j-1)*nx+i
               q1(i,j,2)=q1(i,j,2)+q(i,j,1)*this%sgfx(pos)*dt !! momx
               q1(i,j,3)=q1(i,j,3)+q(i,j,1)*this%sgfy(pos)*dt !! momy
             enddo
           enddo !! do j
         case default
           print *,'rk2_2D: cannot recognize solver'
           stop
       end select
     elseif(eostype .eq. 2) then !!! adiabatic
       if(solverType .eq. 2 .or. solverType .eq. 3) then !! HD solvers
         do j=1,ny
           do i=1,nx
             pos=(j-1)*nx+i
             q1(i,j,2)=q1(i,j,2)+q(i,j,1)*this%sgfx(pos)*dt !! momx
             q1(i,j,3)=q1(i,j,3)+q(i,j,1)*this%sgfy(pos)*dt !! momy
             q1(i,j,4)=q1(i,j,4)+(q(i,j,2)*this%sgfx(pos)+q(i,j,3)*this%sgfy(pos))*dt !!energy
           enddo
         enddo !! do j
       elseif(solverType .eq. 4 .or. solverType .eq. 5) then !!! MHD solvers
         do j=1,ny
           do i=1,nx
             pos=(j-1)*nx+i
             q1(i,j,2)=q1(i,j,2)+q(i,j,1)*this%sgfx(pos)*dt !! momx
             q1(i,j,3)=q1(i,j,3)+q(i,j,1)*this%sgfy(pos)*dt !! momy
             q1(i,j,8)=q1(i,j,8)+(q(i,j,2)*this%sgfx(pos)+q(i,j,3)*this%sgfy(pos))*dt
           enddo
         enddo
       endif !!! if(solverType .eq. 2 or 3 or 4 or 5)
     endif !!! if (eostype .eq. 1 or 2)
   endif !!! end enable_sg: added by hhwang on Sep. 18, 2017



   call this%exchangeBdryMPI(this%q1,this%winq1)
   call this%setBoundary(this%q1)
!!! RK2
   call rieSolver(this,q1,q1,q2,dd=1)
   call rieSolver(this,q1,q2,q2,dd=2)
   call source2D(this,q1,q2)  !!! adding the external source term

   !!!! added on Sep. 18, 2017 by hhwang
   if(this%enable_sg) then
     call this%calcSelfgravity(this%q1)
     dt=this%dt
     if(eostype .eq. 1 .or. eostype .eq. 3) then !! if isothermal HD gas
       select case(solverType)
         case(1,2,4,5,6)
           do j=1,ny
             do i=1,nx
               pos=(j-1)*nx+i
               q2(i,j,2)=q2(i,j,2)+q1(i,j,1)*this%sgfx(pos)*dt !! momx
               q2(i,j,3)=q2(i,j,3)+q1(i,j,1)*this%sgfy(pos)*dt !! momy
             enddo
           enddo !! do j
         case default
           print *,'rk2_2D: cannot recognize solver'
           stop
       end select
     elseif(eostype .eq. 2) then !!! adiabatic
       if(solverType .eq. 2 .or. solverType .eq. 3) then !! HD solvers
         do j=1,ny
           do i=1,nx
             pos=(j-1)*nx+i
             q2(i,j,2)=q2(i,j,2)+q1(i,j,1)*this%sgfx(pos)*dt !! momx
             q2(i,j,3)=q2(i,j,3)+q1(i,j,1)*this%sgfy(pos)*dt !! momy
             q2(i,j,4)=q2(i,j,4)+(q1(i,j,2)*this%sgfx(pos)+q1(i,j,3)*this%sgfy(pos))*dt !!energy
           enddo
         enddo !! do j
       elseif(solverType .eq. 4 .or. solverType .eq. 5) then !!! MHD solvers
         do j=1,ny
           do i=1,nx
             pos=(j-1)*nx+i
             q2(i,j,2)=q2(i,j,2)+q1(i,j,1)*this%sgfx(pos)*dt !! momx
             q2(i,j,3)=q2(i,j,3)+q1(i,j,1)*this%sgfy(pos)*dt !! momy
             q2(i,j,8)=q2(i,j,8)+(q1(i,j,2)*this%sgfx(pos)+q1(i,j,3)*this%sgfy(pos))*dt
           enddo
         enddo
       endif !!! if(solverType .eq. 2 or 3 or 4 or 5)
     endif !!! if (eostype .eq. 1 or 2)
   endif !!! end enable_sg: added by hhwang on Sep. 18, 2017

   call this%exchangeBdryMPI(this%q2,this%winq2)
   call this%setBoundary(this%q2)

   call MPI_ALLREDUCE(this%changeSolver,global_changeSolver,1,MPI_LOGICAL,MPI_LOR,MPI_COMM_WORLD,ierr)
   !!! check the global status of changeSolver
   
   if(global_changeSolver .eqv. .false.) then
      exit
   else
      select case(this%solverType)
        case(5)
          print *,"rk2.f03: myid",myid,"change to use HLL MHD solver"
          solverTypeOriginal=5
          this%solverType=4
          this%changeSolver=.false.
        case(4)
          print *,"rk2.f03: myid=",myid,"HLL MHD solver encounters negative pressure!!"
          stop
      end select
   endif   
enddo

   select case(solverTypeOriginal)
      case(5)
        print *,"rk2.f03: myid=",myid,"change back to HLLD MHD solver"
        this%solverType=solverTypeOriginal
        solverTypeOriginal=0
      case default
   end select   

   q(1:nx,1:ny,:)=0.5d0*(q(1:nx,1:ny,:)+q2(1:nx,1:ny,:))
   call this%exchangeBdryMPI(this%q,this%winq)
   call this%setBoundary(this%q)
   this%t=this%t+this%dt
end subroutine rk2_2D

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! RK2 for 1D AD problems
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine rk2AD_1D(nthis,qn,qn1,qn2,ithis,qi,qi1,qi2)
use gridModule
use riemannSolverModule
implicit none
type(grid)::nthis,ithis
double precision,dimension(1-nthis%nbuf:nthis%nMesh(1)+nthis%nbuf,nthis%nvar)::qn,qn1,qn2
double precision,dimension(1-ithis%nbuf:ithis%nMesh(1)+ithis%nbuf,ithis%nvar)::qi,qi1,qi2
procedure(riemannSolver1D),pointer::rieSolvern,rieSolveri
integer::nx,eosTypen,eosTypei, solverTypen,solverTypei
double precision::t

eosTypen=nthis%eosType
eosTypei=ithis%eosType
solverTypen=nthis%solverType
solverTypei=ithis%solverType
nx=nthis%nMesh(1)

if(eosTypen .eq. 1)then
  select case (solverTypen)
    case (1,2)
      rieSolvern=>solverIso1D
    case default
      print *,"rk2.f03: no suitable 1D solver found!"
      stop
  end select
elseif(eosTypen .eq. 2) then
  select case (solverTypen)
    case (2,3)
      rieSolvern=>solverAdi1D
    case (4,5)
      rieSolvern=>solverAdiMHD1D
    case default
      print *,"rk2.f03: no suitable 1D solver found!"
      stop
  end select
endif

if(eosTypei .eq. 1) then
  select case (solverTypei)
    case (1,2)
      rieSolveri=>solverIso1D
    case default
      print *,"rk2.f03: no suitable 1D solver found!"
      stop
  end select
elseif(eosTypei .eq. 2) then
  select case (solverTypei)
    case (2,3)
      rieSolveri=>solverAdi1D
    case (4,5)
      rieSolveri=>solverAdiMHD1D
    case default
      print *,"rk2.f03: no suitable 1D solver found!"
      stop
  end select
endif

t=nthis%t

call rieSolvern(nthis,qn,qn,qn1,dd=1)
call rieSolveri(ithis,qi,qi,qi1,dd=1)
call evolveAD1D(nthis,ithis,nthis%q1,ithis%q1)
call nthis%exchangeBdryMPI(nthis%q1,nthis%winq1)
call ithis%exchangeBdryMPI(ithis%q1,ithis%winq1)
call nthis%setBoundary(nthis%q1)
call ithis%setBoundary(ithis%q1)

call rieSolvern(nthis,qn1,qn1,qn2,dd=1)
call rieSolveri(ithis,qi1,qi1,qi2,dd=1)
call evolveAD1D(nthis,ithis,nthis%q2,ithis%q2)
call nthis%exchangeBdryMPI(nthis%q2,nthis%winq2)
call ithis%exchangeBdryMPI(ithis%q2,ithis%winq2)
call nthis%setBoundary(nthis%q2)
call ithis%setBoundary(ithis%q2)

qn(1:nx,:)=0.5d0*(qn(1:nx,:)+qn2(1:nx,:))
qi(1:nx,:)=0.5d0*(qi(1:nx,:)+qi2(1:nx,:))
call nthis%exchangeBdryMPI(nthis%q,nthis%winq)
call ithis%exchangeBdryMPI(ithis%q,ithis%winq)
call nthis%setBoundary(nthis%q)
call ithis%setBoundary(ithis%q)
nthis%t=nthis%t+nthis%dt
ithis%t=ithis%t+ithis%dt 

end subroutine rk2AD_1D

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!  RK2 for 1D problems
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine rk2_1D(this,q,q1,q2)
use gridModule
use riemannSolverModule
implicit none
class(grid)::this
double precision,dimension(1-this%nbuf:this%nMesh(1)+this%nbuf,this%nvar)::q,q1,q2
procedure(riemannSolver1D),pointer::rieSolver
integer::nvar,nx,nbuf,ndim,eosType,solverType
double precision::t

eosType=this%eosType
solverType=this%solverType
if(eosType .eq. 1) then
   select case (solverType)
      case (1,2)
         rieSolver=>solverIso1D
      case (4,5)
         rieSolver=>solverIsoMHD1D
      case default
         print *, "rk2.f03: no suitable 1D solver found!"
         stop
   end select
elseif(eosType .eq. 2) then
   select case (solverType)
      case(2,3)
         rieSolver=>solverAdi1D
      case (4,5)
         rieSolver=>solverAdiMHD1D
      case default
         print *, "rk2.f03: no suitable 1D solver found!"
         stop
   end select
elseif(eosType .eq. 3) then
   select case (solverType)
      case (6)
         rieSolver=>solverPoly1D
      case default
         print *, "rk2.f03: no suitable 1D solver found!"
         stop
   end select
endif


nvar=this%nvar
nx=this%nMesh(1)
nbuf=this%nbuf
t=this%t

!!! evolve one time step using RK2
!!! RK1
   call rieSolver(this,q,q,q1,dd=1)
   call this%exchangeBdryMPI(this%q1,this%winq1)
   call this%setBoundary(this%q1)
!!! RK2
   call rieSolver(this,q1,q1,q2,dd=1)
   call this%exchangeBdryMPI(this%q2,this%winq2)
   call this%setBoundary(this%q2)

   q(1:nx,:)=0.5d0*(q(1:nx,:)+q2(1:nx,:))
   call this%exchangeBdryMPI(this%q,this%winq)
   call this%setBoundary(this%q)
   this%t=this%t+this%dt
   
end subroutine rk2_1D
