subroutine calcSG3Dperiodic(this,q,sgDenCmplx,sgfxCmplx,sgfyCmplx,sgfzCmplx,sgfx,sgfy,sgfz)
use gridModule
implicit none
class(grid)::this
double precision,dimension(1-this%nbuf:this%nMesh(1)+this%nbuf,1-this%nbuf:this%nMesh(2)+this%nbuf,&
                           1-this%nbuf:this%nMesh(3)+this%nbuf,this%nvar)::q
complex(C_DOUBLE_COMPLEX),dimension(this%nMesh(1),this%nMesh(2),this%nMesh(3)):: sgDenCmplx,&
                                                                     sgfxCmplx,sgfyCmplx,sgfzCmplx
double precision,dimension(this%nMesh(1),this%nMesh(2),this%nMesh(3))::sgfx,sgfy,sgfz,den
integer::nx_global,ny_global,nz_global
double precision::Lx_global,Ly_global,Lz_global
double precision::ki(this%nMesh_global(1)),kj(this%nMesh_global(2)),kk(this%nMesh_global(3))
double precision::kx(this%nMesh_global(1)),ky(this%nMesh_global(2)),kz(this%nMesh_global(3))
double precision::kxcon(this%nMesh_global(1)),kycon(this%nMesh_global(2)),kzcon(this%nMesh_global(3))
integer::i,j,k
double precision::pi,ksq,dx,dy,dz
complex(C_DOUBLE_COMPLEX)::phifft

pi=4.d0*datan2(1.d0,1.d0)

nx_global=this%nMesh_global(1)
ny_global=this%nMesh_global(2)
nz_global=this%nMesh_global(3)

Lx_global=this%rightBdry_global(1)-this%leftBdry_global(1)
Ly_global=this%rightBdry_global(2)-this%leftBdry_global(2)
Lz_global=this%rightBdry_global(3)-this%leftBdry_global(3)


if(mod(nz_global,nprocs) .ne. 0) then
   print *,"calcSG3Dperiodc: mesh number in z-direction should be divisible to the number of processors..."
   stop
endif


if(mod(nz_global,2) .ne. 0) then
   print *,"calcSG3Dperiodic: mesh number in z-direction needs to be even"
   stop
endif

do i=1,nx_global
  if(i .le. nx_global/2+1) then
    ki(i)=dble(i)
  else
    ki(i)=dble(nx_global-i+2)
  endif
    kx(i)=2.d0*pi/dble(nx_global)*(ki(i)-1.d0)*dble(nx_global)/Lx_global

  if(i .le. nx_global/2) then
    kxcon(i)=kx(i)
  else
    kxcon(i)=-kx(i)
  endif

!!!!!!! issue identified by wgzeng !!!!!!!
!!!!!!! added by hhwang 19.Mar.2019 !!!!!!!
!  if(mod(nx_global,2) .eq. 0) then
!    kxcon(nx_global/2+1) = 0.d0
!  endif
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

enddo

do i=1,ny_global
  if(i .le. ny_global/2+1) then
    kj(i)=dble(i)
  else
    kj(i)=dble(ny_global-i+2)
  endif
    ky(i)=2.d0*pi/dble(ny_global)*(kj(i)-1.d0)*dble(ny_global)/Ly_global

  if(i .le. ny_global/2) then
    kycon(i)=ky(i)
  else
    kycon(i)=-ky(i)
  endif

!!!!!!! issue identified by wgzeng !!!!!!!
!!!!!!! added by hhwang 19.Mar.2019 !!!!!!!
!  if(mod(ny_global,2) .eq. 0) then
!    kycon(ny_global/2+1) = 0.d0
!  endif
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
enddo

do i=1,nz_global
  if(i .le. nz_global/2+1) then
    kk(i)=dble(i)
  else
    kk(i)=dble(nz_global-i+2)
  endif
    kz(i)=2.d0*pi/dble(nz_global)*(kk(i)-1.d0)*dble(nz_global)/Lz_global

  if(i .le. nz_global/2) then
    kzcon(i)=kz(i)
  else
    kzcon(i)=-kz(i)
  endif

!!!!!!! issue identified by wgzeng !!!!!!!
!!!!!!! added by hhwang 19.Mar.2019 !!!!!!!
!  if(mod(nz_global,2) .eq. 0) then
!    kzcon(nz_global/2+1) = 0.d0
!  endif
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
enddo

den=q(1:this%nMesh(1),1:this%nMesh(2),1:this%nMesh(3),1)
sgDenCmplx(1:this%nMesh(1),1:this%nMesh(2),1:this%nMesh(3))=dcmplx(den,0.d0)
call fftw_mpi_execute_dft(this%sgPlanDen,sgDenCmplx,sgDenCmplx)

do k=1,this%nMesh(3)
  do j=1,this%nMesh(2)
    do i=1,this%nMesh(1)
      dx=this%dx(1)%coords(i)
      dy=this%dx(2)%coords(j)
      dz=this%dx(3)%coords(k)     
      ksq=kx(i)**2+ky(j)**2+kz(myid*this%nMesh(3)+k)**2
      phifft=-sgDenCmplx(i,j,k)/ksq
      sgfxCmplx(i,j,k)=phifft*dcmplx(0.d0,-1.d0*kxcon(i))
      sgfyCmplx(i,j,k)=phifft*dcmplx(0.d0,-1.d0*kycon(j))
      sgfzCmplx(i,j,k)=phifft*dcmplx(0.d0,-1.d0*kzcon(k+myid*this%nMesh(3)))
    enddo
  enddo
enddo

sgfxCmplx(1,1,1)=dcmplx(0.d0,0.d0)
sgfyCmplx(1,1,1)=dcmplx(0.d0,0.d0)
sgfzCmplx(1,1,1)=dcmplx(0.d0,0.d0)

call fftw_mpi_execute_dft(this%sgPlanFx,sgfxCmplx,sgfxCmplx)
call fftw_mpi_execute_dft(this%sgPlanFy,sgfyCmplx,sgfyCmplx)
call fftw_mpi_execute_dft(this%sgPlanFz,sgfzCmplx,sgfzCmplx)

sgfx=dreal(sgfxCmplx)/(nx_global*ny_global*nz_global)*4.d0*pi*this%GravConst
sgfy=dreal(sgfyCmplx)/(nx_global*ny_global*nz_global)*4.d0*pi*this%GravConst
sgfz=dreal(sgfzCmplx)/(nx_global*ny_global*nz_global)*4.d0*pi*this%GravConst

end subroutine calcSG3Dperiodic



!!!! The following subroutine is added by zwg for 2D poisson solver with  periodic BCs

subroutine calcSG2Dperiodic(this,q,sgDenCmplx,sgfxCmplx,sgfyCmplx,sgfx,sgfy)
	use gridModule
	implicit none
	class(grid)::this 

	double precision,dimension(1-this%nbuf:this%nMesh(1)+this%nbuf,1-this%nbuf:this%nMesh(2)+this%nbuf,this%nvar)::q
	complex(C_DOUBLE_COMPLEX),dimension(this%nMesh(1),this%nMesh(2)):: sgDenCmplx,sgfxCmplx,sgfyCmplx

	double precision,dimension(this%nMesh(1),this%nMesh(2))::sgfx,sgfy,den
	integer::nx_global,ny_global
	double precision::Lx_global,Ly_global

	double precision::ki(this%nMesh_global(1)),kj(this%nMesh_global(2))
	double precision::kx(this%nMesh_global(1)),ky(this%nMesh_global(2))
	double precision::kxcon(this%nMesh_global(1)),kycon(this%nMesh_global(2))

	integer::i,j
	double precision::pi,ksq,dx,dy

	complex(C_DOUBLE_COMPLEX)::phifft

	pi=4.d0*datan2(1.d0,1.d0)
	nx_global=this%nMesh_global(1)
	ny_global=this%nMesh_global(2)

	Lx_global=this%rightBdry_global(1)-this%leftBdry_global(1)
	Ly_global=this%rightBdry_global(2)-this%leftBdry_global(2)

	if(mod(ny_global,nprocs) .ne. 0) then
	   print *,"calcSG2Dperiodc: mesh number in y-direction should be divisible to the number of processors..."
	   stop
	endif


	if(mod(ny_global,2) .ne. 0) then
	   print *,"calcSG2Dperiodic: mesh number in y-direction needs to be even"
	   stop
	endif



	do i=1,nx_global
	  if(i .le. nx_global/2+1) then
		ki(i)=dble(i)
	  else
		ki(i)=dble(nx_global-i+2)
	  endif
		kx(i)=2.d0*pi/dble(nx_global)*(ki(i)-1.d0)*dble(nx_global)/Lx_global

	  if(i .le. nx_global/2) then
		kxcon(i)=kx(i)
	  else
		kxcon(i)=-kx(i)
	  endif
	!!!!!!! issue identified by wgzeng !!!!!!!
	!!!!!!! added by hhwang 19.Mar.2019 !!!!!!!
	  if(mod(nx_global,2) .eq. 0) then
		kxcon(nx_global/2+1) = 0.d0
	  endif
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	enddo

	do i=1,ny_global
	  if(i .le. ny_global/2+1) then
		kj(i)=dble(i)
	  else
		kj(i)=dble(ny_global-i+2)
	  endif
		ky(i)=2.d0*pi/dble(ny_global)*(kj(i)-1.d0)*dble(ny_global)/Ly_global

	  if(i .le. ny_global/2) then
		kycon(i)=ky(i)
	  else
		kycon(i)=-ky(i)
	  endif

	!!!!!!! issue identified by wgzeng  !!!!!!!
	!!!!!!! added by hhwang 19.Mar.2019 !!!!!!!
	  if(mod(ny_global,2) .eq. 0) then
		kycon(ny_global/2+1) = 0.d0
	  endif
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	enddo



	den=q(1:this%nMesh(1),1:this%nMesh(2),1)
	sgDenCmplx(1:this%nMesh(1),1:this%nMesh(2))=dcmplx(den,0.d0)
	call fftw_mpi_execute_dft(this%sgPlanDen,sgDenCmplx,sgDenCmplx)


	  do j=1,this%nMesh(2)
		do i=1,this%nMesh(1)
		  dx=this%dx(1)%coords(i)
		  dy=this%dx(2)%coords(j)
			
		  !ksq=kx(i)**2+ky(j)**2+kz(myid*this%nMesh(3)+k)**2
		  ksq=kx(i)**2+ky(myid*this%nMesh(2)+j)**2
		  phifft=-sgDenCmplx(i,j)/ksq
		  sgfxCmplx(i,j)=phifft*dcmplx(0.d0,-1.d0*kxcon(i))
		  sgfyCmplx(i,j)=phifft*dcmplx(0.d0,-1.d0*kycon(j+myid*this%nMesh(2)))
		  !sgfzCmplx(i,j,k)=phifft*dcmplx(0.d0,-1.d0*kzcon(k+myid*this%nMesh(3)))
		enddo
	  enddo


	sgfxCmplx(1,1)=dcmplx(0.d0,0.d0)
	sgfyCmplx(1,1)=dcmplx(0.d0,0.d0)


	call fftw_mpi_execute_dft(this%sgPlanFx,sgfxCmplx,sgfxCmplx)
	call fftw_mpi_execute_dft(this%sgPlanFy,sgfyCmplx,sgfyCmplx)


	sgfx=dreal(sgfxCmplx)/(nx_global*ny_global)*4.d0*pi*this%GravConst
	sgfy=dreal(sgfyCmplx)/(nx_global*ny_global)*4.d0*pi*this%GravConst


end subroutine calcSG2Dperiodic 

!!!! The above subroutine is added by zwg for 2D poisson solver with  periodic BCs




subroutine calcSG3D(this,q,sgDensityBuffer,sgDenCmplx,sgfxKernel,sgfyKernel,sgfzKernel,sgfxCmplx,&
                    sgfyCmplx,sgfzCmplx,sgfx,sgfy,sgfz)
use gridModule
implicit none
class(grid)::this
double precision,dimension(1-this%nbuf:this%nMesh(1)+this%nbuf,1-this%nbuf:this%nMesh(2)+this%nbuf,&
                           1-this%nbuf:this%nMesh(3)+this%nbuf,this%nvar)::q
double precision,dimension(this%nMesh(1),this%nMesh(2),2*this%nMesh(3))::sgDensityBuffer
complex(C_DOUBLE_COMPLEX),dimension(2*this%nMesh(1),2*this%nMesh(2),2*this%nMesh(3))::sgfxKernel,sgfyKernel,&
                                 sgfzKernel,sgDenCmplx,sgfxCmplx,sgfyCmplx,sgfzCmplx
double precision,dimension(this%nMesh(1),this%nMesh(2),this%nMesh(3))::sgfx,sgfy,sgfz,den
integer::send_id1,send_id2,datasize,ierr
integer(kind=MPI_ADDRESS_KIND)::offset
integer::i,j,k

datasize=this%nMesh(1)*this%nMesh(2)*this%nMesh(3)
send_id1=int(floor(float(myid)/2.))
sgDenCmplx=dcmplx(0.d0,0.d0)
sgDensityBuffer=0.d0
den=q(1:this%nMesh(1),1:this%nMesh(2),1:this%nMesh(3),1)

call MPI_WIN_FENCE(0,this%winsg,ierr)
  if(mod(myid,2) .eq. 0)then
    offset=0
    call MPI_PUT(den(1,1,1),datasize,MPI_DOUBLE,send_id1,offset,datasize,MPI_DOUBLE,this%winsg,ierr)
  endif
  if(mod(myid,2) .eq. 1)then
    offset=datasize
    call MPI_PUT(den(1,1,1),datasize,MPI_DOUBLE,send_id1,offset,datasize,MPI_DOUBLE,this%winsg,ierr)
  endif
call MPI_WIN_FENCE(0,this%winsg,ierr)

sgDenCmplx(1:this%nMesh(1),1:this%nMesh(2),1:2*this%nMesh(3))=dcmplx(sgDensityBuffer,0.d0)

call fftw_mpi_execute_dft(this%sgPlanDen,sgDenCmplx,sgDenCmplx)

sgfxCmplx=sgfxKernel*sgDenCmplx
sgfyCmplx=sgfyKernel*sgDenCmplx
sgfzCmplx=sgfzKernel*sgDenCmplx

call fftw_mpi_execute_dft(this%sgPlanFx,sgfxCmplx,sgfxCmplx)
call fftw_mpi_execute_dft(this%sgPlanFy,sgfyCmplx,sgfyCmplx)
call fftw_mpi_execute_dft(this%sgPlanFz,sgfzCmplx,sgfzCmplx)

send_id1=2*myid-nprocs
send_id2=2*myid-nprocs+1

!sgDensityBuffer=0.d0
!print *,"calcSG.f03: myid=",myid,"send_id1=",send_id1,"send_id2=",send_id2
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
call MPI_WIN_FENCE(0,this%winsg,ierr)
if(send_id1 .ge. 0) then
  den=dreal(sgfxCmplx(this%nMesh(1)+1:2*this%nMesh(1),this%nMesh(2)+1:2*this%nMesh(2),1:this%nMesh(3)))/ &
      dble(8*this%nMesh_global(1)*this%nMesh_global(2)*this%nMesh_global(3))*this%GravConst
  offset=0
  call MPI_PUT(den(1,1,1),datasize,MPI_DOUBLE,send_id1,offset,datasize,MPI_DOUBLE,this%winsg,ierr)
endif
call MPI_WIN_FENCE(0,this%winsg,ierr)

call MPI_WIN_FENCE(0,this%winsg,ierr)
if(send_id2 .ge. 0) then
  den=dreal(sgfxCmplx(this%nMesh(1)+1:2*this%nMesh(1),this%nMesh(2)+1:2*this%nMesh(2),this%nMesh(3)+1:2*this%nMesh(3)))/ &
       dble(8*this%nMesh_global(1)*this%nMesh_global(2)*this%nMesh_global(3))*this%GravConst
  offset=0
  call MPI_PUT(den(1,1,1),datasize,MPI_DOUBLE,send_id2,offset,datasize,MPI_DOUBLE,this%winsg,ierr)
endif
call MPI_WIN_FENCE(0,this%winsg,ierr)

if(nprocs .gt. 1) then
  sgfx=sgDensityBuffer(1:this%nMesh(1),1:this%nMesh(2),1:this%nMesh(3))
else
  sgfx=dreal(sgfxCmplx(this%nMesh(1)+1:2*this%nMesh(1),this%nMesh(2)+1:2*this%nMesh(2),this%nMesh(3)+1:2*this%nMesh(3)))/ &
        dble(8*this%nMesh_global(1)*this%nMesh_global(2)*this%nMesh_global(3))*this%GravConst
endif
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
call MPI_WIN_FENCE(0,this%winsg,ierr)
if(send_id1 .ge. 0) then
  den=dreal(sgfyCmplx(this%nMesh(1)+1:2*this%nMesh(1),this%nMesh(2)+1:2*this%nMesh(2),1:this%nMesh(3)))/ &
       dble(8*this%nMesh_global(1)*this%nMesh_global(2)*this%nMesh_global(3))*this%GravConst
  offset=0
  call MPI_PUT(den(1,1,1),datasize,MPI_DOUBLE,send_id1,offset,datasize,MPI_DOUBLE,this%winsg,ierr)
endif
call MPI_WIN_FENCE(0,this%winsg,ierr)

call MPI_WIN_FENCE(0,this%winsg,ierr)
if(send_id2 .ge. 0) then
  den=dreal(sgfyCmplx(this%nMesh(1)+1:2*this%nMesh(1),this%nMesh(2)+1:2*this%nMesh(2),this%nMesh(3)+1:2*this%nMesh(3)))/ &
       dble(8*this%nMesh_global(1)*this%nMesh_global(2)*this%nMesh_global(3))*this%GravConst
  offset=0
  call MPI_PUT(den(1,1,1),datasize,MPI_DOUBLE,send_id2,offset,datasize,MPI_DOUBLE,this%winsg,ierr)
endif
call MPI_WIN_FENCE(0,this%winsg,ierr)

if(nprocs .gt. 1) then
  sgfy=sgDensityBuffer(1:this%nMesh(1),1:this%nMesh(2),1:this%nMesh(3))
else
  sgfy=dreal(sgfxCmplx(this%nMesh(1)+1:2*this%nMesh(1),this%nMesh(2)+1:2*this%nMesh(2),this%nMesh(3)+1:2*this%nMesh(3)))/ &
        dble(8*this%nMesh_global(1)*this%nMesh_global(2)*this%nMesh_global(3))*this%GravConst
endif
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
call MPI_WIN_FENCE(0,this%winsg,ierr)
if(send_id1 .ge. 0) then
  den=dreal(sgfzCmplx(this%nMesh(1)+1:2*this%nMesh(1),this%nMesh(2)+1:2*this%nMesh(2),1:this%nMesh(3)))/ &
       dble(8*this%nMesh_global(1)*this%nMesh_global(2)*this%nMesh_global(3))*this%GravConst
  offset=0
  call MPI_PUT(den(1,1,1),datasize,MPI_DOUBLE,send_id1,offset,datasize,MPI_DOUBLE,this%winsg,ierr)
endif
call MPI_WIN_FENCE(0,this%winsg,ierr)

call MPI_WIN_FENCE(0,this%winsg,ierr)
if(send_id2 .ge. 0) then
  den=dreal(sgfzCmplx(this%nMesh(1)+1:2*this%nMesh(1),this%nMesh(2)+1:2*this%nMesh(2),this%nMesh(3)+1:2*this%nMesh(3)))/ &
       dble(8*this%nMesh_global(1)*this%nMesh_global(2)*this%nMesh_global(3))*this%GravConst
  offset=0
  call MPI_PUT(den(1,1,1),datasize,MPI_DOUBLE,send_id2,offset,datasize,MPI_DOUBLE,this%winsg,ierr)
endif
call MPI_WIN_FENCE(0,this%winsg,ierr)

if(nprocs .gt. 1) then
  sgfz=sgDensityBuffer(1:this%nMesh(1),1:this%nMesh(2),1:this%nMesh(3))
else
  sgfz=dreal(sgfxCmplx(this%nMesh(1)+1:2*this%nMesh(1),this%nMesh(2)+1:2*this%nMesh(2),this%nMesh(3)+1:2*this%nMesh(3)))/ &
        dble(8*this%nMesh_global(1)*this%nMesh_global(2)*this%nMesh_global(3))*this%GravConst
endif



!if(myid .eq. 0)then
!    print *,"calcSG2D.f03: myid=",myid
!  do j=1,this%nMesh(2)
    !print *,dreal(sgfxCmplx(this%nMesh(1)+1:2*this%nMesh(1),j))/dble(4*this%nMesh_global(1)*this%nMesh_global(2))
!    print *,sgfx(:,j)
!  enddo
!endif

end subroutine calcSG3D


subroutine calcSG2D(this,q,sgDensityBuffer,sgDenCmplx,sgfxKernel,sgfyKernel,sgfxCmplx,sgfyCmplx,sgfx,sgfy)
use gridModule
implicit none
class(grid)::this
double precision,dimension(1-this%nbuf:this%nMesh(1)+this%nbuf,1-this%nbuf:this%nMesh(2)+this%nbuf,this%nvar)::q
double precision,dimension(this%nMesh(1),2*this%nMesh(2))::sgDensityBuffer
complex(C_DOUBLE_COMPLEX),dimension(2*this%nMesh(1),2*this%nMesh(2))::sgfxKernel,sgfyKernel,sgDenCmplx,sgfxCmplx,sgfyCmplx
double precision,dimension(this%nMesh(1),this%nMesh(2))::sgfx,sgfy,den
integer::send_id1,send_id2,datasize,ierr
integer(kind=MPI_ADDRESS_KIND)::offset
integer::i,j

datasize=this%nMesh(1)*this%nMesh(2)
send_id1=int(floor(float(myid)/2.))
sgDenCmplx=dcmplx(0.d0,0.d0)
sgDensityBuffer=0.d0
den=q(1:this%nMesh(1),1:this%nMesh(2),1)

call MPI_WIN_FENCE(0,this%winsg,ierr)
  if(mod(myid,2) .eq. 0)then
    offset=0
    call MPI_PUT(den(1,1),datasize,MPI_DOUBLE,send_id1,offset,datasize,MPI_DOUBLE,this%winsg,ierr)
  endif
  if(mod(myid,2) .eq. 1)then
    offset=datasize
    call MPI_PUT(den(1,1),datasize,MPI_DOUBLE,send_id1,offset,datasize,MPI_DOUBLE,this%winsg,ierr)
  endif
call MPI_WIN_FENCE(0,this%winsg,ierr)

sgDenCmplx(1:this%nMesh(1),1:2*this%nMesh(2))=dcmplx(sgDensityBuffer,0.d0)

call fftw_mpi_execute_dft(this%sgPlanDen,sgDenCmplx,sgDenCmplx)

sgfxCmplx=sgfxKernel*sgDenCmplx
sgfyCmplx=sgfyKernel*sgDenCmplx

call fftw_mpi_execute_dft(this%sgPlanFx,sgfxCmplx,sgfxCmplx)
call fftw_mpi_execute_dft(this%sgPlanFy,sgfyCmplx,sgfyCmplx)

send_id1=2*myid-nprocs
send_id2=2*myid-nprocs+1

!sgDensityBuffer=0.d0
!print *,"calcSG.f03: myid=",myid,"send_id1=",send_id1,"send_id2=",send_id2
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
call MPI_WIN_FENCE(0,this%winsg,ierr)
if(send_id1 .ge. 0) then
  den=dreal(sgfxCmplx(this%nMesh(1)+1:2*this%nMesh(1),1:this%nMesh(2)))/ &
      dble(4*this%nMesh_global(1)*this%nMesh_global(2))*this%GravConst
  offset=0
  call MPI_PUT(den(1,1),datasize,MPI_DOUBLE,send_id1,offset,datasize,MPI_DOUBLE,this%winsg,ierr)
endif
call MPI_WIN_FENCE(0,this%winsg,ierr)

call MPI_WIN_FENCE(0,this%winsg,ierr)
if(send_id2 .ge. 0) then
  den=dreal(sgfxCmplx(this%nMesh(1)+1:2*this%nMesh(1),this%nMesh(2)+1:2*this%nMesh(2)))/ &
       dble(4*this%nMesh_global(1)*this%nMesh_global(2))*this%GravConst
  offset=0
  call MPI_PUT(den(1,1),datasize,MPI_DOUBLE,send_id2,offset,datasize,MPI_DOUBLE,this%winsg,ierr)
endif
call MPI_WIN_FENCE(0,this%winsg,ierr)

if(nprocs .gt. 1) then
  sgfx=sgDensityBuffer(1:this%nMesh(1),1:this%nMesh(2))
else
  sgfx=dreal(sgfxCmplx(this%nMesh(1)+1:2*this%nMesh(1),this%nMesh(2)+1:2*this%nMesh(2)))/ &
        dble(4*this%nMesh_global(1)*this%nMesh_global(2))*this%GravConst
endif
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
call MPI_WIN_FENCE(0,this%winsg,ierr)
if(send_id1 .ge. 0) then
  den=dreal(sgfyCmplx(this%nMesh(1)+1:2*this%nMesh(1),1:this%nMesh(2)))/ &
       dble(4*this%nMesh_global(1)*this%nMesh_global(2))*this%GravConst
  offset=0
  call MPI_PUT(den(1,1),datasize,MPI_DOUBLE,send_id1,offset,datasize,MPI_DOUBLE,this%winsg,ierr)
  !!! zwg: MPI_Put - Copies data from the origin memory to the target. 
endif
call MPI_WIN_FENCE(0,this%winsg,ierr)

call MPI_WIN_FENCE(0,this%winsg,ierr)
if(send_id2 .ge. 0) then
  den=dreal(sgfyCmplx(this%nMesh(1)+1:2*this%nMesh(1),this%nMesh(2)+1:2*this%nMesh(2)))/ &
       dble(4*this%nMesh_global(1)*this%nMesh_global(2))*this%GravConst
  offset=0
  call MPI_PUT(den(1,1),datasize,MPI_DOUBLE,send_id2,offset,datasize,MPI_DOUBLE,this%winsg,ierr)
endif
call MPI_WIN_FENCE(0,this%winsg,ierr)

if(nprocs .gt. 1) then
  sgfy=sgDensityBuffer(1:this%nMesh(1),1:this%nMesh(2))
else
  sgfy=dreal(sgfxCmplx(this%nMesh(1)+1:2*this%nMesh(1),this%nMesh(2)+1:2*this%nMesh(2)))/ &
        dble(4*this%nMesh_global(1)*this%nMesh_global(2))*this%GravConst
endif



!if(myid .eq. 0)then
!    print *,"calcSG2D.f03: myid=",myid
!  do j=1,this%nMesh(2)
!    print *,dreal(sgfxCmplx(this%nMesh(1)+1:2*this%nMesh(1),j))/dble(4*this%nMesh_global(1)*this%nMesh_global(2))
!    print *,sgfx(:,j)
!  enddo
!endif

end subroutine calcSG2D
