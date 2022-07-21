subroutine sgKernel3d(this,sgfxKernel,sgfyKernel,sgfzKernel)
use mpi
use gridModule
use,intrinsic::iso_c_binding
implicit none
class(grid)::this
complex(C_DOUBLE_COMPLEX),dimension(2*this%nMesh(1),2*this%nMesh(2),2*this%nMesh(3))::sgfxKernel,sgfyKernel,sgfzKernel
double precision::x_ker_loc(1:2*this%nMesh(1))
double precision::y_ker_loc(1:2*this%nMesh(2))
double precision::z_ker_loc(1:2*this%nMesh(3))
double precision::dx,dy,dz,hdx,hdy,hdz,leftBdry_new(3),darea,xl,yl,zl,r_ker_loc
integer(C_INTPTR_T)::i,j,k,nx_ker,ny_ker,nz_ker,nx_ker_global,ny_ker_global,nz_ker_global
type(C_PTR)::planx,plany,planz

!print *,'sgKernel3d.f90: preparing the 3D Kernel...'

if(this%sgBdryType .eq. 0) then
nx_ker=2*this%nMesh(1)
ny_ker=2*this%nMesh(2)
nz_ker=2*this%nMesh(3)
nx_ker_global=2*this%nMesh_global(1)
ny_ker_global=2*this%nMesh_global(2)
nz_ker_global=2*this%nMesh_global(3)
dx=this%dx(1)%coords(1)
dy=this%dx(2)%coords(2)
dz=this%dx(3)%coords(3)

leftBdry_new(1)=2.d0*(this%leftBdry_global(1)-0.5d0*(this%leftBdry_global(1)+this%rightBdry_global(1)))
leftBdry_new(2)=2.d0*(this%leftBdry_global(2)-0.5d0*(this%leftBdry_global(2)+this%rightBdry_global(2)))
leftBdry_new(3)=2.d0*(this%leftBdry_global(3)-0.5d0*(this%leftBdry_global(3)+this%rightBdry_global(3)))

do i=1,nx_ker
  x_ker_loc(i)=leftBdry_new(1)+(dble(i)-0.5d0)*dx
enddo

do i=1,ny_ker
  y_ker_loc(i)=leftBdry_new(2)+(dble(i)-0.5d0)*dy
enddo

do i=1,nz_ker
  z_ker_loc(i)=leftBdry_new(3)+(dble(i)-0.5d0+dble(myid*2*this%nMesh(3)))*dz
enddo

hdx=0.5d0*dx
hdy=0.5d0*dy
hdz=0.5d0*dz
darea=dx*dy*dz

planx=fftw_mpi_plan_dft_3d(nz_ker_global,ny_ker_global,nx_ker_global,sgfxKernel,sgfxKernel,&
                           MPI_COMM_WORLD,FFTW_FORWARD,FFTW_ESTIMATE)
plany=fftw_mpi_plan_dft_3d(nz_ker_global,ny_ker_global,nx_ker_global,sgfyKernel,sgfyKernel,&
                           MPI_COMM_WORLD,FFTW_FORWARD,FFTW_ESTIMATE)
planz=fftw_mpi_plan_dft_3d(nz_ker_global,ny_ker_global,nx_ker_global,sgfzKernel,sgfzKernel,&
                           MPI_COMM_WORLD,FFTW_FORWARD,FFTW_ESTIMATE)

if(myid .eq. 0) then
   print *,"sgKernel3d.f03: Filling kernel matrix for selfgravity ......"
endif

do k=1,nz_ker
  do j=1,ny_ker
    do i=1,nx_ker
      zl=z_ker_loc(k)-hdz
      yl=y_ker_loc(j)-hdy
      xl=x_ker_loc(i)-hdx

      r_ker_loc=dsqrt(zl**2+yl**2+xl**2)
      sgfxKernel(i,j,k)=dcmplx(-xl/r_ker_loc**3*darea,0.d0)
      sgfyKernel(i,j,k)=dcmplx(-yl/r_ker_loc**3*darea,0.d0)
      sgfzKernel(i,j,k)=dcmplx(-zl/r_ker_loc**3*darea,0.d0)

      if(r_ker_loc .lt. 0.5d0*dx) then
        sgfxKernel(i,j,k)=dcmplx(0.d0,0.d0)
        sgfyKernel(i,j,k)=dcmplx(0.d0,0.d0)
        sgfzKernel(i,j,k)=dcmplx(0.d0,0.d0)
      endif
    enddo
  enddo
enddo

if(myid .eq. 0) then
  print *,"sgKernel3d.f03: doing 3D-FFT for kernels......"
endif

call fftw_mpi_execute_dft(planx,sgfxKernel,sgfxKernel)
call fftw_mpi_execute_dft(plany,sgfyKernel,sgfyKernel)
call fftw_mpi_execute_dft(planz,sgfzKernel,sgfzKernel)

if(myid .eq. 0) then
  print *,"sgKernel3d.f03: kernel preparations done......"
endif

call fftw_destroy_plan(planx)
call fftw_destroy_plan(plany)
call fftw_destroy_plan(planz)
endif  !! if this%sgBdryType .eq. 0 (isolated)

if(this%sgBdryType .eq. 1) then
  if(myid .eq. 0) then
    print *,"no need to prepare kernels for periodic selfgravity....."
  endif
endif  !! if this%sgBdryType .eq. 1 (periodic)

end subroutine sgKernel3d


subroutine sgKernel2d(this,sgfxKernel,sgfyKernel)
use mpi
use gridModule
use,intrinsic::iso_c_binding
implicit none
class(grid)::this
complex(C_DOUBLE_COMPLEX),dimension(2*this%nMesh(1),2*this%nMesh(2))::sgfxKernel,sgfyKernel
double precision::x_ker_loc(1:2*this%nMesh(1))
double precision::y_ker_loc(1:2*this%nMesh(2))
double precision::dx,dy,hdx,hdy,leftBdry_new(2),darea,xl,yl,r_ker_loc
integer(C_INTPTR_T)::i,j,nx_ker,ny_ker,nx_ker_global,ny_ker_global
type(C_PTR)::planx,plany

nx_ker=2*this%nMesh(1)
ny_ker=2*this%nMesh(2)
nx_ker_global=2*this%nMesh_global(1)
ny_ker_global=2*this%nMesh_global(2)
dx=this%dx(1)%coords(1)
dy=this%dx(2)%coords(2)
leftBdry_new(1)=2.d0*(this%leftBdry_global(1)-0.5d0*(this%leftBdry_global(1)+this%rightBdry_global(1)))
leftBdry_new(2)=2.d0*(this%leftBdry_global(2)-0.5d0*(this%leftBdry_global(2)+this%rightBdry_global(2)))

do i=1,nx_ker
  x_ker_loc(i)=leftBdry_new(1)+(dble(i)-0.5d0)*dx
enddo

do i=1,ny_ker
  y_ker_loc(i)=leftBdry_new(2)+(dble(i)-0.5d0+dble(myid*2*this%nMesh(2)))*dy
enddo

hdx=0.5d0*dx
hdy=0.5d0*dy
darea=dx*dy

planx=fftw_mpi_plan_dft_2d(ny_ker_global,nx_ker_global,sgfxKernel,sgfxKernel,MPI_COMM_WORLD,FFTW_FORWARD,FFTW_ESTIMATE)
plany=fftw_mpi_plan_dft_2d(ny_ker_global,nx_ker_global,sgfyKernel,sgfyKernel,MPI_COMM_WORLD,FFTW_FORWARD,FFTW_ESTIMATE)

if(myid .eq. 0) then
   print *,"sgKernel2d.f03: Filling kernel matrix for selfgravity ......"
endif

do j=1,ny_ker
  do i=1,nx_ker
    yl=y_ker_loc(j)-hdx
    xl=x_ker_loc(i)-hdx
    r_ker_loc=dsqrt(yl**2+xl**2)
    sgfxKernel(i,j)=dcmplx(-xl/r_ker_loc**3*darea,0.d0)
    sgfyKernel(i,j)=dcmplx(-yl/r_ker_loc**3*darea,0.d0)

    if(r_ker_loc .lt. 0.5d0*dx) then
      sgfxKernel(i,j)=dcmplx(0.d0,0.d0)
      sgfyKernel(i,j)=dcmplx(0.d0,0.d0)
    endif
  enddo
enddo

if(myid .eq. 0) then
  print *,"sgKernel2d.f03: doing 2D-FFT for kernels......"
endif

call fftw_mpi_execute_dft(planx,sgfxKernel,sgfxKernel)
call fftw_mpi_execute_dft(plany,sgfyKernel,sgfyKernel)

if(myid .eq. 0) then
  print *,"sgKernel2d.f03: kernel preparations done......"
endif

call fftw_destroy_plan(planx)
call fftw_destroy_plan(plany)

end subroutine sgKernel2d
