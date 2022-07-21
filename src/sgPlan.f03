subroutine sgPlan3D(this,sgDenCmplx,sgfxCmplx,sgfyCmplx,sgfzCmplx)
use gridModule
use,intrinsic::iso_c_binding
use mpi
implicit none
class(grid)::this
complex(C_DOUBLE_COMPLEX),dimension(2*this%nMesh(1),2*this%nMesh(2),2*this%nMesh(3))::sgDenCmplx,sgfxCmplx,sgfyCmplx,sgfzCmplx
integer(C_INTPTR_T)::nx_global,ny_global,nz_global

if(this%sgBdryType .eq. 0) then
  nx_global=2*this%nMesh_global(1)
  ny_global=2*this%nMesh_global(2)
  nz_global=2*this%nMesh_global(3)
elseif(this%sgBdryType .eq. 1) then
  nx_global=this%nMesh_global(1)
  ny_global=this%nMesh_global(2)
  nz_global=this%nMesh_global(3)
endif

if(myid .eq. 0) then
  print *,"sgPlan3D.f90: creating sgPlan for density"
endif

this%sgPlanDen=fftw_mpi_plan_dft_3d(nz_global,ny_global,nx_global,sgDenCmplx,sgDenCmplx,MPI_COMM_WORLD,FFTW_FORWARD,FFTW_MEASURE)

if(myid .eq. 0) then
  print *,"sgPlan3D.f90: creating sgPlan for sgx"
endif
 this%sgPlanFx=fftw_mpi_plan_dft_3d(nz_global,ny_global,nx_global,sgfxCmplx ,sgfxCmplx ,MPI_COMM_WORLD,FFTW_BACKWARD,FFTW_MEASURE)

if(myid .eq. 0) then
  print *,"sgPlan3D.f90: creating sgPlan for sgy"
endif

 this%sgPlanFy=fftw_mpi_plan_dft_3d(nz_global,ny_global,nx_global,sgfyCmplx ,sgfyCmplx ,MPI_COMM_WORLD,FFTW_BACKWARD,FFTW_MEASURE)

if(myid .eq. 0) then
  print *,"sgPlan3D.f90: creating sgPlan for sgz"
endif

 this%sgPlanFz=fftw_mpi_plan_dft_3d(nz_global,ny_global,nx_global,sgfzCmplx ,sgfzCmplx ,MPI_COMM_WORLD,FFTW_BACKWARD,FFTW_MEASURE)

end subroutine sgPlan3D 


subroutine sgPlan2D(this,sgDenCmplx,sgfxCmplx,sgfyCmplx)
use gridModule
use,intrinsic::iso_c_binding
use mpi
implicit none
class(grid)::this
complex(C_DOUBLE_COMPLEX),dimension(2*this%nMesh(1),2*this%nMesh(2))::sgDenCmplx,sgfxCmplx,sgfyCmplx
integer(C_INTPTR_T)::nx_global,ny_global


if(this%sgBdryType .eq. 0) then  !!!! zwg added for poisson solver with 2D periodic BCs
nx_global=2*this%nMesh_global(1)
ny_global=2*this%nMesh_global(2)
elseif(this%sgBdryType .eq. 1) then  !!!! zwg added for poisson solver with 2D periodic BCs
  !!! the following code is implimented by zwg for poisson solver with 2D periodic BCs, which should be checked before using
  nx_global=this%nMesh_global(1)
  ny_global=this%nMesh_global(2)
  !!! the above code is implimented by zwg for poisson solver with 2D periodic BCs, which should be checked before using
endif !!!! zwg added for poisson solver with 2D periodic BCs

this%sgPlanDen=fftw_mpi_plan_dft_2d(ny_global,nx_global,sgDenCmplx,sgDenCmplx,MPI_COMM_WORLD,FFTW_FORWARD,FFTW_MEASURE)
this%sgPlanFx=fftw_mpi_plan_dft_2d(ny_global,nx_global,sgfxCmplx,sgfxCmplx,MPI_COMM_WORLD,FFTW_BACKWARD,FFTW_MEASURE)
this%sgPlanFy=fftw_mpi_plan_dft_2d(ny_global,nx_global,sgfyCmplx,sgfyCmplx,MPI_COMM_WORLD,FFTW_BACKWARD,FFTW_MEASURE)

end subroutine sgPlan2D
