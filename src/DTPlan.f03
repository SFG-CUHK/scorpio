subroutine DTPlan2D(this,Den_R2C_in,Den_R2C_out,Den_C2R_out, &
        & DTVx_C2R_in,DTVx_C2R_out,DTVy_C2R_in,DTVy_C2R_out)
use gridModule
use,intrinsic::iso_c_binding
use mpi
implicit none
class(grid)::this

  complex(C_DOUBLE_COMPLEX),dimension(this%nMesh(1)/2+1,this%nMesh(2))::Den_R2C_out
  real(C_DOUBLE), dimension(2*(this%nMesh(1)/2+1),this%nMesh(2)) :: Den_R2C_in,Den_C2R_out

  complex(C_DOUBLE_COMPLEX),dimension(this%nMesh(1)/2+1,this%nMesh(2))::DTVx_C2R_in,DTVy_C2R_in   !,DTVz_C2R_in
  real(C_DOUBLE), dimension(2*(this%nMesh(1)/2+1),this%nMesh(2)) :: DTVx_C2R_out,DTVy_C2R_out    !,DTVz_C2R_out

        !call DTPlan2D(this,this%Den_R2C_in,this%Den_R2C_out,this%Den_C2R_out, &
        !& this%DTVx_C2R_in,this%DTVx_C2R_out,this%DTVy_C2R_in,this%DTVy_C2R_out,this%DTVz_C2R_in,this%DTVz_C2R_out)
        !complex(C_DOUBLE_COMPLEX),dimension(2*this%nMesh(1),2*this%nMesh(2))::sgDenCmplx,sgfxCmplx,sgfyCmplx


integer(C_INTPTR_T)::nx_global,ny_global

  nx_global=this%nMesh_global(1)
  ny_global=this%nMesh_global(2)


!this%sgPlanDen=fftw_mpi_plan_dft_2d(ny_global,nx_global,sgDenCmplx,sgDenCmplx,MPI_COMM_WORLD,FFTW_FORWARD,FFTW_MEASURE)
!this%sgPlanFx=fftw_mpi_plan_dft_2d(ny_global,nx_global,sgfxCmplx,sgfxCmplx,MPI_COMM_WORLD,FFTW_BACKWARD,FFTW_MEASURE)
!this%sgPlanFy=fftw_mpi_plan_dft_2d(ny_global,nx_global,sgfyCmplx,sgfyCmplx,MPI_COMM_WORLD,FFTW_BACKWARD,FFTW_MEASURE)

      this%Den_Plan_R2C=fftw_mpi_plan_dft_r2c_2d(ny_global,nx_global, &
                        Den_R2C_in,Den_R2C_out,MPI_COMM_WORLD,FFTW_ESTIMATE)

      this%Den_Plan_C2R=fftw_mpi_plan_dft_c2r_2d(ny_global,nx_global, &
                        Den_R2C_out,Den_C2R_out,MPI_COMM_WORLD,FFTW_ESTIMATE)


      this%DTVx_Plan_C2R=fftw_mpi_plan_dft_c2r_2d(ny_global,nx_global, &
                        DTVx_C2R_in,DTVx_C2R_out,MPI_COMM_WORLD,FFTW_ESTIMATE)

      this%DTVy_Plan_C2R=fftw_mpi_plan_dft_c2r_2d(ny_global,nx_global, &
                        DTVy_C2R_in,DTVy_C2R_out,MPI_COMM_WORLD,FFTW_ESTIMATE)



!!!  type(C_PTR)::DTVx_Plan_C2R,DTVy_Plan_C2R,DTVz_Plan_C2R

end subroutine DTPlan2D





subroutine DTPlan3D(this,Den_R2C_in,Den_R2C_out,Den_C2R_out, &
        & DTVx_C2R_in,DTVx_C2R_out,DTVy_C2R_in,DTVy_C2R_out,DTVz_C2R_in,DTVz_C2R_out)


use gridModule
use,intrinsic::iso_c_binding
use mpi
implicit none
class(grid)::this

  complex(C_DOUBLE_COMPLEX),dimension(this%nMesh(1)/2+1,this%nMesh(2),this%nMesh(3))::Den_R2C_out
  real(C_DOUBLE), dimension(2*(this%nMesh(1)/2+1),this%nMesh(2),this%nMesh(3)) :: Den_R2C_in,Den_C2R_out

  complex(C_DOUBLE_COMPLEX),dimension(this%nMesh(1)/2+1,this%nMesh(2),this%nMesh(3))::DTVx_C2R_in,DTVy_C2R_in,DTVz_C2R_in
  real(C_DOUBLE), dimension(2*(this%nMesh(1)/2+1),this%nMesh(2),this%nMesh(3)) :: DTVx_C2R_out,DTVy_C2R_out,DTVz_C2R_out



integer(C_INTPTR_T)::nx_global,ny_global,nz_global

  nx_global=this%nMesh_global(1)
  ny_global=this%nMesh_global(2)
  nz_global=this%nMesh_global(3)



      this%Den_Plan_R2C=fftw_mpi_plan_dft_r2c_3d(nz_global,ny_global,nx_global, &
                        Den_R2C_in,Den_R2C_out,MPI_COMM_WORLD,FFTW_ESTIMATE)

      this%Den_Plan_C2R=fftw_mpi_plan_dft_c2r_3d(nz_global,ny_global,nx_global, &
                        Den_R2C_out,Den_C2R_out,MPI_COMM_WORLD,FFTW_ESTIMATE)


      this%DTVx_Plan_C2R=fftw_mpi_plan_dft_c2r_3d(nz_global,ny_global,nx_global, &
                        DTVx_C2R_in,DTVx_C2R_out,MPI_COMM_WORLD,FFTW_ESTIMATE)

      this%DTVy_Plan_C2R=fftw_mpi_plan_dft_c2r_3d(nz_global,ny_global,nx_global, &
                        DTVy_C2R_in,DTVy_C2R_out,MPI_COMM_WORLD,FFTW_ESTIMATE)

      this%DTVz_Plan_C2R=fftw_mpi_plan_dft_c2r_3d(nz_global,ny_global,nx_global, &
                        DTVz_C2R_in,DTVz_C2R_out,MPI_COMM_WORLD,FFTW_ESTIMATE)



end subroutine DTPlan3D 



