!!!! The following subroutine is added by WGZeng for 2D and 3D Driving Turbulence Problems

subroutine calcDT2D(this,Den_ini,Den_R2C_in,Den_R2C_out,Den_C2R_out,DTVx_C2R_in,DTVx_C2R_out,DTVy_C2R_in,DTVy_C2R_out,DTVx,DTVy)
use gridModule
use,intrinsic::iso_c_binding
implicit none
class(grid)::this 


     double precision,dimension(this%nMesh(1),this%nMesh(2))::Den_ini,DTVx,DTVy

     complex(C_DOUBLE_COMPLEX),dimension(this%nMesh(1)/2+1,this%nMesh(2)):: Den_R2C_out 
     real(C_DOUBLE),dimension((this%nMesh(1)/2+1)*2,this%nMesh(2)):: Den_R2C_in,Den_C2R_out

     complex(C_DOUBLE_COMPLEX),dimension(this%nMesh(1)/2+1,this%nMesh(2)):: DTVx_C2R_in,DTVy_C2R_in
     real(C_DOUBLE),dimension((this%nMesh(1)/2+1)*2,this%nMesh(2)):: DTVx_C2R_out,DTVy_C2R_out


integer::nx_global,ny_global
double precision::Lx_global,Ly_global

double precision::ki(this%nMesh_global(1)),kj(this%nMesh_global(2))
double precision::kx(this%nMesh_global(1)),ky(this%nMesh_global(2))
double precision::kxcon(this%nMesh_global(1)),kycon(this%nMesh_global(2))

integer::i,j,Nx,Ny
double precision::pi,ksq,dx,dy


    double precision:: sigma
    double precision:: sumkdist=0.0
    double precision:: kkx,kky,km,km1,km2

    double precision::beta,Pxx, Pxy, Pyy
    double precision::A0,Ax,Ay
    double precision::Ax_re,Ax_im,Ay_re,Ay_im

    double precision,external::randn, rand1
    double precision::rand11


    double precision::drivingWN,Energy,zeta



drivingWN=this%drivingWN_DT
Energy=this%Energy_DT
!!!zeta=this%Energy_DT !!! zwg: modifid by lester
zeta=this%zeta_DT

pi=4.d0*datan2(1.d0,1.d0)
nx_global=this%nMesh_global(1)
ny_global=this%nMesh_global(2)

Nx=this%nMesh(1)
Ny=this%nMesh(2)

Lx_global=this%rightBdry_global(1)-this%leftBdry_global(1)
Ly_global=this%rightBdry_global(2)-this%leftBdry_global(2)

if(mod(ny_global,nprocs) .ne. 0) then
   print *,"calcDT2D: mesh number in y-direction should be divisible to the number of processors..."
   stop
endif


if(mod(ny_global,2) .ne. 0) then
   print *,"calcDT2D: mesh number in y-direction needs to be even"
   stop
endif




  !!! zwg : for checking R2C and C2R interface of FFTW
  Den_R2C_in(1:Nx,1:Ny)=Den_ini(1:Nx,1:Ny)
  call fftw_mpi_execute_dft_r2c(this%Den_Plan_R2C,Den_R2C_in,Den_R2C_out)


  if(myid==0)then

    write(*,*)"this Process id isAAAAAAAAAAAAAAAAAAAAAAA: ", myid
    do j=1,Ny
      !do i=1, Nx
      write(*,*) Den_R2C_out(:,j)
      !enddo
   enddo
  
  endif



  call fftw_mpi_execute_dft_c2r(this%Den_Plan_C2R,Den_R2C_out,Den_C2R_out)

  if(myid==0)then

    write(*,*)"this Process id isBBBBBBBBBBBBBBBBBBBBBBB: ", myid
    do j=1,Ny
      !do i=1, Nx
        write(*,*) Den_C2R_out(:,j)/(nx_global*ny_global)
      !enddo
      
   enddo
  
  endif






do i=1,nx_global
  if(i .le. nx_global/2+1) then
    ki(i)=dble(i)
  else
    ki(i)=dble(nx_global-i+2)
  endif
   !kx(i)=(ki(i)-1.d0)		!Lester modified
   kx(i)=2.d0*pi/dble(nx_global)*(ki(i)-1.d0)*dble(nx_global)/Lx_global

  if(i .le. nx_global/2+1) then  !!!zwg modified
    kxcon(i)=kx(i)
  else
    kxcon(i)=-kx(i)
  endif

enddo

do i=1,ny_global
  if(i .le. ny_global/2+1) then
    kj(i)=dble(i)
  else
    kj(i)=dble(ny_global-i+2)
  endif
   !ky(i)=(kj(i)-1.d0)	!Lester modified
    ky(i)=2.d0*pi/dble(ny_global)*(kj(i)-1.d0)*dble(ny_global)/Ly_global

  if(i .le. ny_global/2+1) then  !!!zwg modified
    kycon(i)=ky(i)
  else
    kycon(i)=-ky(i)
  endif
enddo



      do j=1, ny_global
        do i=1,nx_global

          km=sqrt(kxcon(i)**2+kycon(j)**2)
          
	  sumkdist=sumkdist+(km**7)*exp(-8.0*km/drivingWN)	!Lester modified
        !!sumkdist=sumkdist+(km**6)*exp(-8.0*km/drivingWN)

        enddo
     enddo

     sigma=sqrt(2.0*energy/sumkdist)


     do j=1,Ny
       do i=1, Nx/2+1

        kkx=kxcon(i)         !!!kkx=kx(i)  !!!   zwg :  should be checked!!!!
        kky=kycon(myid*Ny+j) !!!kky=ky(myid*Ny+j)  !!!   zwg :  should be checked!!!!

        km=sqrt(kkx**2+kky**2)

        if((i==1) .and.(j==1))then
          DTVx_C2R_in(1,1)=dcmplx(0.d0,0.d0)
          DTVy_C2R_in(1,1)=dcmplx(0.d0,0.d0)
        
       
        else

          beta=(1.0-2.0*zeta)/(km*km)
          Pxx=zeta+beta*kkx*kkx
          Pxy=     beta*kkx*kky
          Pyy=zeta+beta*kky*kky
        
        
	  A0=sigma*sqrt(   (km**7)*exp(-8.0*km/drivingWN)/((3.0*zeta -2.0)*zeta+1.0)  )	!Lester modified

          Ax = A0 * randn()
          Ay = A0 * randn()

        rand11=rand1()
        Ax_re=(Pxx*Ax + Pxy*Ay ) *cos(2*pi*rand11)
        Ax_im=(Pxx*Ax + Pxy*Ay ) *sin(2*pi*rand11)
        !rand11=rand1()
        Ay_re=(Pxy*Ax + Pyy*Ay ) *cos(2*pi*rand11)
        Ay_im=(Pxy*Ax + Pyy*Ay ) *sin(2*pi*rand11)

        DTVx_C2R_in(i,j)=dcmplx(Ax_re,Ax_im)
        DTVy_C2R_in(i,j)=dcmplx(Ay_re,Ay_im)

        endif


       enddo

     enddo


       !DTVx_C2R_in=Den_R2C_out
   call fftw_mpi_execute_dft_c2r(this%DTVx_Plan_C2R,DTVx_C2R_in,DTVx_C2R_out)
       !DTVy_C2R_in=Den_R2C_out
    call fftw_mpi_execute_dft_c2r(this%DTVy_Plan_C2R,DTVy_C2R_in,DTVy_C2R_out)

  if(myid==0)then

    write(*,*)"this Process id isBBBBBBBBBBBBBBBBBBBBBBB11111: ", myid
    do j=1,Ny
      !do i=1, Nx
         write(*,*) DTVx_C2R_out(:,j)/(nx_global*ny_global)
         write(*,*) DTVy_C2R_out(:,j)/(nx_global*ny_global)
      !enddo
      
   enddo
  
  endif


end subroutine calcDT2D

!!!! The above subroutine is added by zwg for 2D Driving Turbulence 




subroutine calcDT3D(this,Den_ini,Den_R2C_in,Den_R2C_out,Den_C2R_out,DTVx_C2R_in, &
                  & DTVx_C2R_out,DTVy_C2R_in,DTVy_C2R_out,DTVz_C2R_in,DTVz_C2R_out,DTVx,DTVy,DTVz)
use gridModule
use,intrinsic::iso_c_binding
implicit none
class(grid)::this


     double precision,dimension(this%nMesh(1),this%nMesh(2),this%nMesh(3))::Den_ini,DTVx,DTVy,DTVz

     complex(C_DOUBLE_COMPLEX),dimension(this%nMesh(1)/2+1,this%nMesh(2),this%nMesh(3)):: Den_R2C_out 
     real(C_DOUBLE),dimension((this%nMesh(1)/2+1)*2,this%nMesh(2),this%nMesh(3)):: Den_R2C_in,Den_C2R_out

     complex(C_DOUBLE_COMPLEX),dimension(this%nMesh(1)/2+1,this%nMesh(2),this%nMesh(3)):: DTVx_C2R_in,DTVy_C2R_in,DTVz_C2R_in
     real(C_DOUBLE),dimension((this%nMesh(1)/2+1)*2,this%nMesh(2),this%nMesh(3)):: DTVx_C2R_out,DTVy_C2R_out,DTVz_C2R_out

integer::nx_global,ny_global,nz_global
double precision::Lx_global,Ly_global,Lz_global
double precision::ki(this%nMesh_global(1)),kj(this%nMesh_global(2)),kk(this%nMesh_global(3))
double precision::kx(this%nMesh_global(1)),ky(this%nMesh_global(2)),kz(this%nMesh_global(3))
double precision::kxcon(this%nMesh_global(1)),kycon(this%nMesh_global(2)),kzcon(this%nMesh_global(3))
integer::i,j,k,Nx,Ny,Nz
double precision::pi,ksq,dx,dy,dz


    double precision:: sigma
    double precision:: sumkdist=0.0
    double precision:: kkx,kky,kkz,km,km1,km2

   !double precision::beta,Pxx, Pxy, Pyy
   !double precision::A0,Ax,Ay
   !double precision::Ax_re,Ax_im,Ay_re,Ay_im

    double precision::beta,Pxx, Pxy,  Pxz, Pyy, Pyz, Pzz
    double precision::A0,Ax,Ay,Az
    double precision::Ax_re,Ax_im,Ay_re,Ay_im,Az_re,Az_im

    double precision,external::randn, rand1
    double precision::rand11


    double precision::drivingWN,Energy,zeta

    drivingWN=this%drivingWN_DT
    Energy=this%Energy_DT
    !!! zeta=this%Energy_DT  
    zeta=this%zeta_DT   !!!! lester modified

    pi=4.d0*datan2(1.d0,1.d0)
    nx_global=this%nMesh_global(1)
    ny_global=this%nMesh_global(2)
    nz_global=this%nMesh_global(3)

    Nx=this%nMesh(1)
    Ny=this%nMesh(2)
    Nz=this%nMesh(3)

    Lx_global=this%rightBdry_global(1)-this%leftBdry_global(1)
    Ly_global=this%rightBdry_global(2)-this%leftBdry_global(2)
    Lz_global=this%rightBdry_global(3)-this%leftBdry_global(3)


if(mod(nz_global,nprocs) .ne. 0) then
   print *,"calcDT3D: mesh number in z-direction should be divisible to the number of processors..."
   stop
endif


if(mod(nz_global,2) .ne. 0) then
   print *,"calcDT3D: mesh number in z-direction needs to be even"
   stop
endif



  !!! zwg : for checking R2C and C2R interface of FFTW
  Den_R2C_in(1:Nx,1:Ny,1:Nz)=Den_ini(1:Nx,1:Ny,1:Nz)
  call fftw_mpi_execute_dft_r2c(this%Den_Plan_R2C,Den_R2C_in,Den_R2C_out)


  if(myid==0)then

     write(*,*)"this Process id isAAAAAAAAAAAAAAAAAAAAAAA: ", myid
    do k=1, Nz
      do j=1,Ny
        
          write(*,*) Den_R2C_out(:,j,k)
        
      enddo
        write(*,*)"=================================================="
   enddo
  
  endif



  call fftw_mpi_execute_dft_c2r(this%Den_Plan_C2R,Den_R2C_out,Den_C2R_out)

  if(myid==0)then

    write(*,*)"this Process id isBBBBBBBBBBBBBBBBBBBBBBB: ", myid
  do k=1, Nz
    do j=1,Ny
      
        write(*,*) Den_C2R_out(:,j,k)/(nx_global*ny_global*nz_global)
     
   enddo
   write(*,*)"=================================================="
  enddo
  
  endif



do i=1,nx_global
  if(i .le. nx_global/2+1) then
    ki(i)=dble(i)
  else
    ki(i)=dble(nx_global-i+2)
  endif

    !kx(i)=(ki(i)-1.d0)		!Lester modified
    kx(i)=2.d0*pi/dble(nx_global)*(ki(i)-1.d0)*dble(nx_global)/Lx_global

  if(i .le. nx_global/2+1) then  !!! zwg modified
    kxcon(i)=kx(i)
  else
    kxcon(i)=-kx(i)
  endif

enddo

do i=1,ny_global
  if(i .le. ny_global/2+1) then
    kj(i)=dble(i)
  else
    kj(i)=dble(ny_global-i+2)
  endif
   !ky(i)=(kj(i)-1.d0)		!Lester modified
    ky(i)=2.d0*pi/dble(ny_global)*(kj(i)-1.d0)*dble(ny_global)/Ly_global

  if(i .le. ny_global/2+1) then  !!! zwg modified
    kycon(i)=ky(i)
  else
    kycon(i)=-ky(i)
  endif
enddo

do i=1,nz_global
  if(i .le. nz_global/2+1) then
    kk(i)=dble(i)
  else
    kk(i)=dble(nz_global-i+2)
  endif
   !kz(i)=(kk(i)-1.d0)		!Lester modified
    kz(i)=2.d0*pi/dble(nz_global)*(kk(i)-1.d0)*dble(nz_global)/Lz_global

  if(i .le. nz_global/2+1) then !!! zwg modified
    kzcon(i)=kz(i)
  else
    kzcon(i)=-kz(i)
  endif
enddo


    do k=1, nz_global
      do j=1, ny_global
        do i=1, nx_global

          km=sqrt(kxcon(i)**2+kycon(j)**2+kzcon(k)**2)
          sumkdist=sumkdist+(km**6)*exp(-8.0*km/drivingWN)

        enddo
     enddo
   enddo

  sigma=sqrt(2.0*energy/sumkdist)




   do k=1,Nz
     do j=1,Ny
       do i=1, Nx/2+1

        kkx=kxcon(i)         !!!kkx=kx(i)        !!!   zwg :  should be checked!!!!
        kky=kycon(j)         !!!kky=ky(j)        !!!   zwg :  should be checked!!!!
        kkz=kzcon(myid*Nz+k) !kkz=kz(myid*Nz+k)  !!!   zwg :  should be checked!!!!

        km=sqrt(kkx**2+kky**2+kkz**2)

        if((i==1) .and.(j==1) .and.(k==1) )then
        DTVx_C2R_in(1,1,1)=dcmplx(0.d0,0.d0)
        DTVy_C2R_in(1,1,1)=dcmplx(0.d0,0.d0)
        DTVz_C2R_in(1,1,1)=dcmplx(0.d0,0.d0)
       
        else

        beta=(1.0-2.0*zeta)/(km*km)




        Pxx=zeta+beta*kkx*kkx
        Pxy=     beta*kkx*kky
        Pxz=     beta*kkx*kkz

        Pyy=zeta+beta*kky*kky
        Pyz=     beta*kky*kkz
        Pzz=zeta+beta*kkz*kkz

        A0=sigma*sqrt( (km**6)*exp(-8.0*km/drivingWN)/((3.0*zeta -2.0)*zeta+1.0)  )

        Ax = A0 * randn()
        Ay = A0 * randn()
        Az = A0 * randn()


        rand11=rand1()
        Ax_re=(Pxx*Ax + Pxy*Ay + Pxz*Az) *cos(2*pi*rand11)
        Ax_im=(Pxx*Ax + Pxy*Ay + Pxz*Az) *sin(2*pi*rand11)
        !rand11=rand1()
        Ay_re=(Pxy*Ax + Pyy*Ay + Pyz*Az) *cos(2*pi*rand11)
        Ay_im=(Pxy*Ax + Pyy*Ay + Pyz*Az) *sin(2*pi*rand11) 
        !rand11=rand1()
        Az_re=(Pxz*Ax + Pyz*Ay + Pzz*Az) *cos(2*pi*rand11)
        Az_im=(Pxz*Ax + Pyz*Ay + Pzz*Az) *sin(2*pi*rand11) 

        DTVx_C2R_in(i,j,k)=dcmplx(Ax_re,Ax_im)
        DTVy_C2R_in(i,j,k)=dcmplx(Ay_re,Ay_im)
        DTVz_C2R_in(i,j,k)=dcmplx(Az_re,Az_im)

        !DTVxCmplx(i,j,k)=dcmplx(Ax_re,Ax_im)
        !DTVyCmplx(i,j,k)=dcmplx(Ay_re,Ay_im)
        !DTVzCmplx(i,j,k)=dcmplx(Az_re,Az_im)

        endif



       enddo
     enddo
   enddo



       DTVx_C2R_in=Den_R2C_out
    call fftw_mpi_execute_dft_c2r(this%DTVx_Plan_C2R,DTVx_C2R_in,DTVx_C2R_out)
       DTVy_C2R_in=Den_R2C_out
    call fftw_mpi_execute_dft_c2r(this%DTVy_Plan_C2R,DTVy_C2R_in,DTVy_C2R_out)
       DTVz_C2R_in=Den_R2C_out
    call fftw_mpi_execute_dft_c2r(this%DTVz_Plan_C2R,DTVz_C2R_in,DTVz_C2R_out)

  if(myid==0)then

    write(*,*)"this Process id isBBBBBBBBBBBBBBBBBBBBBBB11111: ", myid
   do k=1,Nz
    do j=1,Ny
      !do i=1, Nx
         write(*,*) DTVx_C2R_out(:,j,k)/(nx_global*ny_global*nz_global)
         !write(*,*) DTVy_C2R_out(:,j,k)/(nx_global*ny_global*nz_global)
         !write(*,*) DTVz_C2R_out(:,j,k)/(nx_global*ny_global*nz_global)
      !enddo
      
   enddo
  enddo
  
  endif


  if(myid==0)then

    write(*,*)"this Process id isBBBBBBBBBBBBBBBBBBBBBBB22222: ", myid
   do k=1,Nz
    do j=1,Ny
      !do i=1, Nx
         !write(*,*) DTVx_C2R_out(:,j,k)/(nx_global*ny_global*nz_global)
          write(*,*) DTVy_C2R_out(:,j,k)/(nx_global*ny_global*nz_global)
         !write(*,*) DTVz_C2R_out(:,j,k)/(nx_global*ny_global*nz_global)
      !enddo
      
   enddo
  enddo
  
  endif

  if(myid==0)then

    write(*,*)"this Process id isBBBBBBBBBBBBBBBBBBBBBBB33333: ", myid
   do k=1,Nz
    do j=1,Ny
      !do i=1, Nx
         !write(*,*) DTVx_C2R_out(:,j,k)/(nx_global*ny_global*nz_global)
          write(*,*) DTVy_C2R_out(:,j,k)/(nx_global*ny_global*nz_global)
         !write(*,*) DTVz_C2R_out(:,j,k)/(nx_global*ny_global*nz_global)
      !enddo
      
   enddo
  enddo
  
  endif



end subroutine calcDT3D 



!!!! The following subroutine is added by zwg for 2D Driving Turbulence (Modified version)

subroutine calcDT2D_MD(this,DTVx_C2R_in,DTVx_C2R_out,DTVy_C2R_in,DTVy_C2R_out,DTVx,DTVy,q)
use gridModule
use,intrinsic::iso_c_binding
implicit none
class(grid)::this 

!!!zwg added for modification
double precision,dimension(1-this%nbuf:this%nMesh(1)+this%nbuf,1-this%nbuf:this%nMesh(2)+this%nbuf,this%nvar)::q
	
     !double precision,dimension(this%nMesh(1),this%nMesh(2))::Den_ini,DTVx,DTVy
      double precision,dimension(this%nMesh(1),this%nMesh(2))::        DTVx,DTVy
      double precision,dimension(this%nMesh(1),this%nMesh(2))::        BGVx,BGVy,BGrho   !!!!background velocities and density
     !complex(C_DOUBLE_COMPLEX),dimension(this%nMesh(1)/2+1,this%nMesh(2)):: Den_R2C_out 
     !real(C_DOUBLE),dimension((this%nMesh(1)/2+1)*2,this%nMesh(2)):: Den_R2C_in,Den_C2R_out

     complex(C_DOUBLE_COMPLEX),dimension(this%nMesh(1)/2+1,this%nMesh(2)):: DTVx_C2R_in,DTVy_C2R_in
     real(C_DOUBLE),dimension((this%nMesh(1)/2+1)*2,this%nMesh(2)):: DTVx_C2R_out,DTVy_C2R_out



integer::nx_global,ny_global
double precision::Lx_global,Ly_global

double precision::ki(this%nMesh_global(1)),kj(this%nMesh_global(2))

!!! zwg added for modification
integer::kig(this%nMesh_global(1)),kjg(this%nMesh_global(2))

double precision::kx(this%nMesh_global(1)),ky(this%nMesh_global(2))
double precision::kxcon(this%nMesh_global(1)),kycon(this%nMesh_global(2))


integer::i,j,Nx,Ny
double precision::pi,ksq,dx,dy
integer::ierr


    double precision:: sigma
    double precision:: sumkdist=0.0
    double precision:: kkx,kky,km,km1,km2

!!! zwg added for modification
    integer::kkx1,kky1

    double precision::beta,Pxx, Pxy, Pyy
    double precision::A0,Ax,Ay
    double precision::Ax_re,Ax_im,Ay_re,Ay_im

    double precision,external::randn, rand1
    double precision::rand11

    double precision::drivingWN,Energy,zeta
    double precision::netmomx,netmomy	!!!

    double precision::momx_sum1,momx_sum,momy_sum1,momy_sum,rho_sum1,rho_sum	!!!
    double precision::energy1_sum1,energy1_sum, energy2_sum1, energy2_sum, alpha 

    double precision::momx_check,momy_check,momx_check_sum,momy_check_sum

    

   
   

    drivingWN=this%drivingWN_DT
    Energy=this%Energy_DT*this%DTenergyfaction  !!! zwg: modifid by zwg!!!
    !!!zeta=this%Energy_DT !!! zwg: modifid by lester
    zeta=this%zeta_DT
    netmomx=this%netmomx_DT		
    netmomy=this%netmomy_DT		


    
    pi=4.d0*datan2(1.d0,1.d0)
    nx_global=this%nMesh_global(1)
    ny_global=this%nMesh_global(2)

    Nx=this%nMesh(1)
    Ny=this%nMesh(2)

    Lx_global=this%rightBdry_global(1)-this%leftBdry_global(1)
    Ly_global=this%rightBdry_global(2)-this%leftBdry_global(2)

    if(mod(ny_global,nprocs) .ne. 0) then
      print *,"calcDT2D: mesh number in y-direction should be divisible to the number of processors..."
      stop
    endif


    if(mod(ny_global,2) .ne. 0) then
     print *,"calcDT2D: mesh number in y-direction needs to be even"
     stop
    endif





do i=1,nx_global
  if(i .le. nx_global/2+1) then
    ki(i)=dble(i)
  else
    ki(i)=dble(nx_global-i+2)
  endif
    kx(i)=(ki(i)-1.d0)		!Lester modified
   !kx(i)=2.d0*pi/dble(nx_global)*(ki(i)-1.d0)*dble(nx_global)/Lx_global

  if(i .le. nx_global/2+1) then  !!!zwg modified
    kxcon(i)=kx(i)
  else
    kxcon(i)=-kx(i)
  endif
!!!zwg modify
kig(i)=i
enddo

do i=1,ny_global
  if(i .le. ny_global/2+1) then
    kj(i)=dble(i)
  else
    kj(i)=dble(ny_global-i+2)
  endif
    ky(i)=(kj(i)-1.d0)	!Lester modified
   !ky(i)=2.d0*pi/dble(ny_global)*(kj(i)-1.d0)*dble(ny_global)/Ly_global

  if(i .le. ny_global/2+1) then  !!!zwg modified
    kycon(i)=ky(i)
  else
    kycon(i)=-ky(i)
  endif

!!!zwg modify
kjg(i)=i
enddo



      do j=1, ny_global
        do i=1,nx_global

          km=sqrt(kxcon(i)**2+kycon(j)**2)
          
	  sumkdist=sumkdist+(km**7)*exp(-8.0*km/drivingWN)	!Lester modified
        !!sumkdist=sumkdist+(km**6)*exp(-8.0*km/drivingWN)

        enddo
     enddo

     sigma=sqrt(2.0*energy/sumkdist)




     do j=1,Ny
       do i=1, Nx/2+1

        kkx=kxcon(i)         !!!kkx=kx(i)  !!!   zwg :  should be checked!!!!
        kky=kycon(myid*Ny+j) !!!kky=ky(myid*Ny+j)  !!!   zwg :  should be checked!!!!

        !!!zwg added for modification
        kkx1=kig(i)
        kky1=kjg(myid*Ny+j)
        

        km=sqrt(kkx**2+kky**2)

        !!!if((i==1) .and.(j==1))then  !!zwg: a bug should modified as below?
        !!!if((kkx1==1) .and.(kky1==1))then  !!!or 
        if((i==1) .and.(myid*Ny+j==1))then
          DTVx_C2R_in(1,1)=dcmplx(0.d0,0.d0)
          DTVy_C2R_in(1,1)=dcmplx(0.d0,0.d0)
        else

          beta=(1.0-2.0*zeta)/(km*km)
          Pxx=zeta+beta*kkx*kkx
          Pxy=     beta*kkx*kky
          Pyy=zeta+beta*kky*kky
        
        
	  !!!!A0=sigma*sqrt(   (km**7)*exp(-8.0*km/drivingWN)/((3.0*zeta -2.0)*zeta+1.0)  )	!Lester modified
          
          !!!!see  MNRAS 466, 2272¨C2283 (2017), the above is for 3D and the following is for 2D
          A0=sigma*sqrt(   (km**7)*exp(-8.0*km/drivingWN)/((2.0*zeta -2.0)*zeta+1.0)  )	        !zwg  modified

          Ax = A0 * randn()
          Ay = A0 * randn()

          !Ax = A0 !* randn()
          !Ay = A0 !* randn()

        rand11=rand1()

        !rand11=0.0

        Ax_re=(Pxx*Ax + Pxy*Ay ) *cos(2*pi*rand11)
        Ax_im=(Pxx*Ax + Pxy*Ay ) *sin(2*pi*rand11)
        !rand11=rand1()
        Ay_re=(Pxy*Ax + Pyy*Ay ) *cos(2*pi*rand11)
        Ay_im=(Pxy*Ax + Pyy*Ay ) *sin(2*pi*rand11)

        DTVx_C2R_in(i,j)=dcmplx(Ax_re,Ax_im)
        DTVy_C2R_in(i,j)=dcmplx(Ay_re,Ay_im)

        endif


       enddo

     enddo


     
    call fftw_mpi_execute_dft_c2r(this%DTVx_Plan_C2R,DTVx_C2R_in,DTVx_C2R_out)
    call fftw_mpi_execute_dft_c2r(this%DTVy_Plan_C2R,DTVy_C2R_in,DTVy_C2R_out)


    !!!zwg:: need to scale here
    do i=1,Nx
     do j=1,Ny
      DTVx(i,j)=DTVx_C2R_out(i,j)/(nx_global*ny_global)  !!!zwg::need to scale here???
      DTVy(i,j)=DTVy_C2R_out(i,j)/(nx_global*ny_global)  !!!zwg::need to scale here???
     enddo
    enddo


!!!!!!!!!!!!!!Start momentum shifting!!!!!!!!!!
    do i=1,Nx
     do j=1,Ny
      BGrho(i,j)=q(i,j,1)
      BGVx(i,j) =q(i,j,2)/q(i,j,1)
      BGVy(i,j) =q(i,j,3)/q(i,j,1)
     enddo
    enddo

    !double precision::momx_sum1,momx_sum,momy_sum1,momy_sum,rho_sum1,rho_sum	!!!

    momx_sum1=0.0
    momy_sum1=0.0
     rho_sum1=0.0
   
    do i=1,Nx	
      do j=1,Ny
        momx_sum1=momx_sum1+DTVx(i,j)*BGrho(i,j)
        momy_sum1=momy_sum1+DTVy(i,j)*BGrho(i,j)
         rho_sum1= rho_sum1+BGrho(i,j)
      enddo
    enddo

    momx_sum=0.0
    momy_sum=0.0
     rho_sum=0.0

    call mpi_allreduce(momx_sum1,momx_sum,1,MPI_DOUBLE_PRECISION,MPI_sum,MPI_COMM_WORLD,ierr)
    call mpi_allreduce(momy_sum1,momy_sum,1,MPI_DOUBLE_PRECISION,MPI_sum,MPI_COMM_WORLD,ierr)
    call mpi_allreduce(rho_sum1,rho_sum,1,MPI_DOUBLE_PRECISION,MPI_sum,MPI_COMM_WORLD,ierr)

    momx_sum=momx_sum-netmomx
    momy_sum=momy_sum-netmomy

    do i=1,Nx
      do j=1,Ny
        DTVx(i,j)=DTVx(i,j)-momx_sum/rho_sum
        DTVy(i,j)=DTVy(i,j)-momy_sum/rho_sum
      enddo
    enddo


   !momx_check=0.0
   !momy_check=0.0

    !do i=1,Nx
      !do j=1,Ny
        !momx_check=momx_check+DTVx(i,j)*BGrho(i,j)
        !momy_check=momy_check+DTVy(i,j)*BGrho(i,j)
      !enddo
    !enddo


    !momx_check_sum=1.0
    !momy_check_sum=1.0

    !call mpi_allreduce(momx_check,momx_check_sum,1,MPI_DOUBLE_PRECISION,MPI_sum,MPI_COMM_WORLD,ierr)
    !call mpi_allreduce(momy_check,momy_check_sum,1,MPI_DOUBLE_PRECISION,MPI_sum,MPI_COMM_WORLD,ierr)
    !write(*,*)"momx_check_sum= ", momx_check_sum,"momy_check_sum= ", momy_check_sum
    !stop




!!!!!!!!!!!!!!!!Momentum shift finished!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!Start energy shift!!!!!!!!!!!!!!!!!!!!!
    !double precision::energy1_sum1,energy1_sum, energy2_sum1, energy2_sum, alpha 
    energy1_sum1=0
    energy2_sum1=0
    do i=1,Nx
      do j=1,Ny
        energy1_sum1=energy1_sum1+BGrho(i,j)*(DTVx(i,j)*DTVx(i,j)+DTVy(i,j)*DTVy(i,j))
        energy2_sum1=energy2_sum1+BGrho(i,j)*(DTVx(i,j)*BGVx(i,j)+DTVy(i,j)*BGVy(i,j))
      enddo
    enddo

    energy1_sum=0
    energy2_sum=0

    call mpi_allreduce(energy1_sum1,energy1_sum,1,MPI_DOUBLE_PRECISION,MPI_sum,MPI_COMM_WORLD,ierr)
    call mpi_allreduce(energy2_sum1,energy2_sum,1,MPI_DOUBLE_PRECISION,MPI_sum,MPI_COMM_WORLD,ierr)

    
    
    
    if(energy1_sum==0)then
      alpha=0.0d0
    else
      alpha=-(energy2_sum/energy1_sum)+sqrt((energy2_sum/energy1_sum)**2+2*Energy*rho_sum/energy1_sum)
    endif


!write(*,*)alpha
!stop

    do i=1,Nx
      do j=1,Ny
        DTVx(i,j)=alpha*DTVx(i,j)	
        DTVy(i,j)=alpha*DTVy(i,j)
      enddo
    enddo


    do i=1,Nx
      do j=1,Ny
        q(i,j,2)=q(i,j,2)+BGrho(i,j)*DTVx(i,j)
        q(i,j,3)=q(i,j,3)+BGrho(i,j)*DTVy(i,j)
      enddo
    enddo
!!!!!!!!!!!!!!!!Energy shift finished!!!!!!!!!!!!!!!!!!!!!


end subroutine calcDT2D_MD
!!!! The above subroutine is added by zwg for 2D Driving Turbulence (Modified version)



!!!! The following subroutine is added by zwg for 3D Driving Turbulence (Modified version)

subroutine calcDT3D_MD(this,DTVx_C2R_in,DTVx_C2R_out,DTVy_C2R_in,DTVy_C2R_out,DTVz_C2R_in,DTVz_C2R_out,DTVx,DTVy,DTVz,q)

use gridModule
use,intrinsic::iso_c_binding
implicit none
class(grid)::this

     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     !!!zwg added for modification
     double precision,dimension(1-this%nbuf:this%nMesh(1)+this%nbuf,1-this%nbuf:this%nMesh(2)+this%nbuf,&
                           & 1-this%nbuf:this%nMesh(3)+this%nbuf,this%nvar)::q

    !double precision,dimension(this%nMesh(1),this%nMesh(2),this%nMesh(3))::Den_ini,DTVx,DTVy,DTVz
     double precision,dimension(this%nMesh(1),this%nMesh(2),this%nMesh(3))::        DTVx,DTVy,DTVz
     double precision,dimension(this%nMesh(1),this%nMesh(2),this%nMesh(3))::        BGVx,BGVy,BGVz,BGrho   !!!!background velocities and density 


     !complex(C_DOUBLE_COMPLEX),dimension(this%nMesh(1)/2+1,this%nMesh(2),this%nMesh(3)):: Den_R2C_out 
     !real(C_DOUBLE),dimension((this%nMesh(1)/2+1)*2,this%nMesh(2),this%nMesh(3)):: Den_R2C_in,Den_C2R_out

     complex(C_DOUBLE_COMPLEX),dimension(this%nMesh(1)/2+1,this%nMesh(2),this%nMesh(3)):: DTVx_C2R_in,DTVy_C2R_in,DTVz_C2R_in
     real(C_DOUBLE),dimension((this%nMesh(1)/2+1)*2,this%nMesh(2),this%nMesh(3)):: DTVx_C2R_out,DTVy_C2R_out,DTVz_C2R_out

integer::nx_global,ny_global,nz_global
double precision::Lx_global,Ly_global,Lz_global
double precision::ki(this%nMesh_global(1)),kj(this%nMesh_global(2)),kk(this%nMesh_global(3))

!!! zwg added for modification
integer::kig(this%nMesh_global(1)),kjg(this%nMesh_global(2)),kkg(this%nMesh_global(3))

double precision::kx(this%nMesh_global(1)),ky(this%nMesh_global(2)),kz(this%nMesh_global(3))
double precision::kxcon(this%nMesh_global(1)),kycon(this%nMesh_global(2)),kzcon(this%nMesh_global(3))
integer::i,j,k,Nx,Ny,Nz
double precision::pi,ksq,dx,dy,dz
integer::ierr


    double precision:: sigma
    double precision:: sumkdist=0.0
    double precision:: kkx,kky,kkz,km,km1,km2

!!! zwg added for modification
    integer::kkx1,kky1,kkz1  

    double precision::beta,Pxx, Pxy,  Pxz, Pyy, Pyz, Pzz
    double precision::A0,Ax,Ay,Az
    double precision::Ax_re,Ax_im,Ay_re,Ay_im,Az_re,Az_im

    double precision,external::randn, rand1
    double precision::rand11


    double precision::drivingWN,Energy,zeta
    double precision::netmomx,netmomy,netmomz	!!!

    double precision::momx_sum1,momx_sum,momy_sum1,momy_sum,momz_sum1,momz_sum,rho_sum1,rho_sum	!!!
    double precision::energy1_sum1,energy1_sum, energy2_sum1, energy2_sum, alpha 

    double precision::momx_check,momy_check,momz_check,momx_check_sum,momy_check_sum,momz_check_sum



    drivingWN=this%drivingWN_DT
    Energy=this%Energy_DT*this%DTenergyfaction  !!! zwg: modifid by zwg!!!
    !!!zeta=this%Energy_DT !!! zwg: modifid by lester
    zeta=this%zeta_DT

  !write(*,*)this%enable_DT,this%DT_mode
  !write(*,*)"k0= ", this%drivingWN_DT
  !write(*,*)"energy= ", this%Energy_DT*this%DTenergyfaction
  !write(*,*)"zeta= ", this%zeta_DT
  !stop

    netmomx=this%netmomx_DT		
    netmomy=this%netmomy_DT
    netmomz=this%netmomz_DT	



    pi=4.d0*datan2(1.d0,1.d0)
    nx_global=this%nMesh_global(1)
    ny_global=this%nMesh_global(2)
    nz_global=this%nMesh_global(3)

    Nx=this%nMesh(1)
    Ny=this%nMesh(2)
    Nz=this%nMesh(3)

    Lx_global=this%rightBdry_global(1)-this%leftBdry_global(1)
    Ly_global=this%rightBdry_global(2)-this%leftBdry_global(2)
    Lz_global=this%rightBdry_global(3)-this%leftBdry_global(3)


if(mod(nz_global,nprocs) .ne. 0) then
   print *,"calcDT3D: mesh number in z-direction should be divisible to the number of processors..."
   stop
endif


if(mod(nz_global,2) .ne. 0) then
   print *,"calcDT3D: mesh number in z-direction needs to be even"
   stop
endif



do i=1,nx_global
  if(i .le. nx_global/2+1) then
    ki(i)=dble(i)
  else
    ki(i)=dble(nx_global-i+2)
  endif

    kx(i)=(ki(i)-1.d0)		!Lester modified
    !kx(i)=2.d0*pi/dble(nx_global)*(ki(i)-1.d0)*dble(nx_global)/Lx_global

  if(i .le. nx_global/2+1) then  !!! zwg modified
    kxcon(i)=kx(i)
  else
    kxcon(i)=-kx(i)
  endif
 !!!! zwg modify
 kig(i)=i
enddo

do i=1,ny_global
  if(i .le. ny_global/2+1) then
    kj(i)=dble(i)
  else
    kj(i)=dble(ny_global-i+2)
  endif
    ky(i)=(kj(i)-1.d0)		!Lester modified
    !ky(i)=2.d0*pi/dble(ny_global)*(kj(i)-1.d0)*dble(ny_global)/Ly_global

  if(i .le. ny_global/2+1) then  !!! zwg modified
    kycon(i)=ky(i)
  else
    kycon(i)=-ky(i)
  endif
 !!!! zwg modify
kjg(i)=i
enddo

do i=1,nz_global
  if(i .le. nz_global/2+1) then
    kk(i)=dble(i)
  else
    kk(i)=dble(nz_global-i+2)
  endif
   kz(i)=(kk(i)-1.d0)		!Lester modified
   !kz(i)=2.d0*pi/dble(nz_global)*(kk(i)-1.d0)*dble(nz_global)/Lz_global

  if(i .le. nz_global/2+1) then !!! zwg modified
    kzcon(i)=kz(i)
  else
    kzcon(i)=-kz(i)
  endif
 !!!! zwg modify
 kkg(i)=i
enddo

   !drivingWN=2.d0*pi*drivingWN/Lz_global   !!!!zwg added for kz(i)=2.d0*pi/dble(nz_global)*(kk(i)-1.d0)*dble(nz_global)/Lz_global 

    do k=1, nz_global
      do j=1, ny_global
        do i=1, nx_global

          km=sqrt(kxcon(i)**2+kycon(j)**2+kzcon(k)**2)
          sumkdist=sumkdist+(km**6)*exp(-8.0*km/drivingWN)

        enddo
     enddo
   enddo

  sigma=sqrt(2.0*energy/sumkdist)


!write(*,*)"sigma= ",sigma, Nx,Ny,Nz
!call MPI_BARRIER(MPI_COMM_WORLD,ierr)
!stop

   do k=1,Nz
     do j=1,Ny
       do i=1, Nx/2+1

        kkx=kxcon(i)         !!!kkx=kx(i)        !!!   zwg :  should be checked!!!!
        kky=kycon(j)         !!!kky=ky(j)        !!!   zwg :  should be checked!!!!
        kkz=kzcon(myid*Nz+k) !kkz=kz(myid*Nz+k)  !!!   zwg :  should be checked!!!!

        !!!zwg added for modification
        kkx1=kig(i)
        kky1=kjg(j)
        kkz1=kkg(myid*Nz+k)
        !write(*,*)"myid= ",myid, kkx1,kky1,kkz1
        !stop
        
        km=sqrt(kkx**2+kky**2+kkz**2)

        !write(*,*)"km= ", km
        !stop
        
        !if((i==1) .and.(j==1) .and.(k==1) )then  !!!zwg : a bug should modified as following !!! 
        !if((kkx1==1) .and.(kky1==1) .and.(kkz1==1) )then  !!! or
        if((i==1) .and.(j==1) .and.(myid*Nz+k==1) )then
          DTVx_C2R_in(1,1,1)=dcmplx(0.d0,0.d0)
          DTVy_C2R_in(1,1,1)=dcmplx(0.d0,0.d0)
          DTVz_C2R_in(1,1,1)=dcmplx(0.d0,0.d0)
       
        else

        beta=(1.0-2.0*zeta)/(km*km)




        Pxx=zeta+beta*kkx*kkx
        Pxy=     beta*kkx*kky
        Pxz=     beta*kkx*kkz

        Pyy=zeta+beta*kky*kky
        Pyz=     beta*kky*kkz
        Pzz=zeta+beta*kkz*kkz

        A0=sigma*sqrt( (km**6)*exp(-8.0*km/drivingWN)/((3.0*zeta -2.0)*zeta+1.0)  )

        Ax = A0 * randn()
        Ay = A0 * randn()
        Az = A0 * randn()

        !Ax = A0 
        !Ay = A0 
        !Az = A0 

        rand11=rand1()
      
        !rand11=0.0

        Ax_re=(Pxx*Ax + Pxy*Ay + Pxz*Az) *cos(2*pi*rand11)
        Ax_im=(Pxx*Ax + Pxy*Ay + Pxz*Az) *sin(2*pi*rand11)
        !rand11=rand1()
        Ay_re=(Pxy*Ax + Pyy*Ay + Pyz*Az) *cos(2*pi*rand11)
        Ay_im=(Pxy*Ax + Pyy*Ay + Pyz*Az) *sin(2*pi*rand11) 
        !rand11=rand1()
        Az_re=(Pxz*Ax + Pyz*Ay + Pzz*Az) *cos(2*pi*rand11)
        Az_im=(Pxz*Ax + Pyz*Ay + Pzz*Az) *sin(2*pi*rand11) 

        DTVx_C2R_in(i,j,k)=dcmplx(Ax_re,Ax_im)
        DTVy_C2R_in(i,j,k)=dcmplx(Ay_re,Ay_im)
        DTVz_C2R_in(i,j,k)=dcmplx(Az_re,Az_im)

        !DTVxCmplx(i,j,k)=dcmplx(Ax_re,Ax_im)
        !DTVyCmplx(i,j,k)=dcmplx(Ay_re,Ay_im)
        !DTVzCmplx(i,j,k)=dcmplx(Az_re,Az_im)

        endif



       enddo
     enddo
   enddo



       
    call fftw_mpi_execute_dft_c2r(this%DTVx_Plan_C2R,DTVx_C2R_in,DTVx_C2R_out) 
    call fftw_mpi_execute_dft_c2r(this%DTVy_Plan_C2R,DTVy_C2R_in,DTVy_C2R_out)
    call fftw_mpi_execute_dft_c2r(this%DTVz_Plan_C2R,DTVz_C2R_in,DTVz_C2R_out)

    do i=1,Nx
      do j=1,Ny
        do k=1,Nz
           DTVx(i,j,k)=DTVx_C2R_out(i,j,k)/(nx_global*ny_global*nz_global)  !!!zwg::need to scale here???  
           DTVy(i,j,k)=DTVy_C2R_out(i,j,k)/(nx_global*ny_global*nz_global)  !!!zwg::need to scale here???
           DTVz(i,j,k)=DTVz_C2R_out(i,j,k)/(nx_global*ny_global*nz_global)  !!!zwg::need to scale here???
           !write(*,*) DTVx(i,j,k)
        enddo
      enddo
    enddo
!stop
!!!!!!!!!!!!!!Start momentum shifting!!!!!!!!!!

    do i=1,Nx
      do j=1,Ny
        do k=1,Nz
          BGrho(i,j,k)=q(i,j,k,1)           
          BGVx(i,j,k)=q(i,j,k,2)/q(i,j,k,1)    
          BGVy(i,j,k)=q(i,j,k,3)/q(i,j,k,1)
          BGVz(i,j,k)=q(i,j,k,4)/q(i,j,k,1)
        enddo
      enddo
    enddo



    momx_sum1=0.0
    momy_sum1=0.0
    momz_sum1=0.0
     rho_sum1=0.0
   
    do i=1,Nx	
      do j=1,Ny
        do k=1,Nz
        momx_sum1=momx_sum1+DTVx(i,j,k)*BGrho(i,j,k)
        momy_sum1=momy_sum1+DTVy(i,j,k)*BGrho(i,j,k)
        momz_sum1=momz_sum1+DTVz(i,j,k)*BGrho(i,j,k)
         rho_sum1= rho_sum1+BGrho(i,j,k)
         enddo
      enddo
    enddo

    momx_sum=0.0
    momy_sum=0.0
    momz_sum=0.0
     rho_sum=0.0

    call mpi_allreduce(momx_sum1,momx_sum,1,MPI_DOUBLE_PRECISION,MPI_sum,MPI_COMM_WORLD,ierr)
    call mpi_allreduce(momy_sum1,momy_sum,1,MPI_DOUBLE_PRECISION,MPI_sum,MPI_COMM_WORLD,ierr)
    call mpi_allreduce(momz_sum1,momz_sum,1,MPI_DOUBLE_PRECISION,MPI_sum,MPI_COMM_WORLD,ierr)
    call mpi_allreduce(rho_sum1,rho_sum,1,MPI_DOUBLE_PRECISION,MPI_sum,MPI_COMM_WORLD,ierr)

    momx_sum=momx_sum-netmomx
    momy_sum=momy_sum-netmomy
    momz_sum=momz_sum-netmomz

!write(*,*)rho_sum, momx_sum,momy_sum,momz_sum
!stop


    do i=1,Nx	
      do j=1,Ny
        do k=1,Nz

          DTVx(i,j,k)=DTVx(i,j,k)-momx_sum/rho_sum
          DTVy(i,j,k)=DTVy(i,j,k)-momy_sum/rho_sum
          DTVz(i,j,k)=DTVz(i,j,k)-momz_sum/rho_sum

         enddo
      enddo
    enddo



   !momx_check=0.0
   !momy_check=0.0
   !momz_check=0.0

    !do i=1,Nx
      !do j=1,Ny
         !do k=1,Nz
           !momx_check=momx_check+DTVx(i,j,k)*BGrho(i,j,k)
           !momy_check=momy_check+DTVy(i,j,k)*BGrho(i,j,k)
           !momz_check=momz_check+DTVz(i,j,k)*BGrho(i,j,k)
         !enddo
      !enddo
    !enddo


    !momx_check_sum=1.0
    !momy_check_sum=1.0
    !momz_check_sum=1.0

    !call mpi_allreduce(momx_check,momx_check_sum,1,MPI_DOUBLE_PRECISION,MPI_sum,MPI_COMM_WORLD,ierr)
    !call mpi_allreduce(momy_check,momy_check_sum,1,MPI_DOUBLE_PRECISION,MPI_sum,MPI_COMM_WORLD,ierr)
    !call mpi_allreduce(momz_check,momz_check_sum,1,MPI_DOUBLE_PRECISION,MPI_sum,MPI_COMM_WORLD,ierr)
    !write(*,*)"momx_check_sum= ", momx_check_sum, &
    !         & "momy_check_sum= ", momy_check_sum, &
    !         & "momz_check_sum= ", momz_check_sum
    !stop

!!!!!!!!!!!!!!Momentum shifting done!!!!!!!!!!


!!!!!!!!!!!!!!!!Start energy shift!!!!!!!!!!!!!!!!!!!!!
    !double precision::energy1_sum1,energy1_sum, energy2_sum1, energy2_sum, alpha 
    energy1_sum1=0
    energy2_sum1=0
    do i=1,Nx	
      do j=1,Ny
        do k=1,Nz
          energy1_sum1=energy1_sum1+BGrho(i,j,k)*(DTVx(i,j,k)*DTVx(i,j,k)+ &
                                                & DTVy(i,j,k)*DTVy(i,j,k)+ &
                                                & DTVz(i,j,k)*DTVz(i,j,k))

          energy2_sum1=energy2_sum1+BGrho(i,j,k)*(DTVx(i,j,k)*BGVx(i,j,k)+ &
                                                & DTVy(i,j,k)*BGVy(i,j,k)+ &
                                                & DTVz(i,j,k)*BGVz(i,j,k))
        enddo
      enddo
    enddo

    energy1_sum=0
    energy2_sum=0

    call mpi_allreduce(energy1_sum1,energy1_sum,1,MPI_DOUBLE_PRECISION,MPI_sum,MPI_COMM_WORLD,ierr)
    call mpi_allreduce(energy2_sum1,energy2_sum,1,MPI_DOUBLE_PRECISION,MPI_sum,MPI_COMM_WORLD,ierr)

    
    
    
    if(energy1_sum==0)then
      alpha=0.0d0
    else
      alpha=-(energy2_sum/energy1_sum)+sqrt((energy2_sum/energy1_sum)**2+2*Energy*rho_sum/energy1_sum)
    endif

!write(*,*)alpha
!stop


    do i=1,Nx	
      do j=1,Ny
        do k=1,Nz

          DTVx(i,j,k)=alpha*DTVx(i,j,k)	
          DTVy(i,j,k)=alpha*DTVy(i,j,k)
          DTVz(i,j,k)=alpha*DTVz(i,j,k)
          !write(*,*)DTVx(i,j,k)

        enddo
      enddo
    enddo
!stop

    do i=1,Nx	
      do j=1,Ny
        do k=1,Nz
          q(i,j,k,2)=q(i,j,k,2)+BGrho(i,j,k)*DTVx(i,j,k)
          q(i,j,k,3)=q(i,j,k,3)+BGrho(i,j,k)*DTVy(i,j,k)
          q(i,j,k,4)=q(i,j,k,4)+BGrho(i,j,k)*DTVz(i,j,k)

          !!q(i,j,k,4)=q(i,j,k,4)+BGrho(i,j,k)*DTVy(i,j,k)  !!!zwg a bug!!!
        enddo
      enddo
    enddo
!!!!!!!!!!!!!!!!Energy shift done!!!!!!!!!!!!!!!!!!!!!


end subroutine calcDT3D_MD 

!!!! The above subroutine is added by zwg for 3D Driving Turbulence (Modified version)




!-----------------------------------------------------------------------
function randn()
!-----------------------------------------------------------------------
! 
! Returns a random number which is (standard) normal-distributed
! (mean = 0.0, variance = 1.0). This algorithm uses the Box-Muller
! method, which computes two independent normal-distributed numbers
! x,y from two random numbers uniformly distributed on (0,1].
!
!-----------------------------------------------------------------------
   !use real_prec
   !use param
   implicit none

   real*8            :: randn
   real*8            :: u,v,p,q
   integer,save        :: stage = 1
   real*8,save       :: x,y
   real*8:: pi2

   pi2=8.d0*datan2(1.d0,1.d0)

   if (stage==1) then
      ! Random number u must be in (0,1], but Fortran's RNG returns [0,1).
      call random_number(u)
      if (u==0.d0) u=1.d0
      ! Random number v can be from (0,1] or [0,1), it doesn't matter.
      call random_number(v)
      p = sqrt(-2.d0*log(u))
      q = pi2*v
      x = p*cos(q)
      y = p*sin(q)
      randn = x
      stage = 2
   else
      randn = y
      stage = 1
   endif
   return

end function randn


!-----------------------------------------------------------------------
function rand1()
!-----------------------------------------------------------------------
   !use real_prec
   implicit none

   real*8 :: rand1

   call random_number(rand1)
   return

end function rand1











