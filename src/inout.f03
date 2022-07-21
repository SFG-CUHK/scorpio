
subroutine read3d(this,q_loc,nx_global,ny_global,nz_global,nx_loc,ny_loc,nz_loc,inputbufx,&
inputbufy,inputbufz,bufx,bufy,bufz,flnm,dsetname)
use hdf5
use mpi
use gridModule
implicit none
class(grid)::this
integer::i,startx,endx,starty,endy,startz,endz
integer::nx_global,ny_global,nz_global,nx_loc,ny_loc,nz_loc,flag,bufx,inputbufx,bufy,inputbufy,bufz,inputbufz
double precision::q_loc(1-inputbufx:nx_loc+inputbufx,1-inputbufy:ny_loc+inputbufy,1-inputbufz:nz_loc+inputbufz)
double precision, dimension(:,:,:),allocatable::q_temp
character(len=20):: dsetname
character(len=13):: flnm  ! output filename

integer(HID_T)    ::file_id
integer(HID_T)    ::dset_id
integer(HID_T)    ::filespace
integer(HID_T)    ::mspace_id
integer(HID_T)    ::plist_id
integer(HSIZE_T),dimension(3) :: dimsf
integer(HSIZE_T),dimension(3) :: dimsf_loc
integer(HSIZE_T),dimension(3) :: istart, istride, icount, iblock
integer           ::error,ierr


!! create file
call h5open_f(error)
call h5pcreate_f(H5P_FILE_ACCESS_F,plist_id,error)
call h5pset_fapl_mpio_f(plist_id,MPI_COMM_WORLD,MPI_INFO_NULL,error)
if(myid .eq. 0) then
   write(*,*) "reading dataset ",dsetname, " from the file ", flnm
endif
call h5fopen_f(flnm, H5F_ACC_RDWR_F,file_id, error, access_prp=plist_id)
call h5pclose_f(plist_id, error)

!! prepare date to write

dimsf(1) = nx_global+2*bufx
dimsf(2) = ny_global+2*bufy
dimsf(3) = nz_global+2*bufz

if(ny_global .ne. ny_loc .or. nx_global .ne. nx_loc .or. nz_global .ne. nz_loc) then
   dimsf_loc(1) = nx_loc
   dimsf_loc(2) = ny_loc
   dimsf_loc(3) = nz_loc

   istart(1)=0
   istart(2)=0
   istart(3)=0

   do i=0,this%mpicoord(1)-1
     istart(1)=istart(1)+this%nx_list(i)
   enddo
   istart(1)=istart(1)+bufx

   do i=0,this%mpicoord(2)-1
     istart(2)=istart(2)+this%ny_list(i)
   enddo
   istart(2)=istart(2)+bufy

   do i=0,this%mpicoord(3)-1
     istart(3)=istart(3)+this%nz_list(i)
   enddo
   istart(3)=istart(3)+bufz

   istride(1) = 1
   istride(2) = 1
   istride(3) = 1
   icount(1) = 1
   icount(2) = 1
   icount(3) = 1
   iblock(1) = dimsf_loc(1)
   iblock(2) = dimsf_loc(2)
   iblock(3) = dimsf_loc(3)

   startx=1
   endx=nx_loc
   starty=1
   endy=ny_loc
   startz=1
   endz=nz_loc

   if(this%mpicoord(1) .eq. 0) then
      dimsf_loc(1)=dimsf_loc(1)+bufx
      istart(1) = 0
      iblock(1) = dimsf_loc(1)
      startx=1-bufx
   endif

   if(this%mpicoord(2) .eq. 0) then
      dimsf_loc(2)=dimsf_loc(2)+bufy
      istart(2) = 0
      iblock(2) = dimsf_loc(2)
      starty=1-bufy
   endif

   if(this%mpicoord(3) .eq. 0) then
      dimsf_loc(3)=dimsf_loc(3)+bufz
      istart(3) = 0
      iblock(3) = dimsf_loc(3)
      startz=1-bufz
   endif

   if(this%mpicoord(1) .eq. this%dims_mpi(1)-1) then
      dimsf_loc(1)=dimsf_loc(1)+bufx
      iblock(1) = dimsf_loc(1)
      endx=nx_loc+bufx
   endif

   if(this%mpicoord(2) .eq. this%dims_mpi(2)-1) then
      dimsf_loc(2)=dimsf_loc(2)+bufy
      iblock(2) = dimsf_loc(2)
      endy=ny_loc+bufy
   endif

   if(this%mpicoord(3) .eq. this%dims_mpi(3)-1) then
      dimsf_loc(3)=dimsf_loc(3)+bufz
      iblock(3) = dimsf_loc(3)
      endz=nz_loc+bufz
   endif

   allocate(q_temp(startx:endx,starty:endy,startz:endz))

else
   dimsf_loc(1) = dimsf(1)
   dimsf_loc(2) = dimsf(2)
   dimsf_loc(3) = dimsf(3)
   istart(1)=0
   istart(2)=0
   istart(3)=0
   istride(1)=1
   istride(2)=1
   istride(3)=1
   icount(1)=1
   icount(2)=1
   icount(3)=1
   iblock(1)=dimsf_loc(1)
   iblock(2)=dimsf_loc(2)
   iblock(3)=dimsf_loc(3)
   allocate(q_temp(1-bufx:nx_loc+bufx,1-bufy:ny_loc+bufy,1-bufz:nz_loc+bufz))
endif

!! create file/memory space
call h5screate_simple_f(3,dimsf,filespace,error)
call h5screate_simple_f(3,dimsf_loc,mspace_id,error)
!! create dataset
call h5pcreate_f(H5P_DATASET_XFER_F,plist_id,error)
call h5pset_dxpl_mpio_f(plist_id,H5FD_MPIO_COLLECTIVE_F,error)
call h5sselect_hyperslab_f(filespace,H5S_SELECT_SET_F,istart,icount,error,istride,iblock)

call h5dopen_f(file_id,dsetname,dset_id,error)

call h5dread_f(dset_id,H5T_NATIVE_DOUBLE,q_temp,dimsf_loc,error, file_space_id=filespace,mem_space_id=mspace_id,xfer_prp=plist_id)
call h5dclose_f(dset_id,error)
call h5pclose_f(plist_id,error)

call h5sclose_f(filespace,error)
call h5sclose_f(mspace_id,error)
call h5fclose_f(file_id,error)
call h5close_f(error)

if(ny_global .ne. ny_loc .or. nx_global .ne. nx_loc .or. nz_global .ne. nz_loc) then
  q_loc(startx:endx,starty:endy,startz:endz)=q_temp(startx:endx,starty:endy,startz:endz)
else
  q_loc(1-bufx:nx_loc+bufx,1-bufy:ny_loc+bufy,1-bufz:nz_loc+bufz)=q_temp(1-bufx:nx_loc+bufx,1-bufy:ny_loc+bufy,1-bufz:nz_loc+bufz)
endif

deallocate(q_temp)

end subroutine read3d




subroutine read2d(this,q_loc,nx_global,ny_global,nx_loc,ny_loc,inputbufx,inputbufy,bufx,bufy,flnm,dsetname)
use hdf5
use mpi
use gridModule
implicit none
class(grid)::this
integer::i,startx,endx,starty,endy
integer::nx_global,ny_global,nx_loc,ny_loc,flag,bufx,inputbufx,bufy,inputbufy
double precision::q_loc(1-inputbufx:nx_loc+inputbufx,1-inputbufy:ny_loc+inputbufy)
double precision, dimension(:,:),allocatable::q_temp
character(len=20):: dsetname
character(len=13):: flnm  ! output filename

integer(HID_T)    ::file_id
integer(HID_T)    ::dset_id
integer(HID_T)    ::filespace
integer(HID_T)    ::mspace_id
integer(HID_T)    ::plist_id
integer(HSIZE_T),dimension(2) :: dimsf
integer(HSIZE_T),dimension(2) :: dimsf_loc
integer(HSIZE_T),dimension(2) :: istart, istride, icount, iblock
integer           ::error,ierr


!! create file
call h5open_f(error)
call h5pcreate_f(H5P_FILE_ACCESS_F,plist_id,error)
call h5pset_fapl_mpio_f(plist_id,MPI_COMM_WORLD,MPI_INFO_NULL,error)
if(myid .eq. 0) then
   write(*,*) "reading dataset ",dsetname, " from the file ", flnm
endif
call h5fopen_f(flnm, H5F_ACC_RDWR_F,file_id, error, access_prp=plist_id)
call h5pclose_f(plist_id, error)

!! prepare date to write

dimsf(1) = nx_global+2*bufx
dimsf(2) = ny_global+2*bufy

if(ny_global .ne. ny_loc .or. nx_global .ne. nx_loc) then
   dimsf_loc(1) = nx_loc
   dimsf_loc(2) = ny_loc

   istart(1)=0
   istart(2)=0
   do i=0,this%mpicoord(1)-1
     istart(1)=istart(1)+this%nx_list(i)
   enddo
   istart(1)=istart(1)+bufx

   do i=0,this%mpicoord(2)-1
     istart(2)=istart(2)+this%ny_list(i)
   enddo
   istart(2)=istart(2)+bufy

   istride(1) = 1
   istride(2) = 1
   icount(1) = 1
   icount(2) = 1
   iblock(1) = dimsf_loc(1)
   iblock(2) = dimsf_loc(2)

   startx=1
   endx=nx_loc
   starty=1
   endy=ny_loc

   if(this%mpicoord(1) .eq. 0) then
      dimsf_loc(1)=dimsf_loc(1)+bufx
      istart(1) = 0
      iblock(1) = dimsf_loc(1)
      startx=1-bufx
   endif

   if(this%mpicoord(2) .eq. 0) then
      dimsf_loc(2)=dimsf_loc(2)+bufy
      istart(2) = 0
      iblock(2) = dimsf_loc(2)
      starty=1-bufy
   endif

   if(this%mpicoord(1) .eq. this%dims_mpi(1)-1) then
      dimsf_loc(1)=dimsf_loc(1)+bufx
      iblock(1) = dimsf_loc(1)
      endx=nx_loc+bufx
   endif

   if(this%mpicoord(2) .eq. this%dims_mpi(2)-1) then
      dimsf_loc(2)=dimsf_loc(2)+bufy
      iblock(2) = dimsf_loc(2)
      endy=ny_loc+bufy
   endif

   allocate(q_temp(startx:endx,starty:endy))

else
   dimsf_loc(1) = dimsf(1)
   dimsf_loc(2) = dimsf(2)
   istart(1)=0
   istart(2)=0
   istride(1)=1
   istride(2)=1
   icount(1)=1
   icount(2)=1
   iblock(1)=dimsf_loc(1)
   iblock(2)=dimsf_loc(2)
   allocate(q_temp(1-bufx:nx_loc+bufx,1-bufy:ny_loc+bufy))
endif

!! create file/memory space
call h5screate_simple_f(2,dimsf,filespace,error)
call h5screate_simple_f(2,dimsf_loc,mspace_id,error)
!! create dataset
call h5pcreate_f(H5P_DATASET_XFER_F,plist_id,error)
call h5pset_dxpl_mpio_f(plist_id,H5FD_MPIO_COLLECTIVE_F,error)
call h5sselect_hyperslab_f(filespace,H5S_SELECT_SET_F,istart,icount,error,istride,iblock)

call h5dopen_f(file_id,dsetname,dset_id,error)

call h5dread_f(dset_id,H5T_NATIVE_DOUBLE,q_temp,dimsf_loc,error, file_space_id=filespace,mem_space_id=mspace_id,xfer_prp=plist_id)
call h5dclose_f(dset_id,error)
call h5pclose_f(plist_id,error)

call h5sclose_f(filespace,error)
call h5sclose_f(mspace_id,error)
call h5fclose_f(file_id,error)
call h5close_f(error)

if(ny_global .ne. ny_loc .or. nx_global .ne. nx_loc) then
  q_loc(startx:endx,starty:endy)=q_temp(startx:endx,starty:endy)
else
  q_loc(1-bufx:nx_loc+bufx,1-bufy:ny_loc+bufy)=q_temp(1-bufx:nx_loc+bufx,1-bufy:ny_loc+bufy)
endif

deallocate(q_temp)

end subroutine read2d





subroutine read1d(this,dir,q_loc,nx_global,nx_loc,inputbuf,buf,flnm,dsetname)
use gridModule
use hdf5
use mpi
implicit none
class(grid)::this
integer::dir !!!direction
integer::nx_global,nx_loc,flag,buf,inputbuf
double precision::q_loc(1-inputbuf:nx_loc+inputbuf)
double precision, dimension(:),allocatable::q_temp
character(len=20):: dsetname
character(len=13):: flnm  ! output filename

integer(HID_T)    ::file_id
integer(HID_T)    ::dset_id
integer(HID_T)    ::filespace
integer(HID_T)    ::mspace_id
integer(HID_T)    ::plist_id
integer(HSIZE_T),dimension(1) :: dimsf
integer(HSIZE_T),dimension(1) :: dimsf_loc
integer(HSIZE_T),dimension(1) :: istart, istride, icount, iblock
integer           ::error,ierr
integer           ::i,offset


!! create file
call h5open_f(error)
call h5pcreate_f(H5P_FILE_ACCESS_F,plist_id,error)
call h5pset_fapl_mpio_f(plist_id,MPI_COMM_WORLD,MPI_INFO_NULL,error)
if(myid .eq. 0) then
  write(*,*) "reading dataset ",dsetname, " from the file ", flnm
endif
call h5fopen_f(flnm, H5F_ACC_RDWR_F,file_id, error, access_prp=plist_id)
call h5pclose_f(plist_id, error)

!! prepare date to write

dimsf(1) = nx_global+2*buf

if(nx_global .ne. nx_loc) then
   dimsf_loc(1) = nx_loc
   if(nx_loc .eq. nx_global) then
     istart(1) = 0
   else
     istart(1)=0
     do i=0,this%mpicoord(dir)-1
       select case(dir)
         case(1)
           istart(1) = istart(1)+this%nx_list(i)
         case(2)
           istart(1) = istart(1)+this%ny_list(i)
         case(3)
           istart(1) = istart(1)+this%nz_list(i)
         case default
           print *,"inout.f03: extra dimension??"
           stop
       end select
     enddo
     istart(1)=istart(1)+buf
   endif

   istride(1) = 1
   icount(1) = 1
   iblock(1) = dimsf_loc(1)

   if(myid .eq. 0) then
      dimsf_loc(1)=nx_loc+buf
      istart(1) = 0
      iblock = dimsf_loc(1)
      allocate(q_temp(1-buf:nx_loc))
   endif

   if(myid .eq. nprocs-1) then
      dimsf_loc(1)=nx_loc+buf
      iblock(1) = dimsf_loc(1)
      allocate(q_temp(1:nx_loc+buf))
   endif

   if(myid .ne. 0 .and. myid .ne. nprocs-1) then
      allocate(q_temp(1:nx_loc))
   endif

else
   dimsf_loc(1) = dimsf(1)
   istart(1)=0
   istride(1)=1
   icount(1)=1
   iblock=dimsf_loc(1)
   allocate(q_temp(1-buf:nx_loc+buf))
endif

!! create file/memory space
call h5screate_simple_f(1,dimsf,filespace,error)
call h5screate_simple_f(1,dimsf_loc,mspace_id,error)
!! create dataset
call h5pcreate_f(H5P_DATASET_XFER_F,plist_id,error)
call h5pset_dxpl_mpio_f(plist_id,H5FD_MPIO_COLLECTIVE_F,error)
call h5sselect_hyperslab_f(filespace,H5S_SELECT_SET_F,istart,icount,error,istride,iblock)

call h5dopen_f(file_id,dsetname,dset_id,error)


call h5dread_f(dset_id,H5T_NATIVE_DOUBLE,q_temp,dimsf_loc,error, file_space_id=filespace,mem_space_id=mspace_id,xfer_prp=plist_id)
call h5dclose_f(dset_id,error)
call h5pclose_f(plist_id,error)

call h5sclose_f(filespace,error)
call h5sclose_f(mspace_id,error)
call h5fclose_f(file_id,error)
call h5close_f(error)

if(nx_global .ne. nx_loc) then
  if(myid .eq. 0)then
    q_loc(1-buf:nx_loc)=q_temp(1-buf:nx_loc)
  endif
  if(myid .eq. nprocs-1) then
    q_loc(1:nx_loc+buf)=q_temp(1:nx_loc+buf)
  endif
  if(myid .ne. 0 .and. myid .ne. nprocs-1) then
    q_loc(1:nx_loc)=q_temp(1:nx_loc)
  endif
else
  q_loc=q_temp(1-buf:nx_loc+buf)
endif

deallocate(q_temp)

end subroutine read1d








subroutine output1d(this,dir,q_loc,nx_global,nx_loc,inputbuf,buf,flnm,dsetname,flag)
use gridModule
use hdf5
use mpi
implicit none
class(grid)::this
integer::dir !!!direction
integer::nx_global,nx_loc,flag,buf,inputbuf
double precision::q_loc(1-inputbuf:nx_loc+inputbuf)
double precision, dimension(:),allocatable::q_temp
character(len=20):: dsetname
character(len=13):: flnm  ! output filename

integer(HID_T)    ::file_id
integer(HID_T)    ::dset_id
integer(HID_T)    ::filespace
integer(HID_T)    ::mspace_id
integer(HID_T)    ::plist_id
integer(HSIZE_T),dimension(1) :: dimsf
integer(HSIZE_T),dimension(1) :: dimsf_loc
integer(HSIZE_T),dimension(1) :: istart, istride, icount, iblock
integer           ::error,ierr
integer           ::i,offset


!! create file
call h5open_f(error)
call h5pcreate_f(H5P_FILE_ACCESS_F,plist_id,error)
call h5pset_fapl_mpio_f(plist_id,MPI_COMM_WORLD,MPI_INFO_NULL,error)
  if (flag .eq. 1) then
    if(myid .eq. 0) then
      write(*,*) "inout.f03: creating new file:",flnm
      !write(*,*) "writing dataset ",dsetname, " into the file ", flnm
    endif
    call h5fcreate_f(flnm, H5F_ACC_TRUNC_F, file_id, error, access_prp=plist_id)
  else
    if(myid .eq. 0) then
      !write(*,*) "writing dataset ",dsetname, " into the file ", flnm
    endif
    call h5fopen_f(flnm, H5F_ACC_RDWR_F,file_id, error, access_prp=plist_id)
  endif
call h5pclose_f(plist_id, error)

!! prepare date to write

dimsf(1) = nx_global+2*buf

if(nx_global .ne. nx_loc) then
   dimsf_loc(1) = nx_loc
   if(nx_loc .eq. nx_global) then
     istart(1) = 0
   else
     istart(1)=0
     do i=0,this%mpicoord(dir)-1
       select case(dir)
         case(1)
           istart(1) = istart(1)+this%nx_list(i)
         case(2)
           istart(1) = istart(1)+this%ny_list(i)
         case(3)
           istart(1) = istart(1)+this%nz_list(i)
         case default
           print *,"inout.f03: extra dimension??"
           stop
       end select
     enddo
     istart(1)=istart(1)+buf
   endif

   istride(1) = 1
   icount(1) = 1
   iblock(1) = dimsf_loc(1)

   if(myid .eq. 0) then
      dimsf_loc(1)=nx_loc+buf
      istart(1) = 0
      iblock = dimsf_loc(1)
      allocate(q_temp(1-buf:nx_loc))
      q_temp(1-buf:nx_loc)=q_loc(1-buf:nx_loc)
   endif

   if(myid .eq. nprocs-1) then
      dimsf_loc(1)=nx_loc+buf
      iblock(1) = dimsf_loc(1)
      allocate(q_temp(1:nx_loc+buf))
      q_temp(1:nx_loc+buf)=q_loc(1:nx_loc+buf)
   endif

   if(myid .ne. 0 .and. myid .ne. nprocs-1) then
      allocate(q_temp(1:nx_loc))
      q_temp(1:nx_loc)=q_loc(1:nx_loc)
   endif

else
   dimsf_loc(1) = dimsf(1)
   istart(1)=0
   istride(1)=1
   icount(1)=1
   iblock=dimsf_loc(1)
   allocate(q_temp(1-buf:nx_loc+buf))
   q_temp=q_loc(1-buf:nx_loc+buf)
endif

!! create file/memory space
call h5screate_simple_f(1,dimsf,filespace,error)
call h5screate_simple_f(1,dimsf_loc,mspace_id,error)
!! create dataset
call h5pcreate_f(H5P_DATASET_XFER_F,plist_id,error)
call h5pset_dxpl_mpio_f(plist_id,H5FD_MPIO_COLLECTIVE_F,error)
call h5sselect_hyperslab_f(filespace,H5S_SELECT_SET_F,istart,icount,error,istride,iblock)

call h5dcreate_f(file_id,dsetname,H5T_NATIVE_DOUBLE,filespace, dset_id,error)


call h5dwrite_f(dset_id,H5T_NATIVE_DOUBLE,q_temp,dimsf_loc,error, file_space_id=filespace,mem_space_id=mspace_id,xfer_prp=plist_id)

call h5dclose_f(dset_id,error)
call h5pclose_f(plist_id,error)

call h5sclose_f(filespace,error)
call h5sclose_f(mspace_id,error)
call h5fclose_f(file_id,error)
call h5close_f(error)

deallocate(q_temp)

end subroutine output1d





!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine output2d(this,q_loc,nx_global,ny_global,nx_loc,ny_loc,inputbufx,inputbufy,bufx,bufy,flnm,dsetname,flag)
use hdf5
use mpi
use gridModule
implicit none
class(grid)::this
integer::i,startx,endx,starty,endy
integer::nx_global,ny_global,nx_loc,ny_loc,flag,bufx,inputbufx,bufy,inputbufy
double precision::q_loc(1-inputbufx:nx_loc+inputbufx,1-inputbufy:ny_loc+inputbufy)
double precision, dimension(:,:),allocatable::q_temp
character(len=20):: dsetname
character(len=13):: flnm  ! output filename

integer(HID_T)    ::file_id
integer(HID_T)    ::dset_id
integer(HID_T)    ::filespace
integer(HID_T)    ::mspace_id
integer(HID_T)    ::plist_id
integer(HSIZE_T),dimension(2) :: dimsf
integer(HSIZE_T),dimension(2) :: dimsf_loc
integer(HSIZE_T),dimension(2) :: istart, istride, icount, iblock
integer           ::error,ierr


!! create file
call h5open_f(error)
call h5pcreate_f(H5P_FILE_ACCESS_F,plist_id,error)
call h5pset_fapl_mpio_f(plist_id,MPI_COMM_WORLD,MPI_INFO_NULL,error)
  if (flag .eq. 1) then
    if(myid .eq. 0) then
      write(*,*) "inout.f03: creating new file:",flnm
      !write(*,*) "writing dataset ",dsetname, " into the file ", flnm
    endif
    call h5fcreate_f(flnm, H5F_ACC_TRUNC_F, file_id, error, access_prp=plist_id)
  else
    if(myid .eq. 0) then
      !write(*,*) "writing dataset ",dsetname, " into the file ", flnm
    endif
    call h5fopen_f(flnm, H5F_ACC_RDWR_F,file_id, error, access_prp=plist_id)
  endif
call h5pclose_f(plist_id, error)

!! prepare date to write

dimsf(1) = nx_global+2*bufx
dimsf(2) = ny_global+2*bufy

if(ny_global .ne. ny_loc .or. nx_global .ne. nx_loc) then
   dimsf_loc(1) = nx_loc
   dimsf_loc(2) = ny_loc

   istart(1)=0
   istart(2)=0
   do i=0,this%mpicoord(1)-1
     istart(1)=istart(1)+this%nx_list(i)
   enddo
   istart(1)=istart(1)+bufx

   do i=0,this%mpicoord(2)-1
     istart(2)=istart(2)+this%ny_list(i)
   enddo
   istart(2)=istart(2)+bufy

   istride(1) = 1
   istride(2) = 1
   icount(1) = 1
   icount(2) = 1
   iblock(1) = dimsf_loc(1)
   iblock(2) = dimsf_loc(2)

   startx=1
   endx=nx_loc
   starty=1
   endy=ny_loc

   if(this%mpicoord(1) .eq. 0) then
      dimsf_loc(1)=dimsf_loc(1)+bufx
      istart(1) = 0
      iblock(1) = dimsf_loc(1)
      startx=1-bufx
   endif

   if(this%mpicoord(2) .eq. 0) then
      dimsf_loc(2)=dimsf_loc(2)+bufy
      istart(2) = 0
      iblock(2) = dimsf_loc(2)
      starty=1-bufy
   endif

   if(this%mpicoord(1) .eq. this%dims_mpi(1)-1) then
      dimsf_loc(1)=dimsf_loc(1)+bufx
      iblock(1) = dimsf_loc(1)
      endx=nx_loc+bufx
   endif

   if(this%mpicoord(2) .eq. this%dims_mpi(2)-1) then
      dimsf_loc(2)=dimsf_loc(2)+bufy
      iblock(2) = dimsf_loc(2)
      endy=ny_loc+bufy
   endif

   allocate(q_temp(startx:endx,starty:endy))
   q_temp(startx:endx,starty:endy)=q_loc(startx:endx,starty:endy)

else
   dimsf_loc(1) = dimsf(1)
   dimsf_loc(2) = dimsf(2)
   istart(1)=0
   istart(2)=0
   istride(1)=1
   istride(2)=1
   icount(1)=1
   icount(2)=1
   iblock(1)=dimsf_loc(1)
   iblock(2)=dimsf_loc(2)
   allocate(q_temp(1-bufx:nx_loc+bufx,1-bufy:ny_loc+bufy))
   q_temp(1-bufx:nx_loc+bufx,1-bufy:ny_loc+bufy)=q_loc(1-bufx:nx_loc+bufx,1-bufy:ny_loc+bufy)
endif

!! create file/memory space
call h5screate_simple_f(2,dimsf,filespace,error)
call h5screate_simple_f(2,dimsf_loc,mspace_id,error)
!! create dataset
call h5pcreate_f(H5P_DATASET_XFER_F,plist_id,error)
call h5pset_dxpl_mpio_f(plist_id,H5FD_MPIO_COLLECTIVE_F,error)
call h5sselect_hyperslab_f(filespace,H5S_SELECT_SET_F,istart,icount,error,istride,iblock)

call h5dcreate_f(file_id,dsetname,H5T_NATIVE_DOUBLE,filespace, dset_id,error)


call h5dwrite_f(dset_id,H5T_NATIVE_DOUBLE,q_temp,dimsf_loc,error, file_space_id=filespace,mem_space_id=mspace_id,xfer_prp=plist_id)
call h5dclose_f(dset_id,error)
call h5pclose_f(plist_id,error)

call h5sclose_f(filespace,error)
call h5sclose_f(mspace_id,error)
call h5fclose_f(file_id,error)
call h5close_f(error)

deallocate(q_temp)

end subroutine output2d




!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine output3d(this,q_loc,nx_global,ny_global,nz_global,nx_loc,ny_loc,nz_loc,inputbufx,inputbufy,inputbufz,&
& bufx,bufy,bufz,flnm,dsetname,flag)
use hdf5
use mpi
use gridModule
implicit none
class(grid)::this
integer::i,startx,endx,starty,endy,startz,endz
integer::nx_global,ny_global,nx_loc,ny_loc,flag,bufx,inputbufx,bufy,inputbufy
integer::nz_global,nz_loc,bufz,inputbufz
double precision::q_loc(1-inputbufx:nx_loc+inputbufx,1-inputbufy:ny_loc+inputbufy,1-inputbufz:nz_loc+inputbufz)
double precision, dimension(:,:,:),allocatable::q_temp
character(len=20):: dsetname
character(len=13):: flnm  ! output filename

integer(HID_T)    ::file_id
integer(HID_T)    ::dset_id
integer(HID_T)    ::filespace
integer(HID_T)    ::mspace_id
integer(HID_T)    ::plist_id
integer(HSIZE_T),dimension(3) :: dimsf
integer(HSIZE_T),dimension(3) :: dimsf_loc
integer(HSIZE_T),dimension(3) :: istart, istride, icount, iblock
integer           ::error,ierr


!! create file
call h5open_f(error)
call h5pcreate_f(H5P_FILE_ACCESS_F,plist_id,error)
call h5pset_fapl_mpio_f(plist_id,MPI_COMM_WORLD,MPI_INFO_NULL,error)
  if (flag .eq. 1) then
    if(myid .eq. 0) then
      write(*,*) "inout.f03: creating new file:",flnm
      !write(*,*) "writing dataset ",dsetname, " into the file ", flnm
    endif
    call h5fcreate_f(flnm, H5F_ACC_TRUNC_F, file_id, error, access_prp=plist_id)
  else
    if(myid .eq. 0) then
      !write(*,*) "writing dataset ",dsetname, " into the file ", flnm
    endif
    call h5fopen_f(flnm, H5F_ACC_RDWR_F,file_id, error, access_prp=plist_id)
  endif
call h5pclose_f(plist_id, error)

!! prepare date to write

dimsf(1) = nx_global+2*bufx
dimsf(2) = ny_global+2*bufy
dimsf(3) = nz_global+2*bufz

if(ny_global .ne. ny_loc .or. nx_global .ne. nx_loc .or. nz_global .ne. nz_loc) then
   dimsf_loc(1) = nx_loc
   dimsf_loc(2) = ny_loc
   dimsf_loc(3) = nz_loc

   istart(1)=0
   istart(2)=0
   istart(3)=0
   do i=0,this%mpicoord(1)-1
     istart(1)=istart(1)+this%nx_list(i)
   enddo
   istart(1)=istart(1)+bufx

   do i=0,this%mpicoord(2)-1
     istart(2)=istart(2)+this%ny_list(i)
   enddo
   istart(2)=istart(2)+bufy

   do i=0,this%mpicoord(3)-1
     istart(3)=istart(3)+this%nz_list(i)
   enddo
   istart(3)=istart(3)+bufz

   istride(1) = 1
   istride(2) = 1
   istride(3) = 1
   icount(1) = 1
   icount(2) = 1
   icount(3) = 1
   iblock(1) = dimsf_loc(1)
   iblock(2) = dimsf_loc(2)
   iblock(3) = dimsf_loc(3)

   startx=1
   endx=nx_loc
   starty=1
   endy=ny_loc
   startz=1
   endz=nz_loc

   if(this%mpicoord(1) .eq. 0) then
      dimsf_loc(1)=dimsf_loc(1)+bufx
      istart(1) = 0
      iblock(1) = dimsf_loc(1)
      startx=1-bufx
   endif

   if(this%mpicoord(2) .eq. 0) then
      dimsf_loc(2)=dimsf_loc(2)+bufy
      istart(2) = 0
      iblock(2) = dimsf_loc(2)
      starty=1-bufy
   endif

   if(this%mpicoord(3) .eq. 0) then
      dimsf_loc(3)=dimsf_loc(3)+bufz
      istart(3) = 0
      iblock(3) = dimsf_loc(3)
      startz=1-bufz
   endif

   if(this%mpicoord(1) .eq. this%dims_mpi(1)-1) then
      dimsf_loc(1)=dimsf_loc(1)+bufx
      iblock(1) = dimsf_loc(1)
      endx=nx_loc+bufx
   endif

   if(this%mpicoord(2) .eq. this%dims_mpi(2)-1) then
      dimsf_loc(2)=dimsf_loc(2)+bufy
      iblock(2) = dimsf_loc(2)
      endy=ny_loc+bufy
   endif

   if(this%mpicoord(3) .eq. this%dims_mpi(3)-1) then
      dimsf_loc(3)=dimsf_loc(3)+bufz
      iblock(3) = dimsf_loc(3)
      endz=nz_loc+bufz
   endif

   allocate(q_temp(startx:endx,starty:endy,startz:endz))
   q_temp(startx:endx,starty:endy,startz:endz)=q_loc(startx:endx,starty:endy,startz:endz)

else
   dimsf_loc(1) = dimsf(1)
   dimsf_loc(2) = dimsf(2)
   dimsf_loc(3) = dimsf(3)
   istart(1)=0
   istart(2)=0
   istart(3)=0
   istride(1)=1
   istride(2)=1
   istride(3)=1
   icount(1)=1
   icount(2)=1
   icount(3)=1
   iblock(1)=dimsf_loc(1)
   iblock(2)=dimsf_loc(2)
   iblock(3)=dimsf_loc(3)
   allocate(q_temp(1-bufx:nx_loc+bufx,1-bufy:ny_loc+bufy,1-bufz:nz_loc+bufz))
   q_temp(1-bufx:nx_loc+bufx,1-bufy:ny_loc+bufy,1-bufz:nz_loc+bufz)=q_loc(1-bufx:nx_loc+bufx,1-bufy:ny_loc+bufy,1-bufz:nz_loc+bufz)
endif

!! create file/memory space
call h5screate_simple_f(3,dimsf,filespace,error)
call h5screate_simple_f(3,dimsf_loc,mspace_id,error)
!! create dataset
call h5pcreate_f(H5P_DATASET_XFER_F,plist_id,error)
call h5pset_dxpl_mpio_f(plist_id,H5FD_MPIO_COLLECTIVE_F,error)
call h5sselect_hyperslab_f(filespace,H5S_SELECT_SET_F,istart,icount,error,istride,iblock)

call h5dcreate_f(file_id,dsetname,H5T_NATIVE_DOUBLE,filespace, dset_id,error)


call h5dwrite_f(dset_id,H5T_NATIVE_DOUBLE,q_temp,dimsf_loc,error, file_space_id=filespace,mem_space_id=mspace_id,xfer_prp=plist_id)
call h5dclose_f(dset_id,error)
call h5pclose_f(plist_id,error)

call h5sclose_f(filespace,error)
call h5sclose_f(mspace_id,error)
call h5fclose_f(file_id,error)
call h5close_f(error)

deallocate(q_temp)

end subroutine output3d




!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!! coded by zwg to output the results in vtk
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine output2d_MHD_vtk(this,q_tmp)
use gridModule
implicit none
class(grid)::this

double precision,dimension(1-this%nbuf:this%nMesh(1)+this%nbuf,1-this%nbuf:this%nMesh(2)+this%nbuf &
                         &,this%nvar)::q_tmp 

double precision,dimension(1-this%nbuf:this%nMesh(1)+this%nbuf)::xc
double precision,dimension(1-this%nbuf:this%nMesh(2)+this%nbuf)::yc

integer:: filenum,Dirnum,ranknum,gridnum
!integer:: ndim
integer:: i,j,invar
integer:: ib,ie,jb,je,imax,jmax,nx,ny
integer::ierr
integer:: variable(8)

character (len=20):: s_node_num, s_cells, s_imax, s_jmax   !!!for vtk output

character(len=50):: dirname,dirname1
character(len=50)::filename
character(len=50)::rankname
character(len=50)::gridname

logical:: fileDirhave

variable=this%variable
xc=this%xc(1)%coords
yc=this%xc(2)%coords


filenum=this%fnum
gridnum=this%gridID
ranknum=myid

!ndim=this%ndim



write(dirname,'(I6.6)') filenum
write(gridname,'(I3.3)') gridnum

dirname1="g_"//trim(gridname)//"_"//dirname
inquire(file=trim(dirname1),exist=fileDirhave)

if(.not.fileDirhave)then

   if(myid==0)then
    call system("mkdir "//trim(dirname1))
    !write(*,*)"Creat output file direction done !!!"
   endif 
endif
   call MPI_BARRIER(MPI_COMM_WORLD,ierr)

! write(*,*)"here!!",trim(dirname)

write(rankname,'(I6.6)')ranknum
filename="P"//trim(rankname)

open(30,form='formatted',file=trim(dirname1)//'/'//trim(filename)//'.vtk')



        nx=this%nMesh(1)
        ny=this%nMesh(2)

        imax=nx
        jmax=ny

       call i4_to_s_left ( imax*jmax, s_node_num )
       call i4_to_s_left ( (imax-1)*(jmax-1), s_cells )     ! don't use kmax for 2D
       call i4_to_s_left ( imax, s_imax )
       call i4_to_s_left ( jmax, s_jmax ) 


         write ( 30, '(a)' ) '# vtk DataFile Version 2.0 11 22 33 66'
         write ( 30, '(a)' ) 'Scorpio2D'
         write ( 30, '(a)' ) 'ASCII'
         write ( 30, '(a)' ) 'DATASET STRUCTURED_GRID'
         write ( 30, '(a)' ) 'DIMENSIONS ' // (s_imax) // (s_jmax) // '1'  ! kmax = 1 for 2D problems
         write ( 30, '(a)' ) 'POINTS ' // (s_node_num) // 'double'




    ib=1
    ie=imax
    jb=1
    je=jmax

    do j = jb, je
      do i = ib, ie
          write (30, * ) xc(i),yc(j), 0.0
          !write (30, * ) this%xl(1)%coords(i), this%xl(2)%coords(j), 0.0
      end do
    end do


    if(variable(1) .eq. 1) then
       write ( 30, '(a)' ) 'CELL_DATA ' // (s_cells)
       write ( 30, '(a)' ) 'POINT_DATA ' // (s_node_num)

       write ( 30, '(a)' ) 'SCALARS density double '
       write ( 30, '(a)' ) 'LOOKUP_TABLE default '
       write ( 30, * ) ((q_tmp(i,j,1), i=ib,je),j=jb,je)
    endif

    if(variable(2) .eq. 1) then
       write ( 30, '(a)' ) 'CELL_DATA ' // (s_cells)
       write ( 30, '(a)' ) 'POINT_DATA ' // (s_node_num)

       write ( 30, '(a)' ) 'SCALARS vx double '
       write ( 30, '(a)' ) 'LOOKUP_TABLE default '
       write ( 30, * ) ((q_tmp(i,j,2)/q_tmp(i,j,1), i=ib,je),j=jb,je)
    endif

    if(variable(3) .eq. 1) then
       write ( 30, '(a)' ) 'CELL_DATA ' // (s_cells)
       write ( 30, '(a)' ) 'POINT_DATA ' // (s_node_num)

       write ( 30, '(a)' ) 'SCALARS vy double '
       write ( 30, '(a)' ) 'LOOKUP_TABLE default '
       write ( 30, * ) ((q_tmp(i,j,3)/q_tmp(i,j,1), i=ib,je),j=jb,je)
    endif


    if(variable(4) .eq. 1) then
       write ( 30, '(a)' ) 'CELL_DATA ' // (s_cells)
       write ( 30, '(a)' ) 'POINT_DATA ' // (s_node_num)

       write ( 30, '(a)' ) 'SCALARS vz double '
       write ( 30, '(a)' ) 'LOOKUP_TABLE default '
       write ( 30, * ) ((q_tmp(i,j,4), i=ib,je),j=jb,je)
    endif

    if(variable(4) .eq. 1) then
       write ( 30, '(a)' ) 'CELL_DATA ' // (s_cells)
       write ( 30, '(a)' ) 'POINT_DATA ' // (s_node_num)

       write ( 30, '(a)' ) 'SCALARS Bx double '
       write ( 30, '(a)' ) 'LOOKUP_TABLE default '
       write ( 30, * ) ((q_tmp(i,j,5), i=ib,je),j=jb,je)
    endif

    if(variable(4) .eq. 1) then
       write ( 30, '(a)' ) 'CELL_DATA ' // (s_cells)
       write ( 30, '(a)' ) 'POINT_DATA ' // (s_node_num)

       write ( 30, '(a)' ) 'SCALARS By double '
       write ( 30, '(a)' ) 'LOOKUP_TABLE default '
       write ( 30, * ) ((q_tmp(i,j,6), i=ib,je),j=jb,je)
    endif

    if(variable(4) .eq. 1) then
       write ( 30, '(a)' ) 'CELL_DATA ' // (s_cells)
       write ( 30, '(a)' ) 'POINT_DATA ' // (s_node_num)

       write ( 30, '(a)' ) 'SCALARS Bz double '
       write ( 30, '(a)' ) 'LOOKUP_TABLE default '
       write ( 30, * ) ((q_tmp(i,j,7), i=ib,je),j=jb,je)
    endif


    !!!!!!!!!!!!!!!!!!!!!!!!!!zwg!!!!! add other vaiables here!!!

close(30)

end subroutine output2d_MHD_vtk


subroutine output2d_MHD_pvts(this,q)
use gridModule
implicit none
class(grid)::this
double precision,dimension(1-this%nbuf:this%nMesh(1)+this%nbuf,1-this%nbuf:this%nMesh(2)+this%nbuf &
                         &,this%nvar)::q 

double precision,dimension(1-this%nbuf:this%nMesh(1)+this%nbuf)::xl
double precision,dimension(1-this%nbuf:this%nMesh(2)+this%nbuf)::yl
integer,dimension(1:this%ndim)::dims  

integer:: ndim

integer:: gnis,gnie,gnjs,gnje,gnks,gnke
character(len=50)::gnis_str,gnie_str,gnjs_str,gnje_str ,gnks_str,gnke_str


integer:: nis,nie,njs,nje,nks,nke
character(len=50)::nis_str,nie_str,njs_str,nje_str,nks_str,nke_str

integer:: nis1_send,nie1_send,njs1_send,nje1_send,nks1_send,nke1_send
integer:: nis1_rec, nie1_rec, njs1_rec, nje1_rec, nks1_rec,nke1_rec


integer,dimension(:),allocatable:: nis1,nie1,njs1,nje1,nks1,nke1
character(len=50)::nis1_str,nie1_str,njs1_str,nje1_str,nks1_str,nke1_str



integer:: i,j  !k
integer::Ni,Nj !Nk

integer::iproc,Totalprocs

character(len=50)::iproc_str

integer::ierr
integer:: variable(8)

logical:: fileDirhave

integer:: fileid,rankid,gridid
character(len=50)::fileid_str,rankid_str,gridid_str

character(len=50):: dirname,filename_pvts,filename_vts


xl=this%xl(1)%coords
yl=this%xl(2)%coords

variable=this%variable
ndim=this%ndim
dims=this%dims_mpi

Totalprocs=nprocs

!write(*,*)"ndim=",ndim, "dims=",dims, "TotalProcess",nprocs,"nBlockx=",nBlockx,"nBlocky=",nBlocky
!write(*,*)"this%nMesh(1)=",this%nMesh(1),"this%nMesh(2)",this%nMesh(2)
!write(*,*)"myid= ",myid, "this%gnx_list_R=", this%gnx_list_R, "this%gny_list_R=", this%gny_list_R

nie=this%gnx_list_R
nje=this%gny_list_R
nke=0

nis=this%gnx_list_R-this%nMesh(1)
njs=this%gny_list_R-this%nMesh(2)
nks=0

call i4_to_s_left (nis, nis_str )
call i4_to_s_left (njs, njs_str )
call i4_to_s_left (nks, nks_str )

call i4_to_s_left (nie, nie_str )
call i4_to_s_left (nje, nje_str )
call i4_to_s_left (nke, nke_str )


allocate(nis1(0:Totalprocs-1),nie1(0:Totalprocs-1))
allocate(njs1(0:Totalprocs-1),nje1(0:Totalprocs-1))
allocate(nks1(0:Totalprocs-1),nke1(0:Totalprocs-1))


nis1=0
nie1=0

njs1=0
nje1=0

nks1=0
nke1=0




!integer:: nis1_send,nie1_send,njs1_send,nje1_send,nks1_send,nke1_send
!integer:: nis1_rec, nie1_rec, njs1_rec, nje1_rec, nks1_rec,nke1_rec

!do iproc=0, Totalprocs

!if(myid==0 .and. iproc==0) then
   !nis1(iproc)=nis
   !nie1(iproc)=nie

   !njs1(iproc)=njs
   !nje1(iproc)=nje

   !nks1(iproc)=nks
   !nke1(iproc)=nke
!endif

!if(myid>0 .and. myid==iproc) then
   !nis1_send=nis
   !nie1_send=nie

   !njs1_send=njs
   !nje1_send=nje

   !nks1_send=nks
   !nke1_send=nke

!endif

!enddo


!call MPI_GATHER(nis,1,MPI_INT,nis1,1,MPI_INT,0,MPI_COMM_WORLD,ierr)
!call MPI_GATHER(nie,1,MPI_INT,nie1,1,MPI_INT,0,MPI_COMM_WORLD,ierr)

!call MPI_GATHER(njs,1,MPI_INT,njs1,1,MPI_INT,0,MPI_COMM_WORLD,ierr)
!call MPI_GATHER(nje,1,MPI_INT,nje1,1,MPI_INT,0,MPI_COMM_WORLD,ierr)

call MPI_ALLGATHER(nis,1,MPI_INT,nis1,1,MPI_INT,MPI_COMM_WORLD,ierr)
call MPI_ALLGATHER(nie,1,MPI_INT,nie1,1,MPI_INT,MPI_COMM_WORLD,ierr)

call MPI_ALLGATHER(njs,1,MPI_INT,njs1,1,MPI_INT,MPI_COMM_WORLD,ierr)
call MPI_ALLGATHER(nje,1,MPI_INT,nje1,1,MPI_INT,MPI_COMM_WORLD,ierr)

!if(myid==0)then
!write(*,*)nis1
!write(*,*)nie1
!write(*,*)njs1
!write(*,*)nje1
!endif
!stop



!MPI_GATHER(SENDBUF, SENDCOUNT, SENDTYPE, RECVBUF, RECVCOUNT, RECVTYPE, ROOT, COMM, IERROR)
        !<type> SENDBUF(*), RECVBUF(*)
        !INTEGER SENDCOUNT, SENDTYPE, RECVCOUNT, RECVTYPE, ROOT, COMM, IERROR
!call MPI_BARRIER(MPI_COMM_WORLD,ierr)


gnis=0
gnjs=0
gnks=0

gnie=this%nMesh_global(1)
gnje=this%nMesh_global(2)
gnke=0


call i4_to_s_left (gnis, gnis_str )
call i4_to_s_left (gnjs, gnjs_str )
call i4_to_s_left (gnks, gnks_str )

call i4_to_s_left (gnie, gnie_str )
call i4_to_s_left (gnje, gnje_str )
call i4_to_s_left (gnke, gnke_str )


Ni=this%nMesh(1)
Nj=this%nMesh(2)


!write(*,*)"myid= ",myid,"nis=",nis,"njs=",njs,"nie=",nie,"nje=",nje
!write(*,*)"myid= ",myid,"gnie=",gnie,"gnje=",gnje




gridid=this%gridID
fileid=this%fnum
rankid=myid

call i4_to_s_left (gridid, gridid_str )
call i4_to_s_left (fileid, fileid_str )
call i4_to_s_left (rankid, rankid_str )


dirname="MHD_g"//trim(gridid_str)//"_"//trim(fileid_str)

inquire(file=trim(dirname),exist=fileDirhave)

if(.not.fileDirhave)then
   !if(myid==0)then
    call system("mkdir "//trim(dirname))
   !endif  !!! end of if(myid==0)then
endif  ! end of if(.not.fileDirhave)then

!call MPI_BARRIER(MPI_COMM_WORLD,ierr)

!filename_pvts=trim(dirname)//'/'//trim(dirname)//'.pvts'
filename_pvts=trim(dirname)//'/'//trim(dirname)//'_'//trim(rankid_str)//'.pvts'


!if(myid==0)then
    open(10,file=trim(filename_pvts))

          write(10,'(A)')'<?xml version="1.0"?>'
          write(10,'(A)')'<VTKFile type="PStructuredGrid" version="0.1" byte_order="LittleEndian">'
          write(10,'(A)')'  <PStructuredGrid WholeExtent="'//trim(gnis_str)//" "//trim(gnie_str)// &
          & " "//trim(gnjs_str)//" "//trim(gnje_str)//" "//trim(gnks_str)//" "//trim(gnke_str)//'" GhostLevel="#">'

          write(10,'(A)')'    <PCellData>'

          if(variable(1) .eq. 1) then
            write(10,'(A)')'      <PDataArray type="Float32" Name="Density"/>'
          endif

          if(variable(2) .eq. 1) then
            write(10,'(A)')'      <PDataArray type="Float32" Name="Vx"/>'
          endif

          if(variable(3) .eq. 1) then
            write(10,'(A)')'      <PDataArray type="Float32" Name="Vy"/>'
          endif

          if(variable(4) .eq. 1) then
            write(10,'(A)')'      <PDataArray type="Float32" Name="Vz"/>'
          endif

          if(variable(5) .eq. 1) then
            write(10,'(A)')'      <PDataArray type="Float32" Name="Bxl"/>'
          endif

          if(variable(6) .eq. 1) then
            write(10,'(A)')'      <PDataArray type="Float32" Name="Byl"/>'
          endif

          if(variable(7) .eq. 1) then
            write(10,'(A)')'      <PDataArray type="Float32" Name="Bzl"/>'
          endif

          if(variable(8) .eq. 1) then
            !write(10,'(A)')'      <PDataArray type="Float32" Name="Pressure"/>'
            write(10,'(A)')'      <PDataArray type="Float32" Name="Energy"/>'
          endif

          if(variable(5) .eq. 1) then
            write(10,'(A)')'      <PDataArray type="Float32" Name="Bxr"/>'
          endif

          if(variable(6) .eq. 1) then
            write(10,'(A)')'      <PDataArray type="Float32" Name="Byr"/>'
          endif

          if(variable(7) .eq. 1) then
            write(10,'(A)')'      <PDataArray type="Float32" Name="Bzr"/>'
          endif


          write(10,'(A)')'    </PCellData>'
          


          write(10,'(A)')'    <PPoints>'
          write(10,'(A)')'      <PDataArray type="Float32" NumberOfComponents="3"/>'
          write(10,'(A)')'    </PPoints>'


          do iproc=0, Totalprocs-1
               !if(myid==0) then
                 !write(*,*)"_____", iproc
               !endif

               call i4_to_s_left (iproc, iproc_str)

               call i4_to_s_left (nis1(iproc), nis1_str)
               call i4_to_s_left (nie1(iproc), nie1_str)

               call i4_to_s_left (njs1(iproc), njs1_str)
               call i4_to_s_left (nje1(iproc), nje1_str)

               call i4_to_s_left (nks1(iproc), nks1_str)
               call i4_to_s_left (nke1(iproc), nke1_str)

               !call MPI_BARRIER(MPI_COMM_WORLD,ierr)
               write(10,'(A)')'    <Piece Extent="'//trim(nis1_str)//" "//trim(nie1_str)// &
               & " "//trim(njs1_str)//" "//trim(nje1_str)//" "//trim(nks1_str)//" "//trim(nke1_str)// &
               &'"'//' Source="'//trim(iproc_str)//'.vts"/>'      
          enddo

          write(10,'(A)')'  </PStructuredGrid>'
          write(10,'(A)')'</VTKFile>'

    close(10)
!endif  !!! end of if(myid==0)then


   !call MPI_BARRIER(MPI_COMM_WORLD,ierr)

filename_vts=trim(dirname)//'/'//trim(rankid_str)//'.vts'

open(11,file=filename_vts)

         write(11,'(A)')'<?xml version="1.0"?>'
         write(11,'(A)')'<VTKFile type="StructuredGrid" version="0.1" byte_order="LittleEndian">'
         write(11,'(A)')'  <StructuredGrid WholeExtent="'//trim(nis_str)//" "//trim(nie_str)//" " &
                           & //trim(njs_str)//" "//trim(nje_str)//" "//trim(nks_str)//" "//trim(nke_str)//'">'
         write(11,'(A)')'    <Piece Extent="'//trim(nis_str)//" "//trim(nie_str)//" " &
                           & //trim(njs_str)//" "//trim(nje_str)//" "//trim(nks_str)//" "//trim(nke_str)//'">'

         write(11,'(A)')'      <CellData>'

         if(variable(1) .eq. 1) then
         write(11,'(A)')'        <DataArray type="Float32"  Name="Density" format="ascii">'
            do j=1, Nj
               do i=1,Ni
                write(11,*) q(i,j,1)
               enddo
            enddo
         write(11,'(A)')'        </DataArray>'
        endif !! end of if(variable(1) .eq. 1) then

         if(variable(2) .eq. 1) then
         write(11,'(A)')'        <DataArray type="Float32"  Name="Vx" format="ascii">'
            do j=1, Nj
               do i=1,Ni
                write(11,*) q(i,j,2)/q(i,j,1)
               enddo
            enddo
         write(11,'(A)')'        </DataArray>'
        endif !! end of if(variable(2) .eq. 1) then


         if(variable(3) .eq. 1) then
         write(11,'(A)')'        <DataArray type="Float32"  Name="Vy" format="ascii">'
            do j=1, Nj
               do i=1,Ni
                write(11,*) q(i,j,3)/q(i,j,1)
               enddo
            enddo
         write(11,'(A)')'        </DataArray>'
        endif !! end of if(variable(3) .eq. 1) then

         if(variable(4) .eq. 1) then
         write(11,'(A)')'        <DataArray type="Float32"  Name="Vz" format="ascii">'
            do j=1, Nj
               do i=1,Ni
                write(11,*) q(i,j,4)/q(i,j,1)
               enddo
            enddo
         write(11,'(A)')'        </DataArray>'
        endif !! end of if(variable(4) .eq. 1) then

         if(variable(5) .eq. 1) then
         write(11,'(A)')'        <DataArray type="Float32"  Name="Bxl" format="ascii">'
            do j=1, Nj
               do i=1,Ni
                write(11,*) q(i,j,5)
               enddo
            enddo
         write(11,'(A)')'        </DataArray>'
        endif !! end of if(variable(5) .eq. 1) then

        if(variable(6) .eq. 1) then
         write(11,'(A)')'        <DataArray type="Float32"  Name="Byl" format="ascii">'
            do j=1, Nj
               do i=1,Ni
                write(11,*) q(i,j,6)
               enddo
            enddo
         write(11,'(A)')'        </DataArray>'
        endif !! end of if(variable(6) .eq. 1) then

        if(variable(7) .eq. 1) then
         write(11,'(A)')'        <DataArray type="Float32"  Name="Bzl" format="ascii">'
            do j=1, Nj
               do i=1,Ni
                write(11,*) q(i,j,7)
               enddo
            enddo
         write(11,'(A)')'        </DataArray>'
        endif !! end of if(variable(7) .eq. 1) then

        if(variable(8) .eq. 1) then
         !write(11,'(A)')'        <DataArray type="Float32"  Name="Pressure" format="ascii">'
          write(11,'(A)')'        <DataArray type="Float32"  Name="Energy" format="ascii">'
            do j=1, Nj
               do i=1,Ni

                 write(11,*) q(i,j,8)
                 !if(this%eosType=2) then
                 !write(11,*) q(i,j,8)
                 !q(i,j,8)=p/(gam-1.d0)+0.5d0*rho*(vx**2.d0+vy**2.d0+vz**2.d0)+0.5d0*(bxc**2.d0+byc**2.d0+bzc**2.d0)
                 !else
                 !endif
                
               enddo
            enddo
         write(11,'(A)')'        </DataArray>'
        endif !! end of if(variable(8) .eq. 1) then




         if(variable(5) .eq. 1) then
         write(11,'(A)')'        <DataArray type="Float32"  Name="Bxr" format="ascii">'
            do j=1, Nj
               do i=1,Ni
                write(11,*) q(i,j,9)
               enddo
            enddo
         write(11,'(A)')'        </DataArray>'
        endif !! end of if(variable(5) .eq. 1) then

        if(variable(6) .eq. 1) then
         write(11,'(A)')'        <DataArray type="Float32"  Name="Byr" format="ascii">'
            do j=1, Nj
               do i=1,Ni
                write(11,*) q(i,j,10)
               enddo
            enddo
         write(11,'(A)')'        </DataArray>'
        endif !! end of if(variable(6) .eq. 1) then

        if(variable(7) .eq. 1) then
         write(11,'(A)')'        <DataArray type="Float32"  Name="Bzr" format="ascii">'
            do j=1, Nj
               do i=1,Ni
                write(11,*) q(i,j,11)
               enddo
            enddo
         write(11,'(A)')'        </DataArray>'
        endif !! end of if(variable(7) .eq. 1) then




         write(11,'(A)')'      </CellData>'



         write(11,'(A)')'      <Points>'

         !write(*,*)"myid= ",myid,Ni, Nj

         write(11,'(A)')'        <DataArray type="Float32" NumberOfComponents="3" format="ascii">'
            do j=1, Nj+1
               do i=1,Ni+1
                write(11,*) xl(i), yl(j),0.0
               enddo
            enddo
         write(11,'(A)')'        </DataArray>'

         write(11,'(A)')'      </Points>'

         write(11,'(A)')'    </Piece>'
         write(11,'(A)')'  </StructuredGrid>'
         write(11,'(A)')'</VTKFile>'


close(11)

!call MPI_BARRIER(MPI_COMM_WORLD,ierr)
deallocate(nis1,nie1)
deallocate(njs1,nje1)
deallocate(nks1,nke1)

!stop

end subroutine output2d_MHD_pvts



subroutine output2d_HD_pvts(this,q)
use gridModule
implicit none
class(grid)::this
double precision,dimension(1-this%nbuf:this%nMesh(1)+this%nbuf,1-this%nbuf:this%nMesh(2)+this%nbuf &
                         &,this%nvar)::q 

double precision,dimension(1-this%nbuf:this%nMesh(1)+this%nbuf)::xl
double precision,dimension(1-this%nbuf:this%nMesh(2)+this%nbuf)::yl
integer,dimension(1:this%ndim)::dims  

integer:: ndim

integer:: gnis,gnie,gnjs,gnje,gnks,gnke
character(len=50)::gnis_str,gnie_str,gnjs_str,gnje_str ,gnks_str,gnke_str


integer:: nis,nie,njs,nje,nks,nke
character(len=50)::nis_str,nie_str,njs_str,nje_str,nks_str,nke_str

integer:: nis1_send,nie1_send,njs1_send,nje1_send,nks1_send,nke1_send
integer:: nis1_rec, nie1_rec, njs1_rec, nje1_rec, nks1_rec,nke1_rec


integer,dimension(:),allocatable:: nis1,nie1,njs1,nje1,nks1,nke1
character(len=50)::nis1_str,nie1_str,njs1_str,nje1_str,nks1_str,nke1_str



integer:: i,j  !k
integer::Ni,Nj !Nk

integer::iproc,Totalprocs

character(len=50)::iproc_str

integer::ierr
integer:: variable(8)

logical:: fileDirhave

integer:: fileid,rankid,gridid
character(len=50)::fileid_str,rankid_str,gridid_str

character(len=50):: dirname,filename_pvts,filename_vts


xl=this%xl(1)%coords
yl=this%xl(2)%coords

variable=this%variable
ndim=this%ndim
dims=this%dims_mpi

Totalprocs=nprocs

!write(*,*)"ndim=",ndim, "dims=",dims, "TotalProcess",nprocs,"nBlockx=",nBlockx,"nBlocky=",nBlocky
!write(*,*)"this%nMesh(1)=",this%nMesh(1),"this%nMesh(2)",this%nMesh(2)
!write(*,*)"myid= ",myid, "this%gnx_list_R=", this%gnx_list_R, "this%gny_list_R=", this%gny_list_R

nie=this%gnx_list_R
nje=this%gny_list_R
nke=0

nis=this%gnx_list_R-this%nMesh(1)
njs=this%gny_list_R-this%nMesh(2)
nks=0

call i4_to_s_left (nis, nis_str )
call i4_to_s_left (njs, njs_str )
call i4_to_s_left (nks, nks_str )

call i4_to_s_left (nie, nie_str )
call i4_to_s_left (nje, nje_str )
call i4_to_s_left (nke, nke_str )


allocate(nis1(0:Totalprocs-1),nie1(0:Totalprocs-1))
allocate(njs1(0:Totalprocs-1),nje1(0:Totalprocs-1))
allocate(nks1(0:Totalprocs-1),nke1(0:Totalprocs-1))


nis1=0
nie1=0

njs1=0
nje1=0

nks1=0
nke1=0




!integer:: nis1_send,nie1_send,njs1_send,nje1_send,nks1_send,nke1_send
!integer:: nis1_rec, nie1_rec, njs1_rec, nje1_rec, nks1_rec,nke1_rec

!do iproc=0, Totalprocs

!if(myid==0 .and. iproc==0) then
   !nis1(iproc)=nis
   !nie1(iproc)=nie

   !njs1(iproc)=njs
   !nje1(iproc)=nje

   !nks1(iproc)=nks
   !nke1(iproc)=nke
!endif

!if(myid>0 .and. myid==iproc) then
   !nis1_send=nis
   !nie1_send=nie

   !njs1_send=njs
   !nje1_send=nje

   !nks1_send=nks
   !nke1_send=nke

!endif

!enddo


!call MPI_GATHER(nis,1,MPI_INT,nis1,1,MPI_INT,0,MPI_COMM_WORLD,ierr)
!call MPI_GATHER(nie,1,MPI_INT,nie1,1,MPI_INT,0,MPI_COMM_WORLD,ierr)

!call MPI_GATHER(njs,1,MPI_INT,njs1,1,MPI_INT,0,MPI_COMM_WORLD,ierr)
!call MPI_GATHER(nje,1,MPI_INT,nje1,1,MPI_INT,0,MPI_COMM_WORLD,ierr)

call MPI_ALLGATHER(nis,1,MPI_INT,nis1,1,MPI_INT,MPI_COMM_WORLD,ierr)
call MPI_ALLGATHER(nie,1,MPI_INT,nie1,1,MPI_INT,MPI_COMM_WORLD,ierr)

call MPI_ALLGATHER(njs,1,MPI_INT,njs1,1,MPI_INT,MPI_COMM_WORLD,ierr)
call MPI_ALLGATHER(nje,1,MPI_INT,nje1,1,MPI_INT,MPI_COMM_WORLD,ierr)

!if(myid==0)then
!write(*,*)"myid= ", myid
!write(*,*)nis1
!write(*,*)nie1
!write(*,*)njs1
!write(*,*)nje1
!endif
!stop



!MPI_GATHER(SENDBUF, SENDCOUNT, SENDTYPE, RECVBUF, RECVCOUNT, RECVTYPE, ROOT, COMM, IERROR)
        !<type> SENDBUF(*), RECVBUF(*)
        !INTEGER SENDCOUNT, SENDTYPE, RECVCOUNT, RECVTYPE, ROOT, COMM, IERROR
!call MPI_BARRIER(MPI_COMM_WORLD,ierr)


gnis=0
gnjs=0
gnks=0

gnie=this%nMesh_global(1)
gnje=this%nMesh_global(2)
gnke=0


call i4_to_s_left (gnis, gnis_str )
call i4_to_s_left (gnjs, gnjs_str )
call i4_to_s_left (gnks, gnks_str )

call i4_to_s_left (gnie, gnie_str )
call i4_to_s_left (gnje, gnje_str )
call i4_to_s_left (gnke, gnke_str )


Ni=this%nMesh(1)
Nj=this%nMesh(2)


!write(*,*)"myid= ",myid,"nis=",nis,"njs=",njs,"nie=",nie,"nje=",nje
!write(*,*)"myid= ",myid,"gnie=",gnie,"gnje=",gnje




gridid=this%gridID
fileid=this%fnum
rankid=myid

call i4_to_s_left (gridid, gridid_str )
call i4_to_s_left (fileid, fileid_str )
call i4_to_s_left (rankid, rankid_str )


dirname="HD_g"//trim(gridid_str)//"_"//trim(fileid_str)

inquire(file=trim(dirname),exist=fileDirhave)

if(.not.fileDirhave)then
   !if(myid==0)then
    call system("mkdir "//trim(dirname))
   !endif  !!! end of if(myid==0)then
endif  ! end of if(.not.fileDirhave)then

!call MPI_BARRIER(MPI_COMM_WORLD,ierr)

!filename_pvts=trim(dirname)//'/'//trim(dirname)//'.pvts'
filename_pvts=trim(dirname)//'/'//trim(dirname)//'_'//trim(rankid_str)//'.pvts'
!if(myid==0)then
    open(10,file=trim(filename_pvts))

          write(10,'(A)')'<?xml version="1.0"?>'
          write(10,'(A)')'<VTKFile type="PStructuredGrid" version="0.1" byte_order="LittleEndian">'
          write(10,'(A)')'  <PStructuredGrid WholeExtent="'//trim(gnis_str)//" "//trim(gnie_str)// &
          & " "//trim(gnjs_str)//" "//trim(gnje_str)//" "//trim(gnks_str)//" "//trim(gnke_str)//'" GhostLevel="#">'

          write(10,'(A)')'    <PCellData>'

          if(variable(1) .eq. 1) then
            write(10,'(A)')'      <PDataArray type="Float32" Name="Density"/>'
          endif

          if(variable(2) .eq. 1) then
            write(10,'(A)')'      <PDataArray type="Float32" Name="Vx"/>'
          endif

          if(variable(3) .eq. 1) then
            write(10,'(A)')'      <PDataArray type="Float32" Name="Vy"/>'
          endif

          !if(variable(4) .eq. 1) then
            !write(10,'(A)')'      <PDataArray type="Float32" Name="Vz"/>'
          !endif

          !if(variable(5) .eq. 1) then
            !write(10,'(A)')'      <PDataArray type="Float32" Name="Bx"/>'
          !endif

          !if(variable(6) .eq. 1) then
            !write(10,'(A)')'      <PDataArray type="Float32" Name="By"/>'
          !endif

          !if(variable(7) .eq. 1) then
            !write(10,'(A)')'      <PDataArray type="Float32" Name="Bz"/>'
          !endif

          if(variable(8) .eq. 1) then
            !write(10,'(A)')'      <PDataArray type="Float32" Name="Pressure"/>'
            write(10,'(A)')'      <PDataArray type="Float32" Name="Energy"/>'
          endif



          write(10,'(A)')'    </PCellData>'
          


          write(10,'(A)')'    <PPoints>'
          write(10,'(A)')'      <PDataArray type="Float32" NumberOfComponents="3"/>'
          write(10,'(A)')'    </PPoints>'


          do iproc=0, Totalprocs-1
               !if(myid==0) then
                 !write(*,*)"_____", iproc
               !endif

               call i4_to_s_left (iproc, iproc_str)

               call i4_to_s_left (nis1(iproc), nis1_str)
               call i4_to_s_left (nie1(iproc), nie1_str)

               call i4_to_s_left (njs1(iproc), njs1_str)
               call i4_to_s_left (nje1(iproc), nje1_str)

               call i4_to_s_left (nks1(iproc), nks1_str)
               call i4_to_s_left (nke1(iproc), nke1_str)

               !call MPI_BARRIER(MPI_COMM_WORLD,ierr)
               write(10,'(A)')'    <Piece Extent="'//trim(nis1_str)//" "//trim(nie1_str)// &
               & " "//trim(njs1_str)//" "//trim(nje1_str)//" "//trim(nks1_str)//" "//trim(nke1_str)// &
               &'"'//' Source="'//trim(iproc_str)//'.vts"/>'      
          enddo

          write(10,'(A)')'  </PStructuredGrid>'
          write(10,'(A)')'</VTKFile>'

    close(10)
!endif  !!! end of if(myid==0)then


   !call MPI_BARRIER(MPI_COMM_WORLD,ierr)

filename_vts=trim(dirname)//'/'//trim(rankid_str)//'.vts'

open(11,file=filename_vts)

         write(11,'(A)')'<?xml version="1.0"?>'
         write(11,'(A)')'<VTKFile type="StructuredGrid" version="0.1" byte_order="LittleEndian">'
         write(11,'(A)')'  <StructuredGrid WholeExtent="'//trim(nis_str)//" "//trim(nie_str)//" " &
                           & //trim(njs_str)//" "//trim(nje_str)//" "//trim(nks_str)//" "//trim(nke_str)//'">'
         write(11,'(A)')'    <Piece Extent="'//trim(nis_str)//" "//trim(nie_str)//" " &
                           & //trim(njs_str)//" "//trim(nje_str)//" "//trim(nks_str)//" "//trim(nke_str)//'">'

         write(11,'(A)')'      <CellData>'

         if(variable(1) .eq. 1) then
         write(11,'(A)')'        <DataArray type="Float32"  Name="Density" format="ascii">'
            do j=1, Nj
               do i=1,Ni
                write(11,*) q(i,j,1)
               enddo
            enddo
         write(11,'(A)')'        </DataArray>'
        endif !! end of if(variable(1) .eq. 1) then

         if(variable(2) .eq. 1) then
         write(11,'(A)')'        <DataArray type="Float32"  Name="Vx" format="ascii">'
            do j=1, Nj
               do i=1,Ni
                write(11,*) q(i,j,2)/q(i,j,1)
               enddo
            enddo
         write(11,'(A)')'        </DataArray>'
        endif !! end of if(variable(2) .eq. 1) then


         if(variable(3) .eq. 1) then
         write(11,'(A)')'        <DataArray type="Float32"  Name="Vy" format="ascii">'
            do j=1, Nj
               do i=1,Ni
                write(11,*) q(i,j,3)/q(i,j,1)
               enddo
            enddo
         write(11,'(A)')'        </DataArray>'
        endif !! end of if(variable(3) .eq. 1) then

         !if(variable(4) .eq. 1) then
         !write(11,'(A)')'        <DataArray type="Float32"  Name="Vz" format="ascii">'
            !do j=1, Nj
              ! do i=1,Ni
                !write(11,*) q(i,j,4)/q(i,j,1)
               !enddo
            !enddo
         !write(11,'(A)')'        </DataArray>'
        !endif !! end of if(variable(4) .eq. 1) then

         !if(variable(5) .eq. 1) then
         !write(11,'(A)')'        <DataArray type="Float32"  Name="Bx" format="ascii">'
            !do j=1, Nj
               !do i=1,Ni
                !write(11,*) q(i,j,5)
               !enddo
            !enddo
         !write(11,'(A)')'        </DataArray>'
        !endif !! end of if(variable(5) .eq. 1) then

        !if(variable(6) .eq. 1) then
         !write(11,'(A)')'        <DataArray type="Float32"  Name="By" format="ascii">'
            !do j=1, Nj
               !do i=1,Ni
                !write(11,*) q(i,j,6)
              ! enddo
            !enddo
         !write(11,'(A)')'        </DataArray>'
        !endif !! end of if(variable(6) .eq. 1) then

        !if(variable(7) .eq. 1) then
         !write(11,'(A)')'        <DataArray type="Float32"  Name="Bz" format="ascii">'
            !do j=1, Nj
               !do i=1,Ni
                !write(11,*) q(i,j,7)
               !enddo
            !enddo
         !write(11,'(A)')'        </DataArray>'
        !endif !! end of if(variable(7) .eq. 1) then



        if(variable(8) .eq. 1) then   !!!!!! for 2D HD problem, the variable of energy is variable(8) 
         !write(11,'(A)')'        <DataArray type="Float32"  Name="Pressure" format="ascii">'
          write(11,'(A)')'        <DataArray type="Float32"  Name="Energy" format="ascii">'
            do j=1, Nj
               do i=1,Ni

                 write(11,*) q(i,j,4)

                 !if(this%eosType=2) then
                 !write(11,*) q(i,j,4)
                 !q(i,j,8)=p/(gam-1.d0)+0.5d0*rho*(vx**2.d0+vy**2.d0+vz**2.d0)+0.5d0*(bxc**2.d0+byc**2.d0+bzc**2.d0)
                 !else
                 !endif
                
               enddo
            enddo
         write(11,'(A)')'        </DataArray>'
        endif !! end of if(variable(8) .eq. 1) then





         write(11,'(A)')'      </CellData>'



         write(11,'(A)')'      <Points>'

         !write(*,*)"myid= ",myid,Ni, Nj

         write(11,'(A)')'        <DataArray type="Float32" NumberOfComponents="3" format="ascii">'
            do j=1, Nj+1
               do i=1,Ni+1
                write(11,*) xl(i), yl(j),0.0
               enddo
            enddo
         write(11,'(A)')'        </DataArray>'

         write(11,'(A)')'      </Points>'

         write(11,'(A)')'    </Piece>'
         write(11,'(A)')'  </StructuredGrid>'
         write(11,'(A)')'</VTKFile>'




close(11)

!call MPI_BARRIER(MPI_COMM_WORLD,ierr)

deallocate(nis1,nie1)
deallocate(njs1,nje1)
deallocate(nks1,nke1)

!stop

end subroutine output2d_HD_pvts




subroutine output3d_MHD_pvts(this,q)
use gridModule
implicit none
class(grid)::this
double precision,dimension(1-this%nbuf:this%nMesh(1)+this%nbuf,1-this%nbuf:this%nMesh(2)+this%nbuf &
                         &,1-this%nbuf:this%nMesh(3)+this%nbuf,this%nvar)::q 

double precision,dimension(1-this%nbuf:this%nMesh(1)+this%nbuf)::xl
double precision,dimension(1-this%nbuf:this%nMesh(2)+this%nbuf)::yl
double precision,dimension(1-this%nbuf:this%nMesh(3)+this%nbuf)::zl




integer,dimension(1:this%ndim)::dims  

integer:: ndim

integer:: gnis,gnie,gnjs,gnje,gnks,gnke
character(len=50)::gnis_str,gnie_str,gnjs_str,gnje_str ,gnks_str,gnke_str


integer:: nis,nie,njs,nje,nks,nke
character(len=50)::nis_str,nie_str,njs_str,nje_str,nks_str,nke_str

integer:: nis1_send,nie1_send,njs1_send,nje1_send,nks1_send,nke1_send
integer:: nis1_rec, nie1_rec, njs1_rec, nje1_rec, nks1_rec,nke1_rec


integer,dimension(:),allocatable:: nis1,nie1,njs1,nje1,nks1,nke1
character(len=50)::nis1_str,nie1_str,njs1_str,nje1_str,nks1_str,nke1_str



integer:: i,j, k
integer::Ni,Nj,Nk

integer::iproc,Totalprocs

character(len=50)::iproc_str

integer::ierr
integer:: variable(8)

logical:: fileDirhave

integer:: fileid,rankid,gridid
character(len=50)::fileid_str,rankid_str,gridid_str

character(len=50):: dirname,filename_pvts,filename_vts


xl=this%xl(1)%coords
yl=this%xl(2)%coords
zl=this%xl(3)%coords

variable=this%variable
ndim=this%ndim
dims=this%dims_mpi

Totalprocs=nprocs

!write(*,*)"ndim=",ndim, "dims=",dims, "TotalProcess",nprocs,"nBlockx=",nBlockx,"nBlocky=",nBlocky
!write(*,*)"this%nMesh(1)=",this%nMesh(1),"this%nMesh(2)",this%nMesh(2),"this%nMesh(3)",this%nMesh(3)
!write(*,*)"myid= ",myid, "this%gnx_list_R=", this%gnx_list_R, "this%gny_list_R=", this%gny_list_R, & 
!                       & "this%gnz_list_R=", this%gnz_list_R

!stop

nie=this%gnx_list_R
nje=this%gny_list_R
nke=this%gnz_list_R

nis=this%gnx_list_R-this%nMesh(1)
njs=this%gny_list_R-this%nMesh(2)
nks=this%gnz_list_R-this%nMesh(3)

call i4_to_s_left (nis, nis_str )
call i4_to_s_left (njs, njs_str )
call i4_to_s_left (nks, nks_str )

call i4_to_s_left (nie, nie_str )
call i4_to_s_left (nje, nje_str )
call i4_to_s_left (nke, nke_str )


allocate(nis1(0:Totalprocs-1),nie1(0:Totalprocs-1))
allocate(njs1(0:Totalprocs-1),nje1(0:Totalprocs-1))
allocate(nks1(0:Totalprocs-1),nke1(0:Totalprocs-1))


nis1=0
nie1=0

njs1=0
nje1=0

nks1=0
nke1=0



!call MPI_GATHER(nis,1,MPI_INT,nis1,1,MPI_INT,0,MPI_COMM_WORLD,ierr)
!call MPI_GATHER(nie,1,MPI_INT,nie1,1,MPI_INT,0,MPI_COMM_WORLD,ierr)


!call MPI_GATHER(njs,1,MPI_INT,njs1,1,MPI_INT,0,MPI_COMM_WORLD,ierr)
!call MPI_GATHER(nje,1,MPI_INT,nje1,1,MPI_INT,0,MPI_COMM_WORLD,ierr)

!call MPI_GATHER(nks,1,MPI_INT,nks1,1,MPI_INT,0,MPI_COMM_WORLD,ierr)
!call MPI_GATHER(nke,1,MPI_INT,nke1,1,MPI_INT,0,MPI_COMM_WORLD,ierr)


call MPI_ALLGATHER(nis,1,MPI_INT,nis1,1,MPI_INT,MPI_COMM_WORLD,ierr)
call MPI_ALLGATHER(nie,1,MPI_INT,nie1,1,MPI_INT,MPI_COMM_WORLD,ierr)


call MPI_ALLGATHER(njs,1,MPI_INT,njs1,1,MPI_INT,MPI_COMM_WORLD,ierr)
call MPI_ALLGATHER(nje,1,MPI_INT,nje1,1,MPI_INT,MPI_COMM_WORLD,ierr)

call MPI_ALLGATHER(nks,1,MPI_INT,nks1,1,MPI_INT,MPI_COMM_WORLD,ierr)
call MPI_ALLGATHER(nke,1,MPI_INT,nke1,1,MPI_INT,MPI_COMM_WORLD,ierr)

!if(myid==0)then
!write(*,*)nis1
!write(*,*)nie1
!write(*,*)njs1
!write(*,*)nje1
!write(*,*)nks1
!write(*,*)nke1
!endif
!stop



!MPI_GATHER(SENDBUF, SENDCOUNT, SENDTYPE, RECVBUF, RECVCOUNT, RECVTYPE, ROOT, COMM, IERROR)
        !<type> SENDBUF(*), RECVBUF(*)
        !INTEGER SENDCOUNT, SENDTYPE, RECVCOUNT, RECVTYPE, ROOT, COMM, IERROR
!call MPI_BARRIER(MPI_COMM_WORLD,ierr)


gnis=0
gnjs=0
gnks=0

gnie=this%nMesh_global(1)
gnje=this%nMesh_global(2)
gnke=this%nMesh_global(3)


call i4_to_s_left (gnis, gnis_str )
call i4_to_s_left (gnjs, gnjs_str )
call i4_to_s_left (gnks, gnks_str )

call i4_to_s_left (gnie, gnie_str )
call i4_to_s_left (gnje, gnje_str )
call i4_to_s_left (gnke, gnke_str )


Ni=this%nMesh(1)
Nj=this%nMesh(2)
Nk=this%nMesh(3)


!write(*,*)"myid= ",myid,"nis=",nis,"njs=",njs,"nks=",nks,"nie=",nie,"nje=",nje,"nke=",nke
!write(*,*)"myid= ",myid,"gnie=",gnie,"gnje=",gnje,"gnke=",gnke




gridid=this%gridID
fileid=this%fnum
rankid=myid

call i4_to_s_left (gridid, gridid_str )
call i4_to_s_left (fileid, fileid_str )
call i4_to_s_left (rankid, rankid_str )


dirname="MHD_g"//trim(gridid_str)//"_"//trim(fileid_str)

inquire(file=trim(dirname),exist=fileDirhave)

if(.not.fileDirhave)then
   !if(myid==0)then
    call system("mkdir "//trim(dirname))
   !endif  !!! end of if(myid==0)then
endif  ! end of if(.not.fileDirhave)then

!call MPI_BARRIER(MPI_COMM_WORLD,ierr)



!filename_pvts=trim(dirname)//'/'//trim(dirname)//'.pvts'
filename_pvts=trim(dirname)//'/'//trim(dirname)//'_'//trim(rankid_str)//'.pvts'

!if(myid==0)then
    open(10,file=trim(filename_pvts))

          write(10,'(A)')'<?xml version="1.0"?>'
          write(10,'(A)')'<VTKFile type="PStructuredGrid" version="0.1" byte_order="LittleEndian">'
          write(10,'(A)')'  <PStructuredGrid WholeExtent="'//trim(gnis_str)//" "//trim(gnie_str)// &
          & " "//trim(gnjs_str)//" "//trim(gnje_str)//" "//trim(gnks_str)//" "//trim(gnke_str)//'" GhostLevel="#">'

          write(10,'(A)')'    <PCellData>'

          if(variable(1) .eq. 1) then
            write(10,'(A)')'      <PDataArray type="Float32" Name="Density"/>'
          endif

          if(variable(2) .eq. 1) then
            write(10,'(A)')'      <PDataArray type="Float32" Name="Vx"/>'
          endif

          if(variable(3) .eq. 1) then
            write(10,'(A)')'      <PDataArray type="Float32" Name="Vy"/>'
          endif

          if(variable(4) .eq. 1) then
            write(10,'(A)')'      <PDataArray type="Float32" Name="Vz"/>'
          endif

          if(variable(5) .eq. 1) then
            write(10,'(A)')'      <PDataArray type="Float32" Name="Bxl"/>'
          endif

          if(variable(6) .eq. 1) then
            write(10,'(A)')'      <PDataArray type="Float32" Name="Byl"/>'
          endif

          if(variable(7) .eq. 1) then
            write(10,'(A)')'      <PDataArray type="Float32" Name="Bzl"/>'
          endif

          if(variable(8) .eq. 1) then
            !write(10,'(A)')'      <PDataArray type="Float32" Name="Pressure"/>'
            write(10,'(A)')'      <PDataArray type="Float32" Name="Energy"/>'
          endif

          if(variable(5) .eq. 1) then
            write(10,'(A)')'      <PDataArray type="Float32" Name="Bxr"/>'
          endif

          if(variable(6) .eq. 1) then
            write(10,'(A)')'      <PDataArray type="Float32" Name="Byr"/>'
          endif

          if(variable(7) .eq. 1) then
            write(10,'(A)')'      <PDataArray type="Float32" Name="Bzr"/>'
          endif


          write(10,'(A)')'    </PCellData>'
          


          write(10,'(A)')'    <PPoints>'
          write(10,'(A)')'      <PDataArray type="Float32" NumberOfComponents="3"/>'
          write(10,'(A)')'    </PPoints>'


          do iproc=0, Totalprocs-1
               !if(myid==0) then
                 !write(*,*)"_____", iproc
               !endif

               call i4_to_s_left (iproc, iproc_str)

               call i4_to_s_left (nis1(iproc), nis1_str)
               call i4_to_s_left (nie1(iproc), nie1_str)

               call i4_to_s_left (njs1(iproc), njs1_str)
               call i4_to_s_left (nje1(iproc), nje1_str)

               call i4_to_s_left (nks1(iproc), nks1_str)
               call i4_to_s_left (nke1(iproc), nke1_str)

               !call MPI_BARRIER(MPI_COMM_WORLD,ierr)
               write(10,'(A)')'    <Piece Extent="'//trim(nis1_str)//" "//trim(nie1_str)// &
               & " "//trim(njs1_str)//" "//trim(nje1_str)//" "//trim(nks1_str)//" "//trim(nke1_str)// &
               &'"'//' Source="'//trim(iproc_str)//'.vts"/>'      
          enddo

          write(10,'(A)')'  </PStructuredGrid>'
          write(10,'(A)')'</VTKFile>'

    close(10)
!endif  !!! end of if(myid==0)then

   !call MPI_BARRIER(MPI_COMM_WORLD,ierr)

filename_vts=trim(dirname)//'/'//trim(rankid_str)//'.vts'

open(11,file=filename_vts)


         write(11,'(A)')'<?xml version="1.0"?>'
         write(11,'(A)')'<VTKFile type="StructuredGrid" version="0.1" byte_order="LittleEndian">'
         write(11,'(A)')'  <StructuredGrid WholeExtent="'//trim(nis_str)//" "//trim(nie_str)//" " &
                           & //trim(njs_str)//" "//trim(nje_str)//" "//trim(nks_str)//" "//trim(nke_str)//'">'
         write(11,'(A)')'    <Piece Extent="'//trim(nis_str)//" "//trim(nie_str)//" " &
                           & //trim(njs_str)//" "//trim(nje_str)//" "//trim(nks_str)//" "//trim(nke_str)//'">'

         write(11,'(A)')'      <CellData>'

         if(variable(1) .eq. 1) then
         write(11,'(A)')'        <DataArray type="Float32"  Name="Density" format="ascii">'
         do k=1, Nk
            do j=1, Nj
               do i=1,Ni
                write(11,*) q(i,j,k,1)
               enddo
            enddo
         enddo
         write(11,'(A)')'        </DataArray>'
        endif !! end of if(variable(1) .eq. 1) then

         if(variable(2) .eq. 1) then
         write(11,'(A)')'        <DataArray type="Float32"  Name="Vx" format="ascii">'
         do k=1,Nk
            do j=1, Nj
               do i=1,Ni
                write(11,*) q(i,j,k,2)/q(i,j,k,1)
               enddo
            enddo
         enddo
         write(11,'(A)')'        </DataArray>'
        endif !! end of if(variable(2) .eq. 1) then


         if(variable(3) .eq. 1) then
         write(11,'(A)')'        <DataArray type="Float32"  Name="Vy" format="ascii">'
         do k=1,Nk
            do j=1, Nj
               do i=1,Ni
                write(11,*) q(i,j,k,3)/q(i,j,k,1)
               enddo
            enddo
         enddo
         write(11,'(A)')'        </DataArray>'
        endif !! end of if(variable(3) .eq. 1) then

         if(variable(4) .eq. 1) then
         write(11,'(A)')'        <DataArray type="Float32"  Name="Vz" format="ascii">'
         do k=1,Nk
            do j=1, Nj
               do i=1,Ni
                write(11,*) q(i,j,k,4)/q(i,j,k,1)
               enddo
            enddo
         enddo
         write(11,'(A)')'        </DataArray>'
        endif !! end of if(variable(4) .eq. 1) then

         if(variable(5) .eq. 1) then
         write(11,'(A)')'        <DataArray type="Float32"  Name="Bxl" format="ascii">'
          do k=1,Nk
            do j=1, Nj
               do i=1,Ni
                write(11,*) q(i,j,k,5)
               enddo
            enddo
          enddo
         write(11,'(A)')'        </DataArray>'
        endif !! end of if(variable(5) .eq. 1) then

        if(variable(6) .eq. 1) then
         write(11,'(A)')'        <DataArray type="Float32"  Name="Byl" format="ascii">'
          do k=1,Nk
            do j=1, Nj
               do i=1,Ni
                write(11,*) q(i,j,k,6)
               enddo
            enddo
         enddo
         write(11,'(A)')'        </DataArray>'
        endif !! end of if(variable(6) .eq. 1) then

        if(variable(7) .eq. 1) then
         write(11,'(A)')'        <DataArray type="Float32"  Name="Bzl" format="ascii">'
          do k=1,Nk
            do j=1, Nj
               do i=1,Ni
                write(11,*) q(i,j,k,7)
               enddo
            enddo
          enddo
         write(11,'(A)')'        </DataArray>'
        endif !! end of if(variable(7) .eq. 1) then


        if(variable(8) .eq. 1) then
         !write(11,'(A)')'        <DataArray type="Float32"  Name="Pressure" format="ascii">'
          write(11,'(A)')'        <DataArray type="Float32"  Name="Energy" format="ascii">'
          do k=1,Nk
            do j=1, Nj
               do i=1,Ni

                 write(11,*) q(i,j,k,8)
                   !if(this%eosType=2) then
                   !write(11,*) q(i,j,8)
                   !q(i,j,8)=p/(gam-1.d0)+0.5d0*rho*(vx**2.d0+vy**2.d0+vz**2.d0)+0.5d0*(bxc**2.d0+byc**2.d0+bzc**2.d0)
                   !else
                   !endif
                
               enddo
            enddo
         enddo
         write(11,'(A)')'        </DataArray>'
        endif !! end of if(variable(8) .eq. 1) then


         if(variable(5) .eq. 1) then
         write(11,'(A)')'        <DataArray type="Float32"  Name="Bxr" format="ascii">'
          do k=1,Nk
            do j=1, Nj
               do i=1,Ni
                write(11,*) q(i,j,k,9)
               enddo
            enddo
          enddo
         write(11,'(A)')'        </DataArray>'
        endif !! end of if(variable(5) .eq. 1) then

        if(variable(6) .eq. 1) then
         write(11,'(A)')'        <DataArray type="Float32"  Name="Byr" format="ascii">'
          do k=1,Nk
            do j=1, Nj
               do i=1,Ni
                write(11,*) q(i,j,k,10)
               enddo
            enddo
         enddo
         write(11,'(A)')'        </DataArray>'
        endif !! end of if(variable(6) .eq. 1) then

        if(variable(7) .eq. 1) then
         write(11,'(A)')'        <DataArray type="Float32"  Name="Bzr" format="ascii">'
          do k=1,Nk
            do j=1, Nj
               do i=1,Ni
                write(11,*) q(i,j,k,11)
               enddo
            enddo
          enddo
         write(11,'(A)')'        </DataArray>'
        endif !! end of if(variable(7) .eq. 1) then





         write(11,'(A)')'      </CellData>'



         write(11,'(A)')'      <Points>'

         !write(*,*)"myid= ",myid,Ni, Nj,Nk

         write(11,'(A)')'        <DataArray type="Float32" NumberOfComponents="3" format="ascii">'
          do k=1,Nk+1
            do j=1, Nj+1
               do i=1,Ni+1
                write(11,*) xl(i), yl(j),zl(k)
               enddo
            enddo
          enddo
         write(11,'(A)')'        </DataArray>'

         write(11,'(A)')'      </Points>'

         write(11,'(A)')'    </Piece>'
         write(11,'(A)')'  </StructuredGrid>'
         write(11,'(A)')'</VTKFile>'




close(11)

!call MPI_BARRIER(MPI_COMM_WORLD,ierr)

deallocate(nis1,nie1)
deallocate(njs1,nje1)
deallocate(nks1,nke1)


end subroutine output3d_MHD_pvts




subroutine output3d_HD_pvts(this,q)
use gridModule
implicit none
class(grid)::this
double precision,dimension(1-this%nbuf:this%nMesh(1)+this%nbuf,1-this%nbuf:this%nMesh(2)+this%nbuf &
                         &,1-this%nbuf:this%nMesh(3)+this%nbuf,this%nvar)::q 

double precision,dimension(1-this%nbuf:this%nMesh(1)+this%nbuf)::xl
double precision,dimension(1-this%nbuf:this%nMesh(2)+this%nbuf)::yl
double precision,dimension(1-this%nbuf:this%nMesh(3)+this%nbuf)::zl




integer,dimension(1:this%ndim)::dims  

integer:: ndim

integer:: gnis,gnie,gnjs,gnje,gnks,gnke
character(len=50)::gnis_str,gnie_str,gnjs_str,gnje_str ,gnks_str,gnke_str


integer:: nis,nie,njs,nje,nks,nke
character(len=50)::nis_str,nie_str,njs_str,nje_str,nks_str,nke_str

integer:: nis1_send,nie1_send,njs1_send,nje1_send,nks1_send,nke1_send
integer:: nis1_rec, nie1_rec, njs1_rec, nje1_rec, nks1_rec,nke1_rec


integer,dimension(:),allocatable:: nis1,nie1,njs1,nje1,nks1,nke1
character(len=50)::nis1_str,nie1_str,njs1_str,nje1_str,nks1_str,nke1_str



integer:: i,j, k
integer::Ni,Nj,Nk

integer::iproc,Totalprocs

character(len=50)::iproc_str

integer::ierr
integer:: variable(8)

logical:: fileDirhave

integer:: fileid,rankid,gridid
character(len=50)::fileid_str,rankid_str,gridid_str

character(len=50):: dirname,filename_pvts,filename_vts


xl=this%xl(1)%coords
yl=this%xl(2)%coords
zl=this%xl(3)%coords

variable=this%variable
ndim=this%ndim
dims=this%dims_mpi

Totalprocs=nprocs

!write(*,*)"ndim=",ndim, "dims=",dims, "TotalProcess",nprocs,"nBlockx=",nBlockx,"nBlocky=",nBlocky
!write(*,*)"this%nMesh(1)=",this%nMesh(1),"this%nMesh(2)",this%nMesh(2),"this%nMesh(3)",this%nMesh(3)
!write(*,*)"myid= ",myid, "this%gnx_list_R=", this%gnx_list_R, "this%gny_list_R=", this%gny_list_R, & 
!                       & "this%gnz_list_R=", this%gnz_list_R

!stop

nie=this%gnx_list_R
nje=this%gny_list_R
nke=this%gnz_list_R

nis=this%gnx_list_R-this%nMesh(1)
njs=this%gny_list_R-this%nMesh(2)
nks=this%gnz_list_R-this%nMesh(3)

call i4_to_s_left (nis, nis_str )
call i4_to_s_left (njs, njs_str )
call i4_to_s_left (nks, nks_str )

call i4_to_s_left (nie, nie_str )
call i4_to_s_left (nje, nje_str )
call i4_to_s_left (nke, nke_str )


allocate(nis1(0:Totalprocs-1),nie1(0:Totalprocs-1))
allocate(njs1(0:Totalprocs-1),nje1(0:Totalprocs-1))
allocate(nks1(0:Totalprocs-1),nke1(0:Totalprocs-1))


nis1=0
nie1=0

njs1=0
nje1=0

nks1=0
nke1=0



!call MPI_GATHER(nis,1,MPI_INT,nis1,1,MPI_INT,0,MPI_COMM_WORLD,ierr)
!call MPI_GATHER(nie,1,MPI_INT,nie1,1,MPI_INT,0,MPI_COMM_WORLD,ierr)


!call MPI_GATHER(njs,1,MPI_INT,njs1,1,MPI_INT,0,MPI_COMM_WORLD,ierr)
!call MPI_GATHER(nje,1,MPI_INT,nje1,1,MPI_INT,0,MPI_COMM_WORLD,ierr)

!call MPI_GATHER(nks,1,MPI_INT,nks1,1,MPI_INT,0,MPI_COMM_WORLD,ierr)
!call MPI_GATHER(nke,1,MPI_INT,nke1,1,MPI_INT,0,MPI_COMM_WORLD,ierr)

call MPI_ALLGATHER(nis,1,MPI_INT,nis1,1,MPI_INT,MPI_COMM_WORLD,ierr)
call MPI_ALLGATHER(nie,1,MPI_INT,nie1,1,MPI_INT,MPI_COMM_WORLD,ierr)


call MPI_ALLGATHER(njs,1,MPI_INT,njs1,1,MPI_INT,MPI_COMM_WORLD,ierr)
call MPI_ALLGATHER(nje,1,MPI_INT,nje1,1,MPI_INT,MPI_COMM_WORLD,ierr)

call MPI_ALLGATHER(nks,1,MPI_INT,nks1,1,MPI_INT,MPI_COMM_WORLD,ierr)
call MPI_ALLGATHER(nke,1,MPI_INT,nke1,1,MPI_INT,MPI_COMM_WORLD,ierr)

!if(myid==0)then
!write(*,*)nis1
!write(*,*)nie1
!write(*,*)njs1
!write(*,*)nje1
!write(*,*)nks1
!write(*,*)nke1
!endif
!stop



!MPI_GATHER(SENDBUF, SENDCOUNT, SENDTYPE, RECVBUF, RECVCOUNT, RECVTYPE, ROOT, COMM, IERROR)
        !<type> SENDBUF(*), RECVBUF(*)
        !INTEGER SENDCOUNT, SENDTYPE, RECVCOUNT, RECVTYPE, ROOT, COMM, IERROR
!call MPI_BARRIER(MPI_COMM_WORLD,ierr)


gnis=0
gnjs=0
gnks=0

gnie=this%nMesh_global(1)
gnje=this%nMesh_global(2)
gnke=this%nMesh_global(3)


call i4_to_s_left (gnis, gnis_str )
call i4_to_s_left (gnjs, gnjs_str )
call i4_to_s_left (gnks, gnks_str )

call i4_to_s_left (gnie, gnie_str )
call i4_to_s_left (gnje, gnje_str )
call i4_to_s_left (gnke, gnke_str )


Ni=this%nMesh(1)
Nj=this%nMesh(2)
Nk=this%nMesh(3)


!write(*,*)"myid= ",myid,"nis=",nis,"njs=",njs,"nks=",nks,"nie=",nie,"nje=",nje,"nke=",nke
!write(*,*)"myid= ",myid,"gnie=",gnie,"gnje=",gnje,"gnke=",gnke




gridid=this%gridID
fileid=this%fnum
rankid=myid

call i4_to_s_left (gridid, gridid_str )
call i4_to_s_left (fileid, fileid_str )
call i4_to_s_left (rankid, rankid_str )


dirname="HD_g"//trim(gridid_str)//"_"//trim(fileid_str)

inquire(file=trim(dirname),exist=fileDirhave)

if(.not.fileDirhave)then
   !if(myid==0)then
    call system("mkdir "//trim(dirname))
   !endif  !!! end of if(myid==0)then
endif  ! end of if(.not.fileDirhave)then

!call MPI_BARRIER(MPI_COMM_WORLD,ierr)



!filename_pvts=trim(dirname)//'/'//trim(dirname)//'.pvts'
filename_pvts=trim(dirname)//'/'//trim(dirname)//'_'//trim(rankid_str)//'.pvts'
!if(myid==0)then
    open(10,file=trim(filename_pvts))

          write(10,'(A)')'<?xml version="1.0"?>'
          write(10,'(A)')'<VTKFile type="PStructuredGrid" version="0.1" byte_order="LittleEndian">'
          write(10,'(A)')'  <PStructuredGrid WholeExtent="'//trim(gnis_str)//" "//trim(gnie_str)// &
          & " "//trim(gnjs_str)//" "//trim(gnje_str)//" "//trim(gnks_str)//" "//trim(gnke_str)//'" GhostLevel="#">'

          write(10,'(A)')'    <PCellData>'

          if(variable(1) .eq. 1) then
            write(10,'(A)')'      <PDataArray type="Float32" Name="Density"/>'
          endif

          if(variable(2) .eq. 1) then
            write(10,'(A)')'      <PDataArray type="Float32" Name="Vx"/>'
          endif

          if(variable(3) .eq. 1) then
            write(10,'(A)')'      <PDataArray type="Float32" Name="Vy"/>'
          endif

          if(variable(4) .eq. 1) then
            write(10,'(A)')'      <PDataArray type="Float32" Name="Vz"/>'
          endif

          !if(variable(5) .eq. 1) then
            !write(10,'(A)')'      <PDataArray type="Float32" Name="Bx"/>'
          !endif

          !if(variable(6) .eq. 1) then
            !write(10,'(A)')'      <PDataArray type="Float32" Name="By"/>'
          !endif

          !if(variable(7) .eq. 1) then
            !write(10,'(A)')'      <PDataArray type="Float32" Name="Bz"/>'
          !endif

          if(variable(8) .eq. 1) then
            !write(10,'(A)')'      <PDataArray type="Float32" Name="Pressure"/>'
            write(10,'(A)')'      <PDataArray type="Float32" Name="Energy"/>'
          endif



          write(10,'(A)')'    </PCellData>'
          


          write(10,'(A)')'    <PPoints>'
          write(10,'(A)')'      <PDataArray type="Float32" NumberOfComponents="3"/>'
          write(10,'(A)')'    </PPoints>'


          do iproc=0, Totalprocs-1
               !if(myid==0) then
                 !write(*,*)"_____", iproc
               !endif

               call i4_to_s_left (iproc, iproc_str)

               call i4_to_s_left (nis1(iproc), nis1_str)
               call i4_to_s_left (nie1(iproc), nie1_str)

               call i4_to_s_left (njs1(iproc), njs1_str)
               call i4_to_s_left (nje1(iproc), nje1_str)

               call i4_to_s_left (nks1(iproc), nks1_str)
               call i4_to_s_left (nke1(iproc), nke1_str)

               !call MPI_BARRIER(MPI_COMM_WORLD,ierr)
               write(10,'(A)')'    <Piece Extent="'//trim(nis1_str)//" "//trim(nie1_str)// &
               & " "//trim(njs1_str)//" "//trim(nje1_str)//" "//trim(nks1_str)//" "//trim(nke1_str)// &
               &'"'//' Source="'//trim(iproc_str)//'.vts"/>'      
          enddo

          write(10,'(A)')'  </PStructuredGrid>'
          write(10,'(A)')'</VTKFile>'

    close(10)
!endif  !!! end of if(myid==0)then

   !call MPI_BARRIER(MPI_COMM_WORLD,ierr)

filename_vts=trim(dirname)//'/'//trim(rankid_str)//'.vts'

open(11,file=filename_vts)


         write(11,'(A)')'<?xml version="1.0"?>'
         write(11,'(A)')'<VTKFile type="StructuredGrid" version="0.1" byte_order="LittleEndian">'
         write(11,'(A)')'  <StructuredGrid WholeExtent="'//trim(nis_str)//" "//trim(nie_str)//" " &
                           & //trim(njs_str)//" "//trim(nje_str)//" "//trim(nks_str)//" "//trim(nke_str)//'">'
         write(11,'(A)')'    <Piece Extent="'//trim(nis_str)//" "//trim(nie_str)//" " &
                           & //trim(njs_str)//" "//trim(nje_str)//" "//trim(nks_str)//" "//trim(nke_str)//'">'

         write(11,'(A)')'      <CellData>'

         if(variable(1) .eq. 1) then
         write(11,'(A)')'        <DataArray type="Float32"  Name="Density" format="ascii">'
         do k=1, Nk
            do j=1, Nj
               do i=1,Ni
                write(11,*) q(i,j,k,1)
               enddo
            enddo
         enddo
         write(11,'(A)')'        </DataArray>'
        endif !! end of if(variable(1) .eq. 1) then

         if(variable(2) .eq. 1) then
         write(11,'(A)')'        <DataArray type="Float32"  Name="Vx" format="ascii">'
         do k=1,Nk
            do j=1, Nj
               do i=1,Ni
                write(11,*) q(i,j,k,2)/q(i,j,k,1)
               enddo
            enddo
         enddo
         write(11,'(A)')'        </DataArray>'
        endif !! end of if(variable(2) .eq. 1) then


         if(variable(3) .eq. 1) then
         write(11,'(A)')'        <DataArray type="Float32"  Name="Vy" format="ascii">'
         do k=1,Nk
            do j=1, Nj
               do i=1,Ni
                write(11,*) q(i,j,k,3)/q(i,j,k,1)
               enddo
            enddo
         enddo
         write(11,'(A)')'        </DataArray>'
        endif !! end of if(variable(3) .eq. 1) then

         if(variable(4) .eq. 1) then
         write(11,'(A)')'        <DataArray type="Float32"  Name="Vz" format="ascii">'
         do k=1,Nk
            do j=1, Nj
               do i=1,Ni
                write(11,*) q(i,j,k,4)/q(i,j,k,1)
               enddo
            enddo
         enddo
         write(11,'(A)')'        </DataArray>'
        endif !! end of if(variable(4) .eq. 1) then

         !if(variable(5) .eq. 1) then
         !write(11,'(A)')'        <DataArray type="Float32"  Name="Bx" format="ascii">'
          !do k=1,Nk
            !do j=1, Nj
               !do i=1,Ni
                !write(11,*) q(i,j,k,5)
               !enddo
            !enddo
          !enddo
         !write(11,'(A)')'        </DataArray>'
        !endif !! end of if(variable(5) .eq. 1) then

        !if(variable(6) .eq. 1) then
         !write(11,'(A)')'        <DataArray type="Float32"  Name="By" format="ascii">'
          !do k=1,Nk
            !do j=1, Nj
               !do i=1,Ni
                !write(11,*) q(i,j,k,6)
               !enddo
            !enddo
         !enddo
         !write(11,'(A)')'        </DataArray>'
        !endif !! end of if(variable(6) .eq. 1) then

        !if(variable(7) .eq. 1) then
         !write(11,'(A)')'        <DataArray type="Float32"  Name="Bz" format="ascii">'
          !do k=1,Nk
            !do j=1, Nj
               !do i=1,Ni
                !write(11,*) q(i,j,k,7)
               !enddo
            !enddo
          !enddo
         !write(11,'(A)')'        </DataArray>'
        !endif !! end of if(variable(7) .eq. 1) then



        if(variable(8) .eq. 1) then
         !write(11,'(A)')'        <DataArray type="Float32"  Name="Pressure" format="ascii">'
          write(11,'(A)')'        <DataArray type="Float32"  Name="Energy" format="ascii">'
          do k=1,Nk
            do j=1, Nj
               do i=1,Ni

                 write(11,*) q(i,j,k,5)
                   !if(this%eosType=2) then
                   !write(11,*) q(i,j,5)
                   !q(i,j,8)=p/(gam-1.d0)+0.5d0*rho*(vx**2.d0+vy**2.d0+vz**2.d0)+0.5d0*(bxc**2.d0+byc**2.d0+bzc**2.d0)
                   !else
                   !endif
                
               enddo
            enddo
         enddo
         write(11,'(A)')'        </DataArray>'
        endif !! end of if(variable(8) .eq. 1) then





         write(11,'(A)')'      </CellData>'



         write(11,'(A)')'      <Points>'

         !write(*,*)"myid= ",myid,Ni, Nj,Nk

         write(11,'(A)')'        <DataArray type="Float32" NumberOfComponents="3" format="ascii">'
          do k=1,Nk+1
            do j=1, Nj+1
               do i=1,Ni+1
                write(11,*) xl(i), yl(j),zl(k)
               enddo
            enddo
          enddo
         write(11,'(A)')'        </DataArray>'

         write(11,'(A)')'      </Points>'

         write(11,'(A)')'    </Piece>'
         write(11,'(A)')'  </StructuredGrid>'
         write(11,'(A)')'</VTKFile>'




close(11)

!call MPI_BARRIER(MPI_COMM_WORLD,ierr)


deallocate(nis1,nie1)
deallocate(njs1,nje1)
deallocate(nks1,nke1)

end subroutine output3d_HD_pvts








subroutine input3d_MHD_pvts(this,q)
use gridModule
implicit none
class(grid)::this
double precision,dimension(1-this%nbuf:this%nMesh(1)+this%nbuf,1-this%nbuf:this%nMesh(2)+this%nbuf &
                         &,1-this%nbuf:this%nMesh(3)+this%nbuf,this%nvar)::q 

!double precision,dimension(1-this%nbuf:this%nMesh(1)+this%nbuf)::xl
!double precision,dimension(1-this%nbuf:this%nMesh(2)+this%nbuf)::yl
!double precision,dimension(1-this%nbuf:this%nMesh(3)+this%nbuf)::zl




integer,dimension(1:this%ndim)::dims  

integer:: ndim

integer:: gnis,gnie,gnjs,gnje,gnks,gnke
character(len=50)::gnis_str,gnie_str,gnjs_str,gnje_str ,gnks_str,gnke_str


integer:: nis,nie,njs,nje,nks,nke
character(len=50)::nis_str,nie_str,njs_str,nje_str,nks_str,nke_str

integer:: nis1_send,nie1_send,njs1_send,nje1_send,nks1_send,nke1_send
integer:: nis1_rec, nie1_rec, njs1_rec, nje1_rec, nks1_rec,nke1_rec


integer,dimension(:),allocatable:: nis1,nie1,njs1,nje1,nks1,nke1
character(len=50)::nis1_str,nie1_str,njs1_str,nje1_str,nks1_str,nke1_str



integer:: i,j, k
integer::Ni,Nj,Nk

integer::iproc,Totalprocs

character(len=50)::iproc_str

integer::ierr
integer:: variable(8)

logical:: fileDirhave

integer:: fileid,rankid,gridid
character(len=50)::fileid_str,rankid_str,gridid_str

character(len=50):: dirname,filename_pvts,filename_vts


!xl=this%xl(1)%coords
!yl=this%xl(2)%coords
!zl=this%xl(3)%coords

variable=this%variable
ndim=this%ndim
dims=this%dims_mpi

Totalprocs=nprocs

!write(*,*)"ndim=",ndim, "dims=",dims, "TotalProcess",nprocs,"nBlockx=",nBlockx,"nBlocky=",nBlocky
!write(*,*)"this%nMesh(1)=",this%nMesh(1),"this%nMesh(2)",this%nMesh(2),"this%nMesh(3)",this%nMesh(3)
!write(*,*)"myid= ",myid, "this%gnx_list_R=", this%gnx_list_R, "this%gny_list_R=", this%gny_list_R, & 
!                       & "this%gnz_list_R=", this%gnz_list_R

!stop

nie=this%gnx_list_R
nje=this%gny_list_R
nke=this%gnz_list_R

nis=this%gnx_list_R-this%nMesh(1)
njs=this%gny_list_R-this%nMesh(2)
nks=this%gnz_list_R-this%nMesh(3)

call i4_to_s_left (nis, nis_str )
call i4_to_s_left (njs, njs_str )
call i4_to_s_left (nks, nks_str )

call i4_to_s_left (nie, nie_str )
call i4_to_s_left (nje, nje_str )
call i4_to_s_left (nke, nke_str )


allocate(nis1(0:Totalprocs-1),nie1(0:Totalprocs-1))
allocate(njs1(0:Totalprocs-1),nje1(0:Totalprocs-1))
allocate(nks1(0:Totalprocs-1),nke1(0:Totalprocs-1))


nis1=0
nie1=0

njs1=0
nje1=0

nks1=0
nke1=0



!call MPI_GATHER(nis,1,MPI_INT,nis1,1,MPI_INT,0,MPI_COMM_WORLD,ierr)
!call MPI_GATHER(nie,1,MPI_INT,nie1,1,MPI_INT,0,MPI_COMM_WORLD,ierr)


!call MPI_GATHER(njs,1,MPI_INT,njs1,1,MPI_INT,0,MPI_COMM_WORLD,ierr)
!call MPI_GATHER(nje,1,MPI_INT,nje1,1,MPI_INT,0,MPI_COMM_WORLD,ierr)

!call MPI_GATHER(nks,1,MPI_INT,nks1,1,MPI_INT,0,MPI_COMM_WORLD,ierr)
!call MPI_GATHER(nke,1,MPI_INT,nke1,1,MPI_INT,0,MPI_COMM_WORLD,ierr)


call MPI_ALLGATHER(nis,1,MPI_INT,nis1,1,MPI_INT,MPI_COMM_WORLD,ierr)
call MPI_ALLGATHER(nie,1,MPI_INT,nie1,1,MPI_INT,MPI_COMM_WORLD,ierr)


call MPI_ALLGATHER(njs,1,MPI_INT,njs1,1,MPI_INT,MPI_COMM_WORLD,ierr)
call MPI_ALLGATHER(nje,1,MPI_INT,nje1,1,MPI_INT,MPI_COMM_WORLD,ierr)

call MPI_ALLGATHER(nks,1,MPI_INT,nks1,1,MPI_INT,MPI_COMM_WORLD,ierr)
call MPI_ALLGATHER(nke,1,MPI_INT,nke1,1,MPI_INT,MPI_COMM_WORLD,ierr)

!if(myid==0)then
!write(*,*)nis1
!write(*,*)nie1
!write(*,*)njs1
!write(*,*)nje1
!write(*,*)nks1
!write(*,*)nke1
!endif
!stop



!MPI_GATHER(SENDBUF, SENDCOUNT, SENDTYPE, RECVBUF, RECVCOUNT, RECVTYPE, ROOT, COMM, IERROR)
        !<type> SENDBUF(*), RECVBUF(*)
        !INTEGER SENDCOUNT, SENDTYPE, RECVCOUNT, RECVTYPE, ROOT, COMM, IERROR
!call MPI_BARRIER(MPI_COMM_WORLD,ierr)


gnis=0
gnjs=0
gnks=0

gnie=this%nMesh_global(1)
gnje=this%nMesh_global(2)
gnke=this%nMesh_global(3)


call i4_to_s_left (gnis, gnis_str )
call i4_to_s_left (gnjs, gnjs_str )
call i4_to_s_left (gnks, gnks_str )

call i4_to_s_left (gnie, gnie_str )
call i4_to_s_left (gnje, gnje_str )
call i4_to_s_left (gnke, gnke_str )


Ni=this%nMesh(1)
Nj=this%nMesh(2)
Nk=this%nMesh(3)


!write(*,*)"myid= ",myid,"nis=",nis,"njs=",njs,"nks=",nks,"nie=",nie,"nje=",nje,"nke=",nke
!write(*,*)"myid= ",myid,"gnie=",gnie,"gnje=",gnje,"gnke=",gnke




gridid=this%gridID
fileid=this%fnum
rankid=myid

call i4_to_s_left (gridid, gridid_str )
call i4_to_s_left (fileid, fileid_str )
call i4_to_s_left (rankid, rankid_str )


dirname="MHD_g"//trim(gridid_str)//"_"//trim(fileid_str)

inquire(file=trim(dirname),exist=fileDirhave)

if(.not.fileDirhave)then
   write(*,*)"The inputfile dose not exist  !!!!!"
endif  ! end of if(.not.fileDirhave)then


filename_vts=trim(dirname)//'/'//trim(rankid_str)//'.vts'

open(11,file=filename_vts)


         read(11,*)
         read(11,*)
         read(11,*)             
         read(11,*)

         read(11,*)

         read(11,*)
         do k=1, Nk
            do j=1, Nj
               do i=1,Ni
                read(11,*) q(i,j,k,1)
               enddo
            enddo
         enddo
         read(11,*)
       


         read(11,*)
         do k=1,Nk
            do j=1, Nj
               do i=1,Ni
                read(11,*) q(i,j,k,2)
                q(i,j,k,2)=q(i,j,k,1)*q(i,j,k,2)
               enddo
            enddo
         enddo
         read(11,*)
        

         
         read(11,*)
         do k=1,Nk
            do j=1, Nj
               do i=1,Ni
                read(11,*) q(i,j,k,3)
                q(i,j,k,3)=q(i,j,k,1)*q(i,j,k,3)
               enddo
            enddo
         enddo
         read(11,*)
        

         
         read(11,*)
         do k=1,Nk
            do j=1, Nj
               do i=1,Ni
                read(11,*) q(i,j,k,4)
                q(i,j,k,4)=q(i,j,k,1)*q(i,j,k,4)
               enddo
            enddo
         enddo
         read(11,*)
       

        
        read(11,*)
          do k=1,Nk
            do j=1, Nj
               do i=1,Ni
                read(11,*) q(i,j,k,5)
               enddo
            enddo
          enddo
         read(11,*)
        

        
        read(11,*)
          do k=1,Nk
            do j=1, Nj
               do i=1,Ni
               read(11,*) q(i,j,k,6)
               enddo
            enddo
         enddo
         read(11,*)
       

        
         read(11,*)
          do k=1,Nk
            do j=1, Nj
               do i=1,Ni
                read(11,*) q(i,j,k,7)
               enddo
            enddo
          enddo
         read(11,*)
        


        
        
          read(11,*)
          do k=1,Nk
            do j=1, Nj
               do i=1,Ni

                 read(11,*) q(i,j,k,8)

                   !if(this%eosType=2) then
                   !write(11,*) q(i,j,8)
                   !q(i,j,8)=p/(gam-1.d0)+0.5d0*rho*(vx**2.d0+vy**2.d0+vz**2.d0)+0.5d0*(bxc**2.d0+byc**2.d0+bzc**2.d0)
                   !else
                   !endif
                
               enddo
            enddo
         enddo
         read(11,*)


         
         read(11,*)
          do k=1,Nk
            do j=1, Nj
               do i=1,Ni
                read(11,*) q(i,j,k,9)
               enddo
            enddo
          enddo
         read(11,*)
        

       
         read(11,*)
          do k=1,Nk
            do j=1, Nj
               do i=1,Ni
                read(11,*) q(i,j,k,10)
               enddo
            enddo
         enddo
         read(11,*)
        

         read(11,*)
          do k=1,Nk
            do j=1, Nj
               do i=1,Ni
               read(11,*) q(i,j,k,11)
               enddo
            enddo
          enddo
        read(11,*)


        read(11,*)


close(11)

call MPI_BARRIER(MPI_COMM_WORLD,ierr)

deallocate(nis1,nie1)
deallocate(njs1,nje1)
deallocate(nks1,nke1)



end subroutine input3d_MHD_pvts



subroutine input3d_HD_pvts(this,q)
use gridModule
implicit none
class(grid)::this
double precision,dimension(1-this%nbuf:this%nMesh(1)+this%nbuf,1-this%nbuf:this%nMesh(2)+this%nbuf &
                         &,1-this%nbuf:this%nMesh(3)+this%nbuf,this%nvar)::q 

!double precision,dimension(1-this%nbuf:this%nMesh(1)+this%nbuf)::xl
!double precision,dimension(1-this%nbuf:this%nMesh(2)+this%nbuf)::yl
!double precision,dimension(1-this%nbuf:this%nMesh(3)+this%nbuf)::zl




integer,dimension(1:this%ndim)::dims  

integer:: ndim

integer:: gnis,gnie,gnjs,gnje,gnks,gnke
character(len=50)::gnis_str,gnie_str,gnjs_str,gnje_str ,gnks_str,gnke_str


integer:: nis,nie,njs,nje,nks,nke
character(len=50)::nis_str,nie_str,njs_str,nje_str,nks_str,nke_str

integer:: nis1_send,nie1_send,njs1_send,nje1_send,nks1_send,nke1_send
integer:: nis1_rec, nie1_rec, njs1_rec, nje1_rec, nks1_rec,nke1_rec


integer,dimension(:),allocatable:: nis1,nie1,njs1,nje1,nks1,nke1
character(len=50)::nis1_str,nie1_str,njs1_str,nje1_str,nks1_str,nke1_str



integer:: i,j, k
integer::Ni,Nj,Nk

integer::iproc,Totalprocs

character(len=50)::iproc_str

integer::ierr
integer:: variable(8)

logical:: fileDirhave

integer:: fileid,rankid,gridid
character(len=50)::fileid_str,rankid_str,gridid_str

character(len=50):: dirname,filename_pvts,filename_vts


!xl=this%xl(1)%coords
!yl=this%xl(2)%coords
!zl=this%xl(3)%coords

variable=this%variable
ndim=this%ndim
dims=this%dims_mpi

Totalprocs=nprocs

!write(*,*)"ndim=",ndim, "dims=",dims, "TotalProcess",nprocs,"nBlockx=",nBlockx,"nBlocky=",nBlocky
!write(*,*)"this%nMesh(1)=",this%nMesh(1),"this%nMesh(2)",this%nMesh(2),"this%nMesh(3)",this%nMesh(3)
!write(*,*)"myid= ",myid, "this%gnx_list_R=", this%gnx_list_R, "this%gny_list_R=", this%gny_list_R, & 
!                       & "this%gnz_list_R=", this%gnz_list_R

!stop

nie=this%gnx_list_R
nje=this%gny_list_R
nke=this%gnz_list_R

nis=this%gnx_list_R-this%nMesh(1)
njs=this%gny_list_R-this%nMesh(2)
nks=this%gnz_list_R-this%nMesh(3)

call i4_to_s_left (nis, nis_str )
call i4_to_s_left (njs, njs_str )
call i4_to_s_left (nks, nks_str )

call i4_to_s_left (nie, nie_str )
call i4_to_s_left (nje, nje_str )
call i4_to_s_left (nke, nke_str )


allocate(nis1(0:Totalprocs-1),nie1(0:Totalprocs-1))
allocate(njs1(0:Totalprocs-1),nje1(0:Totalprocs-1))
allocate(nks1(0:Totalprocs-1),nke1(0:Totalprocs-1))


nis1=0
nie1=0

njs1=0
nje1=0

nks1=0
nke1=0



!call MPI_GATHER(nis,1,MPI_INT,nis1,1,MPI_INT,0,MPI_COMM_WORLD,ierr)
!call MPI_GATHER(nie,1,MPI_INT,nie1,1,MPI_INT,0,MPI_COMM_WORLD,ierr)


!call MPI_GATHER(njs,1,MPI_INT,njs1,1,MPI_INT,0,MPI_COMM_WORLD,ierr)
!call MPI_GATHER(nje,1,MPI_INT,nje1,1,MPI_INT,0,MPI_COMM_WORLD,ierr)

!call MPI_GATHER(nks,1,MPI_INT,nks1,1,MPI_INT,0,MPI_COMM_WORLD,ierr)
!call MPI_GATHER(nke,1,MPI_INT,nke1,1,MPI_INT,0,MPI_COMM_WORLD,ierr)

call MPI_ALLGATHER(nis,1,MPI_INT,nis1,1,MPI_INT,MPI_COMM_WORLD,ierr)
call MPI_ALLGATHER(nie,1,MPI_INT,nie1,1,MPI_INT,MPI_COMM_WORLD,ierr)


call MPI_ALLGATHER(njs,1,MPI_INT,njs1,1,MPI_INT,MPI_COMM_WORLD,ierr)
call MPI_ALLGATHER(nje,1,MPI_INT,nje1,1,MPI_INT,MPI_COMM_WORLD,ierr)

call MPI_ALLGATHER(nks,1,MPI_INT,nks1,1,MPI_INT,MPI_COMM_WORLD,ierr)
call MPI_ALLGATHER(nke,1,MPI_INT,nke1,1,MPI_INT,MPI_COMM_WORLD,ierr)

!if(myid==0)then
!write(*,*)nis1
!write(*,*)nie1
!write(*,*)njs1
!write(*,*)nje1
!write(*,*)nks1
!write(*,*)nke1
!endif
!stop



!MPI_GATHER(SENDBUF, SENDCOUNT, SENDTYPE, RECVBUF, RECVCOUNT, RECVTYPE, ROOT, COMM, IERROR)
        !<type> SENDBUF(*), RECVBUF(*)
        !INTEGER SENDCOUNT, SENDTYPE, RECVCOUNT, RECVTYPE, ROOT, COMM, IERROR
!call MPI_BARRIER(MPI_COMM_WORLD,ierr)


gnis=0
gnjs=0
gnks=0

gnie=this%nMesh_global(1)
gnje=this%nMesh_global(2)
gnke=this%nMesh_global(3)


call i4_to_s_left (gnis, gnis_str )
call i4_to_s_left (gnjs, gnjs_str )
call i4_to_s_left (gnks, gnks_str )

call i4_to_s_left (gnie, gnie_str )
call i4_to_s_left (gnje, gnje_str )
call i4_to_s_left (gnke, gnke_str )


Ni=this%nMesh(1)
Nj=this%nMesh(2)
Nk=this%nMesh(3)


!write(*,*)"myid= ",myid,"nis=",nis,"njs=",njs,"nks=",nks,"nie=",nie,"nje=",nje,"nke=",nke
!write(*,*)"myid= ",myid,"gnie=",gnie,"gnje=",gnje,"gnke=",gnke




gridid=this%gridID
fileid=this%fnum
rankid=myid

call i4_to_s_left (gridid, gridid_str )
call i4_to_s_left (fileid, fileid_str )
call i4_to_s_left (rankid, rankid_str )


dirname="HD_g"//trim(gridid_str)//"_"//trim(fileid_str)

inquire(file=trim(dirname),exist=fileDirhave)

if(.not.fileDirhave)then
   write(*,*)"The inputfile dose not exist  !!!!!"
endif  ! end of if(.not.fileDirhave)then


filename_vts=trim(dirname)//'/'//trim(rankid_str)//'.vts'

open(11,file=filename_vts)


         read(11,*)
         read(11,*)
         read(11,*)              
         read(11,*)

         read(11,*)

         read(11,*)
         do k=1, Nk
            do j=1, Nj
               do i=1,Ni
                read(11,*) q(i,j,k,1)
               enddo
            enddo
         enddo
        read(11,*)
       

        
         read(11,*)
         do k=1,Nk
            do j=1, Nj
               do i=1,Ni
                read(11,*) q(i,j,k,2)
                 q(i,j,k,2)=q(i,j,k,1)*q(i,j,k,2)
               enddo
            enddo
         enddo
         read(11,*)
      


         
         read(11,*)
         do k=1,Nk
            do j=1, Nj
               do i=1,Ni
                read(11,*) q(i,j,k,3)
                 q(i,j,k,3)=q(i,j,k,1)*q(i,j,k,3)
               enddo
            enddo
         enddo
        read(11,*)
        

       
         read(11,*)
         do k=1,Nk
            do j=1, Nj
               do i=1,Ni
                read(11,*) q(i,j,k,4)
                q(i,j,k,4)=q(i,j,k,1)*q(i,j,k,4)
               enddo
            enddo
         enddo
         read(11,*)
       


       
         
          read(11,*)
          do k=1,Nk
            do j=1, Nj
               do i=1,Ni

                read(11,*) q(i,j,k,5)
                   !if(this%eosType=2) then
                   !read(11,*) q(i,j,5)
                   !q(i,j,8)=p/(gam-1.d0)+0.5d0*rho*(vx**2.d0+vy**2.d0+vz**2.d0)+0.5d0*(bxc**2.d0+byc**2.d0+bzc**2.d0)
                   !else
                   !endif
                
               enddo
            enddo
         enddo
         read(11,*)

         read(11,*)

close(11)

!call MPI_BARRIER(MPI_COMM_WORLD,ierr)


deallocate(nis1,nie1)
deallocate(njs1,nje1)
deallocate(nks1,nke1)

end subroutine input3d_HD_pvts





!*****************************************************************************

        subroutine digit_to_ch ( digit, ch )

        implicit none
        character ch
        integer ( kind = 4 ) digit

        if ( 0 <= digit .and. digit <= 9 ) then
          ch = achar ( digit + 48 )
        else
         ch = '*'
         end if

        return
        end

!*****************************************************************************
       subroutine i4_to_s_left ( i4, s )

        implicit none

        character c
        integer ( kind = 4 ) i
        integer ( kind = 4 ) i4
        integer ( kind = 4 ) idig
        integer ( kind = 4 ) ihi
        integer ( kind = 4 ) ilo
        integer ( kind = 4 ) ipos
        integer ( kind = 4 ) ival
        character ( len = * ) s

         s = ' '
         ilo = 1
         ihi = len ( s )
         if ( ihi <= 0 ) then
          return
         end if

         ival = i4

          if ( ival < 0 ) then
           if ( ihi <= 1 ) then
        s(1:1) = '*'
        return
           end if

         ival = -ival
          s(1:1) = '-'
          ilo = 2

         end if

         ipos = ihi

        do

         idig = mod ( ival, 10 )
         ival = ival / 10

           if ( ipos < ilo ) then
       do i = 1, ihi
        s(i:i) = '*'
       end do
       return
        end if

         call digit_to_ch ( idig, c )

         s(ipos:ipos) = c
         ipos = ipos - 1

        if ( ival == 0 ) then
       exit
        end if

         end do

         s(ilo:ilo+ihi-ipos-1) = s(ipos+1:ihi)
         s(ilo+ihi-ipos:ihi) = ' '

         return
        end

!*****************************************************************************
