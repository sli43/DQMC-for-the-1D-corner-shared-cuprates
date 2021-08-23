module mpi_module
 include 'mpif.h'
 integer nproc, myrank
 integer, parameter :: root = 0
 integer ierr
 integer mpi_ibuff
 double precision, allocatable, dimension(:) :: bin_sgn
 double precision, allocatable, dimension(:) :: bin_one
 double precision, allocatable, dimension(:) :: bin_sgnp
 double precision, allocatable, dimension(:) :: bin_sgnt
 double precision, allocatable, dimension(:) :: bin_buff
 double precision, allocatable, dimension(:) :: bin_buff2 
 double precision, allocatable, dimension(:) :: bin_buff3 
 double complex, allocatable, dimension(:) :: bin_cbuff
 double complex, allocatable, dimension(:) :: bin_cbuff2
contains
 !=============================================================================
 subroutine init_mpi()
 implicit none

 !init mpi_comm_world 
 call mpi_init(ierr)
 if(ierr.ne.0)then
  print*, 'Error calling mpi_init.'
  print*, 'Exited with code: ', ierr
  stop
 endif

 call mpi_comm_rank(mpi_comm_world,myrank,ierr)
 if(ierr.ne.0)then
  print*, 'Error calling mpi_comm_rank.'
  print*, 'Exited with code: ', ierr
  stop
 endif

 call mpi_comm_size(mpi_comm_world,nproc,ierr)
 if(ierr.ne.0)then
  print*, 'Error calling mpi_comm_size.'
  print*, 'Exited with code: ', ierr
  stop
 endif

 if(myrank.eq.root)then
  print*, 'Mpi initialized with nproc = ', nproc, '.'
 endif

 return
 end subroutine init_mpi
 !=============================================================================
 subroutine finalize_mpi()
 implicit none
 call mpi_finalize(ierr)
 if(ierr.ne.0)then
  print*, 'Error calling mpi_finalize.'
  print*, 'Exited with code: ', ierr
  stop
 endif
 return
 end subroutine finalize_mpi
 !=============================================================================
end module mpi_module
