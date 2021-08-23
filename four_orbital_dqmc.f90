include 'random.f90'
include 'linear_algebra.f90'
include 'mpi_module.f90'
!include 'parameters.f90'
include 'fields.f90'
include 'grfun.f90'
include 'measurements.f90'
include 'updates.f90'
include 'io_routines.f90'
include 'inter_file.f90'
!==============================================================================
! The main program is defined here.
!==============================================================================
program main
use parameters
use fields
use random
use mpi_module
use grfun
use measurements
use updates
use io_routines
use inter_file
implicit none
logical warms
integer rand
!integer iseed
integer i, j, ti, bin
integer iwarmbegin,imeasbegin
integer*4 timeArray(3)
double precision sgnup, sgndn, detup, detdn
double precision, allocatable, dimension(:,:) :: Gup, gdn
character data*16,time*16,zone*16,filetime*100
integer,dimension(1:8)::time_value1,time_value2

call init_mpi()
call init_cluster()
call build_kin()

allocate(gup(0:N-1,0:N-1))
allocate(gdn(0:N-1,0:N-1))

iran = -(iseed+myrank)

!init the system
call init_fields()
ti = 0
!if(myrank.eq.root) print*, 'Building spin up Green''s function.'
call getgp(gup,ti,detup,sgnup,spinup)
!if(myrank.eq.root) print*, 'Building spin dn Green''s function.'
call getgp(gdn,ti,detdn,sgndn,spindn)

!==============================================================================
! Step 1. Do the warmup sweeps
!==============================================================================

call zero_counters()
sweeps_per_et_meas = 0
calls_to_sweep_hs = 0
wraps = nwraps

iwarmbegin=1
imeasbegin=1
bin=0

 call inter_readfile(iwarmbegin,imeasbegin,sgnup,sgndn,gup,gdn,detup,detdn,bin)

warms = .true.
do i = iwarmbegin,nwarms
 if(myrank==0.and.mod(i,100).eq.0) print*, i

 if(mod(i,nwarms/10).eq.0)then
  call output_ratios(i,nwarms,warms)
 endif
 !do the woodbury updates
 call single_site_updates(gup,gdn,sgnup,sgndn,warms)
 calls_to_sweep_hs = calls_to_sweep_hs + 1

 if(calls_to_sweep_hs.eq.sweeps_per_block_hs)then
  calls_to_sweep_hs = 0
  if(sweeps_per_block_hs.ne.0)then
   call block_update_hs(gup,gdn,sgnup,sgndn)
  endif
 endif


if(mod(i,stop_numb).eq.0)then
   iwarmbegin=i+1
   call inter_writefile(iwarmbegin,imeasbegin,sgnup,sgndn,gup,gdn,detup,detdn,bin)
   goto 1745
end if
enddo

!==============================================================================
! Step 2. Do the measurement sweeps
!==============================================================================

if(imeasbegin.eq.1) then
  call zero_counters()
  call zero_accumulators()
  bin = 0 !counts the bin number
  sweeps_per_et_meas = 0
end if

warms = .false.
do i = imeasbegin,nsweep
 if(mod(i,nsweep/10).eq.0)then
  call output_ratios(i,nsweep,warms)
  write(*,*) "nwraps = ",nwraps
 endif

 !do the woodbury updates
 call single_site_updates(gup,gdn,sgnup,sgndn,warms)
 calls_to_sweep_hs = calls_to_sweep_hs + 1
 
 !do block updates if needed
 if(calls_to_sweep_hs.eq.sweeps_per_block_hs)then
  calls_to_sweep_hs = 0
  if(sweeps_per_block_hs.ne.0)then
   call block_update_hs(gup,gdn,sgnup,sgndn)
  endif
 endif


 !should we do some measurements?
 if(tauskip.gt.0.and.mod(i,tauskip).eq.0)then
   call meastau(sgnup,sgndn)
 endif

 if(mod(i,nsweep/Nbins).eq.0)then
  bin = bin + 1
  call populate_bins(bin)
  call zero_accumulators()
 endif

if(mod(i+nwarms,stop_numb).eq.0)then
   iwarmbegin=nwarms+1
   imeasbegin=i+1
   call inter_writefile(iwarmbegin,imeasbegin,sgnup,sgndn,gup,gdn,detup,detdn,bin)
   if(i.ne.nsweep) then
       goto 1745
   end if
end if 
enddo

!==============================================================================
! Output the results
!==============================================================================
call output_results()
!==============================================================================
!End the program and exit gracefully
!==============================================================================
1745 continue
call mpi_barrier(mpi_comm_world,ierr)
call deallocate_fields()
call finalize_mpi()
!call date_and_time(data,time,zone,time_value2)
!write(*,*) time_value2
! write(unit=outfile,fmt="('cputime',i6.6,'.dat')") run
! open(file=OUTFILE,unit=81,action='write')
!   write(81,*) "N = ",Nlat
!   write(81,*) "warms = ", nwarms
!   write(81,*) "measure = ",nsweep
!   write(81,*) "number of time slice = ",L
!   write(81,*) "cpu time = ", (time_value2(3)-time_value1(3))*24+(time_value2(5)-time_value1(5))+(time_value2(6)-time_value1(6))/60.0d0
!close(81)

end program main
