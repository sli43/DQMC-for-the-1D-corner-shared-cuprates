module io_routines
 integer, parameter :: stdio = 6
 contains
 subroutine output_ratios(iter,numiter,warms)
 use mpi_module
 use updates, only: accept_hs, reject_hs, & 
                    redo, noredo, zero_counters, & 
                    accept_blk_hs, reject_blk_hs

 use parameters, only: gam
 implicit none
 logical warms
 integer iter, numiter
 integer(KIND=8) accept_tot, reject_tot, redo_tot, noredo_tot
 integer(kind=8) accept_ph_tot, reject_ph_tot
 integer(kind=8) accept_ph_blk_tot, reject_ph_blk_tot
 double precision accept_ratio, reject_ratio, redo_ratio

 call mpi_reduce(accept_hs,accept_tot,1,mpi_integer8,mpi_sum,root,mpi_comm_world,ierr)
 call mpi_reduce(reject_hs,reject_tot,1,mpi_integer8,mpi_sum,root,mpi_comm_world,ierr)
 call mpi_reduce(redo,redo_tot,1,mpi_integer8,mpi_sum,root,mpi_comm_world,ierr)
 call mpi_reduce(noredo,noredo_tot,1,mpi_integer8,mpi_sum,root,mpi_comm_world,ierr)

 if(myrank.eq.root)then
  if(warms)then 
   write(unit=stdio,fmt=98) iter, numiter
  else
   write(unit=stdio,fmt=99) iter, numiter
  endif
  
  if(redo_tot+noredo_tot.ne.0)then
   redo_ratio = dfloat(redo_tot)/dfloat(redo_tot+noredo_tot)
  else
   redo_ratio = 0.0d0
  endif

  if(accept_tot+reject_tot.ne.0)then
   accept_ratio = dfloat(accept_tot)/dfloat(accept_tot+reject_tot)
   reject_ratio = dfloat(reject_tot)/dfloat(accept_tot+reject_tot)
  else
   accept_ratio = 0.0d0
   reject_ratio = 1.0d0
  endif

  print*, 'HS accept/reject ratio = ', accept_ratio, reject_ratio
  print*, 'Current value of gam = ', gam
  print*, 'Redo ratio = ', redo_ratio
 endif

 call mpi_reduce(accept_blk_hs,accept_tot,1,mpi_integer8,mpi_sum,root,mpi_comm_world,ierr)
 call mpi_reduce(reject_blk_hs,reject_tot,1,mpi_integer8,mpi_sum,root,mpi_comm_world,ierr)
 if(myrank.eq.root)then
  if(accept_tot+reject_tot.ne.0)then
   accept_ratio = dfloat(accept_tot)/dfloat(accept_tot+reject_tot)
   reject_ratio = dfloat(reject_tot)/dfloat(accept_tot+reject_tot)
  else
   accept_ratio = 0.0d0
   reject_ratio = 1.0d0
  endif
  print*, 'HS blk accept/reject ratio = ', accept_ratio, reject_ratio
 endif


 call zero_counters()
98 format(' Warm up sweep: ',i8,' of ', i8)
99 format(' Measurement sweep: ',i8,' of ', i8)

 if(myrank.eq.root) print*, ' '

 return
 end subroutine output_ratios

 !======================================================================================
 subroutine output_results()
 use parameters
 use measurements
 use mpi_module
 implicit none
 integer nnn, i, ti,j, orb1, orb2,ii,jj,info
 double precision bin_sgn_buff
 double precision mean, std
 double precision mean1, std1
 double precision mean2, std2
 double precision mean3, std3
 double precision mean4, std4
 double precision mean5, std5
 double precision mean6, std6
 double precision mean7, std7
 double precision mean8, std8
 double precision mean9, std9
 double precision, dimension(1:16) :: buf1, buf2, buf3, buf4
 double precision, dimension(0:Nlat-1,0:L,1:4,1:4)::bin_GF_buff
 double precision, dimension(1:4):: bin_GF_x
 double precision r, k, fk, ek, tau, pi
 double complex, parameter :: eye = (0.0d0,1.0d0)
 double complex, dimension(0:N/2,-L:L) :: gnl0
 double precision, dimension(1:18)::work
 integer, dimension(1:3):: ipiv
 integer, parameter :: nfwrk1 = 60
 integer, parameter :: nfwrk3=63
 integer, parameter :: nfwrk4=64
 integer, parameter :: nfwrk2=62
 integer, parameter :: itag1 = 98
 integer, parameter :: itag2 = 99
 integer, dimension(mpi_status_size) :: istat
 integer, dimension(2) :: req
 integer iproc,kproc 

 pi = 2.0d0*asin(1.0d0)
 if(myrank.eq.root)then
  open(file=OUTFILE,unit=nfwrk1,action='write')
  open(file=outfile3,unit=nfwrk3,action='write')
  open(file=outfile4,unit=nfwrk4,action='write')
  open(file=outfile2,unit=nfwrk2,action='write')

  write(unit=nfwrk1,fmt="(a50)") 'Input parameters:'
  write(unit=nfwrk1,fmt="(a50,i6)") 'nproc = ', nproc
  write(unit=nfwrk1,fmt="(a50,i6)") 'N = ', Nlat
  write(unit=nfwrk1,fmt="(a50,i6)") 'L = ', L
  write(unit=nfwrk1,fmt="(a50,f12.6)") 'Ud = ', Ud
  write(unit=nfwrk1,fmt="(a50,f12.6)") 'Up = ', Up
  write(unit=nfwrk1,fmt="(a50,f12.6)") 'Upd = ', Upd
  write(unit=nfwrk1,fmt="(a50,f12.6)") 'tpxd = ', tpxd
  write(unit=nfwrk1,fmt="(a50,f12.6)") 'tpyd = ', tpyd
  write(unit=nfwrk1,fmt="(a50,f12.6)") 'deltaCu = ', deltaCu
  write(unit=nfwrk1,fmt="(a50,f12.6)") 'deltapx = ', deltaPx
  write(unit=nfwrk1,fmt="(a50,f12.6)") 'deltapy = ', deltapy
  write(unit=nfwrk1,fmt="(a50,f12.6)") 'mu = ', mu
  write(unit=nfwrk1,fmt=*) ' '
  write(unit=nfwrk1,fmt="(a50,f12.6)") 'beta = ', beta
  write(unit=nfwrk1,fmt="(a50,f12.6)") 'delta tau = ', dtau
  write(unit=nfwrk1,fmt="(a50,i10)") 'number of warmup sweeps = ', nwarms
  write(unit=nfwrk1,fmt="(a50,i10)") 'number of measurement sweeps = ', nsweep
  write(unit=nfwrk1,fmt="(a50,i6)") 'orthlen = ', orthlen
  write(unit=nfwrk1,fmt="(a50,i6)") 'starting value of nwraps = ', nwraps
  write(unit=nfwrk1,fmt="(a50,i6)") 'final value of nwraps = ', nwrap
  write(unit=nfwrk1,fmt="(a50,f12.6)") 'difflim = ', difflim
  write(unit=nfwrk1,fmt="(a50,f12.6)") 'errrat = ', errrat
  write(unit=nfwrk1,fmt="(a50,i10)") 'number of bins = ', nbins
  write(unit=nfwrk1,fmt="(a50,i10)") 'number of processor = ', nproc
  write(unit=nfwrk1,fmt="(a50,i10)") 'Block_site_number = ',Block_site_number
  write(unit=nfwrk1,fmt="(a50,i10)") 'Block_spin_number = ',Block_spin_number
  write(unit=nfwrk1,fmt=*) ' '
  write(unit=nfwrk1,fmt="(a50)") 'Average signs'  
 endif

 nnn = nproc*nbins
 allocate(bin_sgn(1:nnn))
 allocate(bin_sgnt(1:nnn))
 allocate(bin_one(1:nnn))
 allocate(bin_buff(1:nnn))
 allocate(bin_buff2(1:nnn))
 allocate(bin_buff3(1:nnn))
 allocate(bin_cbuff(1:nnn))
 allocate(bin_cbuff2(1:nnn))

 !gather the one's array
 call mpi_gather(bone,nbins,MPI_DOUBLE_PRECISION,bin_one,nbins,MPI_DOUBLE_PRECISION,root,mpi_comm_world,ierr)

 !do the signs
 call mpi_gather(bsgnup,nbins,MPI_DOUBLE_PRECISION,bin_sgn,nbins,MPI_DOUBLE_PRECISION,root,mpi_comm_world,ierr)
 if(myrank.eq.root)then 
  call geterr(bin_sgn,bin_one,nnn,mean,std)
  write(unit=nfwrk1,fmt="(a50,f12.6,' +- ',f12.6)") 'Average up sign = ', mean, std
 endif

 call mpi_gather(bsgndn,nbins,MPI_DOUBLE_PRECISION,bin_sgn,nbins,MPI_DOUBLE_PRECISION,root,mpi_comm_world,ierr)
 if(myrank.eq.root)then 
  call geterr(bin_sgn,bin_one,nnn,mean,std)
  write(unit=nfwrk1,fmt="(a50,f12.6,' +- ',f12.6)") 'Average dn sign = ', mean, std
 endif
 
 call mpi_gather(bsgn,nbins,MPI_DOUBLE_PRECISION,bin_sgn,nbins,MPI_DOUBLE_PRECISION,root,mpi_comm_world,ierr)
 if(myrank.eq.root)then 
  call geterr(bin_sgn,bin_one,nnn,mean,std)
  write(unit=nfwrk1,fmt="(a50,f12.6,' +- ',f12.6)") 'Average sign = ', mean, std
 endif


 if(myrank.eq.root)then
  do i = 1,nnn
    if(bin_sgn(i).le.0.05d0)then
      print*, 'warning, average sign in bin ', i, ' is small.'
      print*, '<sign> = ', bin_sgn(i)
     endif
   enddo
 endif


 call mpi_gather(bsgnt,nbins,MPI_DOUBLE_PRECISION,bin_sgnt,nbins,MPI_DOUBLE_PRECISION,root,mpi_comm_world,ierr)
 if(myrank.eq.root)then 
  call geterr(bin_sgn,bin_one,nnn,mean,std)
  write(unit=nfwrk1,fmt="(a50,f12.6,' +- ',f12.6)") 'Average signt = ', mean, std
 endif

if(myrank.eq.root)then
  write(unit=nfwrk1,fmt=*) ' '
  write(unit=nfwrk1,fmt="(a50)") 'Energy results' 
 endif

call mpi_gather(bke,nbins,MPI_DOUBLE_PRECISION,bin_buff,nbins,MPI_DOUBLE_PRECISION,root,mpi_comm_world,ierr)
 if(myrank.eq.root)then 
  call geterr(bin_buff,bin_sgn,nnn,mean,std)
  write(unit=nfwrk1,fmt="(a50,f12.6,' +- ',f12.6)") 'Average kinetic energy = ', mean, std
 endif

 call mpi_gather(bpe1,nbins,MPI_DOUBLE_PRECISION,bin_buff,nbins,MPI_DOUBLE_PRECISION,root,mpi_comm_world,ierr)
 if(myrank.eq.root)then 
  call geterr(bin_buff,bin_sgn,nnn,mean,std)
  write(unit=nfwrk1,fmt="(a50,f12.6,' +- ',f12.6)") 'Average potential energy1 = ', mean, std
 endif

 call mpi_gather(bpe2,nbins,MPI_DOUBLE_PRECISION,bin_buff,nbins,MPI_DOUBLE_PRECISION,root,mpi_comm_world,ierr)
 if(myrank.eq.root)then 
  call geterr(bin_buff,bin_sgn,nnn,mean,std)
  write(unit=nfwrk1,fmt="(a50,f12.6,' +- ',f12.6)") 'Average potential energy2 = ', mean, std
 endif

 call mpi_gather(bpe3,nbins,MPI_DOUBLE_PRECISION,bin_buff,nbins,MPI_DOUBLE_PRECISION,root,mpi_comm_world,ierr)
 if(myrank.eq.root)then 
  call geterr(bin_buff,bin_sgn,nnn,mean,std)
  write(unit=nfwrk1,fmt="(a50,f12.6,' +- ',f12.6)") 'Average potential energy3 = ', mean, std
 endif

 call mpi_gather(bee,nbins,MPI_DOUBLE_PRECISION,bin_buff,nbins,MPI_DOUBLE_PRECISION,root,mpi_comm_world,ierr)
 if(myrank.eq.root)then 
  call geterr(bin_buff,bin_sgn,nnn,mean,std)
  write(unit=nfwrk1,fmt="(a50,f12.6,' +- ',f12.6)") 'Average total energy = ', mean, std
 endif


 if(myrank.eq.root)then
  write(unit=nfwrk1,fmt=*) ' '
  write(unit=nfwrk1,fmt="(a50)") 'Filling results' 
 endif
 call mpi_gather(bn,nbins,MPI_DOUBLE_PRECISION,bin_buff,nbins,MPI_DOUBLE_PRECISION,root,mpi_comm_world,ierr)
 if(myrank.eq.root)then 
  call geterr(bin_buff,bin_sgn,nnn,mean,std)
  write(unit=nfwrk1,fmt="(a50,f12.6,' +- ',f12.6)") 'Average <n> = ', mean, std
 endif

 call mpi_gather(bnup,nbins,MPI_DOUBLE_PRECISION,bin_buff,nbins,MPI_DOUBLE_PRECISION,root,mpi_comm_world,ierr)
 if(myrank.eq.root)then
  call geterr(bin_buff,bin_sgn,nnn,mean,std)
  write(unit=nfwrk1,fmt="(a50,f12.6,' +- ',f12.6)") 'Average <nup> = ', mean, std
 endif

 call mpi_gather(bndn,nbins,MPI_DOUBLE_PRECISION,bin_buff,nbins,MPI_DOUBLE_PRECISION,root,mpi_comm_world,ierr)
 if(myrank.eq.root)then
  call geterr(bin_buff,bin_sgn,nnn,mean,std)
  write(unit=nfwrk1,fmt="(a50,f12.6,' +- ',f12.6)") 'Average <ndn> = ', mean, std
 endif

 call mpi_gather(bn1,nbins,MPI_DOUBLE_PRECISION,bin_buff,nbins,MPI_DOUBLE_PRECISION,root,mpi_comm_world,ierr)
 if(myrank.eq.root)then 
  call geterr(bin_buff,bin_sgn,nnn,mean,std)
  write(unit=nfwrk1,fmt="(a50,f12.6,' +- ',f12.6)") 'Average <n> band 1 = ', mean, std
 endif
 
 call mpi_gather(bnup1,nbins,MPI_DOUBLE_PRECISION,bin_buff,nbins,MPI_DOUBLE_PRECISION,root,mpi_comm_world,ierr)
 if(myrank.eq.root)then 
  call geterr(bin_buff,bin_sgn,nnn,mean,std)
  write(unit=nfwrk1,fmt="(a50,f12.6,' +- ',f12.6)") 'Average <nup> band 1 = ', mean, std
 endif
 call mpi_gather(bndn1,nbins,MPI_DOUBLE_PRECISION,bin_buff,nbins,MPI_DOUBLE_PRECISION,root,mpi_comm_world,ierr)
 if(myrank.eq.root)then 
  call geterr(bin_buff,bin_sgn,nnn,mean,std)
  write(unit=nfwrk1,fmt="(a50,f12.6,' +- ',f12.6)") 'Average <ndn> band 1 = ', mean, std
 endif
 call mpi_gather(bnud1,nbins,MPI_DOUBLE_PRECISION,bin_buff,nbins,MPI_DOUBLE_PRECISION,root,mpi_comm_world,ierr)
 if(myrank.eq.root)then 
  call geterr(bin_buff,bin_sgn,nnn,mean,std)
  write(unit=nfwrk1,fmt="(a50,f12.6,' +- ',f12.6)") 'Average <nup*ndn> band 1 = ', mean, std
 endif

  call mpi_gather(bmz1,nbins,MPI_DOUBLE_PRECISION,bin_buff,nbins,MPI_DOUBLE_PRECISION,root,mpi_comm_world,ierr)
  if(myrank.eq.root)then
   call geterr(bin_buff,bin_sgn,nnn,mean,std)
   write(unit=nfwrk1,fmt="(a50,f12.6,' +- ',f12.6)") 'Average <m_z^2> band 1 = ', mean, std
  endif

 call mpi_gather(bn2,nbins,MPI_DOUBLE_PRECISION,bin_buff,nbins,MPI_DOUBLE_PRECISION,root,mpi_comm_world,ierr)
 if(myrank.eq.root)then 
  call geterr(bin_buff,bin_sgn,nnn,mean,std)
    write(unit=nfwrk1,fmt=*) ' '
  write(unit=nfwrk1,fmt="(a50,f12.6,' +- ',f12.6)") 'Average <n> band 2 = ', mean, std
 endif
 call mpi_gather(bnup2,nbins,MPI_DOUBLE_PRECISION,bin_buff,nbins,MPI_DOUBLE_PRECISION,root,mpi_comm_world,ierr)
 if(myrank.eq.root)then 
  call geterr(bin_buff,bin_sgn,nnn,mean,std)
  write(unit=nfwrk1,fmt="(a50,f12.6,' +- ',f12.6)") 'Average <nup> band 2 = ', mean, std
 endif 
 call mpi_gather(bndn2,nbins,MPI_DOUBLE_PRECISION,bin_buff,nbins,MPI_DOUBLE_PRECISION,root,mpi_comm_world,ierr)
 if(myrank.eq.root)then 
  call geterr(bin_buff,bin_sgn,nnn,mean,std)
  write(unit=nfwrk1,fmt="(a50,f12.6,' +- ',f12.6)") 'Average <ndn> band 2 = ', mean, std
 endif
 call mpi_gather(bnud2,nbins,MPI_DOUBLE_PRECISION,bin_buff,nbins,MPI_DOUBLE_PRECISION,root,mpi_comm_world,ierr)
 if(myrank.eq.root)then 
  call geterr(bin_buff,bin_sgn,nnn,mean,std)
  write(unit=nfwrk1,fmt="(a50,f12.6,' +- ',f12.6)") 'Average <nup*ndn> band 2 = ', mean, std
 endif

 call mpi_gather(bmz2,nbins,MPI_DOUBLE_PRECISION,bin_buff,nbins,MPI_DOUBLE_PRECISION,root,mpi_comm_world,ierr)
  if(myrank.eq.root)then
  call geterr(bin_buff,bin_sgn,nnn,mean,std)
  write(unit=nfwrk1,fmt="(a50,f12.6,' +- ',f12.6)") 'Average <m_z^2> band 2 = ', mean, std
  endif



 call mpi_gather(bn3,nbins,MPI_DOUBLE_PRECISION,bin_buff,nbins,MPI_DOUBLE_PRECISION,root,mpi_comm_world,ierr)
 if(myrank.eq.root)then 
  call geterr(bin_buff,bin_sgn,nnn,mean,std)
    write(unit=nfwrk1,fmt=*) ' '
   write(unit=nfwrk1,fmt="(a50,f12.6,' +- ',f12.6)") 'Average <n> band 3 = ', mean, std
 endif
 call mpi_gather(bnup3,nbins,MPI_DOUBLE_PRECISION,bin_buff,nbins,MPI_DOUBLE_PRECISION,root,mpi_comm_world,ierr)
 if(myrank.eq.root)then 
  call geterr(bin_buff,bin_sgn,nnn,mean,std)
  write(unit=nfwrk1,fmt="(a50,f12.6,' +- ',f12.6)") 'Average <nup> band 3 = ', mean, std
 endif 
 call mpi_gather(bndn3,nbins,MPI_DOUBLE_PRECISION,bin_buff,nbins,MPI_DOUBLE_PRECISION,root,mpi_comm_world,ierr)
 if(myrank.eq.root)then 
  call geterr(bin_buff,bin_sgn,nnn,mean,std)
  write(unit=nfwrk1,fmt="(a50,f12.6,' +- ',f12.6)") 'Average <ndn> band 3 = ', mean, std
 endif
 call mpi_gather(bnud3,nbins,MPI_DOUBLE_PRECISION,bin_buff,nbins,MPI_DOUBLE_PRECISION,root,mpi_comm_world,ierr)
 if(myrank.eq.root)then 
  call geterr(bin_buff,bin_sgn,nnn,mean,std)
  write(unit=nfwrk1,fmt="(a50,f12.6,' +- ',f12.6)") 'Average <nup*ndn> band 3 = ', mean, std
 endif
 call mpi_gather(bmz3,nbins,MPI_DOUBLE_PRECISION,bin_buff,nbins,MPI_DOUBLE_PRECISION,root,mpi_comm_world,ierr)
 if(myrank.eq.root)then 
  call geterr(bin_buff,bin_sgn,nnn,mean,std)
  write(unit=nfwrk1,fmt="(a50,f12.6,' +- ',f12.6)") 'Average <m_z^2> band 3 = ', mean, std
 endif

 call mpi_gather(bn4,nbins,MPI_DOUBLE_PRECISION,bin_buff,nbins,MPI_DOUBLE_PRECISION,root,mpi_comm_world,ierr)
 if(myrank.eq.root)then
   call geterr(bin_buff,bin_sgn,nnn,mean,std)
     write(unit=nfwrk1,fmt=*) ' '
   write(unit=nfwrk1,fmt="(a50,f12.6,' +- ',f12.6)") 'Average <n> band 4 = ', mean, std
endif
call mpi_gather(bnup4,nbins,MPI_DOUBLE_PRECISION,bin_buff,nbins,MPI_DOUBLE_PRECISION,root,mpi_comm_world,ierr)
if(myrank.eq.root)then
 call geterr(bin_buff,bin_sgn,nnn,mean,std)
 write(unit=nfwrk1,fmt="(a50,f12.6,' +- ',f12.6)") 'Average <nup> band 4 = ', mean, std
endif
call mpi_gather(bndn4,nbins,MPI_DOUBLE_PRECISION,bin_buff,nbins,MPI_DOUBLE_PRECISION,root,mpi_comm_world,ierr)
if(myrank.eq.root)then
  call geterr(bin_buff,bin_sgn,nnn,mean,std)
  write(unit=nfwrk1,fmt="(a50,f12.6,' +- ',f12.6)") 'Average <ndn> band 4 = ', mean, std
endif
call mpi_gather(bnud4,nbins,MPI_DOUBLE_PRECISION,bin_buff,nbins,MPI_DOUBLE_PRECISION,root,mpi_comm_world,ierr)
if(myrank.eq.root)then
  call geterr(bin_buff,bin_sgn,nnn,mean,std)
  write(unit=nfwrk1,fmt="(a50,f12.6,' +- ',f12.6)") 'Average <nup*ndn> band 4 = ', mean, std
endif
call mpi_gather(bmz4,nbins,MPI_DOUBLE_PRECISION,bin_buff,nbins,MPI_DOUBLE_PRECISION,root,mpi_comm_world,ierr)
if(myrank.eq.root)then
  call geterr(bin_buff,bin_sgn,nnn,mean,std)
  write(unit=nfwrk1,fmt="(a50,f12.6,' +- ',f12.6)") 'Average <m_z^2> band 4 = ', mean, std
endif



 !output the equal-time green's function
 if(myrank.eq.root)then
  write(unit=nfwrk1,fmt=*) ' '
  write(unit=nfwrk1,fmt=*) ' '
  write(unit=nfwrk1,fmt="(a50)") 'in R space Equal-time Green''s functin:'
 endif 

 do i = 0,Nlat-1
  j = 0
  do orb1 = 1,4
   do orb2 = 1,4
    j=j+1
    call mpi_gather(bgrfun(:,i,orb1,orb2),nbins,MPI_DOUBLE_PRECISION,bin_buff,nbins,MPI_DOUBLE_PRECISION,root,mpi_comm_world,ierr)
    if(myrank.eq.root)then 
     call geterr(bin_buff,bin_sgn,nnn,buf1(j),buf2(j))
    endif
   enddo
  enddo

  if(myrank.eq.root)then 
   write(unit=nfwrk1,fmt="(i3,' ',16(f10.5,' +- ', f10.5))") i, (buf1(j),buf2(j),j=1,16)
  endif 
 enddo

if(myrank.eq.root)then
write(unit=nfwrk1,fmt=*) ' '
write(unit=nfwrk1,fmt=*) ' '
write(unit=nfwrk1,fmt="(a50)") 'in k space Equal-time Green''s functin:'
endif

do i = 0,Nlat-1
  j = 0
  do orb1 = 1,4
  do orb2 = 1,4
    j=j+1
    call mpi_gather(dreal(bgkfun(:,i,orb1,orb2)),nbins,MPI_DOUBLE_PRECISION,&
                 bin_buff,nbins,MPI_DOUBLE_PRECISION,root,mpi_comm_world,ierr)
    if(myrank.eq.root)then
      call geterr(bin_buff,bin_sgn,nnn,buf1(j),buf2(j))
    endif
  enddo
  enddo

  j = 0
  do orb1 = 1,4
  do orb2 = 1,4
    j=j+1
    call mpi_gather(dimag(bgkfun(:,i,orb1,orb2)),nbins,MPI_DOUBLE_PRECISION,&
                   bin_buff,nbins,MPI_DOUBLE_PRECISION,root,mpi_comm_world,ierr)
    if(myrank.eq.root)then
      call geterr(bin_buff,bin_sgn,nnn,buf3(j),buf4(j))
    endif
  enddo
  enddo

if(myrank.eq.root)then
  write(unit=nfwrk1,fmt="(i3,' ',a8 ,16(f10.5,' +- ', f10.5))") i,'Re: ', (buf1(j),buf2(j),j=1,16)
  write(unit=nfwrk1,fmt="(a12 ,16(f10.5,' +- ', f10.5))") 'Im: ', (buf3(j),buf4(j),j=1,16)
endif
enddo







 if(myrank.eq.root)then
  write(unit=nfwrk1,fmt=*) ' '
  write(unit=nfwrk1,fmt="(a50)") 'Spinxx correlation functions in r space:'
 endif 

 do i = 0,Nlat-1
  j = 0
  do orb1 = 1,4
  do orb2 = 1,4
   j=j+1
   call mpi_gather(bspinxx(:,i,orb1,orb2),nbins,MPI_DOUBLE_PRECISION,&
                 bin_buff,nbins,MPI_DOUBLE_PRECISION,root,mpi_comm_world,ierr)
   if(myrank.eq.root)then
     call geterr(bin_buff,bin_sgn,nnn,buf1(j),buf2(j))
   endif
  enddo
  enddo

if(myrank.eq.root)then
write(unit=nfwrk1,fmt="(i3,' ',16(f10.5,' +- ', f10.5))") i, (buf1(j),buf2(j),j=1,16)
endif
enddo


if(myrank.eq.root)then
write(unit=nfwrk1,fmt=*) ' '
write(unit=nfwrk1,fmt=*) ' '
write(unit=nfwrk1,fmt="(a50)") 'Spinxx correlation functions in k space:'
endif

do i = 0,Nlat-1
j = 0
do orb1 = 1,4
do orb2 = 1,4
j=j+1
call mpi_gather(dreal(bspinxxq(:,i,orb1,orb2)),nbins,MPI_DOUBLE_PRECISION,&
                 bin_buff,nbins,MPI_DOUBLE_PRECISION,root,mpi_comm_world,ierr)
if(myrank.eq.root)then
call geterr(bin_buff,bin_sgn,nnn,buf1(j),buf2(j))
endif
enddo
enddo

j = 0
do orb1 = 1,4
do orb2 = 1,4
j=j+1
call mpi_gather(dimag(bspinxxq(:,i,orb1,orb2)),nbins,MPI_DOUBLE_PRECISION,&
                 bin_buff,nbins,MPI_DOUBLE_PRECISION,root,mpi_comm_world,ierr)
if(myrank.eq.root)then
call geterr(bin_buff,bin_sgn,nnn,buf3(j),buf4(j))
endif
enddo
enddo

if(myrank.eq.root)then
write(unit=nfwrk1,fmt="(i3,' ',a8 ,16(f10.5,' +- ', f10.5))") i,'Re: ', (buf1(j),buf2(j),j=1,16)
write(unit=nfwrk1,fmt="(a12 ,16(f10.5,' +- ', f10.5))") 'Im: ', (buf3(j),buf4(j),j=1,16)
endif
enddo


if(myrank.eq.root)then
write(unit=nfwrk1,fmt=*) ' '
write(unit=nfwrk1,fmt="(a50)") 'Spinzz correlation functions in r space:'
endif

do i = 0,Nlat-1
j = 0
do orb1 = 1,4
do orb2 = 1,4
j=j+1
call mpi_gather(bspinzz(:,i,orb1,orb2),nbins,MPI_DOUBLE_PRECISION,&
                 bin_buff,nbins,MPI_DOUBLE_PRECISION,root,mpi_comm_world,ierr)
if(myrank.eq.root)then
call geterr(bin_buff,bin_sgn,nnn,buf1(j),buf2(j))
endif
enddo
enddo

if(myrank.eq.root)then
write(unit=nfwrk1,fmt="(i3,' ',16(f10.5,' +- ', f10.5))") i, (buf1(j),buf2(j),j=1,16)
endif
enddo


if(myrank.eq.root)then
write(unit=nfwrk1,fmt=*) ' '
write(unit=nfwrk1,fmt=*) ' '
write(unit=nfwrk1,fmt="(a50)") 'Spinzz correlation functions in k space:'
endif

do i = 0,Nlat-1
j = 0
do orb1 = 1,4
do orb2 = 1,4
j=j+1
call mpi_gather(dreal(bspinzzq(:,i,orb1,orb2)),nbins,MPI_DOUBLE_PRECISION,&
                bin_buff,nbins,MPI_DOUBLE_PRECISION,root,mpi_comm_world,ierr)
if(myrank.eq.root)then
call geterr(bin_buff,bin_sgn,nnn,buf1(j),buf2(j))
endif
enddo
enddo

j = 0
do orb1 = 1,4
do orb2 = 1,4
j=j+1
call mpi_gather(dimag(bspinzzq(:,i,orb1,orb2)),nbins,MPI_DOUBLE_PRECISION,&
                 bin_buff,nbins,MPI_DOUBLE_PRECISION,root,mpi_comm_world,ierr)
if(myrank.eq.root)then
call geterr(bin_buff,bin_sgn,nnn,buf3(j),buf4(j))
endif
enddo
enddo

if(myrank.eq.root)then
write(unit=nfwrk1,fmt="(i3,' ',a8 ,16(f10.5,' +- ', f10.5))") i,'Re: ', (buf1(j),buf2(j),j=1,16)
write(unit=nfwrk1,fmt="(a12 ,16(f10.5,' +- ', f10.5))") 'Im: ', (buf3(j),buf4(j),j=1,16)
endif
enddo


if(myrank.eq.root)then
write(unit=nfwrk1,fmt=*) ' '
write(unit=nfwrk1,fmt="(a50)") 'Density correlation functions in r space:'
endif

do i = 0,Nlat-1
j = 0
do orb1 = 1,4
do orb2 = 1,4
j=j+1
call mpi_gather(bden(:,i,orb1,orb2),nbins,MPI_DOUBLE_PRECISION,&
              bin_buff,nbins,MPI_DOUBLE_PRECISION,root,mpi_comm_world,ierr)
if(myrank.eq.root)then
call geterr(bin_buff,bin_sgn,nnn,buf1(j),buf2(j))
endif
enddo
enddo

if(myrank.eq.root)then
write(unit=nfwrk1,fmt="(i3,' ',16(f10.5,' +- ', f10.5))") i, (buf1(j),buf2(j),j=1,16)
endif
enddo


if(myrank.eq.root)then
write(unit=nfwrk1,fmt=*) ' '
write(unit=nfwrk1,fmt=*) ' '
write(unit=nfwrk1,fmt="(a50)") 'Density correlation functions in k space:'
endif

do i = 0,Nlat-1
j = 0
do orb1 = 1,4
do orb2 = 1,4
j=j+1
call mpi_gather(dreal(bdenq(:,i,orb1,orb2)),nbins,MPI_DOUBLE_PRECISION,&
                bin_buff,nbins,MPI_DOUBLE_PRECISION,root,mpi_comm_world,ierr)
if(myrank.eq.root)then
call geterr(bin_buff,bin_sgn,nnn,buf1(j),buf2(j))
endif
enddo
enddo

j = 0
do orb1 = 1,4
do orb2 = 1,4
j=j+1
call mpi_gather(dimag(bdenq(:,i,orb1,orb2)),nbins,MPI_DOUBLE_PRECISION,&
           bin_buff,nbins,MPI_DOUBLE_PRECISION,root,mpi_comm_world,ierr)
if(myrank.eq.root)then
call geterr(bin_buff,bin_sgn,nnn,buf3(j),buf4(j))
endif
enddo
enddo

if(myrank.eq.root)then
write(unit=nfwrk1,fmt="(i3,' ',a8 ,16(f10.5,' +- ', f10.5))") i,'Re: ', (buf1(j),buf2(j),j=1,16)
write(unit=nfwrk1,fmt="(a12 ,16(f10.5,' +- ', f10.5))") 'Im: ', (buf3(j),buf4(j),j=1,16)
endif
enddo



 if(myrank.eq.root)then
  write(unit=nfwrk1,fmt=*) ' '
  write(unit=nfwrk1,fmt=*) ' '
  write(unit=nfwrk1,fmt="(a50)") 'G(r, tau) functin:'
 endif 

 do i = 0,Nlat-1 
  do ti = 0,L
   j = 0
   do orb1 = 1,4
    do orb2 = 1,4
     j = j + 1
     call mpi_gather(bgnl(:,i,ti,orb1,orb2),nbins,MPI_DOUBLE_PRECISION,&
             bin_buff,nbins,MPI_DOUBLE_PRECISION,root,mpi_comm_world,ierr)
     if(myrank.eq.root)then 
      call geterr(bin_buff,bin_sgnt,nnn,buf1(j),buf2(j))
     endif
    enddo
   enddo

   if(myrank.eq.root)then 
    write(unit=nfwrk1,fmt="(i3,' ',i3,' ',16(f10.5,' +- ', f10.5))") i, ti, (buf1(j),buf2(j),j=1,16)
   endif
  enddo
  if(myrank.eq.root)then
   write(unit=nfwrk1,fmt=*) ' '
  endif
 enddo


 if(myrank.eq.root)then
  write(unit=nfwrk1,fmt=*) ' '
  write(unit=nfwrk1,fmt=*) ' '
  write(unit=nfwrk1,fmt="(a50)") 'G(q, tau) functin: q=[0,2pi)'
 endif 

 do i = 0,Nlat-1 
  do ti = 0,L
   j = 0
   do orb1 = 1,4
    do orb2 = 1,4
     j = j + 1
     call mpi_gather(real(bgnlq(:,i,ti,orb1,orb2)),nbins,MPI_DOUBLE_PRECISION,&
                 bin_buff,nbins,MPI_DOUBLE_PRECISION,root,mpi_comm_world,ierr)
     if(myrank.eq.root)then 
      call geterr(bin_buff,bin_sgnt,nnn,buf1(j),buf2(j))
     endif
    enddo
   enddo

   j = 0
   do orb1 = 1,4
    do orb2 = 1,4
      j = j + 1
      call mpi_gather(dimag(bgnlq(:,i,ti,orb1,orb2)),nbins,MPI_DOUBLE_PRECISION,&
           bin_buff,nbins,MPI_DOUBLE_PRECISION,root,mpi_comm_world,ierr)
      if(myrank.eq.root)then
        call geterr(bin_buff,bin_sgnt,nnn,buf3(j),buf4(j))
      endif
   enddo
   enddo

   if(myrank.eq.root)then
    write(unit=nfwrk1,fmt="(i3,' ',i3,' ',a8,16(f10.5,' +- ', f10.5))") i,ti,"Re: ", (buf1(j),buf2(j),j=1,16)
    write(unit=nfwrk1,fmt="(a16,16(f10.5,' +- ', f10.5))") "Im: ", (buf3(j),buf4(j),j=1,16)
   endif
  enddo
  if(myrank.eq.root)then
   write(unit=nfwrk1,fmt=*) ' '
  endif
 enddo




if(myrank.eq.root)then
write(unit=nfwrk1,fmt=*) ' '
write(unit=nfwrk1,fmt=*) ' '
write(unit=nfwrk1,fmt="(a50)") 'spinxx(r, tau) functin:'
endif

do i = 0,Nlat-1
do ti = 0,L
j = 0
do orb1 = 1,4
do orb2 = 1,4
j = j + 1
call mpi_gather(bchinl(:,i,ti,orb1,orb2),nbins,MPI_DOUBLE_PRECISION,&
              bin_buff,nbins,MPI_DOUBLE_PRECISION,root,mpi_comm_world,ierr)
if(myrank.eq.root)then
call geterr(bin_buff,bin_sgnt,nnn,buf1(j),buf2(j))
endif
enddo
enddo

if(myrank.eq.root)then
write(unit=nfwrk1,fmt="(i3,' ',i3,' ',16(f10.5,' +- ', f10.5))") i, ti, (buf1(j),buf2(j),j=1,16)
endif
enddo
if(myrank.eq.root)then
write(unit=nfwrk1,fmt=*) ' '
endif
enddo


if(myrank.eq.root)then
write(unit=nfwrk1,fmt=*) ' '
write(unit=nfwrk1,fmt=*) ' '
write(unit=nfwrk1,fmt="(a50)") 'spinxx(q, tau) functin: q=[0,2pi)'
endif

do i = 0,Nlat-1
do ti = 0,L
j = 0
do orb1 = 1,4
do orb2 = 1,4
j = j + 1
call mpi_gather(real(bchinlq(:,i,ti,orb1,orb2)),nbins,MPI_DOUBLE_PRECISION,&
bin_buff,nbins,MPI_DOUBLE_PRECISION,root,mpi_comm_world,ierr)
if(myrank.eq.root)then
call geterr(bin_buff,bin_sgnt,nnn,buf1(j),buf2(j))
endif
enddo
enddo

j = 0
do orb1 = 1,4
do orb2 = 1,4
j = j + 1
call mpi_gather(dimag(bchinlq(:,i,ti,orb1,orb2)),nbins,MPI_DOUBLE_PRECISION,&
bin_buff,nbins,MPI_DOUBLE_PRECISION,root,mpi_comm_world,ierr)
if(myrank.eq.root)then
call geterr(bin_buff,bin_sgnt,nnn,buf3(j),buf4(j))
endif
enddo
enddo

if(myrank.eq.root)then
write(unit=nfwrk1,fmt="(i3,' ',i3,' ',a8,16(f10.5,' +- ', f10.5))") i, ti,"Re: ", (buf1(j),buf2(j),j=1,16)
write(unit=nfwrk1,fmt="(a16,16(f10.5,' +- ', f10.5))") "Im: ", (buf3(j),buf4(j),j=1,16)
endif
enddo
if(myrank.eq.root)then
write(unit=nfwrk1,fmt=*) ' '
endif
enddo


if(myrank.eq.root)then
write(unit=nfwrk1,fmt=*) ' '
write(unit=nfwrk1,fmt=*) ' '
write(unit=nfwrk1,fmt="(a50)") 'spinzz(r, tau) functin:'
endif

do i = 0,Nlat-1
do ti = 0,L
j = 0
do orb1 = 1,4
do orb2 = 1,4
j = j + 1
call mpi_gather(bchinlz(:,i,ti,orb1,orb2),nbins,MPI_DOUBLE_PRECISION,&
          bin_buff,nbins,MPI_DOUBLE_PRECISION,root,mpi_comm_world,ierr)
if(myrank.eq.root)then
call geterr(bin_buff,bin_sgnt,nnn,buf1(j),buf2(j))
endif
enddo
enddo

if(myrank.eq.root)then
write(unit=nfwrk1,fmt="(i3,' ',i3,' ',16(f10.5,' +- ', f10.5))") i, ti, (buf1(j),buf2(j),j=1,16)
endif
enddo
if(myrank.eq.root)then
write(unit=nfwrk1,fmt=*) ' '
endif
enddo


if(myrank.eq.root)then
write(unit=nfwrk1,fmt=*) ' '
write(unit=nfwrk1,fmt=*) ' '
write(unit=nfwrk1,fmt="(a50)") 'spinzz(q, tau) functin: q=[0,2pi)'
endif

do i = 0,Nlat-1
do ti = 0,L
j = 0
do orb1 = 1,4
do orb2 = 1,4
j = j + 1
call mpi_gather(real(bchinlzq(:,i,ti,orb1,orb2)),nbins,MPI_DOUBLE_PRECISION,&
bin_buff,nbins,MPI_DOUBLE_PRECISION,root,mpi_comm_world,ierr)
if(myrank.eq.root)then
call geterr(bin_buff,bin_sgnt,nnn,buf1(j),buf2(j))
endif
enddo
enddo

j = 0
do orb1 = 1,4
do orb2 = 1,4
j = j + 1
call mpi_gather(dimag(bchinlzq(:,i,ti,orb1,orb2)),nbins,MPI_DOUBLE_PRECISION,&
bin_buff,nbins,MPI_DOUBLE_PRECISION,root,mpi_comm_world,ierr)
if(myrank.eq.root)then
call geterr(bin_buff,bin_sgnt,nnn,buf3(j),buf4(j))
endif
enddo
enddo

if(myrank.eq.root)then
write(unit=nfwrk1,fmt="(i3,' ',i3,' ',a8,16(f10.5,' +- ', f10.5))") i, ti,"Re: ", (buf1(j),buf2(j),j=1,16)
write(unit=nfwrk1,fmt="(a16,16(f10.5,' +- ', f10.5))") "Im: ", (buf3(j),buf4(j),j=1,16)
endif
enddo
if(myrank.eq.root)then
write(unit=nfwrk1,fmt=*) ' '
endif
enddo


if(myrank.eq.root)then
write(unit=nfwrk1,fmt=*) ' '
write(unit=nfwrk1,fmt=*) ' '
write(unit=nfwrk1,fmt="(a50)") 'dent(r, tau) functin:'
endif

do i = 0,Nlat-1
do ti = 0,L
j = 0
do orb1 = 1,4
do orb2 = 1,4
j = j + 1
call mpi_gather(bdent(:,i,ti,orb1,orb2),nbins,MPI_DOUBLE_PRECISION,&
       bin_buff,nbins,MPI_DOUBLE_PRECISION,root,mpi_comm_world,ierr)
if(myrank.eq.root)then
call geterr(bin_buff,bin_sgnt,nnn,buf1(j),buf2(j))
endif
enddo
enddo

if(myrank.eq.root)then
write(unit=nfwrk1,fmt="(i3,' ',i3,' ',16(f10.5,' +- ', f10.5))") i, ti, (buf1(j),buf2(j),j=1,16)
endif
enddo
if(myrank.eq.root)then
write(unit=nfwrk1,fmt=*) ' '
endif
enddo


if(myrank.eq.root)then
write(unit=nfwrk1,fmt=*) ' '
write(unit=nfwrk1,fmt=*) ' '
write(unit=nfwrk1,fmt="(a50)") 'dent(q, tau) functin: q=[0,2pi)'
endif

do i = 0,Nlat-1
do ti = 0,L
j = 0
do orb1 = 1,4
do orb2 = 1,4
j = j + 1
call mpi_gather(real(bdentq(:,i,ti,orb1,orb2)),nbins,MPI_DOUBLE_PRECISION,&
bin_buff,nbins,MPI_DOUBLE_PRECISION,root,mpi_comm_world,ierr)
if(myrank.eq.root)then
call geterr(bin_buff,bin_sgnt,nnn,buf1(j),buf2(j))
endif
enddo
enddo

j = 0
do orb1 = 1,4
do orb2 = 1,4
j = j + 1
call mpi_gather(dimag(bdentq(:,i,ti,orb1,orb2)),nbins,MPI_DOUBLE_PRECISION,&
bin_buff,nbins,MPI_DOUBLE_PRECISION,root,mpi_comm_world,ierr)
if(myrank.eq.root)then
call geterr(bin_buff,bin_sgnt,nnn,buf3(j),buf4(j))
endif
enddo
enddo

if(myrank.eq.root)then
write(unit=nfwrk1,fmt="(i3,' ',i3,' ',a8,16(f10.5,' +- ', f10.5))") i, ti,"Re: ", (buf1(j),buf2(j),j=1,16)
write(unit=nfwrk1,fmt="(a16,16(f10.5,' +- ', f10.5))") "Im: ", (buf3(j),buf4(j),j=1,16)
endif
enddo
if(myrank.eq.root)then
write(unit=nfwrk1,fmt=*) ' '
endif
enddo


!========== output for dynamic function ============
   do kproc=1,nbins
     bin_GF_buff(:,:,:,:)=bgnl(kproc,:,:,:,:)
     bin_sgn_buff=bsgnt(kproc)
     call mpi_barrier(MPI_COMM_WORLD,ierr)
       if(myrank.eq.root)then
       do iproc=0,nproc-1
         if(iproc.ne.0) call MPI_RECV(bin_GF_buff,(L+1)*(Nlat)*4*4,&
                                     MPI_DOUBLE_PRECISION,iproc, &
                                     itag1,MPI_COMM_WORLD,istat,ierr)
         if(iproc.ne.0) call MPI_RECV(bin_sgn_buff,1, &
                                     MPI_DOUBLE_PRECISION,iproc, &
                                     itag2,MPI_COMM_WORLD,istat,ierr)
            write(unit=nfwrk3,fmt=*) '    '
            write(unit=nfwrk3,fmt=*) '    '
            write(unit=nfwrk3,fmt=*) ' bin of G(r,tau) from process   ', iproc*nbins+kproc
            write(unit=nfwrk3,fmt=*) ((((bin_GF_buff(i,ti,ii,jj), i=0,Nlat-1),ti=0,L),ii=1,4),jj=1,4),bin_sgn_buff
       end do
         else
        call MPI_ISEND(bin_GF_buff,(L+1)*(Nlat)*4*4,MPI_DOUBLE_PRECISION &
                                   ,0,itag1,MPI_COMM_WORLD,req(1),ierr)
        call MPI_ISEND(bin_sgn_buff,1,MPI_DOUBLE_PRECISION,0,itag2,&
                                            MPI_COMM_WORLD,req(2),ierr)
        call MPI_WAITALL(2,req,istat,ierr)
       endif
   end do


  do kproc=1,nbins
   bin_GF_buff(:,:,:,:)=bchinlz(kproc,:,:,:,:)
   bin_sgn_buff=bsgnt(kproc)
   call mpi_barrier(MPI_COMM_WORLD,ierr)
   if(myrank.eq.root)then
    do iproc=0,nproc-1
     if(iproc.ne.0) call MPI_RECV(bin_GF_buff,(L+1)*(Nlat)*4*4,&
                                   MPI_DOUBLE_PRECISION,iproc, &
                                   itag1,MPI_COMM_WORLD,istat,ierr)
     if(iproc.ne.0) call MPI_RECV(bin_sgn_buff,1, &
                                 MPI_DOUBLE_PRECISION,iproc, &
                                 itag2,MPI_COMM_WORLD,istat,ierr)
     write(unit=nfwrk4,fmt=*) '    '
     write(unit=nfwrk4,fmt=*) '    '
     write(unit=nfwrk4,fmt=*) ' bin of spinzz(r,tau) from process   ', iproc*nbins+kproc
     write(unit=nfwrk4,fmt=*) ((((bin_GF_buff(i,ti,ii,jj), i=0,Nlat-1),ti=0,L),ii=1,4),jj=1,4),bin_sgn_buff
    end do
  else
    call MPI_ISEND(bin_GF_buff,(L+1)*(Nlat)*4*4,MPI_DOUBLE_PRECISION &
                   ,0,itag1,MPI_COMM_WORLD,req(1),ierr)
    call MPI_ISEND(bin_sgn_buff,1,MPI_DOUBLE_PRECISION,0,itag2,&
                   MPI_COMM_WORLD,req(2),ierr)
    call MPI_WAITALL(2,req,istat,ierr)
  endif
 end do





 do kproc=1,nbins
   bin_GF_buff(:,:,:,:)=bdent(kproc,:,:,:,:)
   bin_sgn_buff=bsgnt(kproc)
   call mpi_barrier(MPI_COMM_WORLD,ierr)
   if(myrank.eq.root)then
     do iproc=0,nproc-1
       if(iproc.ne.0) call MPI_RECV(bin_GF_buff,(L+1)*(Nlat)*4*4,&
                                    MPI_DOUBLE_PRECISION,iproc, &
                                    itag1,MPI_COMM_WORLD,istat,ierr)
       if(iproc.ne.0) call MPI_RECV(bin_sgn_buff,1, &
                                    MPI_DOUBLE_PRECISION,iproc, &
                                    itag2,MPI_COMM_WORLD,istat,ierr)
       write(unit=nfwrk2,fmt=*) '    '
       write(unit=nfwrk2,fmt=*) '    '
       write(unit=nfwrk2,fmt=*) ' bin of dent(r,tau) from process   ', iproc*nbins+kproc
       write(unit=nfwrk2,fmt=*) ((((bin_GF_buff(i,ti,ii,jj), i=0,Nlat-1),ti=0,L),ii=1,4),jj=1,4),bin_sgn_buff
    end do
  else
    call MPI_ISEND(bin_GF_buff,(L+1)*(Nlat)*4*4,MPI_DOUBLE_PRECISION &
                   ,0,itag1,MPI_COMM_WORLD,req(1),ierr)
    call MPI_ISEND(bin_sgn_buff,1,MPI_DOUBLE_PRECISION,0,itag2,&
                   MPI_COMM_WORLD,req(2),ierr)
    call MPI_WAITALL(2,req,istat,ierr)
  endif
 end do

 do kproc=1,nbins
    bin_GF_x(1)=bnt1(kproc)
    bin_GF_x(2)=bnt2(kproc)
    bin_GF_x(3)=bnt3(kproc)
    bin_GF_x(4)=bnt4(kproc)
    call mpi_barrier(MPI_COMM_WORLD,ierr)
    if(myrank.eq.root)then
      do iproc=0,nproc-1
        if(iproc.ne.0) call MPI_RECV(bin_GF_x,4,&
                                     MPI_DOUBLE_PRECISION,iproc, &
                                     itag1,MPI_COMM_WORLD,istat,ierr)
        write(unit=nfwrk2,fmt=*) '    '
        write(unit=nfwrk2,fmt=*) '    '
        write(unit=nfwrk2,fmt=*) ' bin of GF_x from process   ', iproc*nbins+kproc
        write(unit=nfwrk2,fmt=*) (bin_GF_x(i), i=1,4)
      end do
    else
       call MPI_ISEND(bin_GF_x,4,MPI_DOUBLE_PRECISION &
                      ,0,itag1,MPI_COMM_WORLD,req(1),ierr)
       call MPI_WAITALL(1,req(1),istat,ierr)
    endif
 enddo





 if(myrank.eq.root)then
  close(unit=nfwrk1)
  close(unit=nfwrk3)
  close(unit=nfwrk4)
  close(unit=nfwrk2)
 endif

 return
 end subroutine output_results
end module io_routines
