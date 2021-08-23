module parameters
 integer, parameter :: nwarms = 50000
 integer, parameter :: nsweep = 10000
 integer, parameter :: stop_numb=2000
 integer, parameter :: seperate_run_time=0
 integer, parameter :: Nlat = 20
 integer, parameter :: N = Nlat*4
 integer, parameter :: Nhis=400  
 integer, parameter :: Block_site_number=1
 integer, parameter :: Block_spin_number=1     
 integer, parameter :: nwraps = 8
 integer wrap,nwrap,wraps
 integer, parameter :: iseed = 24007
 integer, parameter :: run = 1
 integer, parameter :: orthlen = 10
 integer, parameter :: sweeps_per_block_ph = 1
 integer, parameter :: sweeps_per_block_hs = 10
 integer, parameter :: spinup = 1
 integer, parameter :: spindn =-1
 integer, parameter :: left = 1
 integer, parameter :: right = -1
 double precision errrat, difflim 
 double precision, parameter :: beta = 16.0d0
 integer, parameter :: L = 10*beta
 integer, parameter :: tauskip = 50  !tauskip>nsweep/nbin
 integer, parameter :: maxn = 10
 double precision, parameter :: mu=0.0d0
 double precision, parameter :: deltaCu = 0.0d0
 double precision, parameter :: deltapx = 3.0d0
 double precision, parameter :: deltapy = 3.5d0
 double precision, parameter :: tpxd = 1.5d0
 double precision, parameter :: tpyd = 1.8d0
 double precision, parameter :: tpp =  0.75d0
 double precision, parameter :: Ud = 8.0d0
 double precision, parameter :: Up = 4.0d0
 double precision, parameter :: Upd = 1.0d0
 double precision dtau, lambda, lambdap, lambdapp, gam
 double precision, allocatable, dimension(:,:) :: kin, tt,ttk, expk, expki
 character outfile*200, outfile2*200, outfile3*200, outfile4*200
contains
 !=======================================================================
 subroutine init_cluster()
 implicit none
 double precision rtmp
 wrap = 0
 nwrap = nwraps
 gam = 0.5d0 
 dtau = beta/dfloat(L)

 rtmp = dexp(dtau*Ud*0.5d0)
 lambda = dlog(rtmp + dsqrt(rtmp*rtmp - 1.0d0))
 rtmp = dexp(dtau*Up*0.5d0)
 lambdap = dlog(rtmp + dsqrt(rtmp*rtmp - 1.0d0))
 rtmp = dexp(dtau*Upd*0.5d0)
 lambdapp = dlog(rtmp + dsqrt(rtmp*rtmp - 1.0d0))


 errrat = 0.00005d0
 difflim = 0.001d0
 write(unit=outfile,fmt="('DQMC_',i6.6,'.dat')") run
 write(unit=outfile3,fmt="('Green_',i6.6,'.dat')") run
 write(unit=outfile4,fmt="('Spinzz_',i6.6,'.dat')") run
 write(unit=outfile2,fmt="('Charge_',i6.6,'.dat')") run
 print*, trim(outfile)
 return
 end subroutine init_cluster

 !=======================================================================
 !this subroutine constructs the kinetic energy matrix
 !=======================================================================
 subroutine build_kin()
 use linear_algebra, only: exp_symmetric_mat, invertr
 implicit none
 integer i, j, ii, jj
 integer band
 double precision deta, rtmp
 allocate(tt(0:N-1,0:N-1))
 allocate(ttk(0:N-1,0:N-1))
 allocate(kin(0:N-1,0:N-1))
 allocate(expk(0:N-1,0:N-1))
 allocate(expki(0:N-1,0:N-1))

 kin = 0.0d0
 do i = 0,Nlat-1

  !do the on-site energies
  ii = return_idx_for_site(i,0)
  kin(ii,ii) = deltaCu + (Ud + 8.0d0*Upd)*0.5d0 - mu

  ii = return_idx_for_site(i,1)
  kin(ii,ii) = deltapx + (Up + 4.0d0*Upd)*0.5d0 - mu

  ii = return_idx_for_site(i,2)
  kin(ii,ii) = deltapy + (Up + 2.0d0*Upd)*0.5d0 - mu

  ii = return_idx_for_site(i,3)
  kin(ii,ii) = deltapy + (Up + 2.0d0*Upd)*0.5d0 - mu

  ii = return_idx_for_site(i,0)
  jj = return_idx_for_site(i,1)
  kin(ii,jj) = tpxd
  kin(jj,ii) = tpxd

  ii = return_idx_for_site(i,0)
  jj = return_idx_for_site(i+1,1)
  kin(ii,jj) = -tpxd
  kin(jj,ii) = -tpxd

  ii = return_idx_for_site(i,0)
  jj = return_idx_for_site(i,2)
  kin(ii,jj) =  tpyd
  kin(jj,ii) =  tpyd

  ii = return_idx_for_site(i,0)
  jj = return_idx_for_site(i,3)
  kin(ii,jj) = -tpyd
  kin(jj,ii) = -tpyd

  ii = return_idx_for_site(i,1)
  jj = return_idx_for_site(i,2)
  kin(ii,jj) = -tpp
  kin(jj,ii) = -tpp

  ii = return_idx_for_site(i,1)
  jj = return_idx_for_site(i,3)
  kin(ii,jj) = tpp
  kin(jj,ii) = tpp

  ii = return_idx_for_site(i,1)
  jj = return_idx_for_site(i-1,2)
  kin(ii,jj) =  tpp
  kin(jj,ii) =  tpp

  ii = return_idx_for_site(i,1)
  jj = return_idx_for_site(i-1,3)
  kin(ii,jj) = -tpp
  kin(jj,ii) = -tpp


 enddo
 tt = kin  !store the matrix for later

 ttk=kin
 do i=0,Nlat-1
  ii = return_idx_for_site(i,0)
  ttk(ii,ii) = ttk(ii,ii) - (Ud + 6.0d0*Upd)*0.5d0

  ii = return_idx_for_site(i,1)
  ttk(ii,ii) = ttk(ii,ii) - (Up + 2.0d0*Upd)*0.5d0

  ii = return_idx_for_site(i,2)
  ttk(ii,ii) = ttk(ii,ii) - (Up + 2.0d0*Upd)*0.5d0

  ii = return_idx_for_site(i,3)
  ttk(ii,ii) = ttk(ii,ii) - (Up + 2.0d0*Upd)*0.5d0
 enddo


 call exp_symmetric_mat(-dtau*Kin,expk,N)
 call invertr(expk,expki,deta,rtmp,N)

 return
 end subroutine build_kin

 !=======================================================================
 integer function return_idx_for_site(i,orb)
 implicit none
 integer idx, j, i, orb
 j = i
 if(j.ge.Nlat) j = j - Nlat
 if(j.lt.0) j = j + Nlat
 return_idx_for_site = j + orb*Nlat
 return
 end function return_idx_for_site
 !=======================================================================
end module parameters
!==============================================================================
