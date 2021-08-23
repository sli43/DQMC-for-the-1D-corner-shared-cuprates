module updates
 integer(KIND=8) sweeps_per_et_meas
 integer(KIND=8) accept_hs, reject_hs, redo, noredo
 integer(KIND=8) accept_blk_hs, reject_blk_hs
 integer(KIND=8) calls_to_sweep_hs
contains
 !=============================================================================
 subroutine single_site_updates(gup,gdn,sgnup,sgndn,warm)
 use random
 use grfun
 use measurements
 use parameters, only: N, L, Lambda, lambdap, lambdapp, nwrap, &
                       spinup, spindn, left, right, difflim, expk, &
                       errrat, gam, tauskip, wraps
 use fields, only: v_up, v_dn, S
 use linear_algebra 
 implicit none
 logical warm
 integer ii, jj
 integer ti, ti1
 double precision diffup, diffdn, sgnup, sgndn
 double precision sgnup_old, sgndn_old
 double precision detup, detdn, rtmp, accrat
 double precision, dimension(0:N-1) :: bvecup, bvecdn
 double precision, dimension(0:N-1,0:N-1) :: gup, gdn, ogmat

 do ti = 0,L-1
  bvecup(:) = v_up(:,ti)
  bvecdn(:) = v_dn(:,ti)
  call multb( gup,bvecup,-1)
  call multbi(gup,bvecup, 1)
  call multb( gdn,bvecdn,-1)
  call multbi(gdn,bvecdn, 1)  

  !wraps is number of calculations with Sherman-Morrison before recalculating
  !from scratch
  if(wraps.gt.0)then
   wraps = wraps - 1
  else
   wraps = nwrap
   ti1 = mod(ti+1,L)
   call dcopy(N*N,gup,1,ogmat,1)
   call getgp(gup,ti1,detup,sgnup,spinup)
   diffup = diffmat(gup,ogmat,N)
   call dcopy(N*N,gdn,1,ogmat,1)
   call getgp(gdn,ti1,detdn,sgndn,spindn)
   diffdn = diffmat(gdn,ogmat,N)  
 !  print*, diffup, diffdn
   if(diffup.gt.difflim.or.diffdn.gt.difflim)then 
    redo = redo + 1
   else
    noredo = noredo + 1
   endif
  endif

  !now we make the call to a routine to do the woodbury update.
  call swpslice_woodbury(gup,gdn,sgnup,sgndn,accept_hs,reject_hs,ti,gam)
!  call swpslice(gup,gdn,sgnup,sgndn,accept_hs,reject_hs,ti,gam)

  !Here is where code goes to do the measurements.
  if(sweeps_per_et_meas.eq.tauskip.and..not.warm)then
   sweeps_per_et_meas = 0
   call equaltime_measurements(gup,gdn,sgnup,sgndn)
  else
   sweeps_per_et_meas = sweeps_per_et_meas + 1
  endif
 enddo

 !adjust the nwraps if need be
 if((redo.gt.20).and.nwrap.gt.0)then
  rtmp = dfloat(redo)/dfloat(redo+noredo)
  if(rtmp.gt.errrat)then
   nwrap = nwrap - 1
   write(6,*) 'Reducing nwrap to ', nwrap
   redo = 0
   noredo = 0
  endif
 endif 
 
 !recalculate gam
 if(accept_hs+reject_hs.gt.10000)then
  accrat = dfloat(accept_hs)/dfloat(accept_hs+reject_hs)
  if(accrat.gt.0.52d0.or.accrat.lt.0.48d0)then
   gam = gam + accrat - 0.5d0
   gam = dmax1(0.00d0, gam)
   gam = dmin1(1.0d0, gam)
   accept_hs = int(100.0d0*accrat)
   reject_hs = int(100.0d0*(1.0d0-accrat))
  endif
 endif

 return 
 end subroutine single_site_updates
 !=============================================================================
 !=============================================================================
 subroutine swpslice_woodbury(gup,gdn,sgnup,sgndn,accept,reject,ti,gam)
 use parameters, only: N, Nlat, L, lambda,lambdap,lambdapp,spinup,spindn
 use random
 use fields, only: return_v_for_s, S, v_up, v_dn
 use grfun, only: getgp
 use linear_algebra, only: invertr, determinant
 implicit none
 integer(KIND=8) accept, reject
 integer ti, i, ii, ti1, j,jj, i1,j1, a, b, iii, jjj, si,ip
 integer, dimension(0:Nlat-1) :: list
 double precision gam, r, ratio, p, rtmp
 double precision, dimension(1:20) :: Sold, Snew, dS
 double precision, dimension(0:N-1,0:L-1) :: vupold, vdnold
 double precision, dimension(0:N-1,0:N-1) :: gup, gdn, gupold, gdnold
 double precision sgnup, sgndn, detup, detdn, sgnupold, sgndnold
 double precision, dimension(1:4,1:4) :: Rupinv, Rdninv
 double precision, dimension(1:4,1:4) :: Rup, Rdn
 double precision detRup, detRdn, detupold, detdnold
 double precision diiup, dininup, di2ni2nup, di3ni3nup
 double precision diidn, dinindn, di2ni2ndn, di3ni3ndn
 double precision del1, del2, del3, del4
 double precision, dimension(1:4,0:N-1) :: Vup, Vdn
 double precision, dimension(0:N-1,1:4) :: Uup, Udn
 double precision, dimension(0:N-1,0:N-1) :: eye

 eye = 0.0d0
 do i = 0,N-1
  eye(i,i) = 1.0d0
 enddo

 call generate_random_site_table(list,Nlat)
 do ii = 0,Nlat-1
 do si=0,19
  i = list(ii)

  Sold(:) = S(i,ti,:)
  Snew(:) = S(i,ti,:)
 
  j=ran2(iran)*20+1
  if(j>20) print*, "error in field"

  Snew(j) =-Snew(j)
  if(j>16) then
     ip=i-1
	 if(ip<0) ip=ip+Nlat
  else
     ip=i
  endif

  dS = Snew(:)-Sold(:)

  !Decide if we accept the change
  !First calculate the elements of the delta matrix
  if(j<=16) then
     diiup     = exp( lambda*dS(1) + lambdapp*(dS(5) + dS(6) + dS(9) + dS(10) + dS(13) + dS(14)))-1.0d0
     dininup    = exp( lambdap*dS(2) - lambdapp*(dS(5) + dS(7)) )-1.0d0
     di2ni2nup = exp( lambdap*dS(3) - lambdapp*(dS(9) + dS(11)))-1.0d0
     di3ni3nup = exp( lambdap*dS(4) - lambdapp*(dS(13)+dS(15)))-1.0d0

     diidn     = exp(-lambda*dS(1) + lambdapp*(dS(7) + dS(8) + dS(11) + dS(12) + dS(15) + dS(16)))-1.0d0
     dinindn    = exp(-lambdap*dS(2) - lambdapp*(dS(6)+dS(8)))-1.0d0
     di2ni2ndn = exp(-lambdap*dS(3) - lambdapp*(dS(10)+dS(12)))-1.0d0
     di3ni3ndn  = exp(-lambdap*dS(4) - lambdapp*(dS(14)+dS(16)))-1.0d0
  else
  !   print*, "j = ", j
     diiup = exp(lambdapp*(dS(17)+dS(18)))-1.0d0
	 dininup = exp(-lambdapp*(dS(17)+dS(19)))-1.0d0
	 di2ni2nup=0.0d0
	 di3ni3nup=0.0d0
	 
	 diidn = exp(lambdapp*(dS(19)+dS(20)))-1.0d0
	 dinindn = exp(-lambdapp*(dS(18)+dS(20)))-1.0d0
	 di2ni2ndn=0.0d0
	 di3ni3ndn=0.0d0
  endif

 
  !define the V matrix
  Vup(1,:) = -Gup(i,:)
  Vup(2,:) = -Gup(ip+Nlat,:)
  Vup(3,:) = -Gup(i+2*Nlat,:)
  Vup(4,:) = -Gup(i+3*Nlat,:)
  Vup(1,i) = 1.0d0 - Gup(i,i)
  Vup(2,ip+nlat) = 1.0d0 - Gup(ip+nlat,ip+nlat)
  Vup(3,i+2*nlat) = 1.0d0 - Gup(i+2*nlat,i+2*nlat)
  Vup(4,i+3*nlat) = 1.0d0 - Gup(i+3*nlat,i+3*nlat)
  Uup = 0.0d0
  Uup(i,1) = diiup
  Uup(ip+Nlat,2) = dininup
  Uup(i+2*nlat,3) = di2ni2nup
  Uup(i+3*nlat,4) = di3ni3nup

  Vdn(1,:) = -Gdn(i,:)
  Vdn(2,:) = -Gdn(ip+Nlat,:)
  Vdn(3,:) = -Gdn(i+2*Nlat,:)
  Vdn(4,:) = -Gdn(i+3*Nlat,:)
  Vdn(1,i) = 1.0d0 - Gdn(i,i)
  Vdn(2,ip+nlat) = 1.0d0 - Gdn(ip+nlat,ip+nlat)
  Vdn(3,i+2*nlat) = 1.0d0 - Gdn(i+2*nlat,i+2*nlat)
  Vdn(4,i+3*nlat) = 1.0d0 - Gdn(i+3*nlat, i+3*nlat)
  Udn = 0.0d0
  Udn(i,1) = diidn
  Udn(ip+Nlat,2) = dinindn
  Udn(i+2*nlat,3) = di2ni2ndn
  Udn(i+3*nlat,4) = di3ni3ndn

  Rup = matmul(Vup,Uup)
  Rdn = matmul(Vdn,Udn)
  do iii = 1,4
   Rup(iii,iii) = 1.0d0+Rup(iii,iii)
   Rdn(iii,iii) = 1.0d0+Rdn(iii,iii)
  enddo

  call determinant(Rup,detrup,rtmp,4)
  call determinant(Rdn,detrdn,rtmp,4)

  ratio = abs(detRup*detRdn)

  !calculate the probability for acceptance
  if(ratio.le.1.0d0)then
   p = ratio/(1.0d0 + gam*ratio)
  else
   p = ratio/(gam+ratio)
  endif

   !do we accept this probability
  r = ran2(iran)
  if(r.lt.p)then
   accept = accept + 1

   if(detrup.lt.0.0d0) sgnup=-sgnup
   if(detrdn.lt.0.0d0) sgndn=-sgndn
   call invertr(Rup,Rupinv,detup,rtmp,4)
   call invertr(Rdn,Rdninv,detdn,rtmp,4)

   !update the green's function
   gupold = gup 
   gup = matmul(gupold,eye-matmul(Uup,matmul(Rupinv,Vup))) 

   gdnold = gdn
   gdn = matmul(gdnold,eye-matmul(Udn,matmul(Rdninv,Vdn))) 

   !update the fields and the interaction matrix
   S(i,ti,:) = Snew(:)
   v_up(i,ti) = (diiup + 1.0d0)*v_up(i,ti)
   v_up(ip+nlat,ti) = (dininup + 1.0d0)*v_up(ip+Nlat,ti)
   v_up(i+2*nlat,ti) = (di2ni2nup + 1.0d0)*v_up(i+2*Nlat,ti)
   v_up(i+3*nlat,ti) = (di3ni3nup + 1.0d0)*v_up(i+3*Nlat,ti)

   v_dn(i,ti) = (diidn + 1.0d0)*v_dn(i,ti)
   v_dn(ip+nlat,ti) = (dinindn + 1.0d0)*v_dn(ip+Nlat,ti) 
   v_dn(i+2*nlat,ti) = (di2ni2ndn + 1.0d0)*v_dn(i+2*Nlat,ti)
   v_dn(i+3*nlat,ti) = (di3ni3ndn + 1.0d0)*v_dn(i+3*Nlat,ti)

  else
   reject = reject + 1
  endif 
 enddo
 enddo

 return
 end subroutine swpslice_woodbury

 !=============================================================================
 subroutine zero_counters()
 implicit none
 accept_hs = 0
 reject_hs = 0
 redo = 0
 noredo = 0
 accept_blk_hs = 0
 reject_blk_hs = 0
 return
 end subroutine zero_counters
 !=============================================================================
 ! This routine produces a random permutation of the array 0,1,2,3,...N-1
 !=============================================================================
 subroutine generate_random_site_table(list,N)
 use random
 implicit none
 integer i, j, n, k, a, b
 integer, dimension(0:N-1) :: list

 !init the list to an ordered list.
 do i = 0,N-1
  list(i) = i
 enddo

 do i = 1,3*N
  j = int(ran2(iran)*N)
  k = int(ran2(iran)*N)
  a = list(j)
  b = list(k)
  call iswap(a,b)
  list(j) = a
  list(k) = b
 enddo

 return
 end subroutine generate_random_site_table
 !=============================================================================
 ! this subroutine swaps the values of two integers.
 !=============================================================================
 subroutine iswap(a,b)
 implicit none
 integer a,b,tmp
 tmp = b
 b = a
 a = tmp
 return
 end subroutine iswap 
!===========================================================================
  subroutine block_update_HS(gup,gdn,sgnup,sgndn)
 use parameters, only: N, Nlat, L, lambda, lambdap, lambdapp, spinup, spindn,&
                    Block_site_number, Block_spin_number
 use fields, only: S, V_up, V_dn, return_v_for_S
 use grfun
 use linear_algebra
 use random
 implicit none
 logical update
 integer ii, i, ti, tip, j, jj, i1, j1
 double precision r, p, ratio
 double precision sgnup, sgndn, sgnup_old, sgndn_old, detup, detdn
 double precision detup_old, detdn_old
 integer, dimension(0:Nlat-1) :: list
 double precision, dimension(0:N-1,0:N-1) :: gup, gdn, gup_old, gdn_old
 double precision, dimension(0:N-1,0:L-1) :: v_up_old, v_dn_old
 double precision, dimension(0:Nlat-1,0:L-1,1:20) :: S_old

 call generate_random_site_table(list,Nlat)
 do ii = 0,Nlat-1

!pick a site at random
  i = list(ii)
  S_old = S
  v_up_old = v_up
  v_dn_old = v_dn

  tip = 0
  call getgp(gup_old,tip,detup_old,sgnup_old,spinup)
  call getgp(gdn_old,tip,detdn_old,sgndn_old,spindn)

  !pick a field at random
   j = ran2(iran)*20 + 1
  S(i,:,j) = -S(i,:,j)
  call return_v_for_s(S,v_up,v_dn)

  !calculate the new green's function
  call getgp(gup,tip,detup,sgnup,spinup)
  call getgp(gdn,tip,detdn,sgndn,spindn)

  !determine the change in the determinant
  p = abs(exp(detup+detdn-detup_old-detdn_old)) 

  !draw a random number 
  r = ran2(iran)
  if(r.le.p)then
   accept_blk_hs = accept_blk_hs + 1
  else
   reject_blk_hs = reject_blk_hs + 1
   v_up = v_up_old
   v_dn = v_dn_old
   s = s_old
   sgnup = sgnup_old
   sgndn = sgndn_old
   gup = gup_old
   gdn = gdn_old
  endif
 enddo  !ii

 return
 end subroutine block_update_HS
end module updates

