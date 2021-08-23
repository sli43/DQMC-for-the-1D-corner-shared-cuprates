module measurements
 use parameters, only: N, Nlat, L, beta, Nhis
 integer nmeas, nmeast
 integer, parameter :: nbins = 10
 double precision asgn, asgnup, asgndn
 double precision an1, anup1, andn1
 double precision an2, anup2, andn2
 double precision an3, anup3, andn3
 double precision an4, anup4, andn4
 double precision an, anup,andn
 double precision anud1, anud2, anud3,anud4, amz1, amz2, amz3,amz4, asgnt
 double precision, dimension(1:nbins) :: bsgnt
 double precision, dimension(0:Nlat-1,0:L,1:4,1:4) :: agnl
 double precision, dimension(0:Nlat-1,0:L,1:4,1:4) :: achinl
 double precision, dimension(0:Nlat-1,0:L,1:4,1:4) :: achinlz
 double precision, dimension(0:Nlat-1,0:L,1:4,1:4) :: adent
 double precision, dimension(0:Nlat-1,1:4,1:4) :: agrfun
 double precision, dimension(0:Nlat-1,1:4,1:4) :: aspinxx, aspinzz
 double precision, dimension(0:Nlat-1,1:4,1:4) :: aden
 double precision, dimension(1:nbins,0:Nlat-1,0:L,1:4,1:4) :: bgnl
 double complex, dimension(1:nbins,0:Nlat-1,0:L,1:4,1:4) :: bgnlq
 double complex, dimension(1:nbins,0:Nlat-1,1:4,1:4):: qself_ener        ! new add
 double precision, dimension(1:nbins,0:Nlat-1,0:L,1:4,1:4) :: bchinl
 double precision, dimension(1:nbins,0:Nlat-1,0:L,1:4,1:4) :: bchinlz
 double complex, dimension(1:nbins,0:Nlat-1,0:L,1:4,1:4) :: bchinlq, bchinlzq
 double precision, dimension(1:nbins,0:Nlat-1,0:L,1:4,1:4) :: bdent
 double complex, dimension(1:nbins,0:Nlat-1,0:L,1:4,1:4) :: bdentq    ! new add
 double precision, dimension(1:nbins, 0:Nlat-1,1:4,1:4):: quasi_resid   !new add
 double precision, dimension(1:nbins,0:Nlat-1,1:4,1:4) :: bgrfun
 double complex, dimension(1:nbins,0:Nlat-1,1:4,1:4) :: bgkfun
 double precision, dimension(1:nbins,0:Nlat-1,1:4,1:4) :: bspinxx
 double precision, dimension(1:nbins,0:Nlat-1,1:4,1:4) :: bspinzz
 double complex, dimension(1:nbins,0:Nlat-1,1:4,1:4) :: bspinxxq,bspinzzq
 double precision, dimension(1:nbins,0:Nlat-1,1:4,1:4) :: bden
 double complex, dimension(1:nbins,0:Nlat-1,1:4,1:4) :: bdenq
 double precision, dimension(1:nbins) :: bone
 double precision, dimension(1:nbins) :: bn, bn1, bn2, bn3, bn4
 double precision, dimension(1:nbins) :: bmz1, bmz2, bmz3, bmz4
 double precision, dimension(1:nbins) :: bnud1, bnud2, bnud3, bnud4
 double precision, dimension(1:nbins) :: bnup, bnup1, bnup2, bnup3, bnup4
 double precision, dimension(1:nbins) :: bndn, bndn1, bndn2, bndn3, bndn4
 double precision, dimension(1:nbins) :: bsgn, bsgnup, bsgndn
 double precision ake, ape1, ape2,ape3, aee
 double precision, dimension(1:nbins) :: bke, bpe1, bpe2,bpe3,bee
 double precision ant1,ant2,ant3,ant4
 double precision, dimension(1:nbins):: bnt1, bnt2, bnt3, bnt4
contains
 subroutine meastau(sgnup,sgndn)
 use parameters, only: N, Nlat, L, spinup, spindn
 use grfun, only: makegt
 implicit none
 integer ti, tip, i, j, ii, jj, k, k1, orb1, orb2
 double precision sgnup, sgndn, sgn, rtmp, pi
 double precision nt1, nt2, nt3, nt4
 double precision, dimension(0:N-1,0:N-1,0:L) :: G
 double precision, dimension(0:Nlat-1,0:L,1:4,1:4) :: gnl
 double precision, dimension(0:Nlat-1,0:L,1:4,1:4) :: chinl,chinlz
 double complex, dimension(0:Nlat-1,0:L,1:4,1:4) :: chinlq,chinlzq
 double precision, dimension(0:Nlat-1,0:L,1:4,1:4) :: dent
 double precision, dimension(0:N-1,0:N-1) :: gttup, gttdn
 double precision, dimension(0:N-1,0:N-1) :: gt0up, gt0dn
 double precision, dimension(0:N-1,0:N-1) :: g0tup, g0tdn
 double precision, dimension(0:N-1,0:N-1) :: g00up, g00dn, dummymat

 pi = 2.0d0*asin(1.0d0)
 G = 0.0d0
 chinl = 0.0d0
 chinlz= 0.0d0
 dent = 0.0d0

 ti = 0
 call makegt(ti,ti,g00up,dummymat,spinup)
 call makegt(ti,ti,g00dn,dummymat,spindn)

 do ti = 0,L
  tip = 0
  call makegt(mod(ti,L),tip,gt0up,g0tup,spinup)  
  call makegt(mod(ti,L),tip,gt0dn,g0tdn,spindn) 
  call makegt(mod(ti,L),mod(ti,L),gttup,dummymat,spinup)  
  call makegt(mod(ti,L),mod(ti,L),gttdn,dummymat,spindn) 
  do i = 0,N-1
   do j = 0,N-1
    g(i,j,ti) = g(i,j,ti) + 0.5d0*(gt0up(i,j)+gt0dn(i,j)) 
   enddo
  enddo
 
  !get the dynamic spin and charge susceptibilities
  do k = 0,Nlat-1
    do i=0,Nlat-1
       j=mod(i+k,Nlat)     
      do orb1 = 1,4
        ii = i + nlat*(orb1-1)
      do orb2 = 1,4
      jj = j + nlat*(orb2-1)

      chinl(k,ti,orb1,orb2) = chinl(k,ti,orb1,orb2) - g0tup(jj,ii)*gt0dn(ii,jj) &
                              - g0tdn(jj,ii)*gt0up(ii,jj)

      chinlz(k,ti,orb1,orb2) = chinlz(k,ti,orb1,orb2)-(g00up(ii,ii)-g00dn(ii,ii))*(gttup(jj,jj)-gttdn(jj,jj)) &
                                - g0tup(jj,ii)*gt0up(ii,jj) - g0tdn(jj,ii)*gt0dn(ii,jj)

      dent(k,ti,orb1,orb2) = dent(k,ti,orb1,orb2) &
                             - (2.0d0-gttup(ii,ii)-gttdn(ii,ii))*(2.0d0-g00up(jj,jj)-g00dn(jj,jj)) &
                             - gt0up(ii,jj)*g0tup(jj,ii) - gt0dn(ii,jj)*g0tdn(jj,ii)
      enddo !orb1
     enddo  !orb2
   end do !i
  end do  !k
enddo


 !build the proper green's function
 gnl = 0.0d0
 do k = 0,Nlat-1
  do i=0,Nlat-1
    j=mod(i+k,Nlat)
   do orb1 = 0,3
    ii = i + orb1*Nlat
    do orb2 = 0,3
     jj = j + orb2*Nlat
     do ti = 0,L
      gnl(k,ti,orb1+1,orb2+1) = gnl(k,ti,orb1+1,orb2+1)+g(ii,jj,ti)
     enddo
    enddo
   enddo
 end do
 enddo

chinl=chinl/dble(Nlat)
chinlz=chinlz/dble(nlat)
dent=dent/dble(nlat)
gnl=gnl/dble(nlat)


 gnl(:,L,:,:) =-gnl(:,0,:,:)
 gnl(0,L,1,1) =1.0d0 - gnl(0,0,1,1)
 gnl(0,L,2,2) =1.0d0 - gnl(0,0,2,2)
 gnl(0,L,3,3) =1.0d0 - gnl(0,0,3,3)
 gnl(0,L,4,4) =1.0d0 - gnl(0,0,4,4)


 nt1=0.0d0
 nt2=0.0d0
 nt3=0.0d0
 nt4=0.0d0
 nt1=(1.0d0-gnl(0,0,1,1))*2.0d0
 nt2=(1.0d0-gnl(0,0,2,2))*2.0d0
 nt3=(1.0d0-gnl(0,0,3,3))*2.0d0
 nt4=(1.0d0-gnl(0,0,4,4))*2.0d0



 !store things in the accumulator
 nmeast = nmeast + 1
 sgn = sgnup*sgndn
 asgnt = asgnt + sgn 
 agnl = agnl + gnl*sgn
 achinl = achinl + chinl*sgn
 achinlz = achinlz + chinlz*sgn
 adent = adent + dent*sgn

 ant1=ant1+nt1*sgn
 ant2=ant2+nt2*sgn
 ant3=ant3+nt3*sgn
 ant4=ant4+nt4*sgn

 return
 end subroutine meastau

 !=============================================================================
 ! Subroutine equaltime_measurements
 !=============================================================================
 subroutine equaltime_measurements(gup,gdn,sgnup,sgndn)
 use parameters, only: N, Nlat, dtau, L, ttk, Ud, Up, Upd
 implicit none
 integer i, j, ii, jj,jjj, k, k1,ti,tip,iA,iB,jA,jB
 integer orb1, orb2
 double precision sgn, sgnup, sgndn, nup1, nup2, n1, n2, ndn1, ndn2
 double precision n3, nup3, ndn3, nud3, n4, nup4, ndn4, nud4
 double precision nud1, nud2, norm, mz1, mz2, mz3, mz4, norm2
 double precision, dimension(0:Nlat-1,1:4,1:4) :: spinzz, spinxx
 double precision, dimension(0:Nlat-1,1:4,1:4) :: den
 double precision, dimension(0:N-1,0:N-1) :: gup, gdn
 double precision, dimension(0:Nlat-1,1:4,1:4) :: grfun
 double precision xxtemp, zztemp, rtmp
 double precision pe1, pe2,pe3, ke, ee


 norm = 1.0d0/dble(Nlat)
 norm2 = 1.0d0/dble(Nlat*L)
 ndn1 = 0.0d0
 nup1 = 0.0d0
 nud1 = 0.0d0
 n1 = 0.0d0
 ndn2 = 0.0d0
 nup2 = 0.0d0
 nud2 = 0.0d0
 n2 = 0.0d0
 ndn3 = 0.0d0
 nup3 = 0.0d0
 nud3 = 0.0d0
 n3 = 0.0d0
 ndn4=0.0d0
 nup4=0.0d0
 nud4=0.0d0
 n4=0.0d0
 mz1 = 0.0d0
 mz2 = 0.0d0
 mz3 = 0.0d0
 mz4=0.0d0


 ke=0.0d0
 pe1=0.0d0
 pe2=0.0d0
 pe3=0.0d0
  
 sgn = sgnup*sgndn
! measure kinetic energy and potential
 do i=0,N-1
   do j=0, N-1
     if(i.ne.j) then
          ke=ke-gup(i,j)*ttk(i,j)-gdn(i,j)*ttk(i,j)
     else
          ke = ke + (1.0d0-gup(i,j))*ttk(i,j)+(1.0d0-gdn(i,j))*ttk(i,j)
     end if
   end do
 end do

do i=0,Nlat-1
  pe1=pe1+Ud*(1.0d0-gup(i,i))*(1.0d0-gdn(i,i))
  pe2=pe2+Up*(1.0d0-gup(i+nlat,i+nlat))*(1.0d0-gdn(i+nlat,i+nlat))+&
      Up*(1.0d0-gup(i+2*nlat,i+2*nlat))*(1.0d0-gdn(i+2*nlat,i+2*nlat))+&
      Up*(1.0d0-gup(i+3*nlat,i+3*nlat))*(1.0d0-gdn(i+3*nlat,i+3*nlat))
  pe3=pe3+Upd*(2.0d0-gup(i,i)-gdn(i,i))*(6.0d0-gup(i+nlat,i+nlat)-gdn(i+nlat,i+nlat)&
             - gup(i+2*nlat,i+2*nlat) - gdn(i+2*nlat,i+2*nlat) - gup(i+3*nlat,i+3*nlat) &
             - gdn(i+3*nlat,i+3*nlat))
end do

ee=ke+pe1+pe2+pe3


 do i = 0,Nlat-1
  j = i + Nlat
  jj= i + 2*nlat
  jjj=i+3*Nlat

  n1 = n1 + (1.0d0 - Gup(i,i)) + (1.0d0 - Gdn(i,i))
  nup1 = nup1 + (1.0d0 - Gup(i,i))
  ndn1 = ndn1 + (1.0d0 - Gdn(i,i))
  nud1 = nud1 + (1.0d0 - Gup(i,i))*(1.0d0 - Gdn(i,i))
  mz1 = ((1.0d0-gup(i,i)) - (1.0d0 - gdn(i,i)))**2.0d0

  n2 = n2 + (1.0d0 - Gup(j,j)) + (1.0d0 - Gdn(j,j))
  nup2 = nup2 + (1.0d0 - Gup(j,j))
  ndn2 = ndn2 + (1.0d0 - Gdn(j,j))
  nud2 = nud2 + (1.0d0 - Gup(j,j))*(1.0d0 - Gdn(j,j))
  mz2 = ((1.0d0-gup(j,j)) - (1.0d0 - gdn(j,j)))**2.0d0

  n3 = n3 + (1.0d0 - Gup(jj,jj)) + (1.0d0 - Gdn(jj,jj))
  nup3 = nup3 + (1.0d0 - Gup(jj,jj))
  ndn3 = ndn3 + (1.0d0 - Gdn(jj,jj))
  nud3 = nud3 + (1.0d0 - Gup(jj,jj))*(1.0d0 - Gdn(jj,jj))
  mz3 = ((1.0d0-gup(jj,jj)) - (1.0d0 - gdn(jj,jj)))**2.0d0

  n4 = n4 + (1.0d0 - Gup(jjj,jjj)) + (1.0d0 - Gdn(jjj,jjj))
  nup4 = nup4 + (1.0d0 - Gup(jjj,jjj))
  ndn4 = ndn4 + (1.0d0 - Gdn(jjj,jjj))
  nud4 = nud4 + (1.0d0 - Gup(jjj,jjj))*(1.0d0 - Gdn(jjj,jjj))
  mz4 = ((1.0d0-gup(jjj,jjj)) - (1.0d0 - gdn(jjj,jjj)))**2.0d0

 enddo


 den=0.0d0
 spinxx=0.0d0
 spinzz=0.0d0
 do k=0,Nlat-1
  do i=0,Nlat-1
    j=mod(i+k,Nlat)

  do orb1=1,4
   do orb2=1,4
      ii=i+nlat*(orb1-1)
      jj=j+nlat*(orb2-1)

      if(ii.eq.jj) then
        spinxx(k,orb1,orb2) = spinxx(k,orb1,orb2)+ 2.0d0-gup(ii,jj)-gdn(ii,jj)-&
                            2.0d0*(1.0d0-gup(ii,jj))*(1.0d0-gdn(ii,jj))

        spinzz(k,orb1,orb2) = spinzz(k,orb1,orb2) +2.0d0-gup(ii,jj)-gdn(ii,jj)-&
                              2.0d0*(1.0d0-gup(ii,jj))*(1.0d0-gdn(ii,jj))

        den(k,orb1,orb2) = den(k,orb1,orb2)+2.0d0-gup(ii,ii)-gdn(ii,ii) &
                           +2.0d0*(1.0d0-gup(ii,ii))*(1.0d0-gdn(ii,ii))
     else
       spinxx(k,orb1,orb2) = spinxx(k,orb1,orb2)-gup(ii,jj)*gdn(jj,ii)-gup(jj,ii)*gdn(ii,jj)

       spinzz(k,orb1,orb2) = spinzz(k,orb1,orb2)+gup(ii,ii)*gup(jj,jj)+gdn(ii,ii)*gdn(jj,jj)&
                             -gdn(ii,ii)*gup(jj,jj)-gup(ii,ii)*gdn(jj,jj)-gup(ii,jj)*gup(jj,ii)&
                             -gdn(ii,jj)*gdn(jj,ii)

       den(k,orb1,orb2) = den(k,orb1,orb2) +(2.0d0-gup(ii,ii)-gdn(ii,ii))*&
                         (2.0d0-gup(jj,jj)-gdn(jj,jj)) &
                         -gup(jj,ii)*gup(ii,jj)-gdn(jj,ii)*gdn(ii,jj)
      endif

   enddo
  enddo

 enddo
end do

spinxx=spinxx/4.0d0/dble(nlat)
spinzz=spinzz/4.0d0/dble(nlat)
den=den/dble(nlat)



 !measure the equaltime green's function
 grfun(:,:,:) = 0.0d0
 do k = 0,Nlat-1
  do i=0,Nlat-1
    j=mod(i+k,Nlat)
   do orb1 = 0,3
    ii = i + nlat*orb1
    do orb2 = 0,3
     jj = j + nlat*orb2
     grfun(k,orb1+1,orb2+1) = grfun(k,orb1+1,orb2+1) + (Gup(ii,jj)+Gdn(ii,jj))/2.0d0
    enddo
   enddo
 enddo
end do

grfun=grfun/dble(nlat)


 nmeas = nmeas + 1

 !Store all the values in the accumulators 
 asgn = asgn + sgn
 asgnup = asgnup + sgnup
 asgndn = asgndn + sgndn

 agrfun = agrfun + grfun*sgn
 aspinxx = aspinxx + spinxx*sgn
 aspinzz = aspinzz + spinzz*sgn
 aden = aden + den*sgn

 ake=ake+ke*sgn*norm
 ape1=ape1+pe1*sgn*norm
 ape2=ape2+pe2*sgn*norm
 ape3=ape3+pe3*sgn*norm
 aee=aee+ee*sgn*norm


 !filling for band 1
 an = an + (n1+n2+n3+n4)*sgn*norm
 anup=anup+(nup1+nup2+nup3+nup4)*sgn*norm
 andn=andn+(ndn1+ndn2+ndn3+ndn4)*sgn*norm
 an1 = an1 + n1*sgn*norm
 anup1 = anup1 + nup1*sgn*norm
 andn1 = andn1 + ndn1*sgn*norm
 anud1 = anud1 + nud1*sgn*norm
 amz1 = amz1 + mz1*sgn*norm

 !fillings for band 2
 an2 = an2 + n2*sgn*norm
 anup2 = anup2 + nup2*sgn*norm
 andn2 = andn2 + ndn2*sgn*norm
 anud2 = anud2 + nud2*sgn*norm
 amz2 = amz2 + mz2*sgn*norm
 !fillings for band 3
 an3 = an3 + n3*sgn*norm
 anup3 = anup3 + nup3*sgn*norm
 andn3 = andn3 + ndn3*sgn*norm
 anud3 = anud3 + nud3*sgn*norm
 amz3 = amz3 + mz3*sgn*norm

!fillings for band 4
 an4 = an4 + n4*sgn*norm
 anup4 = anup4 + nup4*sgn*norm
 andn4 = andn4 + ndn4*sgn*norm
 anud4 = anud4 + nud4*sgn*norm
 amz4 = amz4 + mz4*sgn*norm

 return
 end subroutine equaltime_measurements

 !=============================================================================
 subroutine zero_accumulators()
 implicit none
 nmeas = 0
 nmeast = 0
 asgnt = 0.0d0
 asgn = 0.0d0
 asgnup = 0.0d0
 asgndn = 0.0d0
 an = 0.0d0
 anup=0.0d0
 andn=0.0d0
 an1 = 0.0d0
 anup1 = 0.0d0
 andn1 = 0.0d0
 anud1 = 0.0d0
 an2 = 0.0d0
 anup2 = 0.0d0
 andn2 = 0.0d0
 anud2 = 0.0d0
 an3 = 0.0d0
 anup3 = 0.0d0
 andn3 = 0.0d0
 anud3 = 0.0d0
 an4 = 0.0d0
 anup4 = 0.0d0
 andn4 = 0.0d0
 anud4 = 0.0d0
 amz1 = 0.0d0
 amz2 = 0.0d0
 amz3 = 0.0d0
 amz4=0.0d0

 agrfun = 0.0d0
 aspinxx = 0.0d0
 aspinzz = 0.0d0
 aden = 0.0d0
 agnl = 0.0d0
 achinl = 0.0d0
 achinlz = 0.0d0
 adent = 0.0d0
 ake=0.0d0
 ape1=0.0d0
 ape2=0.0d0
 ape3=0.0d0
 aee=0.0d0

 ant1=0.0d0
 ant2=0.0d0
 ant3=0.0d0
 ant4=0.0d0

 return
 end subroutine zero_accumulators


 subroutine populate_bins(k)
 use parameters, only: N, Nlat, L,dtau,beta,tt
 use mpi_module
 implicit none
 integer k
 integer i, ti, ii, jj
 integer j, orb1, orb2
 double precision pi
 integer info
 double complex, dimension(1:24)::work
 double complex, dimension(0:Nlat-1,1:4,1:4)::exilen
 double complex, dimension(1:4,1:4):: orimatrix,invermatrix
 double complex, dimension(1:4,1:4)::identity
 integer, dimension(1:4):: ipiv

 if(myrank.eq.root) print*, 'Populating bin ',k
 bone(k) = 1.0d0
 if(nmeas.eq.0) nmeas = 1
 bsgn(k) = asgn/dfloat(nmeas)
 bsgnup(k) = asgnup/dfloat(nmeas)
 bsgndn(k) = asgndn/dfloat(nmeas)

 bke(k)=ake/dfloat(nmeas)
 bpe1(k)=ape1/dfloat(nmeas)
 bpe2(k)=ape2/dfloat(nmeas)
 bpe3(k)=ape3/dfloat(nmeas)
 bee(k)=aee/dfloat(nmeas)

 bn(k) = an/dfloat(nmeas)
 bnup(k)=anup/dfloat(nmeas)
 bndn(k)=andn/dfloat(nmeas)

 bnud1(k) = anud1/dfloat(nmeas)
 bn1(k) = an1/dfloat(nmeas)
 bnup1(k) = anup1/dfloat(nmeas)
 bndn1(k) = andn1/dfloat(nmeas)

 bnud2(k) = anud2/dfloat(nmeas)
 bn2(k) = an2/dfloat(nmeas)
 bnup2(k) = anup2/dfloat(nmeas)
 bndn2(k) = andn2/dfloat(nmeas)

 bnud3(k) = anud3/dfloat(nmeas)
 bn3(k) = an3/dfloat(nmeas)
 bnup3(k) = anup3/dfloat(nmeas)
 bndn3(k) = andn3/dfloat(nmeas)

 bnud4(k) = anud4/dfloat(nmeas)
 bn4(k) = an4/dfloat(nmeas)
 bnup4(k) = anup4/dfloat(nmeas)
 bndn4(k) = andn4/dfloat(nmeas)

 bmz1(k) = amz1/dfloat(nmeas)
 bmz2(k) = amz2/dfloat(nmeas)
 bmz3(k) = amz3/dfloat(nmeas)
 bmz4(k) = amz4/dfloat(nmeas)


 bgrfun(k,:,:,:) = agrfun(:,:,:)/dfloat(nmeas)
 bspinxx(k,:,:,:) = aspinxx(:,:,:)/dfloat(nmeas)
 bspinzz(k,:,:,:) = aspinzz(:,:,:)/dfloat(nmeas)
 bden(k,:,:,:) = aden(:,:,:)/dfloat(nmeas)

 call ftntok(bspinxx(k,:,:,:),bspinxxq(k,:,:,:),Nlat)
 call ftntok(bspinzz(k,:,:,:),bspinzzq(k,:,:,:),Nlat)
 call ftntok(bden(k,:,:,:),bdenq(k,:,:,:),Nlat)
 call ftntok(bgrfun(k,:,:,:),bgkfun(k,:,:,:),Nlat)

 bsgnt(k) = asgnt/dfloat(nmeast)
 bgnl(k,:,:,:,:) = agnl(:,:,:,:)/dfloat(nmeast)
 bchinl(k,:,:,:,:) = achinl(:,:,:,:)/dfloat(nmeast)
 bchinlz(k,:,:,:,:) = achinlz(:,:,:,:)/dfloat(nmeast)
 bdent(k,:,:,:,:) = adent(:,:,:,:)/dfloat(nmeast)


 do ti=0,L
  call ftntok(bgnl(k,:,ti,:,:),bgnlq(k,:,ti,:,:),Nlat)
  call ftntok(bdent(k,:,ti,:,:),bdentq(k,:,ti,:,:),Nlat)
  call ftntok(bchinl(k,:,ti,:,:),bchinlq(k,:,ti,:,:),Nlat)
  call ftntok(bchinlz(k,:,ti,:,:),bchinlzq(k,:,ti,:,:),Nlat)
end do !ti

bnt1(k)=ant1/dble(nmeast)
bnt2(k)=ant2/dble(nmeast)
bnt3(k)=ant3/dble(nmeast)
bnt4(k)=ant4/dble(nmeast)



!call ftntok(ttn,exilen,Nlat)
!identity=dcmplx(0.0d0,0.0d0)
!do i=1,4
!  identity(i,i)=dcmplx(0.0d0,1.0d0)
!end do

!pi=2.0d0*asin(1.0d0)
!orimatrix=dcmplx(0.0d0,0.0d0)
!invermatrix=dcmplx(0.0d0,0.0d0)

! call ftltow(bgnlq(k,:,:,:,:),qself_ener(k,:,:,:))

!do i=0,Nlat-1
!orimatrix=qself_ener(k,i,:,:)
!call ZGETRF(4,4,orimatrix,4,IPIV,info)
!call ZGETRI(4,orimatrix,4,IPIV,work,24,info)

!qself_ener(k,i,:,:)=pi/beta*identity(:,:) - exilen(i,:,:) + orimatrix(:,:)
!end do


 return
 end subroutine populate_bins
 !============================================================================

 !rebinned jackknife estimation. 10 bins from each process are added to form a
 !data point. does not work for single processor jobs as in that case there would be
 !only one data point. probably at least four processors should be used to get
 !meaningful error estimates. this method does not underestimate the error since
 !separate processes are completely independent (assuming the rng works well)

 subroutine geterr_new(dat, sgn, N, mean, std)
 implicit none
 integer N
 double precision mean, std
 double precision, dimension(1:N) :: dat, sgn
 double precision, dimension(1:10,1:N/10) :: dat_blocked, sgn_blocked
 double precision, dimension(1:N/10) :: dat_avgs, sgn_avgs !=10*avg per process
 double precision, dimension(1:N/10) :: x_j
 double precision sum_dat, sum_sgn, F

 dat_blocked = reshape(dat, shape(dat_blocked))
 sgn_blocked = reshape(sgn, shape(sgn_blocked))
 dat_avgs = sum(dat_blocked, DIM=1) !divide by 10 to get actual averages
 sgn_avgs = sum(sgn_blocked, DIM=1)

 sum_dat = sum(dat_avgs)
 sum_sgn = sum(sgn_avgs)

 x_j = (sum_dat - dat_avgs)/(sum_sgn - sgn_avgs) !jackknife averages
 mean = sum(x_j)/dfloat(N/10)
 std = sqrt(sum((x_j-mean)**2) * dfloat((N/10)-1)/dfloat(N/10))
 !correct for O(1/N) bias
 mean = (dfloat(N/10) * sum_dat / sum_sgn) - (dfloat((N/10)-1) * mean)

 !call getF(dat_blocked, 10, N/10, F)
 !if (F.gt.2.1) std = -std

 return
 end subroutine geterr_new
 !====================================================================
 subroutine geterr(dat,sgn,N,mean,std)
 implicit none
 integer N
 double precision, parameter :: tiny = 1e-8
 double precision mean, std
 double precision, dimension(1:N) :: dat, sgn, X
 integer i

 !S. johnston, 29 Aug 2013. I inserted an if statement to set sgn = tiny 
 !if it is less than tiny.  This is to avoid an IEEE
 !exception for any division by zero.  The bin should only be zero if no
 !measurements were done. 
 do i=1,N
  if(sgn(i).le.tiny)then
   sgn(i) = tiny
  endif
 enddo
 X = dat/sgn
 mean = sum(x)/dfloat(N)
 std = sqrt(sum((x-mean)**2)/dfloat(n-1))
 return
 end subroutine geterr
!=====================================================================
 subroutine geterr1(dat,N,mean)
 implicit none
 integer N
 double precision, parameter :: tiny = 1e-8
 double precision mean
 double precision, dimension(1:N) :: dat
 integer i,j

 mean=0
 mean=sum(dat(:))

 return
 end subroutine geterr1
!======================================================================
! This is a fourier transform routine lifted directly from Richard's 
! oringal code. 
!
! This routine transforms from n to k.  Ie from real space to momentum
! space. 
!======================================================================
      subroutine ftntok(gn,gk,N)
      implicit none
      integer N
      integer nkx, nx, lx, nn, nqx
      double precision kx, pi, arg 
      double precision, dimension(0:N-1,1:4,1:4) :: gn
      double complex, dimension(0:N-1,1:4,1:4):: gk
      double complex eye

      !constants
      pi = 2.0d0*asin(1.0d0)
      gk = dcmplx(0.0d0,0.0d0)
      eye=dcmplx(0.0d0,1.0d0)
      nn = N*N
      !loop over the allowed values of K
      ! the momentum changes from 0 to 2pi
      do nkx = 0,N-1
         kx = 2.0d0*pi*dfloat(nkx)/dfloat(N)
        !loop over position.  Recall that the input is the translational invarient
        !version of G(n,tau)
        do nx = 0,N-1
           arg = dfloat(nx)*kx
           gk(nkx,1,1)=gk(nkx,1,1)+exp(eye*arg)*gn(nx,1,1)
           gk(nkx,2,2)=gk(nkx,2,2)+exp(eye*arg)*gn(nx,2,2)
           gk(nkx,3,3)=gk(nkx,3,3)+exp(eye*arg)*gn(nx,3,3)
           gk(nkx,4,4)=gk(nkx,4,4)+exp(eye*arg)*gn(nx,4,4)

           gk(nkx,1,3)=gk(nkx,1,3)+exp(eye*arg)*gn(nx,1,3)
           gk(nkx,3,1)=gk(nkx,3,1)+exp(eye*arg)*gn(nx,3,1)

           gk(nkx,1,4)=gk(nkx,1,4)+exp(eye*arg)*gn(nx,1,4)
           gk(nkx,4,1)=gk(nkx,4,1)+exp(eye*arg)*gn(nx,4,1)

           gk(nkx,3,4)=gk(nkx,3,4)+exp(eye*arg)*gn(nx,3,4)
           gk(nkx,4,3)=gk(nkx,4,3)+exp(eye*arg)*gn(nx,4,3)


           arg= dfloat(nx)*kx - kx*0.5d0
           gk(nkx,1,2)=gk(nkx,1,2)+exp(eye*arg)*gn(nx,1,2)
           gk(nkx,3,2)=gk(nkx,3,2)+exp(eye*arg)*gn(nx,3,2)
           gk(nkx,4,2)=gk(nkx,4,2)+exp(eye*arg)*gn(nx,4,2)

           arg= dfloat(nx)*kx + kx*0.5d0
           gk(nkx,2,1)=gk(nkx,2,1)+exp(eye*arg)*gn(nx,2,1)
           gk(nkx,2,3)=gk(nkx,2,3)+exp(eye*arg)*gn(nx,2,3)
           gk(nkx,2,4)=gk(nkx,2,4)+exp(eye*arg)*gn(nx,2,4)



        end do !nx
     end do !nkx
       
      return
 end subroutine ftntok

!====================================================================================
      subroutine ftltow(gt,gw)
      use parameters, only: Nlat, L, dtau, beta
      implicit none
      integer, parameter :: nsteps = 10000
      integer nn, i, j, nx, ny, ti, ii, maxn
      double precision pi, wn, tau, M, B, arg
      double precision ddtau, xx, im, re
      double precision Reyp1, Reypn, Imyp1, Imypn
      double precision, dimension(0:L) :: Rey, Imy, Rey2, Imy2, X
      double precision, dimension(0:Nlat-1,0:L,1:4,1:4) :: gt
      double complex, dimension(0:Nlat-1,1:4,1:4) :: gw
      double complex cfac1, cfac2, f0, f1, f2
      integer io1,io2

      ddtau = beta/dfloat(nsteps*L)
      gw = dcmplx(0.0d0, 0.0d0)
      pi = 2.0d0*asin(1.0d0)

      do io1=1,4
      do io2=1,4
      do nx = 0,Nlat-1
        
        do nn = 0,0
         wn = dfloat(2*nn+1)*pi/beta

         Rey = 0.0d0
         Imy = 0.0d0

         do ti = 0,L-1
          tau = dfloat(ti)*dtau
          X(ti) = tau

          if(ti.lt.L)then
           cfac2 = dcmplx(gt(nx,ti,io1,io2), 0.0d0)
          else
           if(io1.eq.io2) then
             cfac2 = dcmplx(1.0d0 - gt(nx,0,io1,io2), 0.0d0)
           else 
             cfac2=dcmplx(-gt(nx,0,io1,io2),0.0d0)
           end if
          endif

          if(abs(dreal(cfac2)).gt.10e-12) Rey(ti) = dreal(cfac2)
          if(abs(dimag(cfac2)).gt.10e-12) Imy(ti) = dimag(cfac2)
         enddo

         call spline(X,Rey,L+1,2d30,2d30,Rey2)
         call spline(X,Imy,L+1,2d30,2d30,Imy2)

         do ti = 1,nsteps*L-2,2
          tau = ddtau*dfloat(ti-1)
          cfac1 = dcmplx(cos(wn*tau),sin(wn*tau))
          call splint(x,rey,rey2,L+1,tau,Re)
          call splint(x,imy,imy2,L+1,tau,Im)
          cfac2 = dcmplx(Re,Im)
          f0 = cfac2*cfac1

          tau = ddtau*dfloat(ti)
          cfac1 = dcmplx(cos(wn*tau),sin(wn*tau))
          call splint(x,rey,rey2,L+1,tau,Re)
          call splint(x,imy,imy2,L+1,tau,Im)
          cfac2 = dcmplx(Re,Im)
          f1 = cfac2*cfac1

          tau = ddtau*dfloat(ti+1)
          cfac1 = dcmplx(cos(wn*tau),sin(wn*tau))
          call splint(x,rey,rey2,L+1,tau,Re)
          call splint(x,imy,imy2,L+1,tau,Im)
          cfac2 = dcmplx(Re,Im)
          f2 = cfac2*cfac1

          gw(nx,io1,io2) = gw(nx,io1,io2)+ddtau*(f0 + 4.0d0*f1 + f2)/3.0d0
         enddo
        enddo
      enddo
      end do
      end do
      return
      end subroutine ftltow
!=====================================================================
       SUBROUTINE SPLINE(X,Y,N,YP1,YPN,Y2)
       implicit none
       integer nmax,n
       PARAMETER (NMAX=400)
       double precision x(n),y(n),y2(n),u(nmax)
       double precision yp1,ypn,p,qn,sig,un
       integer i,k
       IF (YP1.GT..99D30) THEN
         Y2(1)=0.d0
         U(1)=0.d0
       ELSE
         Y2(1)=-0.5d0
         U(1)=(3.d0/(X(2)-X(1)))*((Y(2)-Y(1))/(X(2)-X(1))-YP1)
       ENDIF
       DO I=2,N-1
         SIG=(X(I)-X(I-1))/(X(I+1)-X(I-1))
         P=SIG*Y2(I-1)+2.d0
         Y2(I)=(SIG-1.d0)/P
         U(I)=(6.d0*((Y(I+1)-Y(I))/(X(I+1)-X(I))-(Y(I)-Y(I-1))&
              /(X(I)-X(I-1)))/(X(I+1)-X(I-1))-SIG*U(I-1))/P
       enddo
       
       IF (YPN.GT..99D30) THEN
         QN=0.d0
         UN=0.d0
       ELSE
         QN=0.5d0
         UN=(3.d0/(X(N)-X(N-1)))*(YPN-(Y(N)-Y(N-1))/(X(N)-X(N-1)))
       ENDIF
       Y2(N)=(UN-QN*U(N-1))/(QN*Y2(N-1)+1.d0)
       DO K=N-1,1,-1
              Y2(K)=Y2(K)*Y2(K+1)+U(K)
       enddo
       RETURN
       END
! ********************************************************************
      SUBROUTINE SPLINT(XA,YA,Y2A,N,X,Y)
      implicit none
      integer n
      double precision XA(N),YA(N),Y2A(N)
      integer klo,khi,k
      double precision x,y,h,a,b
       KLO=1
       KHI=N 
       do while(KHI-KLO.GT.1)
         K=(KHI+KLO)/2
         IF(XA(K).GT.X)THEN
           KHI=K
         ELSE
           KLO=K
         ENDIF
       enddo
       
       H=XA(KHI)-XA(KLO)
       IF (H.EQ.0.d0) then
         write(6,*) 'Bad XA input.'
         stop
       endif
       A=(XA(KHI)-X)/H
       B=(X-XA(KLO))/H
       Y=A*YA(KLO)+B*YA(KHI)+&
          ((A**3-A)*Y2A(KLO)+(B**3-B)*Y2A(KHI))*(H**2)/6.d0
       RETURN
       END    

end module measurements
