module inter_file
contains
subroutine inter_writefile(iwarmbegin,imeasbegin,sgnup,sgndn,gup,gdn,detup,detdn,bin)
use parameters
use fields
use random
use mpi_module
use grfun
use measurements
use updates
use io_routines
implicit none
double precision sgnup, sgndn, detup, detdn
integer iwarmbegin,imeasbegin
double precision, dimension(0:N-1,0:N-1) :: Gup, gdn
integer i,j,ti,bin,k,i1,i2,j1,j2, ix, iy
integer, parameter::fwrte=582
character  inter_output_file*100

 write(unit=inter_output_file,fmt="('inter_run_',i7.7,'+rank_',i3.3,'+runtime_',i3.3,'.dat')") 1000000+run,100+myrank,&
                                                                                            100+seperate_run_time+1
 open(file=inter_output_file,unit=fwrte,action='write')

!--------------------------------------------------------------------------------------------
!                parameters part
write(unit=fwrte,fmt=*) nwrap,wraps
write(unit=fwrte,fmt=*) lambda,lambdap,lambdapp,gam
!--------------------------------------------------------------------------------------------
!               random part
write(unit=fwrte,fmt=*) iran,label_ia,label_ic,label_iy,label_j,label_m
write(unit=fwrte,fmt=*) (label_IR(i1),i1=1,97)
write(unit=fwrte,fmt=*) label_iff
!--------------------------------------------------------------------------------------------
!               main part
write(unit=fwrte,fmt=*) iwarmbegin,imeasbegin,bin
write(unit=fwrte,fmt=*) sgnup,sgndn,detup,detdn
do ix=0,N-1
  write(unit=fwrte,fmt=*) (gup(ix,iy),iy=0,N-1)
end do

do ix=0,N-1
  write(unit=fwrte,fmt=*) (gdn(ix,iy),iy=0,N-1)
end do
!----------------------------------------------------------------------------------------------
!               fields part
do ix=0,Nlat-1
 do ti=0,L-1
   write(unit=fwrte,fmt=*) (S(ix,ti,i1),i1=1,20)
 enddo
end do

do ix=0,N-1
   write(unit=fwrte,fmt=*) (v_up(ix,ti),ti=0,L-1)
end do

do ix=0,N-1
   write(unit=fwrte,fmt=*) (V_dn(ix,ti),ti=0,L-1)
end do

!----------------------------------------------------------------------------------------------
!                updates
 write(unit=fwrte,fmt=*) sweeps_per_et_meas,&
                        accept_hs, reject_hs, redo, noredo,&
                        calls_to_sweep_hs,&
                        accept_blk_hs, reject_blk_hs
!-----------------------------------------------------------------------------------------------
!               measurements
 write(unit=fwrte,fmt=*) nmeas, nmeast
 write(unit=fwrte,fmt=*) asgn, asgnup, asgndn, &
                         an1, anup1, andn1, an2, anup2, andn2,&
                         an3, anup3, andn3, an4, anup4, andn4,&
                         anud1, anud2, anud3, anud4, amz1, amz2,&
                         amz3, amz4, asgnt, an, anup, andn, ake,&
                        ape1, ape2, ape3, aee, ant1, ant2, ant3, ant4

 write(unit=fwrte,fmt=*) (bsgn(i2),i2=1,nbins)
 write(unit=fwrte,fmt=*) (bsgnup(i2),i2=1,nbins)
 write(unit=fwrte,fmt=*) (bsgndn(i2),i2=1,nbins)
 write(unit=fwrte,fmt=*) (bn1(i2),i2=1,nbins)
 write(unit=fwrte,fmt=*) (bnup1(i2),i2=1,nbins)
 write(unit=fwrte,fmt=*) (bndn1(i2),i2=1,nbins)
 write(unit=fwrte,fmt=*) (bn2(i2),i2=1,nbins)
 write(unit=fwrte,fmt=*) (bnup2(i2),i2=1,nbins)
 write(unit=fwrte,fmt=*) (bndn2(i2),i2=1,nbins)
 write(unit=fwrte,fmt=*) (bn3(i2),i2=1,nbins)
 write(unit=fwrte,fmt=*) (bnup3(i2),i2=1,nbins)
 write(unit=fwrte,fmt=*) (bndn3(i2),i2=1,nbins)
 write(unit=fwrte,fmt=*) (bn4(i2),i2=1,nbins)
 write(unit=fwrte,fmt=*) (bnup4(i2),i2=1,nbins)
 write(unit=fwrte,fmt=*) (bndn4(i2),i2=1,nbins)
 write(unit=fwrte,fmt=*) (bnud1(i2),i2=1,nbins)
 write(unit=fwrte,fmt=*) (bnud2(i2),i2=1,nbins)
 write(unit=fwrte,fmt=*) (bnud3(i2),i2=1,nbins)
 write(unit=fwrte,fmt=*) (bnud4(i2),i2=1,nbins)
 write(unit=fwrte,fmt=*) (bmz1(i2),i2=1,nbins)
 write(unit=fwrte,fmt=*) (bmz2(i2),i2=1,nbins)
 write(unit=fwrte,fmt=*) (bmz3(i2),i2=1,nbins)
 write(unit=fwrte,fmt=*) (bmz4(i2),i2=1,nbins)
 write(unit=fwrte,fmt=*) (bsgnt(i2),i2=1,nbins)
 write(unit=fwrte,fmt=*) (bn(i2),i2=1,nbins)
 write(unit=fwrte,fmt=*) (bnup(i2),i2=1,nbins)
 write(unit=fwrte,fmt=*) (bndn(i2),i2=1,nbins)
 write(unit=fwrte,fmt=*) (bsgn(i2),i2=1,nbins)
 write(unit=fwrte,fmt=*) (bke(i2),i2=1,nbins)
 write(unit=fwrte,fmt=*) (bpe1(i2),i2=1,nbins)
 write(unit=fwrte,fmt=*) (bpe2(i2),i2=1,nbins)
 write(unit=fwrte,fmt=*) (bpe3(i2),i2=1,nbins)
 write(unit=fwrte,fmt=*) (bee(i2),i2=1,nbins)
 write(unit=fwrte,fmt=*) (bone(i2),i2=1,nbins)
 write(unit=fwrte,fmt=*) (bnt1(i2),i2=1,nbins)
 write(unit=fwrte,fmt=*) (bnt2(i2),i2=1,nbins)
 write(unit=fwrte,fmt=*) (bnt3(i2),i2=1,nbins)
 write(unit=fwrte,fmt=*) (bnt4(i2),i2=1,nbins)


 do i=0,Nlat-1
  do ti=0,L
      write(unit=fwrte,fmt=*) ((agnl(i,ti,j1,j2),j1=1,4),j2=1,4)
  enddo
 enddo

  do i=0,Nlat-1
  do ti=0,L
      write(unit=fwrte,fmt=*) ((achinl(i,ti,j1,j2),j1=1,4),j2=1,4)
  enddo
 enddo

 do i=0,Nlat-1
  do ti=0,L
      write(unit=fwrte,fmt=*) ((achinlz(i,ti,j1,j2),j1=1,4),j2=1,4)
  enddo
 enddo

 do i=0,Nlat-1
  do ti=0,L
      write(unit=fwrte,fmt=*) ((adent(i,ti,j1,j2),j1=1,4),j2=1,4)
  enddo
 enddo

 do i=0,Nlat-1
  do ti=0,L
    write(unit=fwrte,fmt=*) (((bgnl(k,i,ti,j1,j2),j1=1,4),j2=1,4),k=1,nbins)
  enddo
 enddo

 do i=0,Nlat-1
  do ti=0,L
    write(unit=fwrte,fmt=*) (((bchinl(k,i,ti,j1,j2),j1=1,4),j2=1,4),k=1,nbins)
  enddo
 enddo

 do i=0,Nlat-1
  do ti=0,L
    write(unit=fwrte,fmt=*) (((bchinlz(k,i,ti,j1,j2),j1=1,4),j2=1,4),k=1,nbins)
  enddo
 enddo

 do i=0,Nlat-1
  do ti=0,L
    write(unit=fwrte,fmt=*) (((bdent(k,i,ti,j1,j2),j1=1,4),j2=1,4),k=1,nbins)
  enddo
 enddo

 do i=0,Nlat-1
  do ti=0,L
   do j1=1,4
    do j2=1,4
    write(unit=fwrte,fmt=*) (real(bgnlq(k,i,ti,j1,j2)),dimag(bgnlq(k,i,ti,j1,j2)),k=1,nbins)
    enddo
   enddo
  enddo
 enddo

 do i=0,Nlat-1
  do ti=0,L
   do j1=1,4
    do j2=1,4
    write(unit=fwrte,fmt=*) (real(bchinlq(k,i,ti,j1,j2)),dimag(bchinlq(k,i,ti,j1,j2)),k=1,nbins)
    enddo
   enddo
  enddo
 enddo

 do i=0,Nlat-1
  do ti=0,L
   do j1=1,4
    do j2=1,4
    write(unit=fwrte,fmt=*) (real(bchinlzq(k,i,ti,j1,j2)),dimag(bchinlzq(k,i,ti,j1,j2)),k=1,nbins)
    enddo
   enddo
  enddo
 enddo

 do i=0,Nlat-1
  do ti=0,L
   do j1=1,4
    do j2=1,4
    write(unit=fwrte,fmt=*) (real(bdentq(k,i,ti,j1,j2)),dimag(bdentq(k,i,ti,j1,j2)),k=1,nbins)
    enddo
   enddo
  enddo
 enddo

do i=0,Nlat-1
    write(unit=fwrte,fmt=*) ((agrfun(i,j1,j2),j1=1,4),j2=1,4)
enddo

do i=0,Nlat-1
    write(unit=fwrte,fmt=*) ((aspinxx(i,j1,j2),j1=1,4),j2=1,4)
enddo

do i=0,Nlat-1
    write(unit=fwrte,fmt=*) ((aspinzz(i,j1,j2),j1=1,4),j2=1,4)
enddo

do i=0,Nlat-1
    write(unit=fwrte,fmt=*) ((aden(i,j1,j2),j1=1,4),j2=1,4)
enddo

do i=0,Nlat-1
    write(unit=fwrte,fmt=*) (((bgrfun(k,i,j1,j2),j1=1,4),j2=1,4),k=1,nbins)
enddo

do i=0,Nlat-1
    write(unit=fwrte,fmt=*) (((bspinxx(k,i,j1,j2),j1=1,4),j2=1,4),k=1,nbins)
enddo

do i=0,Nlat-1
    write(unit=fwrte,fmt=*) (((bspinzz(k,i,j1,j2),j1=1,4),j2=1,4),k=1,nbins)
enddo

do i=0,Nlat-1
    write(unit=fwrte,fmt=*) (((bden(k,i,j1,j2),j1=1,4),j2=1,4),k=1,nbins)
enddo


do i=0,Nlat-1
 do j1=1,4
  do j2=1,4
   write(unit=fwrte,fmt=*) (real(bgkfun(k,i,j1,j2)), dimag(bgkfun(k,i,j1,j2)),k=1,nbins)
  enddo
 enddo
enddo

do i=0,Nlat-1
 do j1=1,4
  do j2=1,4
   write(unit=fwrte,fmt=*) (real(bspinxxq(k,i,j1,j2)), dimag(bspinxxq(k,i,j1,j2)),k=1,nbins)
  enddo
 enddo
enddo

do i=0,Nlat-1
 do j1=1,4
  do j2=1,4
   write(unit=fwrte,fmt=*) (real(bspinzzq(k,i,j1,j2)), dimag(bspinzzq(k,i,j1,j2)),k=1,nbins)
  enddo
 enddo
enddo

do i=0,Nlat-1
 do j1=1,4
  do j2=1,4
   write(unit=fwrte,fmt=*) (real(bdenq(k,i,j1,j2)), dimag(bdenq(k,i,j1,j2)),k=1,nbins)
  enddo
 enddo
enddo



close(fwrte)
return
end subroutine inter_writefile
!=======================================================================================
subroutine inter_readfile(iwarmbegin,imeasbegin,sgnup,sgndn,gup,gdn,detup,detdn,bin)
use parameters
use fields
use random
use mpi_module
use grfun
use measurements
use updates
use io_routines
implicit none
double precision sgnup, sgndn, detup, detdn
integer iwarmbegin,imeasbegin
double precision, dimension(0:N-1,0:N-1) :: Gup, gdn
double precision, dimension(1:nbins):: rra, rrb
integer i,j,ti,bin,k,i1,i2,j1,j2, ix, iy
integer, parameter::fwrte=581
character  inter_output_file*100,rcmpt*52

write(unit=inter_output_file,fmt="('inter_run_',i7.7,'+rank_',i3.3,'+runtime_',i3.3,'.dat')") 1000000+run,100+myrank,&
                                                                                             100+seperate_run_time
 open(file=inter_output_file,unit=fwrte,action='read')

!--------------------------------------------------------------------------------------------
!                parameters part
read(unit=fwrte,fmt=*) nwrap,wraps
read(unit=fwrte,fmt=*) lambda,lambdap,lambdapp,gam
!--------------------------------------------------------------------------------------------
!               random part
read(unit=fwrte,fmt=*) iran,label_ia,label_ic,label_iy,label_j,label_m
read(unit=fwrte,fmt=*) (label_IR(i1),i1=1,97)
read(unit=fwrte,fmt=*) label_iff
!--------------------------------------------------------------------------------------------
!               main part
read(unit=fwrte,fmt=*) iwarmbegin,imeasbegin,bin
read(unit=fwrte,fmt=*) sgnup,sgndn,detup,detdn
do ix=0,N-1
  read(unit=fwrte,fmt=*) (gup(ix,iy),iy=0,N-1)
end do

do ix=0,N-1
  read(unit=fwrte,fmt=*) (gdn(ix,iy),iy=0,N-1)
end do
!----------------------------------------------------------------------------------------------
!               fields part
do ix=0,Nlat-1
 do ti=0,L-1
   read(unit=fwrte,fmt=*) (S(ix,ti,i1),i1=1,20)
 enddo
end do

do ix=0,N-1
   read(unit=fwrte,fmt=*) (v_up(ix,ti),ti=0,L-1)
end do

do ix=0,N-1
   read(unit=fwrte,fmt=*) (V_dn(ix,ti),ti=0,L-1)
end do

!----------------------------------------------------------------------------------------------
!                updates
 read(unit=fwrte,fmt=*) sweeps_per_et_meas,&
                        accept_hs, reject_hs, redo, noredo,&
                        calls_to_sweep_hs,&
                        accept_blk_hs, reject_blk_hs
!-----------------------------------------------------------------------------------------------
!               measurements
 read(unit=fwrte,fmt=*) nmeas, nmeast
 read(unit=fwrte,fmt=*) asgn, asgnup, asgndn, &
                         an1, anup1, andn1, an2, anup2, andn2,&
                         an3, anup3, andn3, an4, anup4, andn4,&
                         anud1, anud2, anud3, anud4, amz1, amz2,&
                         amz3, amz4, asgnt, an, anup, andn, ake,&
                        ape1, ape2, ape3, aee, ant1, ant2, ant3, ant4

 read(unit=fwrte,fmt=*) (bsgn(i2),i2=1,nbins)
 read(unit=fwrte,fmt=*) (bsgnup(i2),i2=1,nbins)
 read(unit=fwrte,fmt=*) (bsgndn(i2),i2=1,nbins)
 read(unit=fwrte,fmt=*) (bn1(i2),i2=1,nbins)
 read(unit=fwrte,fmt=*) (bnup1(i2),i2=1,nbins)
 read(unit=fwrte,fmt=*) (bndn1(i2),i2=1,nbins)
 read(unit=fwrte,fmt=*) (bn2(i2),i2=1,nbins)
 read(unit=fwrte,fmt=*) (bnup2(i2),i2=1,nbins)
 read(unit=fwrte,fmt=*) (bndn2(i2),i2=1,nbins)
 read(unit=fwrte,fmt=*) (bn3(i2),i2=1,nbins)
 read(unit=fwrte,fmt=*) (bnup3(i2),i2=1,nbins)
 read(unit=fwrte,fmt=*) (bndn3(i2),i2=1,nbins)
 read(unit=fwrte,fmt=*) (bn4(i2),i2=1,nbins)
 read(unit=fwrte,fmt=*) (bnup4(i2),i2=1,nbins)
 read(unit=fwrte,fmt=*) (bndn4(i2),i2=1,nbins)
 read(unit=fwrte,fmt=*) (bnud1(i2),i2=1,nbins)
 read(unit=fwrte,fmt=*) (bnud2(i2),i2=1,nbins)
 read(unit=fwrte,fmt=*) (bnud3(i2),i2=1,nbins)
 read(unit=fwrte,fmt=*) (bnud4(i2),i2=1,nbins)
 read(unit=fwrte,fmt=*) (bmz1(i2),i2=1,nbins)
 read(unit=fwrte,fmt=*) (bmz2(i2),i2=1,nbins)
 read(unit=fwrte,fmt=*) (bmz3(i2),i2=1,nbins)
 read(unit=fwrte,fmt=*) (bmz4(i2),i2=1,nbins)
 read(unit=fwrte,fmt=*) (bsgnt(i2),i2=1,nbins)
 read(unit=fwrte,fmt=*) (bn(i2),i2=1,nbins)
 read(unit=fwrte,fmt=*) (bnup(i2),i2=1,nbins)
 read(unit=fwrte,fmt=*) (bndn(i2),i2=1,nbins)
 read(unit=fwrte,fmt=*) (bsgn(i2),i2=1,nbins)
 read(unit=fwrte,fmt=*) (bke(i2),i2=1,nbins)
 read(unit=fwrte,fmt=*) (bpe1(i2),i2=1,nbins)
 read(unit=fwrte,fmt=*) (bpe2(i2),i2=1,nbins)
 read(unit=fwrte,fmt=*) (bpe3(i2),i2=1,nbins)
 read(unit=fwrte,fmt=*) (bee(i2),i2=1,nbins)
 read(unit=fwrte,fmt=*) (bone(i2),i2=1,nbins)
 read(unit=fwrte,fmt=*) (bnt1(i2),i2=1,nbins)
 read(unit=fwrte,fmt=*) (bnt2(i2),i2=1,nbins)
 read(unit=fwrte,fmt=*) (bnt3(i2),i2=1,nbins)
 read(unit=fwrte,fmt=*) (bnt4(i2),i2=1,nbins)


 do i=0,Nlat-1
  do ti=0,L
      read(unit=fwrte,fmt=*) ((agnl(i,ti,j1,j2),j1=1,4),j2=1,4)
  enddo
 enddo

  do i=0,Nlat-1
  do ti=0,L
      read(unit=fwrte,fmt=*) ((achinl(i,ti,j1,j2),j1=1,4),j2=1,4)
  enddo
 enddo

 do i=0,Nlat-1
  do ti=0,L
      read(unit=fwrte,fmt=*) ((achinlz(i,ti,j1,j2),j1=1,4),j2=1,4)
  enddo
 enddo

 do i=0,Nlat-1
  do ti=0,L
      read(unit=fwrte,fmt=*) ((adent(i,ti,j1,j2),j1=1,4),j2=1,4)
  enddo
 enddo

 do i=0,Nlat-1
  do ti=0,L
    read(unit=fwrte,fmt=*) (((bgnl(k,i,ti,j1,j2),j1=1,4),j2=1,4),k=1,nbins)
  enddo
 enddo

 do i=0,Nlat-1
  do ti=0,L
    read(unit=fwrte,fmt=*) (((bchinl(k,i,ti,j1,j2),j1=1,4),j2=1,4),k=1,nbins)
  enddo
 enddo

 do i=0,Nlat-1
  do ti=0,L
    read(unit=fwrte,fmt=*) (((bchinlz(k,i,ti,j1,j2),j1=1,4),j2=1,4),k=1,nbins)
  enddo
 enddo

 do i=0,Nlat-1
  do ti=0,L
    read(unit=fwrte,fmt=*) (((bdent(k,i,ti,j1,j2),j1=1,4),j2=1,4),k=1,nbins)
  enddo
 enddo

 do i=0,Nlat-1
  do ti=0,L
   do j1=1,4
    do j2=1,4
      read(unit=fwrte,fmt=*) ( rra(k), rrb(k), k=1,nbins)
        bgnlq(:,i,ti,j1,j2)=dcmplx(rra(:),rrb(:)) 
    enddo
   enddo
  enddo
 enddo

 do i=0,Nlat-1
  do ti=0,L
   do j1=1,4
    do j2=1,4
       read(unit=fwrte,fmt=*) ( rra(k), rrb(k), k=1,nbins)
       bchinlq(:,i,ti,j1,j2) = dcmplx(rra(:),rrb(:))
    enddo
   enddo
  enddo
 enddo

 do i=0,Nlat-1
  do ti=0,L
   do j1=1,4
    do j2=1,4
      read(unit=fwrte,fmt=*) ( rra(k), rrb(k), k=1,nbins)
      bchinlzq(:,i,ti,j1,j2)=dcmplx(rra(:),rrb(:))
    enddo
   enddo
  enddo
 enddo

 do i=0,Nlat-1
  do ti=0,L
   do j1=1,4
    do j2=1,4
       read(unit=fwrte,fmt=*) ( rra(k), rrb(k), k=1,nbins)
       bdentq(:,i,ti,j1,j2)=dcmplx(rra(:),rrb(:))
    enddo
   enddo
  enddo
 enddo

do i=0,Nlat-1
    read(unit=fwrte,fmt=*) ((agrfun(i,j1,j2),j1=1,4),j2=1,4)
enddo

do i=0,Nlat-1
    read(unit=fwrte,fmt=*) ((aspinxx(i,j1,j2),j1=1,4),j2=1,4)
enddo

do i=0,Nlat-1
    read(unit=fwrte,fmt=*) ((aspinzz(i,j1,j2),j1=1,4),j2=1,4)
enddo

do i=0,Nlat-1
    read(unit=fwrte,fmt=*) ((aden(i,j1,j2),j1=1,4),j2=1,4)
enddo

do i=0,Nlat-1
    read(unit=fwrte,fmt=*) (((bgrfun(k,i,j1,j2),j1=1,4),j2=1,4),k=1,nbins)
enddo

do i=0,Nlat-1
    read(unit=fwrte,fmt=*) (((bspinxx(k,i,j1,j2),j1=1,4),j2=1,4),k=1,nbins)
enddo

do i=0,Nlat-1
    read(unit=fwrte,fmt=*) (((bspinzz(k,i,j1,j2),j1=1,4),j2=1,4),k=1,nbins)
enddo

do i=0,Nlat-1
    read(unit=fwrte,fmt=*) (((bden(k,i,j1,j2),j1=1,4),j2=1,4),k=1,nbins)
enddo


do i=0,Nlat-1
 do j1=1,4
  do j2=1,4
     read(unit=fwrte,fmt=*) ( rra(k), rrb(k), k=1,nbins)
     bgkfun(:,i,j1,j2)=dcmplx(rra(:),rrb(:))
  enddo
 enddo
enddo

do i=0,Nlat-1
 do j1=1,4
  do j2=1,4
    read(unit=fwrte,fmt=*) ( rra(k), rrb(k), k=1,nbins)
    bspinxxq(:,i,j1,j2)=dcmplx(rra(:),rrb(:))
  enddo
 enddo
enddo

do i=0,Nlat-1
 do j1=1,4
  do j2=1,4
     read(unit=fwrte,fmt=*) ( rra(k), rrb(k), k=1,nbins)
     bspinzzq(:,i,j1,j2)=dcmplx(rra(:),rrb(:))
  enddo
 enddo
enddo

do i=0,Nlat-1
 do j1=1,4
  do j2=1,4
   read(unit=fwrte,fmt=*) ( rra(k), rrb(k), k=1,nbins)
   bdenq(:,i,j1,j2)=dcmplx(rra(:),rrb(:))
  enddo
 enddo
enddo


close(fwrte)
return
end subroutine inter_readfile

end module inter_file
