module fields
use parameters,only: Nlat, N, L
implicit none
 double precision, dimension(0:Nlat-1,0:L-1,1:20) :: S
 double precision, dimension(0:N-1,0:L-1) :: v_up
 double precision, dimension(0:N-1,0:L-1) :: v_dn
contains
 subroutine init_fields()
 use random
 use parameters, only: N, L, Nlat
 implicit none
 double precision r
 integer ti, i, j

 do i = 0,nlat-1
  do ti = 0,L-1
   do j = 1,20
    r = ran2(iran)-0.5d0
    if(r.le.0)then
     S(i,ti,j) = -1.0d0
    else
     S(i,ti,j) = 1.0d0
    endif
   enddo
  enddo
 enddo
 call return_v_for_s(S,v_up,v_dn)




 return
 end subroutine init_fields
 !=============================================================================
 subroutine return_v_for_S(S,v_up,v_dn)
 use parameters, only: N, L, Nlat, lambda, lambdap, lambdapp
 implicit none
 integer i, ti, i1
 double precision, dimension(0:Nlat-1,0:L-1,1:20) :: S
 double precision, dimension(0:N-1,0:L-1) :: v_up, v_dn
 !init the exp_V matrix to zero.
 v_up = 0.0d0
 v_dn = 0.0d0
 do i = 0,Nlat-1
  i1=i+1
  if(i1>=Nlat) i1=i1-Nlat
  do ti = 0,L-1
   !orbital 1
   v_up(i,ti) = dexp( S(i,ti,1)*lambda &
              + lambdapp*S(i,ti,5) + lambdapp*S(i,ti,6) &
              + lambdapp*S(i,ti,9) + lambdapp*S(i,ti,10)&
              + lambdapp*S(i,ti,13) + lambdapp*S(i,ti,14)&
              + 1.0*lambdapp*S(i,ti,17) + 1.0*lambdapp*S(i,ti,18) )
   v_dn(i,ti) = dexp(-S(i,ti,1)*lambda &
              + lambdapp*S(i,ti,7) + lambdapp*S(i,ti,8) &
              + lambdapp*S(i,ti,11) + lambdapp*S(i,ti,12) &
              + lambdapp*S(i,ti,15) + lambdapp*S(i,ti,16) &
              + 1.0*lambdapp*S(i,ti,19) + 1.0*lambdapp*S(i,ti,20))
   !orbital 2
   v_up(i+Nlat,ti) = dexp( S(i,ti,2)*lambdap &
                   - lambdapp*S(i,ti,5) - lambdapp*S(i,ti,7) &
                   - 1.0*lambdapp*S(i1,ti,17)- 1.0*lambdapp*S(i1,ti,19))
   v_dn(i+Nlat,ti) = dexp(-S(i,ti,2)*lambdap &
                   - lambdapp*S(i,ti,6) - lambdapp*S(i,ti,8)&
                   - 1.0*lambdapp*S(i1,ti,18)- 1.0*lambdapp*S(i1,ti,20))
   !orbital 3
   v_up(i+2*Nlat,ti) = dexp( S(i,ti,3)*lambdap &
                     - lambdapp*S(i,ti,9) - lambdapp*S(i,ti,11))

   v_dn(i+2*Nlat,ti) = dexp(-S(i,ti,3)*lambdap &
                     - lambdapp*S(i,ti,10) - lambdapp*S(i,ti,12) )
   !orbital 4
   v_up(i+3*Nlat,ti) = dexp( S(i,ti,4)*lambdap &
                       - lambdapp*S(i,ti,13) - lambdapp*S(i,ti,15))

   v_dn(i+3*Nlat,ti) = dexp(-S(i,ti,4)*lambdap &
                        - lambdapp*S(i,ti,14) - lambdapp*S(i,ti,16) )


  enddo
 enddo
 return
 end subroutine return_v_for_s
 !=============================================================================

 subroutine deallocate_fields()
 implicit none

 return
 end subroutine deallocate_fields

end module fields
