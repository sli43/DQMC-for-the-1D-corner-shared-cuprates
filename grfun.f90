module grfun
contains
 !=============================================================================
 ! Subroutine getgp
 ! This subroutine constructs the equaltime green's function on slice (ti-1),
 ! where ti is the input integer. It returns the following:
 ! 1) The green's function for the spin specified by "spin"
 ! 2) the logarithm of the det(inv(G)) (in logdet)
 ! 3) the sign of the determinant of inv(G).
 !=============================================================================
 subroutine getgp(G,ti,logdet,sgndet,spin)
 use fields, only: v_up, v_dn
 use linear_algebra, only: identity, udr_decomp, invertr, determinant
 use parameters, only: N, L, spinup, spindn, left, right, expk, orthlen
 implicit none
 integer jj, row, nl, spin, ti, slice
 integer, dimension(0:N-1) :: jpvt
 double precision, dimension(0:N-1,0:N-1) :: G, I, U, R, Rold, rmati
 double precision, dimension(0:N-1,0:N-1) :: stor, stori, temp
 double precision, dimension(0:N-1) :: bvec, tmpvec
 double precision sgndet, det, sgn,rtmp, logdet
 
 G = 0.0d0
 sgndet = 1.0d0
 logdet = 0.0d0
 call identity(I,N)
 U = I
 rold = I
 bvec = 1.0d0

 !set the time slice starting at ti and loop around to ti-1
 do nl = 0,L-1
  slice = mod(nl+ti,L)
  !we for a product of B(L) matrices where B(L) is exp(K)exp(V(l))
  do jj = 0,N-1
   if(spin.eq.spinup)then
    tmpvec(jj) = V_up(jj,slice)
   else
    tmpvec(jj) = V_dn(jj,slice)
   endif
  enddo
  call multb(U,tmpvec,-1)

  !once we have multiplied the first "orthlen" terms in the product we have to
  !perform a UDR decomposition. 
  if(mod(nl+1,orthlen).eq.0.or.(nl.eq.L-1))then
   !We start by rescaling the columns of U by the elements of bvec.  
   !bvec is the matrix D in the UDR notation. (D = Identity on the first pass).     
   !Note that on the first pass U is also the identity.  On later iterations U is the
   !previous orthogonal matrix.       
   do jj = 0,N-1
    call dscal(N,bvec(jj),U(0,jj),1)
   enddo

   !At this point we have B(L-m)B(L-m+1)...B(L)*U*D.  
   !On the first pass U = D = I.  We now form a UDR decomposition
   !of this product. 
   call UDR_Decomp(U,bvec,R,JPVT,N)
   !now apply the jpvt to the R matrix from the previous iteration and place it
   !in a temporary matrix "store"
   do row = 0,N-1
    stor(row,:) = Rold(JPVT(row)-1,:)
   enddo
   !we now multiply the new R matrix by the old one to complete the product
   call dtrmm('l','u','n','u',N,N,1.0d0,R,N,stor,N)
   call dcopy(N*N,stor,1,Rold,1)
  endif
  !nwo move on to the next group of B(l) matrices
 enddo !nl = 0,L-1
 !once we have excited this loop we have U*R where R is the
 !produce of each of the R's from the UDR decompositions.  U will hold the 
 !product of the orthogonal matrices.  U will also contain the output
 !of the last UDR decomposition and so it containts omat after the 
 !final UDR decomposition that is done.

 !We now can put together the Green's function using equation (20) 
 !of reference 1.  
 ! (20)   inv(G) = I + U*D*R
 !               = U(UR + D)R
 !               = UDR   <- which is inverted in the
 !               last step.
 !Calculate the inverse of the R matrix then store it back into R.
 call invertr(Rold,rmati,det,sgn,N)
 if(sgn.eq.0.0d0)then
  print*, "Warning (in getgp): Inverstion of a non-invertable matrix."
 endif
 if(sgn.lt.0.0d0) sgndet = -sgndet
 logdet = det
 R = rmati

 !now we calculate inv(U)*inv(R).  Since U is orthogonal we only need to 
 !multiply by the transpose of U.
 call DGEMM('T','N',N,N,N,1.0d0,U,N,R,N,0.0d0,stor,N)

 !add the elements of the b-vector to the diagonal of stor
 do jj = 0,N-1
  stor(jj,jj) = stor(jj,jj) + bvec(jj)
 enddo
 !now we invert stor and calculate its determinant. Place the inverse back into
 !stor.
 call invertr(stor,stori,det,sgn,N)
 stor = stori
 if(sgn.lt.0.0d0) sgndet = -sgndet
 logdet = logdet + det

 !multiply inv(R)*inv(stor) and place the result in tmp
 call DGEMM('n','n',N,N,N,1.0d0,R,N,stor,N,0.0d0,temp,N)
 !multiply temp*inv(Q) and store in G
 call DGEMM('n','t',N,N,N,1.0d0,temp,N,U,N,0.0d0,G,N)
 !finally, get the determinant of U to determine if we need to flip the sign 
 !of the determinant

 call determinant(U,det,rtmp,N)
 if(rtmp.lt.0.0d0) sgndet = -sgndet
 logdet = logdet + log(abs(det))

 return
 end subroutine getgp

 !=====================================================================
 ! subroutine multb(A,V,dir)
 ! This routine mutliplies a matrix A by B(l) where B(l) is a 
 ! diagonal matrix used to define the interactions.
 ! Direction = Left means   B*A
 !           = right means A*B 
 !=====================================================================
 subroutine multb(A,V,dir)
 use parameters, only: N, left, right, expk
 implicit none
 integer dir, i
 double precision, dimension(0:N-1) :: V
 double precision, dimension(0:N-1,0:N-1) :: A, store
 if(dir.eq.-1)then
  !copy the input matrix to store
  call dcopy(N*N,A,1,store,1)
  !multiply the input array by exp(k)
  call dgemm('n','n',N,N,N,1.0d0,expk,n,store,n,0.0d0,A,n)
  !V is a diagonal matrix so we simply rescale the rows
  do i = 0,N-1
   call dscal(N,V(i),A(i,0),N)
  enddo

 elseif(dir.eq.1)then
  !same as above but in the opposite directions
  !rescale the columns
  do i = 0,N-1
   call dscal(N,V(i),A(0,i),1)
  enddo
  !copy the input matrix to store
  call dcopy(N*N,A,1,store,1)
  !now multiply by exp(k) on the right left side and store the result in A.
  call dgemm('n','n',N,N,N,1.0d0,store,n,expk,n,0.0d0,A,n)
 else
  print*, 'Invalid Direction specified in multb.'
  stop
 endif

 return
 end subroutine multb

 !=====================================================================
 ! subroutine multbi(A,V,dir)
 ! This routine mutliplies a matrix A by inv(B(l)) where B(l) is a 
 ! diagonal matrix used to define the interactions.  B(l) is given by  
 ! a column vector V.
 ! Direction = Left means   B*A
 !           = right means A*B 
 !=====================================================================
 subroutine multbi(A,V,dir)
 use parameters, only: N, left, right, expki
 implicit none
 integer dir, i
 double precision, dimension(0:N-1) :: V
 double precision, dimension(0:N-1,0:N-1) :: A, store
 if(dir.eq.1)then
  !copy the input matrix to store
  call dcopy(N*N,A,1,store,1)
  !multiply the input array by exp(k)
  call dgemm('n','n',N,N,N,1.0d0,store,n,expki,n,0.0d0,A,n)
  !V is a diagonal matrix so we simply rescale the rows
  do i = 0,N-1
   call dscal(N,1.0d0/V(i),A(0,i),1)
  enddo

 elseif(dir.eq.-1)then
  !same as above but in the opposite directions
  !rescale the columns
  do i = 0,N-1
   call dscal(N,1.0d0/V(i),A(i,0),N)
  enddo
  !copy the input matrix to store
  call dcopy(N*N,A,1,store,1)
  !now multiply by exp(k) on the right left side and store the result in A.
  call dgemm('n','n',N,N,N,1.0d0,expki,n,store,n,0.0d0,A,n)
 else
  print*, 'Invalid Direction specified in multb.'
  stop
 endif
 return
 end subroutine multbi
 !=======================================================================
 subroutine get_seqB(ti,tip,U,D,R,spin)
 use fields, only: V_up, V_dn
 use linear_algebra, only: UDR_Decomp, identity
 use parameters, only: N, L, spinup, spindn, left, right, expk, orthlen, wraps, nwrap
 implicit none
 integer jj,ti, tip, spin, counter, slice
 integer, dimension(0:N-1) :: jpvt
 double precision, dimension(0:N-1,0:N-1) :: U, R, rold, stor
 double precision, dimension(0:N-1) :: D, vec 

 call identity(U,N)
 call identity(Rold,N)
 D = 1.0d0
 slice = ti
 counter = 0
 do while(slice.ne.mod(tip+1,L))
  counter = counter + 1
  if(spin.eq.spinup)then
   vec(:) = V_up(:,slice)
  else
   vec(:) = V_dn(:,slice)
  endif
  call multb(U,vec,-1) 

  !increment the slice and wrap back to zero if need be
  slice = slice + 1
  if(slice.eq.L) slice = 0
  if(mod(counter,orthlen).eq.0.or.slice.eq.mod(tip+1,L))then
   do jj = 0,N-1
    call dscal(N,D(jj),U(0,jj),1)
   enddo
   call UDR_Decomp(U,D,R,JPVT,N)
   do jj = 0,N-1
    stor(jj,:) = Rold(JPVT(jj)-1,:)
   enddo
   !we now multiply the new R matrix by the old one to complete the product
   call dtrmm('l','u','n','u',N,N,1.0d0,R,N,stor,N)
   call dcopy(N*N,stor,1,Rold,1)
  endif
 enddo
 R = Rold  
 return
 end subroutine get_seqB
 !=======================================================================
 subroutine get_seqBi(ti,tip,U,D,R,spin)
 use fields, only: V_up, V_dn
 use linear_algebra, only: udr_decomp, identity
 use parameters, only: N, L, spinup, spindn, left, right, expk, orthlen,&
                        wraps, nwrap
 implicit none
 integer jj,ti, tip, spin,  exit_value, counter, slice
 integer, dimension(0:N-1) :: jpvt
 double precision, dimension(0:N-1,0:N-1) :: U, R, rold, stor
 double precision, dimension(0:N-1) :: D, vec 
 
 call identity(U,N)
 call identity(Rold,N)
 D = 1.0d0
 slice = ti
 counter = 0
 exit_value = tip-1
 if(exit_value.lt.0) exit_value = exit_value + L

 do while(slice.ne.exit_value)
  counter = counter + 1
  if(spin.eq.spinup)then
   vec(:) = V_up(:,slice)
  else
   vec(:) = V_dn(:,slice)
  endif
  call multbi(U,vec,-1)

  !increment the slice and wrap back to zero if need be
  slice = slice - 1
  if(slice.lt.0) slice = L-1

  if(mod(counter,orthlen).eq.0.or.slice.eq.exit_value)then
   do jj = 0,N-1
    call dscal(N,D(jj),U(0,jj),1)
   enddo

   call UDR_Decomp(U,D,R,JPVT,N)
   do jj = 0,N-1
    stor(jj,:) = Rold(JPVT(jj)-1,:)
   enddo
   !we now multiply the new R matrix by the old one to complete the product
   call dtrmm('l','u','n','u',N,N,1.0d0,R,N,stor,N)
   call dcopy(N*N,stor,1,Rold,1)
  endif
 enddo
 R = Rold
 return
 end subroutine get_seqBi
 

 !=======================================================================
 ! This next routine calculates the unequal time green's function for 
 ! slices ti to ti1 i.e. 
 !   (G(\tau,\tau')) where \tau = \delta\tau*ti
 !                         \tau' = \delta\tau*tip
 ! This green's fucntion is placed into gt0.  The routine also calculates 
 ! G(\tau',\tau) and store the corresponding result in g0t. 
 ! 
 ! the Green's functions are defined as follows. 
 ! G(ti,tip) = G(i,j) where i,j denote the blocks for ti and tip.
 !
 !     G_ii = inv(I+B_iB_{i-1}...B_1B_L...B_{i+1})   
 !     G_ij = -Gup_ii(B_iB_{i-1}...B_1B_L...B_{j+1})  for j = i...L  
 !     G_ij =  Gup_ii(B_iB_{i-1}...B_1B_L...B_{j+1})  for j = 1...i-1  
 !
 ! Recall that the B(ti) matricies are quite stiff at low temperature and/or
 ! with strong interactions. In order to evaluate the product of B matricies 
 ! we need to play some games in order to maintain stability during 
 ! the calculation. (This procedure was adabpted from Richard's quest source 
 ! code.) 
 !
 ! Write G_ij = (+-) inv(I + A_1*A_2)*A1
 !            = (+-) inv(inv(A_1)+A_2)
 ! where A_1 = B_{i}...B{j+1}
 !       A_2 = B_{j}...B{i+1}
 !
 ! step 1. Perform UDR decomp on inv(A_1) and A_2
 !         inv(A_1) = U_1D_1R_1
 !             A_2  = U_2D_2R_2
 ! step 2. Decompose D_1 = \bar{D}_1*\hat{D}_1
 !                   D_2 = \bar{D}_2*\hat{D}_2
 !         where bar{D}_1)i,i) = max(1,D_1(i,i))
 !               hat{D}_2(i,i) = min(1,D_1(i,i))  
 ! step 3. 
 !   C = hat{D}_2*T_2*inv(T_1)*inv(bar{D}_1)
 !     + inv(bar{D}_2)*inv(U_2)*U_1*hat{D}_1
 ! step 4. 
 !    G = inv(T_1)*inv(barD_1)*inv(C)*inv(barD_2)*inv(U_2)
 !      = inv(T_1)*inv(D_2T_2inv(T_1)+inv(U_2)U_1D_1)*inv(U_2)
 !      = inv(U_1D_1T_1+U_2D_2T_2)
 !      = inv(inv(A_1)+A_2)
 !      = inv(I+A_1A_2)A_1
 !      = inv(I+B_{i}...B_1*B_l...B_{i-1})B_i...B_{j+1}
 ! similar tricks apply for G_ji
 !=======================================================================
 subroutine makegt(ti,ti1,gt0,g0t,spin)
 use fields, only: V_up, V_dn
 use linear_algebra, only: UDR_Decomp, identity
 use parameters, only: N, L, spinup, spindn, left, right, expk, orthlen,&
                       wraps, nwrap
 implicit none
 integer spin, ta, tb, info, i, i1, i2, ti1, ti, lwork
 integer, dimension(0:N-1) :: jpvt, jpvt2
 double precision det, sgndet
 double precision, parameter :: ONE = 1.0d0
 double precision, parameter :: ZERO = 0.0d0
 double precision, dimension(0:N-1) :: work
 double precision, dimension(0:N-1,0:N-1) :: U1, U2, R1, R2, W1, W2
 double precision, dimension(0:N-1) :: D1, D2, bar1i, bar2i, hat1, hat2
 double precision, dimension(0:N-1,0:N-1) :: gt0,g0t

 !If ti = ti1 then we want the equal time green's function.  In this case we
 !simply call getgp for slice ti+1 (recall that getgp returns the green's
 !function on the slice ti-1 for input ti.  (This is done because of how I
 !perform the wrapping in the sweep routines.)
 if(ti.eq.ti1)then
  i = mod(ti+1,L)
  call getgp(Gt0,i,det,sgndet,spin)
  do i1 = 0,N-1
   do i2 = 0,N-1
    if(i1.eq.i2)then
     g0t(i1,i2) =1.0d0-gt0(i1,i2)
    else
     g0t(i1,i2) =-gt0(i1,i2)
    endif
   enddo
  enddo
 else
  !if ti doesn't equal ti1 then we have some work to do. 
  lwork = N

  g0t = 0.0d0
  gt0 = 0.0d0

  !Step 1 for A_1
  !Form the product of B(i) matricies
  ta = mod(ti+1,L)
  tb = mod(ti1,L)
  call get_seqB(ta,tb,U2,D2,R2,spin)

  !step 1 form inv(A_1)
  ta = mod(ti,L)
  tb = mod(ti1+1,L)
  call get_seqBi(ta,tb,U1,D1,R1,spin)

  !step 2
  do i = 0,N-1
   bar1i(i) = ONE/max(one,abs(D1(i)))
   hat1(i) = D1(i)*bar1i(i)
   bar2i(i) = ONE/max(ONE,abs(D2(i)))
   hat2(i) = D2(i)*bar2i(i)
  enddo

  !steps 3 and 4 are completed with a series of lapack calls
  ! STEP 3. Compute C = hatD_2*T_2*inv(T_1)*inv(barD_1)+
  !                     inv(barD_2)*inv(U_2)*U_1*hatD_1
  !Copy R1 to storage array W2, we will need it later
  call dcopy(N*N,R1,1,W2,1)
  !set W1 = transpose(R2)
  W1 = transpose(R2)

  !W_1 = inv(W2')*W_1 = inv(R_1')*R_2'
  call dgetrf(n,n,W2,n,jpvt,info)
  call dgetrs('T',N,N,W2,N,jpvt,W1,N,info)

  !T2 = transpose(W1) = transpose(inv(R1')*R2') = R2*inv(R1)
  R2 = transpose(W1)

  !U1 = G_ij = U2'*U1
  call dgemm('T', 'N', n, n, n, ONE, U2, n, U1, n, ZERO, Gt0, n)
  call dcopy(n*n,Gt0,1,U1,1)

  ! U_1 = bar{D}_2*U_2'*U1*hat{D}_1
  call scalerow(N,Gt0,bar2i)
  call scalecol(N,gt0,hat1)

  ! W_1 = hat{D}_2*T_2(inv(T_1))*bar{D_1}
  call dcopy(N*N,R2,1,W1,1)
  call scalerow(N,W1,hat2)
  call scalecol(n,W1,bar1i)

  !W_1 = W_1 + G_ij
  call daxpy(N*N,one,gt0,1,W1,1)

  ! Step 4. compute inv(T_1)*inv(barD_1)*inv(C)*inv(barD_2)*inv(U_2)
  Gt0 = transpose(U2)
  call scalerow(N,Gt0,bar2i)
  call dgesv(N,N,W1,n,jpvt2,Gt0,N,info)

  call scalerow(N,gt0,bar1i)
  call dgetrs('N',N,N,W2,N,jpvt,gt0,n,info) !G_ij is now stored in gt0

  !now we need to get g0t (G_ji) by repeating the above...
  do i = 0,N-1
   !check for division by zero
   if(D1(i).eq.ZERO)then
    print*, "Error in makegt: D1(i) = 0.0, i = ", i
    stop
   endif
   if(D2(i).eq.ZERO)then
    print*, "Error in makegt: D2(i) = 0.0, i = ", i
    stop
   endif

   !form \bar{D}_1, \bar{D}_2, \hat{D}_1 and \hat{D}_2
   D1(i) = one/D1(i)
   bar1i(i) = one/max(ONE,D1(i))
   hat1(i) = D1(i)*bar1i(i)
   D2(i) = one/D2(i)
   bar2i(i) = one/max(ONE,D2(i))
   hat2(i) = D2(i)*bar2i(i)
  enddo

  !some stuff is already set up for us from constructing G_ij
  !We currently have R2 = R2*inv(R1) and U1 = inv(U2)*U1
  !therefore we just need to invert the two
  call dgetrf(n,n,R2,n,jpvt,info)
  call dgetri(N,R2,n,jpvt,work,lwork,info)
  !W1 = inv(inv(U2)*U1) = inv(U1)*U2
  W1 = Transpose(U1)

  !compute inv(bar{D}_1)*R1*inv(R2)*hat{D}_2
  call scalerow(N,R2,bar1i)
  call scalecol(N,R2,hat2)
  !compute hat{D}_1*inv(U1)*U2*inv(bar{D}_2)
  call scalerow(N,W1,hat1)
  call scalecol(N,W1,bar2i)

  !W1 = W1 + R2
  call daxpy(N*N,ONE,R2,1,W1,1)

  !Compute U_2*inv(barD_2)*inv(...)*inv(barD_1)*R_1
  call scalerow(N,R1,bar1i)
  call dgesv(n,n,W1,n,jpvt,R1,N,info)
  !inv(bar(D)_2)*inv(...)*inv(bar(D)_1)*R1
  call scalerow(N,R1,bar2i)
  !store previous result
  call dcopy(N*N,g0t,1,W2,1)
  !multiply -U2
  call dgemm('N','N',n,n,n,-ONE,U2,N,R1,N,ZERO,G0t,N)
  g0t=-g0t  
 endif
 

 return
 end subroutine makegt
 !=======================================================================
 subroutine ScaleCol(N,A,D)
 implicit none
 integer n, i
 double precision, dimension(0:N-1,0:N-1) :: A
 double precision, dimension(0:N-1) :: D
 do i = 0,N-1
  call dscal(n,D(i),A(0,i),1)
 enddo
 return
 end subroutine ScaleCol
 !=======================================================================
 !=======================================================================
 subroutine ScaleRow(N,A,D)
 implicit none
 integer n, i
 double precision, dimension(0:N-1,0:N-1) :: A
 double precision, dimension(0:N-1) :: D
 do i = 0,N-1
  call dscal(n,D(i),A(i,0),N)
 enddo
 return
 end subroutine ScaleRow 

end module grfun
