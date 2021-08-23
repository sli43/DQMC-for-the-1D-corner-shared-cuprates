module linear_algebra
contains
 !======================================================================
 double precision function diffmat(A,B,N)
 implicit none
 integer N, i, j
 double precision, dimension(1:N,1:N) :: A, B
 !S. Johnston, recoded to remove underflow error.
 diffmat = 1e-10
 do i = 1,N
  do j = 1,N
   diffmat = max(diffmat,abs(A(i,j)-B(i,j)))
  enddo 
 enddo
 return
 end function diffmat
 !======================================================================
 ! Subroutine UDR_DECOMP(U,D,R,JPVT,N) 
 ! This subroutine takes an input matrix U and determines its UDR-
 ! decomposition.  
 !======================================================================
 subroutine UDR_Decomp(U,D,R,JPVT,N)
 implicit none
 integer N, info, i, j
 integer Lwork 
 integer, dimension(0:N-1) :: jpvt
 double precision, dimension(1:3*N+1) :: WORK
 double precision, dimension(0:N-1) :: tau, D
 double precision, dimension(0:N-1,0:N-1) :: A, U, R, tmp
 LWORK = 3*N+1
 !clear the pivot array
 jpvt = 0
 !call DGEQP3 to form the UDR decomposition
 call DGEQP3(N,N,U,N,jpvt,tau,work,lwork,info)
 if(info.ne.0)then
  print*, 'Error: in routine UDR_Decomp, calling dgeqp3.'
  stop
 endif
 !We now build the upper triangular matrix R.  
 !Copy A into R so we can work with it. 
 R = U
 !now zero out the parts of R we don't need. Ie the Lower triangular part.
 do i = 0,N-1
  do j = i+1,N-1
   R(j,i) = 0.0d0
  enddo
  D(i) = R(i,i)
  if(D(i).eq.0.0d0) D(i) = 1.0d0
  R(i,:) = R(i,:)/D(i)
 enddo
 !finally, find the orthogonal matrix U using a BLAS Call
 call DORGQR(N,N,N,U,N,tau,WORK,Lwork,info)
 if(info.ne.0)then
  print*, 'Error: in routine UDR_Decomp, calling dorgqr.'
  stop
 endif
 return
 end subroutine UDR_Decomp
 !======================================================================
  
 !======================================================================
 !  subtrouine exp_symmetric_mat(A,exp_A,N)
 !
 ! This routine evaluates the exponential of a NxN matrix.  This is
 ! done using a series of lapack calls to diagonalize the matrix 
 ! first.  
 !
 ! NOTE: This routine assumes that the matrix is symmetric!
 !======================================================================
 subroutine exp_symmetric_mat(A,exp_A,N)
 implicit none
 integer i, N, INFO
 integer LWORK 
 double precision det, sgn
 double precision, dimension(1:3*N-1) :: WORK
 double precision, dimension(1:N,1:N) :: A, exp_A, tmp
 double precision, dimension(1:N,1:N) :: VL,VR,P, invP, D
 double precision, dimension(1:N) :: WR,WI

 LWORK = 3*N-1
 !First, perform a lapack call to diagonalize A.
 !P is the transformation matrix, D is a diagonal matrix with 
 !the eigenvalues of A.
 tmp = A
 call DSYEV('V','U',N,A,N,WR,WORK,LWORK,INFO)
 !Now build D, and caluculate Pinv
 P = A
 call invertr(P,invP,det,sgn,N)
 D = 0.0d0
 do i = 1,N
  D(i,i) = dexp(WR(i))
 enddo
 exp_A = matmul(P,matmul(D,invP))
 A = tmp
 return
 end subroutine exp_symmetric_mat
 !======================================================================
 ! subroutine determinant returns (det(A))) and sign(det(A))
 !======================================================================
 subroutine determinant(A,deta,sgn,N)
 implicit none
 integer i, N, info
 integer, dimension(1:N) :: IPIV
 double precision detA, sgn
 double precision, dimension(1:N,1:N) :: A, Ainv
 double precision, dimension(1:N) :: Work

 Ainv = A
 call dgetrf(N,N,Ainv,N,ipiv,info)
 detA = 1.0d0
 sgn = 1.0d0
 do i = 1,N
  detA = detA*Ainv(i,i)
  if(IPIV(i).ne.i)then
   detA = -detA
  endif
 enddo
 if(deta.lt.0.0d0) sgn = -sgn

 return
 end subroutine determinant
 !======================================================================
 ! Subroutine invertr
 ! This routine inverts a real NxN matrix.
 !======================================================================
 subroutine invertr(A,Ainv,deta,sgn,N)
 implicit none
 integer I, N, info
 integer, dimension(1:N) :: ipiv
 double precision deta, sgn
 double precision, dimension(1:N,1:N) :: A, Ainv
 double precision, dimension(1:N) :: WORK

 Ainv = A
 !Make a function call to dgetrf to perform an LU decomposition.
 call DGETRF(N,N,Ainv,N,IPIV,INFO)

 !now calculate the determinant from the LU decomposition using the 
 !property that det(A*B) = det(A)*det(B).  A is lower triangular with 1's on its 
 !diagonal.
 sgn = 1.0d0
 detA = 0.0d0
 do i = 1,N
  if(IPIV(i).ne.i) sgn = -sgn
  sgn = sgn*Ainv(i,i)/abs(Ainv(i,i))
  deta = DetA + log(abs(Ainv(i,i)))
 enddo

 if(sgn.eq.0.0d0)then
  Ainv = 0.0d0
  print*, 'Warning: Attempted to invert a non-invertable matrix'
 else
  !Now call dgetri to perform the matrix inversion.
  call DGETRI(N,Ainv,N,IPIV,WORK,N,INFO)  
 endif
 return
 end subroutine invertr
 !======================================================================
 !======================================================================
 subroutine identity(A,N)
 implicit none
 integer N, i
 double precision, dimension(1:N,1:N) :: A
 A = 0.0d0
 do i = 1,N
  A(i,i) = 1.0d0
 enddo
 return
 end subroutine identity
 !======================================================================
end module linear_algebra
