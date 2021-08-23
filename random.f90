module random
 double precision lastr
 integer iran
 integer,dimension(1:97):: label_IR
 DATA label_IFF /0/
 integer label_ia,label_ic,label_iy,label_j,label_m
contains
 DOUBLE PRECISION FUNCTION RAN2(idum)
 implicit none
 double precision rtmp, rm
 integer ia,ic,iff,ir,iy,j,m,idum
 save
 PARAMETER (M=714025,IA=1366,IC=150889,RM=1.4005112d-6)
 dimension IR(97)
 DATA IFF /0/

 ir=label_ir
 iff=label_iff
 iy=label_iy
 j=label_j

 IF(IDUM.LT.0.OR.IFF.EQ.0)THEN
  IFF=1
  IDUM=MOD(IC-IDUM,M)
  DO 11 J=1,97
   IDUM=MOD(IA*IDUM+IC,M)
   IR(J)=IDUM
11 CONTINUE
  IDUM=MOD(IA*IDUM+IC,M)
  IY=IDUM
 ENDIF
 J=1+(97*IY)/M
 IF(J.GT.97.OR.J.LT.1)PAUSE
 IY=IR(J)
 RAN2=IY*RM
 IDUM=MOD(IA*IDUM+IC,M)
 IR(J)=IDUM
 label_ir=ir
 label_iff=iff
 label_ia=ia
 label_ic=ic
 label_iy=iy
 label_j=j
 label_m=m
 RETURN
 END function ran2
end module random
