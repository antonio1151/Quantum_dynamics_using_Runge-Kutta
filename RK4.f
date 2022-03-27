       SUBROUTINE RK(DT,F0,F,N,T,RPAR)
       IMPLICIT NONE
       DOUBLE PRECISION RPAR(10),HBAR,DT,T
       INTEGER I,N,J
       COMPLEX*8 F(N),F0(N),F1(N),F2(N),F3(N),F4(N)
       COMPLEX*8 F10(N),F20(N),F30(N),F40(N),H2(N,N)
       COMPLEX*8 :: eye
       COMPLEX*8 H(N,N)
        HBAR=0.6582D0   !PLANCK'S CONSTANT IN eV*fs 
       eye = (0.0d0,1.0d0)                     !IMAGINARY NUMBER



       CALL HH(N,T,H,RPAR)
       H2=H
       F1=-eye*MATMUL(H2,F0)/HBAR

       F20=0.0D0
          DO J=1,N
           F20(J)=F0(J)+(F1(J)*DT)/2.0
          END DO

       CALL HH(N,T+DT/2.0,H,RPAR)

       H2=H
       F2=-eye*MATMUL(H2,F20)/HBAR

       F30=0.0D0
          DO J=1,N
           F30(J)=F0(J)+(F2(J)*DT)/2.0
          END DO

       CALL HH(N,T+DT/2.0,H,RPAR)

       H2=H
       F3=-eye*MATMUL(H2,F30)/HBAR

       F40=0.0D0
          DO J=1,N
           F40(J)=F0(J)+(F3(J)*DT)
          END DO


       CALL HH(N,T+DT,H,RPAR)
       H2=H
       F4=-eye*MATMUL(H2,F40)/HBAR

       F=0.0D0
          DO J=1,N
            F(J)=F0(J)+(DT/6.0)*(F1(J)+F4(J)+2.0*(F3(J)+F2(J)))
          END DO


       RETURN
       END 
