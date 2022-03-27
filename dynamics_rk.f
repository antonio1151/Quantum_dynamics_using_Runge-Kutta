       PROGRAM RK DYNAMICS

       IMPLICIT NONE
       COMPLEX*8,ALLOCATABLE,DIMENSION(:)::C_INIT,C_NOW
       INTEGER I,J,N,IDT,JOUT
       COMPLEX*8 EYE
       DOUBLE PRECISION TMAX,TIME_STEP,TOUT,T,NORM
       DOUBLE PRECISION RPAR(6),R,RI
       DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:)::PSI0,PSI
C=======================================E====================================
C----CONSTANTS
       EYE=(0.0D0,1.0D0)  !IMAGINARY NUMBER i
C------PARAMETRS FOR THE DYNAMICS
       TIME_STEP=0.01D0   !TIME STEP FOR THE DYNAMICS IN fs
       TOUT=0.10D0         !TIME INTERVAL IN fs FOR PRINTING OUTPUT 
       TMAX=1000.0D0         !TOTAL TIME OF THE DYNAMICS
C===========================================================================
C-----PARAMETER OF THE SYSTEM
       N=2              !NUMBER OF COEFFICIENTS
C------PULSE PARAMETERS

       RPAR(1)=1.0D0   !FREQUENCY OF THE PULSE IN eV
       RPAR(2)=300.0D0 !TIME TO TURN THE PULSE IN fs
       RPAR(3)=100.0D0   !SOFTNESS   OF THE PULSE IN fs 
       RPAR(4)=0.1D0    !FIELD AMPLITUDE IN V/A

C============================================================================
C---HAMILTONIAN
        RPAR(5)=1.0D0    !DIAGONAL TERMS FOR THE HAMILTONIAN IN eV 
        RPAR(6)=1.0D0      !DIPOLE MOMENT IN ANGSTROMS
C=============================================================================
C----SETTING DIMENSIONS
      ALLOCATE(C_INIT(N),C_NOW(N))
C=============================================================================
C----INITIAL CONDITION

        OPEN(UNIT=10,FILE='INITIAL_STATE.dat',STATUS='UNKNOWN')


        DO I=1,N
          READ(10,*)R,RI
         C_INIT(I)=CMPLX(R,RI)
        END DO
C===============================================================================
C----PROPAGATION
      IDT=INT(TMAX/TIME_STEP)
      JOUT=INT(TOUT/TIME_STEP)
      T=0.0D0

       OPEN(UNIT=12,FILE='DYNAMICS.dat',STATUS='UNKNOWN')   

       WRITE(12,*)T,(REAL(C_INIT(I)*CONJG(C_INIT(I))),I=1,N) 


       DO I=1,IDT
          CALL RK(TIME_STEP,C_INIT,C_NOW,N,T,RPAR)
          C_INIT=C_NOW
          
          NORM=0.0D0
          DO J=1,N
             NORM=NORM+C_INIT(J)*CONJG(C_INIT(J))
          END DO
          DO J=1,N
             C_INIT(J)=C_INIT(J)/SQRT(NORM)
          END DO

          IF(MOD(I,JOUT).EQ.0) THEN
            WRITE(12,*)T,(REAL(C_INIT(J)*CONJG(C_INIT(J))),J=1,N)
          END IF

       T=T+TIME_STEP
       END DO


      END PROGRAM

C================================================================================
C---HAMILTONIAN
       SUBROUTINE HH(N,T,H,RPAR)
       IMPLICIT NONE
       INTEGER I,N
       DOUBLE PRECISION RPAR(10),T,HBAR,FIELD
       COMPLEX*8 H(N,N)
       HBAR=0.6582D0   !PLANCK'S CONSTANT IN eV*fs 



       IF(T.LE.RPAR(2)) THEN
         FIELD=RPAR(4)*EXP(-((T-RPAR(2))/RPAR(3))**2.0)*
     .      COS(RPAR(1)*(T-RPAR(2))/HBAR)
       ELSE IF((T.GT.RPAR(2)).and.(t.le.600d0)) THEN
         FIELD=RPAR(4)*COS(RPAR(1)*(T-RPAR(2))/HBAR)
       END IF


       DO I=1,N/2
          H(I,I)=-RPAR(5)/2.0
       END DO
       DO I=1+N/2,N
          H(I,I)=RPAR(5)/2.0
       END DO
       DO I=1,N-1
         H(I,I+1)=RPAR(6)*FIELD
         H(I+1,I)=RPAR(6)*FIELD
       END DO


       RETURN
       END SUBROUTINE










