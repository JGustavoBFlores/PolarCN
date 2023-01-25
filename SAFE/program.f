      PROGRAM POLAR
      INCLUDE 'commons'
      COMPLEX*16 PSI(N),PHI(N)
      COMPLEX*16 AMX(N,N),BMX(N,N)
      COMPLEX*16 AM2(N,N)
      CHARACTER*11 NAME
C The following are needed for lapack to work:
      EXTERNAL ZGETRF
      EXTERNAL ZGETRI
      COMPLEX*16 WORK(N)
      INTEGER IPIV(N),INFO,LWORK

C First we must initialize our wave function.
C  this must be a solution to the polar Schrodinger equation.
C  so it must be continuous along the cyclical run.
C Considering the radius to be constant, the solution to the polar SE
C  would be |Psi> = EXP[I*k*x]/SQRT[2Pi]; where I is the imaginary  unit i,
C  and k is the unique quantum number of our system.

      CALL INITIALIZE_PSI(N,PSI)
C NEXT SECTIONS:

C PLACEHOLDER for the potential function. for now we will consider a free particle.

C Calculate matrixes A and B
C   to make the last and firs elements depend on each other, 
C   its as easy as to consider in A and B the (N,1) and (1,N) element
C   to be non-zero. This is the step to have the polar CN Matrixes.
      CALL AB_MATRIXES(AMX,BMX)

      AM2=AMX

      CALL ZGETRF(N,N,AM2,N,IPIV,INFO)
      CALL ZGETRI(N,AM2,N,IPIV,WORK,N,INFO)

C Now AM2 holds the inverse of AMX. 
C To calculate the Crank Nicholson Matrix C=AM2*BMX we must do a product:
C As an extra step to save memory, lets save this new C matrix in AMX
      AMX=0.0D0
      DO I=1,N
      DO J=1,N
      DO K=1,N
       AMX(I,J)=AMX(I,J)+AM2(I,K)*BMX(J,K)
      END DO
      END DO
      END DO  

C     OPEN(UNIT=29,FILE='Matrix')
C     WRITE(29,*) AMX
C     CLOSE(29)
C     STOP
C     PRINT*, (AMX(1,J),J=1,3)
C     STOP
C To evolve |Psi> a unit of dT in time we must multiply it by C
      DO K=1,1000 
      PHI=0.0D0
      DO I=1,N
      DO J=1,N
       PHI(I)=PHI(I)+AMX(I,J)*PSI(J)
      END DO
      END DO

C Now we must save this into a file to print later
C because we can't really see imaginary parts of a wave function
C lets just save the probability density for now:
      WRITE(NAME,"('file',i4.4)")K
      OPEN(UNIT=20,FILE=NAME)
      DO I=1,N
      WRITE(20,*) COS((2*Pi)*I/N),SIN(I*(2*Pi)/N)
     &,dble(PHI(I)*CONJG(PHI(I)))
      END DO
      CLOSE(20)
      IF(I.EQ.(I/10.0)*10)CALL SYSTEM('mv file* EVOLUTION/.')
      PSI=PHI
      END DO
      
      END PROGRAM

      SUBROUTINE AB_MATRIXES(AMX,BMX)
      INCLUDE 'commons'
      COMPLEX*16 AMX(N,N),BMX(N,N)

C     PRINT*, dt,R,dX,dt/(R*dX)**2
      AMX=0.0D0
      
      DO I=1,N
C First we will define the diagonal elements
       AMX(I,I)=COMPLEX(2.0E0,2.0E0*dt/((r*dX)**2))
      END DO
C Now the off diagonals:
      DO I=1,N-1
       AMX(I,I+1)=COMPLEX(0.0E0,-(dt/((r*dX)**2)))
      END DO
      DO I=2,N
       AMX(I,I-1)=COMPLEX(0.0E0,-(dt/((r*dX)**2)))
      END DO

C Now we must put down the special elements for this to be a polar
C   Crank Nicholson Code.
       AMX(1,N)=COMPLEX(0.0E0,-(dt/((r*dX)**2)))
       AMX(N,1)=COMPLEX(0.0E0,-(dt/((r*dX)**2)))
 
C We now have our A matrix, and to define B is as easy as to conjugate A:

       BMX=CONJG(AMX)
      END SUBROUTINE

C The following subroutine sets a value for the complex wave function 
C   |Psi>
      SUBROUTINE INITIALIZE_PSI(N,PSI)
      IMPLICIT REAL*16 (A-H,O-Z)
      COMPLEX*16 PSI(N)
      PARAMETER(Pi=4*ATAN(1.0E0))
      PSI=(0.0E0,0.0E0)
      PSI(1)=COMPLEX(1.0E0,0.0E0)
C       DO I=1,90   
C      PSI(I)=COMPLEX(COS((2*I*Pi/N)),0.0E0)
C      IF(abs(dble(PSI(I))).LT.0.1E-10) PSI(I)=COMPLEX(0.0d0,0.0D0)
C      END DO
C      DO I=270,360
C      PSI(I)=COMPLEX(COS((2*I*Pi/N)),0.0E0)
C      IF(abs(dble(PSI(I))).LT.0.1E-10) PSI(I)=COMPLEX(0.0d0,0.0D0)
C      END DO

C     PSI(1)=(0.0E0,1.0E0)
      OPEN(UNIT=10,FILE='PSI')
      DO I=1,N
       WRITE(10,*) COS(2*Pi*I/N),SIN(2*Pi*I/N),psi(i)
      END DO
      OPEN(UNIT=11,FILE='PSI2')
      DO I=1,N
       WRITE(11,*) COS(2*Pi*I/N),SIN(2*Pi*I/N),
     &DBLE(PSI(i)*CONJG(PSI(I)))
      END DO
      CLOSE(11)
C     STOP
      END SUBROUTINE
