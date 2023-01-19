      PROGRAM POLAR
      IMPLICIT REAL*8 (A-H,O-Z)
      PARAMETER (dT=1.2D0,dX=2.0D0,R=1.0D0,m=1.0-16)
      PARAMETER (N=360)
      COMPLEX PSI(N)
      
C First we must initialize our wave function.
C  this must be a solution to the polar Schrodinger equation.
C  so it must be continuous along the cyclical run.
C Considering the radius to be constant, the solution to the polar SE
C  would be |Psi> = EXP[I*k*x]/SQRT[2Pi]; where I is the imaginary unit i,
C  and k is the unique quantum number of our system.

      CALL INITIALIZE_PSI(N,PSI)
      
C NEXT SECTIONS:

C Calculate matrixes A and B
C Calculate matrix D=A^-1
C Calculate Crank Nicholson Matrix C=D*B
C Multiply |Psi> by C (Remember to make the first and last elements 
C  consider each other as neighbors
C Write the potential initializer subroutine so it can be time dependent
      
      END PROGRAM

      SUBROUTINE INITIALIZE_PSI(N,PSI)
      IMPLICIT REAL*8 (A-H,O-Z)
      COMPLEX PSI(N)
      PARAMETER(Pi=4*ATAN(1.0D0))
      PSI=(0.0D0,0.0D0)
      DO I=1,N
       PSI(I)=(1/SQRT(Pi))*SIN(I*2*Pi/N)
       IF(abs(dble(PSI(I))).LT.0.1D-10) PSI(I)=(0.0d0,0.0D0)         
      END DO

      OPEN(UNIT=10,FILE='PSI')
      DO I=1,N
       WRITE(10,*) i,psi(i)
      END DO
      CLOSE(10)
      END SUBROUTINE
