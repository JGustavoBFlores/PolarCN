      PROGRAM POLAR
      INCLUDE 'commons'
      COMPLEX*16 PSI(N)
      COMPLEX*16 AMX(N,N),BMX(N,N)
      COMPLEX*16 AM2(N,N),BM2(N,N)
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
      BM2=BMX

      CALL ZGETRF(N,N,AM2,N,IPIV,INFO)
      CALL ZGETRI(N,AM2,N,IPIV,WORK,N,INFO)
C Test of the 
      DO I=1,N
      HOLD=0.0D0
      DO J=1,N
       HOLD=HOLD+AM2(I,J)*AMX(J,I)
      END DO
       PRINT*, HOLD
      END DO

C Calculate matrix D=A^-1
C Calculate Crank Nicholson Matrix C=D*B
C Multiply |Psi> by C 
C Write the potential initializer subroutine so it can be time dependent
      
      END PROGRAM

      SUBROUTINE AB_MATRIXES(AMX,BMX)
      INCLUDE 'commons'
      COMPLEX*16 AMX(N,N),BMX(N,N)
      
      AMX=0.0D0
      
      DO I=1,N
C First we will define the diagonal elements
       AMX(I,I)=COMPLEX(2.0D0,(hb/am)*dt/((r*dX)**2))
      END DO
C Now the off diagonals:
      DO I=1,N-1
       AMX(I,I+1)=COMPLEX(0.0D0,-((hb/(2*am))/((r*dX)**2)))
      END DO
      DO I=2,N
       AMX(I,I-1)=COMPLEX(0.0D0,-((hb/(2*am))/((r*dX)**2)))
      END DO

C Now we must put down the special elements for this to be a polar
C   Crank Nicholson Code.
       AMX(1,N)=COMPLEX(0.0D0,-((hb/(2*am))/((r*dX)**2)))
       AMX(N,1)=COMPLEX(0.0D0,-((hb/(2*am))/((r*dX)**2)))

C We now have our A matrix, and to define B is as easy as to conjugate A:

       BMX=CONJG(AMX)
      END SUBROUTINE

C The following subroutine sets a value for the complex wave function 
C   |Psi>
      SUBROUTINE INITIALIZE_PSI(N,PSI)
      IMPLICIT REAL*16 (A-H,O-Z)
      COMPLEX*16 PSI(N)
      PARAMETER(Pi=4*ATAN(1.0D0))
      PSI=(0.0D0,0.0D0)
      DO I=1,N
       PSI(I)=COMPLEX((1/SQRT(Pi))*SIN(I*2*Pi/N),0.0D0)
       IF(abs(dble(PSI(I))).LT.0.1D-10) PSI(I)=COMPLEX(0.0d0,0.0D0)
      END DO

      OPEN(UNIT=10,FILE='PSI')
      DO I=1,N
       WRITE(10,*) i,psi(i)
      END DO
      CLOSE(10)
      END SUBROUTINE
