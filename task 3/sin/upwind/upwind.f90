PROGRAM Pr
IMPLICIT NONE
INTEGER, parameter:: IO = 12 ! input-output unit        
INTEGER NX, NT, I, J, m
REAL(8),ALLOCATABLE :: U(:), UN(:), X(:), UI(:)
REAL(8) L, h, CFL, dt, t, Time
REAL(8) C0, C1, pi, k, C

WRITE(*,*) 'Read input file' 
OPEN(IO,FILE='Input.txt')
READ(IO,*) L
READ(IO,*) m
READ(IO,*) C0, C1
READ(IO,*) NX
READ(IO,*) NT
READ(IO,*) CFL
READ(IO,*) C
CLOSE(IO)      
      
ALLOCATE(U(0:NX+1), UN(0:NX+1), X(0:NX+1), UI(0:NX+1))

pi=3.14159265359d0
k = m*pi/L
h = L/(NX-1)
dt =  (CFL*h)/C
Time = dt*NT

WRITE(*,*) 'L=', L, 'h=', h, 'NX=', NX
WRITE(*,*) 'CFL=', CFL, 'dt=', dt, 'Time=', Time, 'NT=', NT

X(0)=-h
DO I=1, NX+1
      X(I)=X(I-1)+h
END DO
      
U(:)=0.0
UN(:)=0.0

t=0.0d0

DO I=0, NX+1
      U(I)=sin(m*pi*X(I))
      UI(I)=sin(m*pi*X(I))
END DO

!-------------------------  Solve equation ------------------
DO J=1,NT
	DO I=1, NX
 		IF (U(I)>=0) THEN
			UN(I)=U(I)-(CFL/2)*(U(I)**2-U(I-1)**2)
        ELSE
            UN(I)=U(I)-(CFL/2)*(U(I+1)**2-U(I)**2)
        END IF 
	END DO
	U=UN
END DO

!-------------------------  Output results ------------------
OPEN(IO,FILE='Res_I.dat')
WRITE(IO,*) 'VARIABLES = "X", "U" '
WRITE(IO, *) "ZONE I=", NX
DO I = 1, NX
	WRITE(IO,*) X(I), UI(I)
END DO

OPEN(IO,FILE='Res.dat')
WRITE(IO,*) 'VARIABLES = "X", "U" '
WRITE(IO, *) "ZONE I=", NX
DO I = 1, NX
	WRITE(IO,*) X(I), U(I)
END DO

CLOSE(IO)      
            
END PROGRAM