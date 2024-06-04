PROGRAM pr1
IMPLICIT NONE

INTEGER, parameter:: IO = 12
INTEGER NX, NT1, NT2, I, J, n, m
REAL, ALLOCATABLE :: U(:), UN(:), X(:)
REAL L, h, VNM, dt1, dt2, t, Time, thetha
REAL C0, C1, a1, a2, lambda1, lambda2, rho1, rho2, c_1, c_2, TL, TR

WRITE(*,*) 'Read input file' 
OPEN(IO,FILE='Input.txt')
READ(IO,*) L
READ(IO,*) NX
READ(IO,*) Time
READ(IO,*) VNM
READ(IO,*) C0, C1
READ(IO,*) lambda1
READ(IO,*) lambda2
READ(IO,*) rho1
READ(IO,*) rho2
READ(IO,*) c_1
READ(IO,*) c_2
READ(IO,*) TL
READ(IO,*) TR
READ(IO,*) thetha
CLOSE(IO)

ALLOCATE(U(1:NX), UN(1:NX),X(1:NX))      

h = 2*L/(NX-1)
a1 = lambda1/(rho1*c_1)
a2 = lambda2/(rho2*c_2)
dt1 = (VNM*h**2)/a1
dt2 = (VNM*h**2)/a2
dt1 = min(dt1, dt2)
dt2 = dt1
NT1 = Time/dt1
NT2 = Time/dt2

WRITE(*,*) 'L=',L, 'h=', h, 'NX=', NX
WRITE(*,*) 'VNM=', VNM, 'dt1=', dt1, 'dt2=', dt2,'Time=', Time, 'NT1=', NT1, 'NT2=', NT2, 'TL=', TL, 'TR=', TR
!pause

X(1)=-L
DO I=2, NX/2
    X(I)=X(I-1)+h
END DO
X(NX/2+1)=0.0
DO I=NX/2+2, NX-1
    X(I)=X(I-1)+h
END DO
X(NX)=L

U=0.0
UN=0.0

!НАЧАЛЬНЫЕ УСЛОВИЯ
DO I = 1, NX/2
	U(I)=TL
END DO

DO I = NX/2+2, NX
	U(I)=TR
END DO

U(NX/2+1)=(lambda1*U(NX/2) + lambda2*U(NX/2+2)) / (lambda1+lambda2)

!-------------------------  Solve equation ------------------
DO m = 1,MAX(NT1,NT2)
	U(NX/2+1)=(lambda1*U(NX/2) + lambda2*U(NX/2+2)) / (lambda1+lambda2)
	UN(NX/2+1) = U(NX/2+1)

	DO I = 2, NX/2
		UN(I)=U(I)+VNM*(U(I+1)-2*U(I)+U(I-1))
	END DO

    DO I =  NX-1, NX/2+2, -1
		UN(I)=U(I)+VNM*(U(I+1)-2*U(I)+U(I-1))
	END DO

	UN(NX)=TR
	UN(1)=TL

	U=UN

END DO


OPEN(IO,FILE='Res_explicit.dat')

WRITE(IO,*) 'VARIABLES = "X", "U" '
WRITE(IO, *) "ZONE I=", NX
DO I = 1, NX
	WRITE(IO,*) X(I), UN(I)
END DO

close(IO)
    print*, '1'
End program