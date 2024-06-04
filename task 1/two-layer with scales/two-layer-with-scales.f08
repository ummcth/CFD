PROGRAM pr1
IMPLICIT NONE

INTEGER, parameter:: IO = 12
INTEGER NX, NT1, NT2, I, J, n, m
REAL, ALLOCATABLE :: U(:), UN(:), X(:), A(:), B(:), C(:), D(:), alpha(:), beta(:)
REAL L, h, VNM, dt1, dt2, t_0, Time, thetha, K
REAL C0, C1, a_1, a_2, lambda1, lambda2, rho1, rho2, c_1, c_2, TL, TR

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

ALLOCATE(U(1:NX), UN(1:NX), X(1:NX))
ALLOCATE(A(1:NX), B(1:NX), C(1:NX), D(1:NX))
ALLOCATE(alpha(3:NX-1), beta(3:NX-1))

h = 2*L/(NX-1)
a_1 = lambda1/(rho1*c_1)
a_2 = lambda2/(rho2*c_2)
dt1 = (VNM*h**2)/a_1
dt2 = (VNM*h**2)/a_2
dt1 = min(dt1, dt2)
dt2 = dt1
NT1 = Time/dt1
NT2 = Time/dt2

WRITE(*,*) 'L=',L, 'h=', h, 'NX=', NX
WRITE(*,*) 'VNM=', VNM, 'dt1=', dt1, 'dt2=', dt2,'Time=', Time, 'NT1=', NT1, 'NT2=', NT2, 'TL=', TL, 'TR=', TR
!pause

X(1)=-L
DO I=2, NX-1
    X(I)=X(I-1)+h
END DO
X(NX)=L

U=0.0

!НАЧАЛЬНЫЕ УСЛОВИЯ
DO I = 2, NX/2
	U(I)=TL
END DO
DO I = NX/2+2, NX-1
	U(I)=TR
END DO


!ГРАНИЧНЫЕ УСЛОВИЯ
U(NX)=TR
U(1)=TL

UN = U

DO J=1, MAX(NT1,NT2)
	U(NX/2+1)=(lambda1*U(NX/2) + lambda2*U(NX/2+2)) / (lambda1+lambda2)
	!U(NX/2+1) = (U(NX/2) + U(NX/2+2)) / 2
	UN(NX/2+1) = U(NX/2+1)

	!КОЭФФИЦИЕНТЫ ПРОГОНКИ
	
	!ЛЕВАЯ ЧАСТЬ
	
	A(2) = 0.0
	B(2) = 1/dt1 + 2*thetha*a_1/h**2
	C(2) = -thetha*a_1/h**2
	D(2) = U(2)/dt1 + (((1 - thetha)*a_1)/h**2)*(U(1)-2*U(2)+U(3)) + (thetha * a_1 / h**2)*U(1)
	
	A(NX/2) = -thetha*a_1/h**2
	B(NX/2) = 1/dt1 + 2*thetha*a_1/h**2
	C(NX/2) = 0.0
	D(NX/2) = U(NX/2)/dt1 + (((1 - thetha)*a_1)/h**2)*(U(NX/2-1)-2*U(NX/2)+U(NX/2 + 1)) + (thetha * a_1 / h**2)*U(NX/2 + 1)

	DO I = 3, NX/2 - 1
	A(I) = -thetha*a_1/h**2
	B(I) = 1/dt1 + 2*thetha*a_1/h**2
	C(I) = -thetha*a_1/h**2
	D(I) = U(I)/dt1 + (((1 - thetha)*a_1)/h**2)*(U(I-1)-2*U(I)+U(I+1))
	END DO

	!ПРЯМОЙ ХОД
	alpha(3) = - C(2)/B(2)
	beta(3) = D(2)/B(2)

	DO I = 4, NX/2
		K = (B(I-1)+A(I-1)*alpha(I-1))
		alpha(I) = - C(I-1)/K
		beta(I) = (D(I-1)-A(I-1)*beta(I-1))/K
	END DO
	
	UN(NX/2) = (D(NX/2)-A(NX/2)*beta(NX/2))/(B(NX/2)+A(NX/2)*alpha(NX/2))

	!ОБРАТНЫЙ ХОД
	DO I = NX/2 - 1, 2, -1
		UN(I) = alpha(I+1)*UN(I+1)+beta(I+1)
	END DO
	
	!ПРАВАЯ ЧАСТЬ
	
	A(NX/2 + 2) = 0
	B(NX/2 + 2) = 1/dt2 + 2*thetha*a_2/h**2
	C(NX/2 + 2) = -thetha*a_2/h**2
	D(NX/2 + 2) = U(NX/2+2)/dt2 + (((1 - thetha)*a_2)/h**2)*(U(NX/2+1)-2*U(NX/2+2)+U(NX/2+3)) + (thetha * a_2 / h**2)*U(NX/2 + 1)
	
	A(NX-1) = -thetha*a_2/h**2
	B(NX-1) = 1/dt2 + 2*thetha*a_2/h**2
	C(NX-1) = 0.0
	D(NX-1) = U(NX-1)/dt2 + (((1 - thetha)*a_2)/h**2)*(U(NX-2)-2*U(NX-1)+U(NX)) + (thetha * a_2 / h**2)*U(NX)
	
	DO I = NX/2+3, NX-2
		A(I) = -thetha*a_2/h**2
		B(I) = 1/dt2 + 2*thetha*a_2/h**2
		C(I) = -thetha*a_2/h**2
		D(I) = U(I)/dt2 + (((1 - thetha)*a_2)/h**2)*(U(I-1)-2*U(I)+U(I+1))
	END DO
	
	!ПРЯМОЙ ХОД
	
	alpha(NX/2+3) = - C(NX/2+2)/B(NX/2+2)
	beta(NX/2+3) = D(NX/2+2)/B(NX/2+2)
	
	DO I = NX/2+4, NX-1
		K = (B(I-1)+A(I-1)*alpha(I-1))
		alpha(I) = - C(I-1)/K
		beta(I) = (D(I-1)-A(I-1)*beta(I-1))/K
	END DO
	
	UN(NX-1) = (D(NX-1)-A(NX-1)*beta(NX-1))/(B(NX-1)+A(NX-1)*alpha(NX-1))
	
	!ОБРАТНЫЙ ХОД
	DO I = NX-2, NX/2+2, -1
		UN(I) = alpha(I+1)*UN(I+1)+beta(I+1)
	END DO
	
	U = UN
END DO

OPEN(IO,FILE='Res_scales.plt')
!---------------------------Output results--------------------                       
WRITE(IO,*) 'VARIABLES = "X", "U" '
WRITE(IO, *) "ZONE I=", NX
DO I = 1, NX
	WRITE(IO,*) X(I), U(I)
END DO

close(IO)
    print*, '1'
END PROGRAM