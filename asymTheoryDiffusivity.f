PROGRAM asymTheoryDiffusivity

!  Purpose:
!    This program calculates the diffusivities by solving
!	 the asymptotic theory equation from Aharon and Shaw
! 	 using Newton's method with damping.
!  Record of revisions:
!      Date       Programmer          Description of change
!      ====       ==========          =====================
!     07/19/17    C. Vang             Original code
!
IMPLICIT NONE

REAL :: dc_sq, K, do_measured, p, yo, dc, y_ofc, percent_increase, &
tau, H_0minus, H_1minus, H2_minus, h_match, phi, epsilon, h_0minus, &
h1_minus

dc_sq = 3.570341984
K 	  = 0.997561297
do_measured = 3.341
p     = 1
yo    = 0.05
dc    = SQRT(dc_sq)




END PROGRAM asymTheoryDiffusivity



!subroutine to calculate boilng point temperature at a specified pressure, Psat
!using Antoinne Equation
SUBROUTINE TsatAntoinne(Psat, A,B,C,D,Tmin,Tmax,Tsat)
	IMPLICIT NONE
	REAL, INTENT(IN) :: Psat, A, B, C, D, Tmin, Tmax
	REAL, INTENT(OUT) :: Tsat
	REAL :: To, FTOL, XTOL, eps, kmax, DX, &
	Fx, EuNormFx, lambda, DxKbar, T, Fxbar, Tbar, &
	EunormDxKbar, EunormDxK, Told, Dfx, DxK, temp, theta
	INTEGER :: l, lmax, k

	T = (Tmax + Tmin)/2.	!initial guess for Tsat
	eps = 1e-12		!error tolerance for norms
	FTOL = eps		!max norm tolerance	for func iterats
	XTOL = eps		!max norm tolerance for x iterates
	kmax = 100		!max iteration
	lmax = 15		!min damping factor
	k    = 0		!initialize iteration counter

	DX = 1			!initialize ||xk-xk-1||

	Fx  = (A*LOG(T)) + ((B/T)) + C + (D*T**2) - (LOG(Psat))
	Dfx = (A/T) - (B/(T**2)) + (2*D*T)

	!Compute Tsat using Newton's Method with damping
	IF (Fx == eps) THEN
		WRITE(*,*) "Stopped. Initial guess for Antoinne Eqn Produces F(x)=0 exactly."
		WRITE(*,40) T 
		40 FORMAT("Initial guess: " ES14.6)
	ELSE
		EuNormFx = SQRT( Fx**2 )  !initialize || fxk-fxk-1 ||
		DO
			DxK = -Fx/Dfx
			!Monotiniticty Test 
			Monotinicity : DO l=0,lmax
				IF (k == 0) THEN
					lambda = 1/(2**l)
				ELSE
					temp = 2*lambda/(2**l)
					lambda = MIN(1.0, temp)
				END IF
				!Solve for DeltaXkBar
				Tbar = T + lambda*DxK
				Fxbar = A*LOG(Tbar) + (B/Tbar) + C + D*Tbar**2 - LOG(Psat)
				DxKbar = -Fxbar / Dfx

				theta = 1. - (lambda/2.)
				EunormDxKbar = SQRT( DxKbar**2 )
				EunormDxK    = SQRT( DxK**2 )
				IF (EunormDxKbar <= theta*EunormDxK) EXIT
			END DO Monotinicity

			Told = T
			T = T + lambda*DxK
			DX = SQRT( (T - Told)**2 )

			Fx  = A*LOG(T) + (B/T) + C + D*T**2 - LOG(Psat)
			Dfx = (A/T) - B/(T**2) + 2*D*T

			EuNormFx = SQRT( (Fx)**2 )
			IF ( (EuNormFx <= FTOL) .OR. (DX <= XTOL) .OR. (k > kmax) ) EXIT
			k = k + 1
		END DO
	END IF 

	IF ( (T < Tmax).AND.(T > Tmin) ) THEN
		Tsat = T
	ELSE 
		WRITE(*,*) "Tsat computed is out of range of validity in Antoinnne Eqn!!!0"
	END IF
END SUBROUTINE TsatAntoinne


SUBROUTINE linspace(varArray, lowLimit, upLimit, numPts) 
	IMPLICIT NONE

	REAL, INTENT(OUT), DIMENSION(numPts) :: varArray
	REAL, INTENT(IN) :: lowLimit, upLimit
	INTEGER, INTENT(IN) :: numPts
	INTEGER :: i
	REAL :: intervalSize

	intervalSize = (upLimit - lowLimit) / numPts
	varArray(1) = lowLimit
	DO i = 2, numPts-1
		varArray(i) = varArray(i-1) + intervalSize
	END DO
	varArray(1) = lowLimit
	varArray(numPts) = upLimit
END SUBROUTINE linspace


