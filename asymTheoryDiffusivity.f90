PROGRAM asymTheoryDiffusivity

!  Purpose:
!    This program calculates the diffusivities by solving
!	 the asymptotic theory equation from Aharon and Shaw
! 	 using Newton's method with damping.
!  Record of revisions:
!      Date       Programmer          Description of change
!      ====       ==========          =====================
!     07/20/17    C. Vang             Original code
!
IMPLICIT NONE

! INTEGER, PARAMTER :: DP = SELECTED_REAL_KIND(14)
REAL(8) :: dc_sq, K, do_measured, p, yo, d_c, y_ofc, percent_increase, &
d_o, tau, LHS, epsilon, eps_ig, D, err_tol, EuNormFx, ep_min, ep_max
REAL(8), DIMENSION(6) :: fc_props
REAL(8), DIMENSION(500) :: fvalues, ep_vals
CHARACTER(len=11) :: filename
CHARACTER(len=100) :: junk
INTEGER :: i, status, N

!Get the file name and echo back to user
! WRITE(*,*) "Enter the name of the file with experimental measurements: "
! READ(*,10) filename
! 10 FORMAT(A30)

filename = "fcprops.txt"

WRITE(*,20) filename
20 FORMAT(' ','The input file name is: ' A)

!Open the file and check for errors on Open 
OPEN(UNIT=1, FILE=filename, STATUS='OLD', ACTION='READ', IOSTAT=status)

openif : IF (status == 0) THEN
	WRITE(*,*) 'OPENING OF FILE IS SUCCESSFUL!'

	readjunkLines: DO i=1,4 			!read lines 1-4
	  READ(1,*,IOSTAT=status) junk
	END DO readjunkLines

	readData: DO i=1,6 					!read lines 5-9
	  READ(1,*,IOSTAT=status) fc_props(i)
	END DO readData

	p 			= fc_props(1)  !chamber pressure [atm]
	dc_sq 		= fc_props(2)  !drop diameter @ onset of fc squared [mm^2]
	K 			= fc_props(3)  !burning rate constant prior to onset of fc
	do_measured = fc_props(4)  !initial drop diameter measured [mm]
	yo 			= fc_props(5)  !initial mass frac of low volitiliy component
	err_tol 	= fc_props(6)  !error tolerrance for newton's method

	CLOSE(UNIT=1)
ELSE

	WRITE(*,*) 'Open of file NOT sucessful!'

ENDIF openif

d_c    = SQRT(dc_sq)   !drop diameter at onset of flame contraction [mm]
y_ofc = 1. 			   !mass fraction of low volatility comp at onset of fc

IF (p == 1) THEN
	percent_increase = 0.055
ELSE
	percent_increase = 0.086
ENDIF

d_o = do_measured * (1.0 + percent_increase) !do accounting for droplet swelling
tau = LOG(d_o / d_c)
LHS = y_ofc / yo

!prompt user to enter initial guess to start Newton iteration
WRITE(*,*) "Enter an initial guess for epsilon (0.01-0.1): "
READ(*,*) eps_ig

EuNormFx = 1.  !initialize function residual (as sanity check)

DO WHILE ( EuNormFx > err_tol )

	CALL NewtonSolve(tau, LHS, eps_ig, epsilon, err_tol, EuNormFx)
	WRITE(*,30) epsilon
	30 FORMAT(" The value of epsilon is....................." ES14.6)

	IF ( EuNormFx > err_tol) THEN 
		WRITE(*,*) "!! INITIAL GUESS IS TOO FAR AWAY. TRY DIFFERENT INITIAL GUESS !!"
		WRITE(*,*) "Enter an initial guess for epsilon (0.01-0.1): "
		READ(*,*) eps_ig	
	END IF	

END DO

D = epsilon * K/ 8.0
WRITE(*,50) D
50 FORMAT(" The effective liquid diffusivity is........." ES14.6, " mm^2/s")
WRITE(*,52) D*(1.0/1000.)**2.
52 FORMAT(" The effective liquid diffusivity is........." ES14.6, " m^2/s")

!plot function around epsilon as a visual check for multiple solutions nearby
ep_min = epsilon / 4.
ep_max = epsilon*10
N      = 500

CALL linspace(ep_vals, ep_min, ep_max, N)

OPEN(UNIT=10,FILE='asymTheoryData')
DO i=1,N 
	CALL Fx_eval(tau, ep_vals(i), LHS, fvalues(i) )  
	WRITE(10,*) ep_vals(i), fvalues(i)
END DO
CLOSE(UNIT = 10)

CALL SYSTEM('gnuplot -p data_plot.plt')

END PROGRAM asymTheoryDiffusivity



SUBROUTINE NewtonSolve(tau, LHS, eps_ig, epsilon, err_tol, EuNormFx)
	IMPLICIT NONE 
	REAL(8), INTENT(IN) :: tau, LHS, eps_ig, err_tol
	REAL(8), INTENT(OUT) :: epsilon, EuNormFx
	REAL(8) :: FTOL, XTOL, DX, &
	lambda, DxKbar, temp, eps_bar, Fxbar, &
	theta, EunormDxKbar, EunormDxK, eps_old, DxK

	REAL(8) :: Fx, Dfx, eps
	
	REAL(8), PARAMETER :: PI = 3.1415927

	INTEGER :: l, lmax, kmax, k 

	eps = eps_ig 		!set eps equal to initial guess
	FTOL = err_tol		!max norm tolerance	for func iterats
	XTOL = err_tol		!max norm tolerance for x iterates
	kmax = 100			!max iteration
	lmax = 15			!min damping factor

	CALL Fx_eval(tau, eps, LHS, Fx)  			
	CALL Dfx_eval(tau, eps, Dfx)

	DX = 1.0			!initialize norm between successive iterates ||x^k-x^(k-1)||

	!Compute eps using Newton's Method with damping. Information on this
	!particular method as well as some notation consistent with
	!the variables used here can be found in Stoer and Bulirsch
	IF (Fx == eps) THEN
		WRITE(*,*) "Stopped. Initial guess produces F(x)=0 exactly."
		WRITE(*,60) eps
		60 FORMAT("' ',Initial guess: " ES14.6)
	ELSE
		EuNormFx = SQRT( Fx**2 )  !initialize || fx^k - fx^(k-1) ||

		newtonIterates : DO k = 0, kmax     
			DxK = -Fx/Dfx

			!Monotiniticty Test 
			Monotinicity : DO l=0,lmax
				IF (k == 0) THEN
					lambda = 1./(2.**l)
				ELSE
					temp = 2.*lambda/(2.**l)
					lambda = MIN(1.0, temp)
				END IF

				!Solve for DeltaXkBar
				eps_bar = eps + lambda*DxK
				CALL Fx_eval(tau, eps_bar, LHS, Fxbar)
				DxKbar = -Fxbar / Dfx

				theta = 1. - (lambda/2.)
				EunormDxKbar = SQRT( DxKbar**2 )
				EunormDxK    = SQRT( DxK**2 )
				IF (EunormDxKbar <= theta*EunormDxK) EXIT
			END DO Monotinicity

			eps_old = eps 
			eps = eps + lambda*DxK 

			DX = SQRT( (eps - eps_old)**2 )
			CALL Fx_eval(tau, eps, LHS, Fx)  			
			CALL Dfx_eval(tau, eps, Dfx)
			EuNormFx = SQRT( (Fx)**2 )	

			WRITE(*,70) k, EuNormFx
			70 FORMAT(' ',"iteration count:" I3, "   function residual:" ES14.6)			
			
			! IF ( (EuNormFx <= FTOL) .OR. (DX <= XTOL) .OR. (k > kmax) ) EXIT
			IF ( (EuNormFx <= FTOL) .OR. (k > kmax) ) EXIT

		END DO newtonIterates
	END IF 
	epsilon = eps
END SUBROUTINE NewtonSolve

SUBROUTINE Fx_eval(tau, eps, LHS, Fx)
	IMPLICIT NONE
	REAL(8), PARAMETER :: PI = 3.1415927
	REAL(8), INTENT(IN) :: tau, eps, LHS
	REAL(8), INTENT(OUT) :: Fx
	REAL(8) :: H0minus, H1minus, H2minus, h_0minus, h_1minusIntegral, &
	h_1minus, h_match, Dfx

	!define parts of asymptotic equation that do not depend on epsilon
	H0minus = ( EXP(3.*tau) - 1. ) / 3.
	H1minus = ( EXP(3.*tau) + 2. ) / 3.
	H2minus = (22. / 9.) +  ( 14.-24.*tau )*EXP(3*tau) / 9.

	!define parts of asymptotic equation that depend on epsilon
	h_0minus = ((tau/eps)/2.) * ( 1. + ERF( SQRT(tau/eps)/2. ) ) &
				+ ERF( SQRT(tau/eps)/2. ) &
				+ SQRT( (tau/eps)/PI ) * EXP( -(tau/eps) / 4. )
	!NOTE: h_1minusIntegral is the alernate form of the integral in
	!      the variable h_1minus, as given by Mathematica
	h_1minusIntegral = (1./2.) * SQRT(PI) * SQRT(tau/eps) * ( -4. + (3.*tau/eps) ) &
				- (1./4.)*EXP( (tau/eps) / 4. ) * PI & 
				* ( 8. + ( tau*(2.*eps + 3.*tau) ) / eps**2. ) &
				* ERFC( (tau/eps) / 2. )
	h_1minus = 4. + (tau/eps) + ( 3.*(tau/eps)**2. / 2. ) &
				- 2. * EXP( -(tau/eps) / 4.  ) &
				+ ( EXP( -(tau/eps) / 4. )/ PI ) * h_1minusIntegral
	h_match = 1.0 + (tau/eps) + eps * ( 4. + (tau/eps) + (3./2.)*(tau**2. / eps**2.) )
	Fx  = 1.0 + h_0minus + eps * h_1minus + (H0minus/eps) + &
		H1minus + eps * H2minus - h_match - LHS
END SUBROUTINE Fx_eval


SUBROUTINE Dfx_eval(tau, eps, Dfx)
	IMPLICIT NONE 
	REAL(8), PARAMETER :: PI = 3.1415927
	REAL(8), INTENT(IN) :: tau, eps 
	REAL(8), INTENT(OUT) :: Dfx

	!be careful not to change anything here!
	Dfx = EXP( -tau/(4*eps) ) * ( -2. * ( -EXP(tau/(4*eps))*SQRT(PI)*(6.+44*eps**2 + 9*tau) &
		+ 9*eps*( SQRT(PI)*(4*eps + tau) + SQRT(tau/eps)*(2. + 4*eps + 3*tau ) ) &
		+ 2*EXP( (tau/4.)*( 12.+(1./eps) ) )*SQRT(PI)*( 3. + 2*eps**2 * (-7. + 12*tau) ) ) & 
		- 18*EXP ( tau/(4*eps) )*SQRT(PI)*tau*ERF( SQRT(tau/eps) / 2. ) & 
		- 9*EXP( tau/(4*eps) )*SQRT(PI)*( 8*eps**2 - 3*tau**2 )*ERFC( SQRT(tau/eps)/2. )  ) 
END SUBROUTINE Dfx_eval


SUBROUTINE linspace(varArray, lowLimit, upLimit, numPts) 
	IMPLICIT NONE

	REAL(8), INTENT(OUT), DIMENSION(numPts) :: varArray
	REAL(8), INTENT(IN) :: lowLimit, upLimit
	INTEGER, INTENT(IN) :: numPts
	INTEGER :: i
	REAL(8) :: intervalSize

	intervalSize = (upLimit - lowLimit) / numPts
	varArray(1) = lowLimit
	DO i = 2, numPts-1
		varArray(i) = varArray(i-1) + intervalSize
	END DO
	varArray(numPts) = upLimit
END SUBROUTINE linspace



