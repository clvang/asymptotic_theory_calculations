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
REAL(KIND=8) :: dc_sq, K, do_measured, p, yo, d_c, y_ofc, percent_increase, &
d_o, tau, LHS, epsilon, D, err_tol, EuNormFx, ep_min, ep_max, dfdo, epsilon_bisect
REAL(KIND=8), DIMENSION(6) :: fc_props
REAL(KIND=8), DIMENSION(50) :: ig_vector
REAL(KIND=8), DIMENSION(700) :: fvalues, ep_vals
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

!!!!!!!!!!!!!!!!!!!!!!!!!! SOLUTION USING NEWTON'S METHOD !!!!!!!!!!!!!!!!!!!!!!!!!!
ep_min = 0.01 
ep_max = 0.5
N = 50
!create array filled with intial guesses
CALL linspace(ig_vector, ep_min, ep_max, N)
EuNormFx = 1.  !initialize function residual (as sanity check)
i = 1
DO WHILE ( EuNormFx > err_tol )
	WRITE(*,22) ig_vector(i)
	22 FORMAT(' ','============ Initial guess: ' ES10.3,'=============')

	CALL NewtonSolve(tau, LHS, ig_vector(i), epsilon, err_tol, EuNormFx)

	IF ( i > N ) WRITE(*,*) "!! NUMBER OF INITIAL GUESSES EXHAUSTED !!"
	i = i + 1
END DO

WRITE(*,*) " ======================Newton Iteration Results========================"
WRITE(*,30) epsilon
30 FORMAT(" || The value of epsilon is....................." ES14.6, "        ||")
D = epsilon * K/ 8.0
WRITE(*,50) D
50 FORMAT(" || The effective liquid diffusivity is........." ES14.6, " mm^2/s ||")
WRITE(*,52) D*(1.0/1000.)**2.
52 FORMAT(" || The effective liquid diffusivity is........." ES14.6, " m^2/s  ||")
WRITE(*,*) " ====================================================================="

!plot function around epsilon as a visual check for multiple solutions nearby
ep_min = epsilon / 4.
ep_max = epsilon*100
N      = 700
CALL linspace(ep_vals, ep_min, ep_max, N)
OPEN(UNIT=10,FILE='asymTheoryData')
DO i=1,N 
	CALL Fx_eval(tau, ep_vals(i), LHS, fvalues(i) )  
	WRITE(10,*) ep_vals(i), fvalues(i)
END DO
CLOSE(UNIT = 10)
!call gnuplot through shell to plot results
CALL SYSTEM('gnuplot -p data_plot.plt')

!!!!!!!!!!!!!!!!!!!!!!!!!! SOLUTION USING BISECTION METHOD !!!!!!!!!!!!!!!!!!!!!!!!!!
CALL bisection(tau, LHS, err_tol, epsilon_bisect)
WRITE(*,*) " ==================== Bisection Iteration Results ====================="
WRITE(*,54) epsilon_bisect
54 FORMAT(" || The value of epsilon is....................." ES14.6, "        ||")
D = epsilon_bisect * K/ 8.0
WRITE(*,56) D
56 FORMAT(" || The effective liquid diffusivity is........." ES14.6, " mm^2/s ||")
WRITE(*,58) D*(1.0/1000.)**2.
58 FORMAT(" || The effective liquid diffusivity is........." ES14.6, " m^2/s  ||")
WRITE(*,*) " ====================================================================="

!calculate uncertainties in eff diffusivity using TSM
!calculate df_do
! CALL partialF_partial_do(d_o, dfdo, d_c, yo, y_ofc, err_tol)


END PROGRAM asymTheoryDiffusivity



SUBROUTINE partialF_partial_do(d_o, dfdo_out, d_c, yo, y_ofc, err_tol)
	IMPLICIT NONE
	REAL(KIND = 8), INTENT(IN) :: d_o, d_c, yo, y_ofc, err_tol
	REAL(KIND = 8), INTENT(OUT) :: dfdo_out
	REAL(KIND = 8), DIMENSION(4) :: delta_do, do_1, do_2, tau_1, tau_2, df_do
	REAL(KIND = 8) :: p_incr, successive_norm
	INTEGER :: i, N 

	p_incr = 0.01
	N = 4
	DO i=1,N 
		delta_do(i) = p_incr / ( i*2.0 )
		do_1(i) = d_o + delta_do(i) 
		do_2(i) = d_o - delta_do(i)
		tau_1(i) = LOG( do_1(i) / d_c )
		tau_2(i) = LOG( do_2(i) / d_c )

		CALL partialF_partial_x( do_1(i), do_2(i), d_c, yo, yo, delta_do(i), y_ofc, df_do(i), err_tol )
	END DO 

	DO i=2,N
		successive_norm = ABS( df_do(i) - df_do(i-1) )
		WRITE(*,10) successive_norm
		10 FORMAT(' ','Successive Norm: ', ES14.6)
	END DO
	
	dfdo_out = df_do(N) 
	WRITE(*,20) dfdo_out
	20 FORMAT(' ','df_do approximately: ', ES14.6)

END SUBROUTINE partialF_partial_do


SUBROUTINE partialF_partial_x(do_1, do_2, d_c, yo_1, yo_2, delta_x, y_ofc, df_dx, err_tol)
	IMPLICIT NONE 
	REAL(KIND = 8), INTENT(IN) :: do_1, do_2, d_c, yo_1, yo_2, &
	delta_x, y_ofc, err_tol
	REAL(KIND = 8), INTENT(OUT) :: df_dx
	! REAL(KIND = 8) :: tau_1, tau_2, epsilon_1, epsilon_2, &
	! EuNormFx_1, EuNormFx_2, LHS_1, LHS_2, ep_min, ep_max
	REAL(KIND = 8) :: tau_1, tau_2, epsilon_1, epsilon_2, &
	LHS_1, LHS_2
	REAL(KIND=8), DIMENSION(100) :: ig_vector
	INTEGER :: N, i

	tau_1 = LOG(do_1 / d_c)
	tau_2 = LOG(do_2 / d_c)
	LHS_1 = y_ofc / do_1
	LHS_2 = y_ofc / do_2

	! ep_min = 0.001 
	! ep_max = 0.5
	! N = 100
	! !create array filled with intial guesses
	! CALL linspace(ig_vector, ep_min, ep_max, N)

	! DO i=1,N 
	! 	CALL NewtonSolve(tau_1, LHS_1, ig_vector(i), epsilon_1, err_tol, EuNormFx_1)
	! END DO 

	! DO i=1,N 
	! 	CALL NewtonSolve(tau_2, LHS_2, ig_vector(i), epsilon_2, err_tol, EuNormFx_2)
	! END DO 

	CALL bisection(tau_1, LHS_1, err_tol, epsilon_1)
	CALL bisection(tau_2, LHS_2, err_tol, epsilon_2)

	df_dx = (epsilon_2 - epsilon_1) / (2. * delta_x)		

END SUBROUTINE partialF_partial_x

SUBROUTINE bisection(tau, LHS, err_tol, c)
	IMPLICIT NONE
	REAL(KIND=8), INTENT(IN) :: tau, LHS, err_tol
	REAL(KIND=8), INTENT(OUT) :: c
	REAL(KIND=8) :: a_o, b_o, f_c, f_ao, f_bo, b_previous
	INTEGER :: maxit, i 

	a_o = 0.001
	b_o = 0.5
	maxit = 300
	DO i=1,maxit
		c = 0.5*(a_o + b_o)
		CALL Fx_eval(tau, c, LHS, f_c)  		
		CALL Fx_eval(tau, a_o, LHS, f_ao)
		CALL Fx_eval(tau, b_o, LHS, f_bo) 

		IF ( ABS(f_c) < err_tol ) THEN
			EXIT 
		ELSE IF ( f_ao * f_bo < 0 ) THEN
			b_previous = b_o 
			b_o = c 
		ELSE 
			a_o = c 
			b_o = b_previous
		END IF  

		WRITE(*,20) i, ABS(f_c)
		20 FORMAT(' ',"iteration count:" I3, "   function residual:" ES14.6)	

		IF ( ABS(f_c) < err_tol )  EXIT  

	END DO 
END SUBROUTINE bisection


SUBROUTINE NewtonSolve(tau, LHS, eps_ig, epsilon, err_tol, EuNormFx)
	IMPLICIT NONE 
	REAL(KIND=8), INTENT(IN) :: tau, LHS, eps_ig, err_tol
	REAL(KIND=8), INTENT(OUT) :: epsilon, EuNormFx
	REAL(KIND=8) :: FTOL, XTOL, DX, &
	lambda, DxKbar, temp, eps_bar, Fxbar, &
	theta, EunormDxKbar, EunormDxK, eps_old, DxK

	REAL(KIND=8) :: Fx, Dfx, eps
	
	REAL(KIND=8), PARAMETER :: PI = 3.1415927

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
	REAL(KIND=8), PARAMETER :: PI = 3.1415927
	REAL(KIND=8), INTENT(IN) :: tau, eps, LHS
	REAL(KIND=8), INTENT(OUT) :: Fx
	REAL(KIND=8) :: H0minus, H1minus, H2minus, h_0minus, h_1minusIntegral, &
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
	REAL(KIND=8), PARAMETER :: PI = 3.1415927
	REAL(KIND=8), INTENT(IN) :: tau, eps 
	REAL(KIND=8), INTENT(OUT) :: Dfx

	!be careful not to change anything here!
	Dfx = EXP( -tau/(4*eps) ) * ( -2. * ( -EXP(tau/(4*eps))*SQRT(PI)*(6.+44*eps**2 + 9*tau) &
		+ 9*eps*( SQRT(PI)*(4*eps + tau) + SQRT(tau/eps)*(2. + 4*eps + 3*tau ) ) &
		+ 2*EXP( (tau/4.)*( 12.+(1./eps) ) )*SQRT(PI)*( 3. + 2*eps**2 * (-7. + 12*tau) ) ) & 
		- 18*EXP ( tau/(4*eps) )*SQRT(PI)*tau*ERF( SQRT(tau/eps) / 2. ) & 
		- 9*EXP( tau/(4*eps) )*SQRT(PI)*( 8*eps**2 - 3*tau**2 )*ERFC( SQRT(tau/eps)/2. )  ) 
END SUBROUTINE Dfx_eval


SUBROUTINE linspace(varArray, lowLimit, upLimit, numPts) 
	IMPLICIT NONE

	REAL(KIND=8), INTENT(OUT), DIMENSION(numPts) :: varArray
	REAL(KIND=8), INTENT(IN) :: lowLimit, upLimit
	INTEGER, INTENT(IN) :: numPts
	INTEGER :: i
	REAL(KIND=8) :: intervalSize

	intervalSize = (upLimit - lowLimit) / numPts
	varArray(1) = lowLimit
	DO i = 2, numPts-1
		varArray(i) = varArray(i-1) + intervalSize
	END DO
	varArray(numPts) = upLimit
END SUBROUTINE linspace



