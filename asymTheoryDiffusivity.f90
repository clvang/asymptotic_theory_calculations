PROGRAM asymTheoryDiffusivity

    !  Purpose:
    !    This program calculates the diffusivities by solving
    !	 the asymptotic theory equation from Aharon and Shaw
    ! 	 using the Bisection Method.
    !  Record of revisions:
    !      Date       Programmer          Description of change
    !      ====       ==========          =====================
    !     07/20/17    C. Vang             Original code
    !
    IMPLICIT NONE

    REAL(KIND=8) :: dc_sq, K, do_measured, p, yo, d_c, y_ofc, percent_increase, &
        d_o, tau, LHS, D, err_tol, ep_min, ep_max, dfdo, epsilon_bisect, dfdc, &
        dfdY, Uk_sq, Udo_sq, Udc_sq, UYo_sq, Ueps, Udiff
    REAL(KIND=8), DIMENSION(10) :: fc_props
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
20  FORMAT(' ','The input file name is: ' A)

    !Open the file and check for errors on Open 
    OPEN(UNIT=1, FILE=filename, STATUS='OLD', ACTION='READ', IOSTAT=status)

    openif : IF (status == 0) THEN
        WRITE(*,*) 'OPENING OF FILE IS SUCCESSFUL!'

        readjunkLines: DO i=1,4 			!read lines 1-4
            READ(1,*,IOSTAT=status) junk
        END DO readjunkLines

        readData: DO i=1,10 					!read lines 5-9
            READ(1,*,IOSTAT=status) fc_props(i)
        END DO readData

        p 			= fc_props(1)  !chamber pressure [atm]
        dc_sq 		= fc_props(2)  !drop diameter @ onset of fc squared [mm^2]
        K 			= fc_props(3)  !burning rate constant prior to onset of fc
        do_measured = fc_props(4)  !initial drop diameter measured [mm]
        yo 			= fc_props(5)  !initial mass frac of low volitiliy component
        err_tol 	= fc_props(6)  !error tolerrance for newton's method
        Uk_sq       = fc_props(7)**2
        Udo_sq      = fc_props(8)**2
        Udc_sq      = fc_props(9)**2
        UYo_sq      = (fc_props(10)*yo)**2
        CLOSE(UNIT=1)
    ELSE
        WRITE(*,*) 'Open of file NOT sucessful!'
    ENDIF openif

    d_c    = SQRT(dc_sq)   !drop diameter at onset of flame contraction [mm]
    y_ofc = 1.0 			   !mass fraction of low volatility comp at onset of fc

    IF (p == 1) THEN
        percent_increase = 0.055
    ELSE
        percent_increase = 0.086
    ENDIF
    d_o = do_measured * (1.0 + percent_increase) !do accounting for droplet swelling
    tau = LOG(d_o / d_c)
    LHS = y_ofc / yo

    !!!!!!!!!!!!!!!!!!!!!!!!!! SOLUTION USING BISECTION METHOD !!!!!!!!!!!!!!!!!!!!!!!!!!
    WRITE(*,*) '************** BEGIN CALCULATING EPSILON ***************'    
    CALL bisection(tau, LHS, err_tol, epsilon_bisect)
    WRITE(*,*) '************** END CALCULATING EPSILON ***************'       

    !plot function around epsilon as a visual check for multiple solutions nearby
    ep_min = (epsilon_bisect )/ 4.0  
    ep_max = epsilon_bisect*100     
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

    !!!!!!!!!!!!!!!!!!!! CALCULATE UNCERTAINTIES IN D USING TSM !!!!!!!!!!!!!!!!!!!!!!
    !calculate df_do
    CALL partialF_partial_do(d_o, dfdo, d_c, yo, y_ofc, err_tol)
    CALL partialF_partial_dc(d_o, dfdc, d_c, yo, y_ofc, err_tol)  
    CALL partialF_partial_dY(d_o, dfdY, d_c, yo, y_ofc, err_tol)  

    WRITE(*,*) " ==================== Bisection Iteration Results ====================="
    WRITE(*,54) epsilon_bisect
54  FORMAT(" || The value of epsilon is....................." ES14.6, "        ||")
    D = epsilon_bisect * K/ 8.0
    WRITE(*,56) D
56  FORMAT(" || The effective liquid diffusivity is........." ES14.6, " mm^2/s ||")
    WRITE(*,58) D*(1.0/1000.)**2.
58  FORMAT(" || The effective liquid diffusivity is........." ES14.6, " m^2/s  ||")
    WRITE(*,*) " ====================================================================="  

    CALL uncertainty_diffusivity(dfdo, Udo_sq, dfdc, Udc_sq, dfdY, UYo_sq, D, Uk_sq, K, epsilon_bisect, Ueps, Udiff)

    WRITE(*,*) " ------------------------ Equation Sensitivies ------------------------"
    WRITE(*,60) dfdo
60  FORMAT(" || Sensitivity in d_o (df_do) is:.............." ES14.6, "        ||")
    WRITE(*,70) dfdc
70  FORMAT(" || Sensitivity in d_c (df_dc) is:.............." ES14.6, "        ||")
    WRITE(*,80) dfdY
80  FORMAT(" || Sensitivity in y_o (df_yo) is:.............." ES14.6, "        ||")
    WRITE(*,90) Ueps
90  FORMAT(" || Uncertainty in epsilon (U_eps) is:.........." ES14.6, "        ||")
    WRITE(*,100) Udiff
100 FORMAT(" || Uncertainty in D (U_diff) is:..............." ES14.6, "        ||")
    WRITE(*,*) " ----------------------------------------------------------------------"      

END PROGRAM asymTheoryDiffusivity

SUBROUTINE uncertainty_diffusivity(dfdo, Udo_sq, dfdc, Udc_sq, dfdY, UYo_sq, D, Uk_sq, K, epsilon, Ueps, Udiff)
    IMPLICIT NONE 
    REAL(KIND=8), INTENT(IN) :: dfdo, Udo_sq, dfdc, Udc_sq, dfdY, UYo_sq, D, Uk_sq, K, epsilon 
    REAL(KIND=8), INTENT(OUT) :: Ueps, Udiff 

    Ueps = SQRT( (dfdo**2)*(Udo_sq) + (dfdc**2)*(Udc_sq) + (dfdY**2)*(UYo_sq) )

    Udiff = SQRT( ( Uk_sq/(K*K) ) + ( Ueps*Ueps/(epsilon*epsilon) ) ) * D
END SUBROUTINE uncertainty_diffusivity 


SUBROUTINE partialF_partial_do(d_o, dfdo_out, d_c, yo, y_ofc, err_tol)
    IMPLICIT NONE
    REAL(KIND=8), INTENT(IN) :: d_o, d_c, yo, y_ofc, err_tol
    REAL(KIND=8), INTENT(OUT) :: dfdo_out
    REAL(KIND=8), DIMENSION(8) :: delta_do, do_1, do_2, df_do !tau_1, tau_2, 
    REAL(KIND=8) :: successive_norm
    INTEGER :: i, N 

    WRITE(*,*) '************** BEGIN CALCULATING DF_DO ***************'

    N = 8
    ! TODO: add restrction on do_1 and do_2 s.t. do/d_c > 0 iff tau > 0
    ! do_max = d_c 
    ! do_min = d_o / 4.0

    delta_do(1) = 0.01*d_o
    DO i=1,N 
        delta_do(i+1) = 0.01*d_o / ( (i+1)*2.0 )
        do_1(i) = d_o + delta_do(i) 
        do_2(i) = d_o - delta_do(i)

        CALL partialF_partial_x( do_1(i), do_2(i), d_c, d_c, yo, yo, delta_do(i), y_ofc, df_do(i), err_tol )
    END DO 

    DO i=2,N
        successive_norm = ABS( df_do(i) - df_do(i-1) )
        WRITE(*,10) successive_norm
10      FORMAT(' ','Successive Norm between df_do Values: ', ES14.6)
    END DO
    dfdo_out = df_do(N) 
    WRITE(*,*) '************** END CALCULATING DF_DO ***************'

END SUBROUTINE partialF_partial_do

SUBROUTINE partialF_partial_dc(d_o, dfdc_out, d_c, yo, y_ofc, err_tol)
    IMPLICIT NONE
    REAL(KIND=8), INTENT(IN) :: d_o, d_c, yo, y_ofc, err_tol
    REAL(KIND=8), INTENT(OUT) :: dfdc_out
    REAL(KIND=8), DIMENSION(8) :: delta_dc, dc_1, dc_2, df_dc !tau_1, tau_2, 
    REAL(KIND=8) :: successive_norm
    INTEGER :: i, N 

    WRITE(*,*) '************** BEGIN CALCULATING DF_DC ***************'
    N = 8

    delta_dc(1) = 0.01*d_c
    DO i=1,N 
        delta_dc(i+1) = 0.01*d_c / ( (i+1)*2.0 )
        dc_1(i) = d_c + delta_dc(i) 
        dc_2(i) = d_c - delta_dc(i)

        CALL partialF_partial_x( d_o, d_o, dc_1(i), dc_2(i), yo, yo, delta_dc(i), y_ofc, df_dc(i), err_tol )
    END DO 

    DO i=2,N
        successive_norm = ABS( df_dc(i) - df_dc(i-1) )
        WRITE(*,10) successive_norm
10      FORMAT(' ','Successive Norm between df_dc Values: ', ES14.6)
    END DO
	
    dfdc_out = df_dc(N) 

    WRITE(*,*) '************** END CALCULATING DF_DC ***************'

END SUBROUTINE partialF_partial_dc


SUBROUTINE partialF_partial_dY(d_o, dfdY_out, d_c, yo, y_ofc, err_tol)
    IMPLICIT NONE
    REAL(KIND=8), INTENT(IN) :: d_o, d_c, yo, y_ofc, err_tol
    REAL(KIND=8), INTENT(OUT) :: dfdY_out
    REAL(KIND=8), DIMENSION(8) :: delta_dY, Y_1, Y_2, df_dY !tau_1, tau_2, 
    REAL(KIND=8) :: successive_norm
    INTEGER :: i, N 

    WRITE(*,*) '************** BEGIN CALCULATING DF_DY ***************'
    N = 8

    delta_dY(1) = 0.01*yo 
    DO i=1,N 
        delta_dY(i+1) = 0.01*yo / ( (i+1)*2.0 )
        Y_1(i) = yo + delta_dY(i) 
        Y_2(i) = yo - delta_dY(i)

        CALL partialF_partial_x( d_o, d_o, d_c, d_c, Y_1(i), Y_2(i), delta_dY(i), y_ofc, df_dY(i), err_tol )
    END DO 

    DO i=2,N
        successive_norm = ABS( df_dY(i) - df_dY(i-1) )
        WRITE(*,10) successive_norm
10      FORMAT(' ','Successive Norms between df_dY Values: ', ES14.6)
    END DO
	
    dfdY_out = df_dY(N) 

    WRITE(*,*) '************** END CALCULATING DF_DY ***************'

END SUBROUTINE partialF_partial_dY


SUBROUTINE partialF_partial_x(do_1, do_2, dc_1, dc_2, yo_1, yo_2, delta_x, y_ofc, df_dx, err_tol)
    IMPLICIT NONE 
    REAL(KIND=8), INTENT(IN) :: do_1, do_2, dc_1, dc_2, yo_1, yo_2, &
        delta_x, y_ofc, err_tol
    REAL(KIND=8), INTENT(OUT) :: df_dx
    REAL(KIND=8) :: tau_1, tau_2, epsilon_1, epsilon_2, &
        LHS_1, LHS_2
    INTEGER :: N, i

    tau_1 = LOG(do_1 / dc_1)  
    tau_2 = LOG(do_2 / dc_2)  
    LHS_1 = y_ofc / yo_1
    LHS_2 = y_ofc / yo_2

    WRITE(*,10) tau_1, LHS_1
10  FORMAT('(upper bound values) tau_1 =',ES14.6, ',  LHS_1 =', ES14.6)   
    WRITE(*,20) tau_2, LHS_2
20  FORMAT('(lower bound values) tau_2 =',ES14.6, ',  LHS_2 =', ES14.6)

    CALL bisection(tau_1, LHS_1, err_tol, epsilon_1)
    CALL bisection(tau_2, LHS_2, err_tol, epsilon_2)

    WRITE(*,30) epsilon_1, epsilon_2, delta_x
30	FORMAT('eps_1:', ES14.6, '  eps_2:', ES14.6, '  delta_dx:', ES14.6)

    !central difference formula to calculate derivative
    df_dx = (epsilon_1 - epsilon_2) / (2.0 * delta_x)		

END SUBROUTINE partialF_partial_x

SUBROUTINE bisection(tau, LHS, err_tol, p)
    IMPLICIT NONE
    REAL(KIND=8), INTENT(IN) :: tau, LHS, err_tol
    REAL(KIND=8), INTENT(OUT) :: p
    REAL(KIND=8) :: a, b, f_p, f_a, f_b 
    INTEGER :: maxit, i 

    a = 1.0E-02
    b = 1.0
    CALL Fx_eval(tau, a, LHS, f_a)
    CALL Fx_eval(tau, b, LHS, f_b) 
    IF ( f_a * f_b > 0 ) THEN 
        WRITE(*,10) a, b 
10      FORMAT('!!!!! ROOT IS NOT BOUNDED IN: [', ES8.2, ',', ES8.2 '] !!!!!!!!!!!!')  
    ELSE 
        maxit = 300
        DO i=1,maxit
            p = a + (b - a)/2.0
            CALL Fx_eval(tau, p, LHS, f_p)  

            IF ( ( ABS(f_p) < err_tol ) .OR. ( (b-a)/2. < err_tol ) ) EXIT 

            IF ( f_a*f_p > 0 ) THEN 
                a = p 
                f_a = f_p 
            ELSE 
                b = p 
            END IF 

            WRITE(*,20) i, ABS(f_p)
20          FORMAT(' ',"iteration count:" I3,  "   f(epsilon):" ES14.6)	
        END DO 
    END IF 
END SUBROUTINE bisection

SUBROUTINE Fx_eval(tau, eps, LHS, Fx)
    IMPLICIT NONE
    REAL(KIND=8), PARAMETER :: PI = ATAN(1.0)*4.0 
    REAL(KIND=8), INTENT(IN) :: tau, eps, LHS
    REAL(KIND=8), INTENT(OUT) :: Fx
    REAL(KIND=8) :: H0minus, H1minus, H2minus, h_0minus, h_1minusIntegral, &
        h_1minus, h_match, Dfx, phi


    IF (tau < 0 .OR. eps < 0) THEN 
        WRITE(*,10) tau 
        WRITE(*,12) eps 
10      FORMAT(' ','!!!!!!!!!!!!! ERROR :: tau =' ES14.6, " < 0 !!!!!!!!!!!!!!!!!!")
12      FORMAT(' ','!!!!!!!!!!!!! ERROR :: tau =' ES14.6, " < 0 !!!!!!!!!!!!!!!!!!")
        WRITE(*,*) '!!!!!!!! CANNOT EVALUATE SQUARE ROOTS OF VALUES <0 !!!!!!!!!!!!!!'	
    END IF 

    phi = tau/eps
    !define parts of asymptotic equation that do not depend on epsilon
    H0minus = ( EXP(3.*tau) - 1. ) / 3.
    H1minus = ( EXP(3.*tau) + 2. ) / 3.
    H2minus = (22. / 9.) +  ( 14.-24.*tau )*EXP(3*tau) / 9.

    !define parts of asymptotic equation that depend on epsilon
    h_0minus = (phi/2.) * ( 1. + ERF( SQRT(phi)/2. ) ) &
        + ERF( SQRT(phi)/2. ) &
        + SQRT( phi/PI ) * EXP( -phi / 4. )
    !NOTE: h_1minusIntegral is the alernate form of the integral in
    !      the variable h_1minus, as given by Mathematica
    h_1minusIntegral = ( EXP(-phi/4.0) * SQRT(phi) * (3.0*phi - 4.0) )/ ( 2.0*SQRT(PI) ) &
    	- ( 1.0/4.0 ) * ( 8.0 + phi*(2.0 + 3.0*phi) ) * ERFC( SQRT(phi)/2.0 )

    h_1minus = 4. + phi + ( 3.*(phi**2.) / 2. ) &
        - 2. * EXP( -phi / 4.  ) 

    h_match = 1.0 + phi + eps * ( 4. + phi + (3./2.)*(phi**2) )

    Fx  = 1.0 + h_0minus + eps*h_1minus + eps*h_1minusIntegral &
    	+ (H0minus/eps) + H1minus + eps*H2minus - h_match - LHS

END SUBROUTINE Fx_eval

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



