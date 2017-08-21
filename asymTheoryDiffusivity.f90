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
        d_o, tau, LHS, D, err_tol, ep_min, ep_max, dfdo, epsilon_bisect
    REAL(KIND=8), DIMENSION(6) :: fc_props
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
    CALL bisection(tau, LHS, err_tol, epsilon_bisect)
    WRITE(*,*) " ==================== Bisection Iteration Results ====================="
    WRITE(*,54) epsilon_bisect
54  FORMAT(" || The value of epsilon is....................." ES14.6, "        ||")
    D = epsilon_bisect * K/ 8.0
    WRITE(*,56) D
56  FORMAT(" || The effective liquid diffusivity is........." ES14.6, " mm^2/s ||")
    WRITE(*,58) D*(1.0/1000.)**2.
58  FORMAT(" || The effective liquid diffusivity is........." ES14.6, " m^2/s  ||")
    WRITE(*,*) " ====================================================================="

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

END PROGRAM asymTheoryDiffusivity

SUBROUTINE partialF_partial_do(d_o, dfdo_out, d_c, yo, y_ofc, err_tol)
    IMPLICIT NONE
    REAL(KIND = 8), INTENT(IN) :: d_o, d_c, yo, y_ofc, err_tol
    REAL(KIND = 8), INTENT(OUT) :: dfdo_out
    REAL(KIND = 8), DIMENSION(8) :: delta_do, do_1, do_2, df_do !tau_1, tau_2, 
    REAL(KIND = 8) :: p_incr, successive_norm
    INTEGER :: i, N 

    WRITE(*,*) '************** BEGIN CALCULATING DF_DO ***************'

    p_incr = 0.01
    N = 8
    ! TODO: add restrction on do_1 and do_2 s.t. do/d_c > 0 iff tau > 0
    ! do_max = d_c 
    ! do_min = d_o / 4.0
    ! delta_do = (do_max - do_min)

    delta_do(1) = 0.01*d_o 
    DO i=1,N 
        delta_do(i+1) = p_incr / ( (i+1)*2.0 )
        do_1(i) = d_o + delta_do(i) 
        do_2(i) = d_o - delta_do(i)

        CALL partialF_partial_x( do_1(i), do_2(i), d_c, yo, yo, delta_do(i), y_ofc, df_do(i), err_tol )
    END DO 

    DO i=2,N
        successive_norm = ABS( df_do(i) - df_do(i-1) )
        WRITE(*,10) successive_norm
10      FORMAT(' ','Successive Norm: ', ES14.6)
    END DO
	
    dfdo_out = df_do(N-1) 
    WRITE(*,20) dfdo_out
20  FORMAT(' ','df_do approximately: ', ES14.6)

    WRITE(*,*) '************** END CALCULATING DF_DO ***************'

END SUBROUTINE partialF_partial_do


SUBROUTINE partialF_partial_x(do_1, do_2, d_c, yo_1, yo_2, delta_x, y_ofc, df_dx, err_tol)
    IMPLICIT NONE 
    REAL(KIND = 8), INTENT(IN) :: do_1, do_2, d_c, yo_1, yo_2, &
        delta_x, y_ofc, err_tol
    REAL(KIND = 8), INTENT(OUT) :: df_dx
    REAL(KIND = 8) :: tau_1, tau_2, epsilon_1, epsilon_2, &
        LHS_1, LHS_2
    INTEGER :: N, i

    tau_1 = LOG(do_1 / d_c)  
    tau_2 = LOG(do_2 / d_c)  
    LHS_1 = y_ofc / yo_1
    LHS_2 = y_ofc / yo_2


    WRITE(*,10) tau_1, LHS_1
10  FORMAT('(upper bound values) tau_1 =',ES14.6, ',  LHS_1 =', ES14.6)   
    CALL bisection(tau_1, LHS_1, err_tol, epsilon_1)

    WRITE(*,20) tau_2, LHS_2
20  FORMAT('(lower bound values) tau_2 =',ES14.6, ',  LHS_2 =', ES14.6)
    CALL bisection(tau_2, LHS_2, err_tol, epsilon_2)

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
    REAL(KIND=8), PARAMETER :: PI = 3.1415927
    REAL(KIND=8), INTENT(IN) :: tau, eps, LHS
    REAL(KIND=8), INTENT(OUT) :: Fx

    IF (tau < 0 .OR. eps < 0) THEN 
        WRITE(*,10) tau 
        WRITE(*,12) eps 
10      FORMAT(' ','!!!!!!!!!!!!! ERROR :: tau =' ES14.6, " < 0 !!!!!!!!!!!!!!!!!!")
12      FORMAT(' ','!!!!!!!!!!!!! ERROR :: tau =' ES14.6, " < 0 !!!!!!!!!!!!!!!!!!")
        WRITE(*,*) '!!!!!!!! CANNOT EVALUATE SQUARE ROOTS OF VALUES <0 !!!!!!!!!!!!!!'	
    END IF 

    ! A possible form of the asymptotic theory equation
    ! as output by Mathematica (See liquid_diff.Rev03)
    Fx = (1.0/36.0)*( 2.0 * ( -6.0 + 12.0*eps + 44.0*eps**2 - 9.0*tau )/eps &
        - 18.0*EXP(-tau/(4.0*eps)) * ( 4.0*SQRT(PI)*eps + (-2.0 + 4.0*eps - 3.0*tau)*SQRT(tau/eps) ) / SQRT(PI) & 
        - 4.0*EXP(3.0*tau) * ( -3.0 - 3.0*eps + 2.0*eps**2 * (-7.0 + 12.0*tau) ) / eps  &
        + 18.0*( 2*eps + tau )*ERF( SQRT(tau/eps) / 2.0 ) / eps &
        - 9.0 * ( 8.0*eps**2 + 2.0*eps*tau + 3.0*tau**2 )*ERFC( SQRT(tau/eps) / 2.0 ) / eps - 36.0*LHS )

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



