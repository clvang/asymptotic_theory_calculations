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