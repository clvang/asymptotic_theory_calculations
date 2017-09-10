SUBROUTINE bisection(tau, LHS, err_tol, p)
! subroutine to solve asymptotic equation using
! bisection method.  I've tried using Newton's Method
! with damping - this doesn't work so well since convergence
! was very sensitive to initial guess.

    IMPLICIT NONE
    REAL(KIND=8), INTENT(IN) :: tau, LHS, err_tol
    REAL(KIND=8), INTENT(OUT) :: p
    REAL(KIND=8) :: a, b, f_p, f_a, f_b 
    INTEGER :: maxit, i 

    a = 1.0E-10
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
        END DO 
            ! WRITE(*,20) i, ABS(f_p)
! 20          FORMAT(' ',"iteration count:" I3,  "   f(epsilon):" ES14.6)         
    END IF 

END SUBROUTINE bisection