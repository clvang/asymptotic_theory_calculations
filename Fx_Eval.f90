SUBROUTINE Fx_eval(tau, eps, LHS, Fx)
!   This subroutine defines the asympotic theory equation
!   from Aharon and Shaw and is used in the program
!   for functional evaluations.

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