SUBROUTINE uncertainty_diffusivity(dfdo, Udo_sq, dfdc, Udc_sq, dfdY, UYo_sq, D, Uk_sq, K, epsilon, Ueps, Udiff)
    IMPLICIT NONE 
    REAL(KIND=8), INTENT(IN) :: dfdo, Udo_sq, dfdc, Udc_sq, dfdY, UYo_sq, D, Uk_sq, K, epsilon 
    REAL(KIND=8), INTENT(OUT) :: Ueps, Udiff 

    Ueps = SQRT( (dfdo**2)*(Udo_sq) + (dfdc**2)*(Udc_sq) + (dfdY**2)*(UYo_sq) )

    Udiff = SQRT( ( Uk_sq/(K*K) ) + ( Ueps*Ueps/(epsilon*epsilon) ) ) * D
END SUBROUTINE uncertainty_diffusivity 