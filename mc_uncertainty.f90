
SUBROUTINE mc_uncertainty( tau_vector, LHS_vector, K_vector, err_tol, N, D_vector )

  IMPLICIT NONE

  INTEGER (KIND=4), INTENT(IN) :: N
  REAL(KIND = 8), DIMENSION(N), INTENT(IN) :: tau_vector, LHS_vector, K_vector
  REAL (KIND=8), INTENT(IN) :: err_tol  
  REAL (KIND=8), DIMENSION(N), INTENT(OUT) :: D_vector

  INTEGER ( KIND = 4 ) :: i
  REAL (KIND=8), DIMENSION(N) :: eps_mc

  DO i = 1, N
    CALL bisection( tau_vector(i), &
            LHS_vector(i), &
            err_tol, &
            eps_mc(i) )    
  END DO


  D_vector = K_vector * eps_mc / 8.0


END SUBROUTINE mc_uncertainty 


