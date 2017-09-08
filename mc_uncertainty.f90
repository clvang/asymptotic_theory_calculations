
SUBROUTINE mc_uncertainty( dc_mu, dc_sigma, do_mu, do_sigma, yo_mu, & 
						yo_sigma, y_ofc, err_tol, &
						dc_mc_out, do_mc_out, yo_mc_out, tau_mc_out, LHS_mc_out, eps_mc, & 
						NR, K_mu, dK_sigma, K_mc_out, D_mc )
  IMPLICIT NONE

  REAL(KIND = 8), INTENT(IN) :: dc_mu, dc_sigma, do_mu, do_sigma, yo_mu, & 
  								yo_sigma, y_ofc, err_tol, K_mu, dK_sigma
  INTEGER (KIND=4), INTENT(IN) :: NR
  INTEGER ( KIND = 4 ) :: i
  REAL ( KIND = 4 ) ::  dc_mu_k4, dc_sigma_k4, do_mu_k4, do_sigma_k4, &
  						yo_mu_k4, yo_sigma_k4, dc_mc_temp, yofc_k4, K_mu_k4, dK_sigma_k4
  REAL ( KIND = 4 ) :: r4_normal_ab
  INTEGER ( KIND = 4 ) :: seed
  REAL (KIND=4), DIMENSION(NR) :: dc_mc, do_mc, &
  													yo_mc, tau_mc, LHS_mc, K_mc
  REAL (KIND=8), DIMENSION(NR), INTENT(OUT) :: eps_mc, dc_mc_out, do_mc_out, &
                                                yo_mc_out, tau_mc_out, LHS_mc_out, &
                                                K_mc_out, D_mc

  WRITE ( *, '(a)' ) ' '
  WRITE ( *, '(a)' ) '  Generating random values from normal distribution '
  WRITE ( *, '(a)' ) '  with mean MU and standard deviation SIGMA.'
  WRITE ( *, '(a)' ) ' '  

  !convert intent(in) variables from KIND8 to KIND4
  dc_mu_k4 = REAL(dc_mu,KIND(dc_mu_k4))   			
  dc_sigma_k4 = REAL(dc_sigma, KIND(dc_sigma_k4) )   

  do_mu_k4 = REAL(do_mu,KIND(do_mu_k4))
  do_sigma_k4 = REAL(do_sigma, KIND(do_sigma_k4) )

  yo_mu_k4 = REAL(yo_mu, KIND(yo_mu_k4))
  yo_sigma_k4 = REAL(yo_sigma, KIND(yo_sigma_k4))

  K_mu_k4 = REAL(K_mu, KIND(K_mu_k4))
  dK_sigma_k4 = REAL(dK_sigma, KIND(dK_sigma_k4))

  yofc_k4 = REAL( y_ofc, KIND(yofc_k4) )
  seed = 5

  DO i = 1, NR
    do_mc(i) = r4_normal_ab ( do_mu_k4, do_sigma_k4, seed )
    yo_mc(i) = r4_normal_ab ( yo_mu_k4, yo_sigma_k4, seed )
    K_mc(i) = r4_normal_ab ( K_mu_k4, dK_sigma_k4, seed )
    ! WRITE ( *, '(2x,i8,2x,f14.8,2x, f14.8)' ) i, do_mc(i), yo_mc(i)
  END DO
  LHS_mc = yofc_k4 / yo_mc

  DO i = 1,NR
  	dc_mc_temp = r4_normal_ab ( dc_mu_k4, dc_sigma_k4, seed )

  	DO WHILE ( dc_mc_temp > do_mc(i) )
  		dc_mc_temp = r4_normal_ab ( dc_mu_k4, dc_sigma_k4, seed )
  		WRITE(*,*) dc_mc_temp
  	END DO
  	dc_mc(i) = dc_mc_temp
  	tau_mc(i) = LOG( do_mc(i) / dc_mc(i) )
    ! WRITE ( *, '(2x,i8,2x,f14.8)' ) i, tau_mc(i)
  END DO

  !convert KIND4 vectors to KIND8 for output
  dc_mc_out = REAL(dc_mc,KIND(dc_mu))
  do_mc_out = REAL(do_mc,KIND(do_mu))
  yo_mc_out = REAL(yo_mc,KIND(yo_mu))
  tau_mc_out = REAL(tau_mc,KIND(do_mu))
  LHS_mc_out = REAL(LHS_mc,KIND(y_ofc))
  K_mc_out = REAL(K_mc,KIND(K_mu))

  DO i = 1, NR
    CALL bisection( tau_mc_out(i), &
            LHS_mc_out(i), &
            err_tol, &
            eps_mc(i) )    
    ! WRITE(*,'(2x,i8,2x,ES14.6)') i, eps_mc(i)
  END DO


  D_mc = K_mc_out * eps_mc / 8.0
  ! DO i = 1,NR 
  !   D_mc(i) = K_mc_out(i)*eps_mc(i) / 8.0
  ! END DO 

END SUBROUTINE mc_uncertainty 


