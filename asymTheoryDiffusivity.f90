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
        dfdY, Uk_sq, Udo_sq, Udc_sq, UYo_sq, Ueps, Udiff, Udiff_95
    INTEGER :: sol_id
    REAL(KIND=8), DIMENSION(11) :: fc_props
    REAL(KIND=8), DIMENSION(700) :: fvalues, ep_vals
    CHARACTER(len=11) :: filename
    CHARACTER(len=100) :: junk
    INTEGER :: i, status, N
    REAL(KIND=8), DIMENSION(1000000) :: dc_mc, do_mc, yo_mc, tau_mc, LHS_mc, & 
                                        eps_mc, K_mc, D_mc

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

        readData: DO i=1,11 					!read lines 5-9
            READ(1,*,IOSTAT=status) fc_props(i)
        END DO readData

        p 			= fc_props(1)    ! chamber pressure [atm]
        dc_sq 		= fc_props(2)    ! drop diameter @ onset of fc squared [mm^2]
        K 			= fc_props(3)    ! burning rate constant prior to onset of fc [mm^2/s]
        do_measured = fc_props(4)    ! initial drop diameter measured [mm]
        yo 			= fc_props(5)    ! initial mass frac of low volitiliy component
        err_tol 	= fc_props(6)    ! error tolerrance for bisection method
        Uk_sq       = fc_props(7)**2 ! uncertainty in K squared (U_k ^2) [mm^2/s]^2
        Udo_sq      = fc_props(8)**2 ! uncertainty in do squared (U_do^2) [mm^2]
        Udc_sq      = fc_props(9)**2 ! uncertainty in dc squared (U_dc^2) [mm^2]
        UYo_sq      = (fc_props(10)*yo)**2 !uncertainty in Yo squared
        sol_id      = fc_props(11)   ! solvent id (1)-Heptane (2)-Propanol
        CLOSE(UNIT=1)
    ELSE
        WRITE(*,*) 'Open of file NOT sucessful!'
    ENDIF openif

    d_c    = SQRT(dc_sq)   ! drop diameter at onset of flame contraction [mm]
    y_ofc = 1.0 		   ! mass fraction of low volatility comp at onset of fc

    ! correct d_o for heptane-hexadecane droplets
    IF (sol_id == 1) THEN
        IF (p == 1) THEN
            percent_increase = 0.055  !so far these numbers are only valid for hep-hex exp's
        ELSE
            percent_increase = 0.086
        ENDIF
    ENDIF

    ! correct d_o for propanol-glycerol droplets
    IF (sol_id == 2) THEN
        IF (p == 1) THEN
            percent_increase = 0.03  
        ELSE
            percent_increase = 0.05
        ENDIF
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

    CALL uncertainty_diffusivity(dfdo, Udo_sq, dfdc, Udc_sq, dfdY, UYo_sq, D, Uk_sq, &
                                K, epsilon_bisect, Ueps, Udiff, Udiff_95)

    WRITE(*,*) " ------------------------ Equation Sensitivies ------------------------"
    WRITE(*,60) dfdo
60  FORMAT(" || Sensitivity in d_o (df_do) is:.............." ES14.6, "        ||")
    WRITE(*,70) dfdc
70  FORMAT(" || Sensitivity in d_c (df_dc) is:.............." ES14.6, "        ||")
    WRITE(*,80) dfdY
80  FORMAT(" || Sensitivity in y_o (df_yo) is:.............." ES14.6, "        ||")
    WRITE(*,90) Ueps
90  FORMAT(" || Uncertainty in epsilon (U_eps) is:.........." ES14.6, "        ||")
!     WRITE(*,100) Udiff
! 100 FORMAT(" || Uncertainty in D (U_diff) is:..............." ES14.6, " mm^2/s ||")
    WRITE(*,110) Udiff_95
110 FORMAT(" || Uncertainty in D (95% conf interval):......." ES14.6, " mm^2/s ||")
    WRITE(*,112) Udiff_95*(1.0/1000.)**2.
112 FORMAT(" || Uncertainty in D (95% conf interval):......." ES14.6, "  m^2/s ||")
    WRITE(*,120) (Udiff_95 / D) * 100
120 FORMAT(" || Relative error in D (95% conf interval):...." ES14.6, " %      ||")
    WRITE(*,*) " ----------------------------------------------------------------------"      

! SUBROUTINE mc_uncertainty( dc_mu, dc_sigma, do_mu, do_sigma, yo_mu, & 
!                         yo_sigma, y_ofc, err_tol, &
!                         dc_mc_out, do_mc_out, yo_mc_out, tau_mc_out, LHS_mc_out, eps_mc, & 
!                         NR, K_mu, dK_sigma, K_mc_out )

    CALL mc_uncertainty(d_c, SQRT(Udc_sq), d_o, SQRT(Udo_sq), yo, & 
                        SQRT(UYo_sq), y_ofc, err_tol, &
                        dc_mc, do_mc, yo_mc, tau_mc, LHS_mc, eps_mc, 1000000, & 
                        K, SQRT(Uk_sq), K_mc, D_mc )

    ! OPEN(UNIT=20,FILE='eps_mc_valuest')
    ! WRITE(20,*) "dc_mc, do_mc,  yo_mc, tau_mc, LHS_mc, eps_mc, K_mc, D_mc"
    ! DO i=1,1000000 
    !     WRITE(20,*) dc_mc(i), do_mc(i), yo_mc(i), tau_mc(i), LHS_mc(i), eps_mc(i), K_mc(i), D_mc(i)
    ! END DO
    ! CLOSE(UNIT = 20)    

    WRITE(*,*) ".......PROGRAM COMPLETED RUNNING......."
END PROGRAM asymTheoryDiffusivity

