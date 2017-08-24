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