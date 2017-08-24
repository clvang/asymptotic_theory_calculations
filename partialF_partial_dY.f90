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