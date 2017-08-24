SUBROUTINE partialF_partial_do(d_o, dfdo_out, d_c, yo, y_ofc, err_tol)
    IMPLICIT NONE
    REAL(KIND=8), INTENT(IN) :: d_o, d_c, yo, y_ofc, err_tol
    REAL(KIND=8), INTENT(OUT) :: dfdo_out
    REAL(KIND=8), DIMENSION(8) :: delta_do, do_1, do_2, df_do !tau_1, tau_2, 
    REAL(KIND=8) :: successive_norm
    INTEGER :: i, N 

    WRITE(*,*) '************** BEGIN CALCULATING DF_DO ***************'

    N = 8
    ! TODO: add restrction on do_1 and do_2 s.t. do/d_c > 0 iff tau > 0
    ! do_max = d_c 
    ! do_min = d_o / 4.0

    delta_do(1) = 0.01*d_o
    DO i=1,N 
        delta_do(i+1) = 0.01*d_o / ( (i+1)*2.0 )
        do_1(i) = d_o + delta_do(i) 
        do_2(i) = d_o - delta_do(i)

        CALL partialF_partial_x( do_1(i), do_2(i), d_c, d_c, yo, yo, delta_do(i), y_ofc, df_do(i), err_tol )
    END DO 

    DO i=2,N
        successive_norm = ABS( df_do(i) - df_do(i-1) )
        WRITE(*,10) successive_norm
10      FORMAT(' ','Successive Norm between df_do Values: ', ES14.6)
    END DO
    dfdo_out = df_do(N) 
    WRITE(*,*) '************** END CALCULATING DF_DO ***************'

END SUBROUTINE partialF_partial_do