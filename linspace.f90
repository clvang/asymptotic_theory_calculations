
SUBROUTINE linspace(varArray, lowLimit, upLimit, numPts) 
    IMPLICIT NONE

    REAL(KIND=8), INTENT(OUT), DIMENSION(numPts) :: varArray
    REAL(KIND=8), INTENT(IN) :: lowLimit, upLimit
    INTEGER, INTENT(IN) :: numPts
    INTEGER :: i
    REAL(KIND=8) :: intervalSize

    intervalSize = (upLimit - lowLimit) / numPts
    varArray(1) = lowLimit
    DO i = 2, numPts-1
        varArray(i) = varArray(i-1) + intervalSize
    END DO
    varArray(numPts) = upLimit
END SUBROUTINE linspace