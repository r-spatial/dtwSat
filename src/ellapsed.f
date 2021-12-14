C   Computation ellapsed time in days
C
C     TD - time difference in days
      SUBROUTINE ellapsed(TD)
      DOUBLE PRECISION TD, HPC
      PARAMETER(PC=366.0)
      HPC = PC/2
C     Compute ellapsed time difference 
      TD = SQRT(TD * TD)
C     Correct ellapsed time with year cycle
      IF (TD.GT.HPC) THEN
         TD = PC - TD
      ENDIF
      TD = ABS(TD)
      RETURN
      END

