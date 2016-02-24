CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                             C
C   (c) Victor Maus <vwmaus1@gmail.com>                       C
C       Institute for Geoinformatics (IFGI)                   C
C       University of Muenster (WWU), Germany                 C
C                                                             C
C       Earth System Science Center (CCST)                    C
C       National Institute for Space Research (INPE), Brazil  C
C                                                             C
C                                                             C
C                   Computation allapsed time - 2016-01-22    C
C                                                             C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C     TM - Time difference matrix 
C     N  - Number of rows in CM
C     M  - Number of columns CM
C     PC - Cycle length in the unity of the measurements 
      SUBROUTINE g(TM, N, M, PC)
C     I/O Variables       
      INTEGER N, M
      DOUBLE PRECISION TM(N,M), PC
C     Internals
      DOUBLE PRECISION HPC
      INTEGER I, J
      PARAMETER(TWO=2.0)
      HPC = PC/2
C     Compute ellapsed time matrix 
      DO 30 J = 1, M 
         DO 20 I = 1, N
            IF (TM(I,J).GT.HPC) THEN
               TM(I,J) = ABS(PC - TM(I,J))
            ENDIF
   20    CONTINUE
   30 CONTINUE
      END
      