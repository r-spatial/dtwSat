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
C     XM - matrix with the time series (N,D)
C     YM - matrix with the temporal profile (M,D)
C     N  - Number of rows in CM, DM, and VM - time series
C     M  - Number of columns CM, DM, and VM - temporal profile
C     D  - Number of spectral dimensions including time in XM and YM
C     I  - Single point in the time series to calculate the local distance
C     J  - Single point in the temporal profile to calculate the local distance
      REAL FUNCTION distance(YM, XM, N, M, D, I, J)
      INTEGER N, M, D, I, J, K 
      DOUBLE PRECISION XM(M,D), YM(N,D), PC, TD, BD, CD, HPC
      PARAMETER(ZERO=0,ONE=1)
      REAL NAN, INF
      NAN  = ZERO 
      NAN  = NAN / NAN
      PC = 366
      HPC = PC/2
      distance = NAN
C     Compute ellapsed time difference 
      TD = SQRT((YM(I,1) - XM(J,1)) * (YM(I,1) - XM(J,1)))
      IF (TD.GT.HPC) THEN
         TD = ABS(PC - TD)
      ENDIF
      CD = ZERO
      DO 30 K = 2, D
         BD = YM(I,K) - XM(J,K)
         CD = CD + (BD * BD)
   30 CONTINUE
C      WRITE (*,*) 'The value of X J', J, 'is', XM(J,D)
      distance = SQRT(CD) + 1.0 / (1.0 + EXP(-0.1 * (TD - 50.0)))
      RETURN
      END

