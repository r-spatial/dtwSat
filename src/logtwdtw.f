C   Compute TWDTW distance using logistic weight
C
C     XM - matrix with the time series (N,D)
C     YM - matrix with the temporal profile (M,D)
C     N  - Number of rows in CM, DM, and VM - time series
C     M  - Number of columns CM, DM, and VM - temporal profile
C     D  - Number of spectral dimensions including time in XM and YM
C     I  - Single point in the time series to calculate the local distance
C     J  - Single point in the temporal profile to calculate the local distance
C     A  - Time-Weight parameter alpha
C     B  - Time-Weight parameter beta
      REAL FUNCTION distance(YM, XM, N, M, D, I, J, TW, TD)
      INTEGER N, M, D, I, J, K 
      DOUBLE PRECISION XM(M,D), YM(N,D), TD, BD, CD, TW(2)
      PARAMETER(ZERO=0)
      REAL NAN
      NAN  = ZERO 
      NAN  = NAN / NAN
      distance = NAN
C      TD = YM(I,1) - XM(J,1)
C      CALL ellapsed(TD)
      CD = ZERO
      DO 30 K = 2, D
         BD = YM(I,K) - XM(J,K)
         CD = CD + (BD * BD)
   30 CONTINUE
C      WRITE (*,*) 'The value of X J', J, 'is', XM(J,D)
      distance = SQRT(CD) + 1.0 / (1.0 + EXP(TW(1) * (TD - TW(2))))
      RETURN
      END

