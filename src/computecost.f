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
C   Efficient computation of DTW cost matrix  - 2015-10-16    C
C                                                             C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C     CM - Input local cost and output cumulative cost matrix
C     DM - Direction matrix 
C     SM - Matrix of step patterns 
C     N  - Number of rows in CM and DM 
C     M  - Number of columns CM and DM 
C     NS - Number of rows in SM 
      SUBROUTINE computecost(CM, DM, SM, N, M, NS)
C  800 FORMAT(I4,I4,I4,I4,F8.0,F8.0)
C     I/O Variables       
      INTEGER N, M, NS, SM(NS,4), DM(N,M)
      DOUBLE PRECISION CM(N,M)
C     Internals
      DOUBLE PRECISION W, LM(N,M), CP(NS), VMIN
      INTEGER I, J, IL, JL, K, PK, KMIN, ZERO, ONE
      PARAMETER(ZERO=0,ONE=1)
      REAL NAN, INF
      NAN  = ZERO 
      NAN  = NAN / NAN
      INF  = HUGE(ZERO)
C     Initialize cost and direction matrices 
      DO 30 J = 1, M 
         DO 20 I = 1, N 
            LM(I,J) = CM(I,J)
            IF (I.NE.ONE.AND.J.NE.ONE) THEN
              CM(I,J) = NAN
            ENDIF
   20    CONTINUE
   30 CONTINUE
C     Initialize the firt row/col in cost and direction matrices 
      DO 21 I = 2, N 
         CM(I,1) = CM(I-1,1) + CM(I,1)
         DM(I,1) = 3
   21 CONTINUE
      DO 31 J = 2, M 
         CM(1,J) = CM(1,J-1) + CM(1,J)
         DM(1,J) = 2
   31 CONTINUE
C     Compute cumulative cost matrix
      DO 32 J = 1, M
         DO 22 I = 1, N
            IF ( .NOT.ISNAN(CM(I,J)) ) THEN
                GOTO 22 
            ENDIF
C           Initialize list of step cost 
            DO 10 K = 1, NS
               CP(K) = NAN
   10       CONTINUE  
            DO 11 K = 1, NS
               PK = SM(K,1)
               IL = I - SM(K,2) 
               JL = J - SM(K,3)  
               IF ((IL.GT.ZERO).AND.(JL.GT.ZERO)) THEN
                  W = SM(K,4)
                  IF (W.EQ.-ONE) THEN
                    CP(PK) = CM(IL,JL)
                  ELSE
                    CP(PK) = CP(PK) + LM(IL,JL)*W
                  ENDIF
               ENDIF
   11       CONTINUE
            KMIN = -ONE
            VMIN =  INF
            DO 12 K = 1, NS
               PK = SM(K,1)
               IF (.NOT.ISNAN(CP(PK)).AND.CP(PK).LT.VMIN) THEN
                  KMIN = PK
                  VMIN = CP(PK)
               ENDIF
   12       CONTINUE
            IF (KMIN.GT.-ONE) THEN 
               CM(I,J) = VMIN
               DM(I,J) = KMIN
            ENDIF
   22    CONTINUE
   32 CONTINUE
      END