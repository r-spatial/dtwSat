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
C   This function was adpted from the C function 'computeCM'  C
C   implemented in the R package 'dtw' by Toni Giorgino.      C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C     XM - matrix with the time series (N,D)
C     YM - matrix with the temporal profile (M,D)
C     CM - Output cumulative cost matrix
C     DM - Direction matrix 
C     VM - Starting points matrix 
C     SM - Matrix of step patterns 
C     N  - Number of rows in CM, DM, and VM - time series
C     M  - Number of columns CM, DM, and VM - temporal profile
C     D  - Number of spectral dimensions including time in XM and YM
C     NS - Number of rows in SM 
      SUBROUTINE computecostfast(XM, YM, CM, DM, VM, SM, N, M, D, NS)
C     I/O Variables       
      INTEGER N, M, D, NS, SM(NS,4), DM(N+1,M), VM(N+1,M)
      DOUBLE PRECISION XM(M,D), YM(N,D), CM(N+1,M)
C     Internals
      DOUBLE PRECISION W, CP(NS), VMIN
      INTEGER I, J, IL(NS), JL(NS), K, PK, KMIN, ZERO, ONE
      PARAMETER(ZERO=0,ONE=1)
      REAL NAN, INF
      NAN  = ZERO 
      NAN  = NAN / NAN
      INF  = HUGE(ZERO)
      VM(1,1) = 1
C     Initialize the firt row and col of the matrices 
      DO 21 I = 2, N+1 
         CM(I,1) = CM(I-1,1) + distance(YM, XM, N, M, D, I-1, 1)
C         WRITE (*,*) 'The distance ',I,',1 is ', CM(I,1)
         DM(I,1) = 3
         VM(I,1) = 1
   21 CONTINUE
      DO 31 J = 2, M 
         CM(2,J) = CM(2,J-1) + distance(YM, XM, N, M, D, 1, J)
C         WRITE (*,*) 'The distance 2,',J,' is ', CM(2,J)
         DM(1,J) = 2
         VM(1,J) = J
   31 CONTINUE
C     Compute cumulative cost matrix
      DO 32 J = 2, M
         DO 22 I = 2, N+1
C           Calculate local distance 
            CM(I,J) = distance(YM, XM, N, M, D, I-1, J)
C           Initialize list of step cost 
            DO 10 K = 1, NS
               CP(K) = NAN
   10       CONTINUE  
            DO 11 K = 1, NS
               PK = SM(K,1)
               IL(K) = I - SM(K,2) 
               JL(K) = J - SM(K,3)  
               IF ((IL(K).GT.ZERO).AND.(JL(K).GT.ZERO)) THEN
                  W = SM(K,4)
                  IF (W.EQ.-ONE) THEN
                    CP(PK) = CM(IL(K),JL(K))
                  ELSE
                    CP(PK) = CP(PK) + CM(IL(K),JL(K))*W
                  ENDIF
               ENDIF
   11       CONTINUE
            KMIN = -ONE
            VMIN =  INF
            ILMIN = -ONE
            JLMIN = -ONE
            DO 12 K = 1, NS
               PK = SM(K,1)
               IF (CP(PK).EQ.CP(PK).AND.CP(PK).LT.VMIN) THEN
                  KMIN = PK
                  VMIN = CP(PK)
                  ILMIN = IL(K)
                  JLMIN = JL(K)
               ENDIF
   12       CONTINUE
            IF (KMIN.GT.-ONE) THEN 
               CM(I,J) = VMIN
               DM(I,J) = KMIN
               VM(I,J) = VM(ILMIN, JLMIN)
            ENDIF
   22    CONTINUE
   32 CONTINUE
   99 CONTINUE
      END
