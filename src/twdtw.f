C   Computation of TWDTW cost matrix
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
C     TW - Time-Weight parameters alpha and beta 
C     LB - Constrain TWDTW calculation to band given by TW(2)
      SUBROUTINE twdtw(XM, YM, CM, DM, VM, SM, N, M, D, NS, TW, LB, JB)
C     I/O Variables
      INTEGER N, M, D, NS, SM(NS,4), DM(N+1,M), VM(N+1,M), JB(N)
      DOUBLE PRECISION XM(M,D), YM(N,D), CM(N+1,M), TW(2)
      LOGICAL LB
C     Internals
      DOUBLE PRECISION W, CP(NS), VMIN, A, B, TD
      INTEGER I, J, IL(NS), JL(NS), K, PK, KMIN, ZERO, ONE, JM
      PARAMETER(ZERO=0,ONE=1)
      DOUBLE PRECISION NAN, INF
      NAN  = 0.0
      NAN  = NAN / NAN
      INF  = HUGE(0.0)
      IML  = 1
      VM(1,1) = 1

C     Initialize the firt row and col of the matrices 
      DO 21 I = 2, N+1
         TD = YM(I-1,1) - XM(1,1)
         CALL ellapsed(TD)
C         IF (TD.GT.TW(2)) THEN
C           CM(I,1) = INF
C         ELSE
           CM(I,1) = CM(I-1,1) + 
     &               distance(YM, XM, N, M, D, I-1, 1, TW, TD)
C         ENDIF
C         WRITE (*,*) 'The distance ',I,',1 is ', CM(I,1)
         DM(I,1) = 3
         VM(I,1) = 1
   21 CONTINUE
      DO 31 J = 2, M 
         TD = YM(2,1) - XM(J,1)
         CALL ellapsed(TD)
C         IF (TD.GT.TW(2)) THEN
C           CM(2,J) = INF
C         ELSE
           CM(2,J) = CM(2,J-1) + 
     &               distance(YM, XM, N, M, D, 1, J, TW, TD)
C         ENDIF
C         WRITE (*,*) 'The distance 2,',J,' is ', CM(2,J)
         DM(1,J) = 2
         VM(1,J) = J
   31 CONTINUE
C     Compute cumulative cost matrix
      J = 2
      DO 32 WHILE ( J .LE. M )
         I = 2
         DO 22 WHILE ( I .LE. N+1 ) 
C         PRINT *, "J: ", J, "I: ", I
C           Calculate local distance 
C           # the call takes I-1 because local matrix has an additional row at the begning
            TD = YM(I-1,1) - XM(J,1)
            CALL ellapsed(TD)
            IF (LB.AND.(TD.GT.TW(2))) THEN
C              print *, "I: ", I, "TD: ", TD, " -- TW: ", TW(2)
              CM(I,J) = INF
              DM(I,J) = -ONE
              VM(I,J) = ZERO
              GOTO 44
            ELSE
              CM(I,J) = distance(YM, XM, N, M, D, I-1, J, TW, TD) 
            ENDIF
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
   44       CONTINUE
            I = I + 1
   22    CONTINUE
         J = J + 1
   32 CONTINUE
   99 CONTINUE
      J = 1
      K = ZERO
C      PRINT *, "DONE: LOOP 1"
      DO 69 WHILE ( J .LE. M )
         IF (VM(N+1,J).NE.ZERO) THEN
            IF (K.EQ.ZERO) THEN
               K = 1
               JB(K) = J
               JM = VM(N+1,J)
               GOTO 68
            ENDIF
            IF (VM(N+1,J).NE.JM) THEN
               K = K + 1
               JB(K) = J
               JM = VM(N+1,J)
               GOTO 68
            ENDIF
C            PRINT *, J, "JB:",JB(k),"-", CM(N+1,J),"-",CM(N+1,JB(K))
            IF (CM(N+1,J).LT.CM(N+1,JB(K))) THEN
               JB(K) = J
               GOTO 68
            ENDIF
         ENDIF
   68    CONTINUE
         J = J + 1
   69 CONTINUE
C      PRINT *, "XM", XM
C      PRINT *, "YM", YM
C      PRINT *, "CM", CM
C      PRINT *, "DM", DM
C      PRINT *, "VM", VM
C      PRINT *, "SM", SM
C      PRINT *, "N", N
C      PRINT *, "M", M
C      PRINT *, "D", D
C      PRINT *, "NS", NS
C      PRINT *, "TW", TW
C      PRINT *, "LB", LB
C      PRINT *, "JB", JB
      END
