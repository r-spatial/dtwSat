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
C   Efficient computation of DTW trace back - 2015-10-17      C
C                                                             C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C     DM   - Direction matrix 
C     SM   - Matrix of step patterns 
C     JMIN - Positions o the minimum points 
C     IND1 - Alignment indices in the pattern
C     IND2 - Alignment indices in the template
C     POS  - Starting points in IND1 and IND2 
C     N    - Number of rows in DM 
C     M    - Number of columns DM 
C     NS   - Number of rows in SM 
C     NJ   - Number of minimum points 
C     AL   - Length of IND1 and IND2
      SUBROUTINE tracepath(DM,SM,JMIN,IND1,IND2,POS,N,M,NS,NJ,AL)
C  800 FORMAT(I4,I4,I4,I4,I4,I4,I4,I4)
C     I/O Variables       
      INTEGER N, M, NS, NJ, AL
      INTEGER IND1(AL),IND2(AL),POS(NJ+1),SM(NS,4),DM(N,M),JMIN(NJ)
C     Internals
      INTEGER K, PK, I, J, P, PI, PJ, S, IS(NS), STEPS(NS,4*NS)
      INTEGER ZERO, ONE
      PARAMETER(ZERO=0,ONE=1)
C     Initialize steps 
      DO 10 PK = 1, SM(NS,1)
         STEPS(PK,1) = ZERO
         IS(PK) = 2
   10 CONTINUE
C     Get steps direction
      DO 11 K = 1, NS 
         PK = SM(K,1)
         PI = SM(K,2)
         PJ = SM(K,3)
         IF (PI.EQ.ZERO.AND.PJ.EQ.ZERO) THEN 
            GO TO 11
         ENDIF 
         STEPS(PK,1) = STEPS(PK,1) + ONE
         STEPS(PK,IS(PK)) = PI
         STEPS(PK,IS(PK)+1) = PJ
         IS(PK) = IS(PK) + 2
   11 CONTINUE
C     Trace back      
      P = ONE
      POS(1) = ZERO
      DO 30 JM = 1, NJ
         I = N 
         J = JMIN(JM) 
         S = DM(I,J) 
   20    IF (I.NE.ONE.AND.P.LT.AL) THEN 
            IND1(P) = I 
            IND2(P) = J 
            NS = STEPS(S,1)
            DO 12 K = 1, NS
               P = P + 1
               IND1(P) = I - STEPS(S,2*K)
               IND2(P) = J - STEPS(S,2*K+1)
               I = IND1(P)
               J = IND2(P)
   12       CONTINUE
            S = DM(I,J)
            GO TO 20
         ENDIF
         POS(JM+1) = P
         P = P + 1
   30 CONTINUE
      END
