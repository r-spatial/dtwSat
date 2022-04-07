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
C   Find best matches TWDTW - 2016-03-26                      C
C                                                             C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C     XM - Two columns matrix with 'from' and 'to' dates as integers 
C     AM - Matrix (P-1 x L) with classification intervals (P-1) and possible classes (L) 
C     DM - Vector length K DTW distance for each alignment 
C     DP - Vector length P-1 with classification intervals 
C     X  - Vector length K with alignments index 
C     IM - Matrix to return best matches, intervals, and class labels 
C     A  - Vector with alignment index
C     K  - Number of alignments 
C     P  - Number of dates defineing classification intervals 
C     L  - Number of classes 
C     OV - Minimum temporal overlap 
      SUBROUTINE bestmatches(XM,AM,DM,DP,X,IM,DB,A,K,P,L,OV)
C     I/O Variables 
      INTEGER K, P, L, XM(K,2), X(K), DP(P), IM(P-1,4), A(K)
      DOUBLE PRECISION AM(P-1,L), DD, DM(K), OV, DB(P-1)
C     Internals
      DOUBLE PRECISION R 
      INTEGER I, J, IL, B1, B2, D1, D2
C     For all time intervals 
      DO 30 J = 1, P-1
         B1 = DP(J)
         B2 = DP(J+1)
         DD = AM(J,1)
         
C     For all TWDTW matches
         DO 20 I = 1, K
C            print *, "I: ", I
            D1 = XM(I,1)
            D2 = XM(I,2)
            IL = X(I)
            IF ((D2.LT.B1).OR.(D1.GT.B2)) THEN
C                print *, "D1: ", D1, "D2: ", D2
                GOTO 20 
            ENDIF
            IF (D1.LT.B1) THEN
                D1 = B1
            ENDIF
            IF (B2.LT.D2) THEN
                D2 = B2
            ENDIF
            R = REAL(D2 - D1) / REAL(B2 - B1)
            IF( .NOT.(OV.LE.R.AND.R.LE.(2-OV)) ) THEN
C                print *, "R: ", R
                GOTO 20 
            ENDIF
            IF( DM(I).GE.AM(J,IL) ) THEN
C                print *, "DM(I): ", DM(I)
                GOTO 20 
            ENDIF
            AM(J,IL) = DM(I)
            IF( DM(I).LT.DD ) THEN
                DD = DM(I)
                IM(J,1) = IL
                IM(J,2) = I
                IM(J,3) = A(I)
C                print *, "IM: ", IM
                DB(J) = DD
            ENDIF
   20    CONTINUE
   30 CONTINUE
      END
      
      
      
      
      
