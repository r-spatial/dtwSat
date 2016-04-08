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
C     CM - Input local cost and output cumulative cost matrix
C     DM - Direction matrix 
C     VM - Starting points matrix 
C     SM - Matrix of step patterns 
C     N  - Number of rows in CM, DM, and VM 
C     M  - Number of columns CM, DM, and VM 
C     NS - Number of rows in SM 
      SUBROUTINE bestmatches(XM, AM, DM, DP, X, IM, A, K, P, L, OV)
C  800 FORMAT('i: ',I5,'   i: ',I5,'   VALUE: ',F8.4,'   VALUE: ',F10.4)
C  801 FORMAT('Here: ',I5)
C     I/O Variables       
      INTEGER K, P, L, XM(K,2), X(K), DP(P), IM(P-1,3), A(K)
      DOUBLE PRECISION AM(P-1,L), DD, DM(K), OV
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
            D1 = XM(I,1)
            D2 = XM(I,2)
            IL = X(I)
            IF ((D2.LT.B1).OR.(D1.GT.B2)) THEN
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
                GOTO 20 
            ENDIF
            IF( DM(I).GE.AM(J,IL) ) THEN
                GOTO 20 
            ENDIF
            AM(J,IL) = DM(I)
            IF( DM(I).LT.DD ) THEN
                IM(J,1) = IL
                IM(J,2) = I
                IM(J,3) = A(I)
                DD = DM(I)
            ENDIF
   20    CONTINUE
   30 CONTINUE
      END
      
      
      
      
      
