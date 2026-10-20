
C --------------------------------------------------------------------
      SUBROUTINE SUBALTERNATIVPOSONE(X1,Y1,Z1,X2,Y2,Z2,X3,Y3,Z3,
     >                   XALTERNATIVE,YALTERNATIVE,ZALTERNATIVE,
     >                   BONDVECTOR,THISMOVE,QSUMSQUARE)
C --------------------------------------------------------------------
C
C  SUBROUTINE TO SELECT _THE_ ALTERNATIVE POSITION FOR THE SELECTED BEAD
C  (X2,Y2,Z2). THE NEW POSITION (XALTERNATIVE,YALTERNATIVE,ZALTERNATIVE)
C  IS DETERMINED UNIQUELY!

      INTEGER  I
      INTEGER  X1,Y1,Z1, X2,Y2,Z2, X3,Y3,Z3
D     INTEGER  QSUM
      INTEGER  QSUMSQUARE
      INTEGER  VECTORX,VECTORY,VECTORZ
      INTEGER  XALTERNATIVE,YALTERNATIVE,ZALTERNATIVE
      INTEGER  XALTERNATIVE1,YALTERNATIVE1,ZALTERNATIVE1
      INTEGER  XALTERNATIVE2,YALTERNATIVE2,ZALTERNATIVE2
      INTEGER  THISMOVE
      INTEGER  BONDVECTOR(12,3)

      THISMOVE = 0

      IF( QSUMSQUARE .EQ. 2 ) THEN

D       QSUM = (X3-X1) + (Y3-Y1) + (Z3-Z1)
D       IF(( QSUM .NE. 2 ) .AND. ( QSUM .NE. -2 ) ) THEN
D         WRITE(4,*) 'PROBLEM WITH MOVES QSUMSQUARE = 2 / QSUM = 0'
D       ENDIF

        IF( X3 .EQ. X1 ) THEN
          XALTERNATIVE1 = X2
          YALTERNATIVE1 = Y2 + 1
          ZALTERNATIVE1 = Z2 - 1
          XALTERNATIVE2 = X2
          YALTERNATIVE2 = Y2 - 1
          ZALTERNATIVE2 = Z2 + 1
          THISMOVE = 7
D         IF( Y3-Y1 .EQ. 1 ) THEN
D           IF( Z3-Z1 .NE. 1) THEN
D             WRITE(4,*) 'PROBLEM WITH MOVE 7'
D           END IF
D         ELSE IF( Y3-Y1 .EQ. -1) THEN
D           IF( Z3-Z1 .NE. -1) THEN
D             WRITE(4,*) 'PROBLEM WITH MOVE 7'
D           END IF
D         ELSE
D           WRITE(4,*) 'PROBLEM WITH MOVE 7'
D         END IF
        ELSE IF( Y3 .EQ. Y1 ) THEN
          XALTERNATIVE1 = X2 - 1
          YALTERNATIVE1 = Y2
          ZALTERNATIVE1 = Z2 + 1
          XALTERNATIVE2 = X2 + 1
          YALTERNATIVE2 = Y2
          ZALTERNATIVE2 = Z2 - 1
          THISMOVE = 8
D         IF( X3-X1 .EQ. 1 ) THEN
D           IF( Z3-Z1 .NE. 1) THEN
D             WRITE(4,*) 'PROBLEM WITH MOVE 8'
D           END IF
D         ELSE IF( X3-X1 .EQ. -1) THEN
D           IF( Z3-Z1 .NE. -1) THEN
D             WRITE(4,*) 'PROBLEM WITH MOVE 8'
D           END IF
D         ELSE
D           WRITE(4,*) 'PROBLEM WITH MOVE 8'
D         END IF
        ELSE
          XALTERNATIVE1 = X2 - 1
          YALTERNATIVE1 = Y2 + 1
          ZALTERNATIVE1 = Z2
          XALTERNATIVE2 = X2 + 1
          YALTERNATIVE2 = Y2 - 1
          ZALTERNATIVE2 = Z2
          THISMOVE = 9
D         IF( Z3 .NE. Z1 ) THEN
D           WRITE(4,*) 'PROBLEM WITH QSUMSQUARE = 2 / QSUM = 2 ; -2'
D         ENDIF
D         IF( X3-X1 .EQ. 1 ) THEN
D           IF( Y3-Y1 .NE. 1) THEN
D             WRITE(4,*) 'PROBLEM WITH MOVE 9'
D           END IF
D         ELSE IF( X3-X1 .EQ. -1) THEN
D           IF( Y3-Y1 .NE. -1) THEN
D             WRITE(4,*) 'PROBLEM WITH MOVE 9'
D           END IF
D         ELSE
D           WRITE(4,*) 'PROBLEM WITH MOVE 9'
D         END IF
        END IF

      ELSE IF( QSUMSQUARE .EQ. 5 ) THEN

        IF( ABS(X3-X1) .EQ. 1 ) THEN
          XALTERNATIVE1 = X2 + 1
          YALTERNATIVE1 = Y2
          ZALTERNATIVE1 = Z2
          XALTERNATIVE2 = X2 - 1
          YALTERNATIVE2 = Y2
          ZALTERNATIVE2 = Z2
          THISMOVE = 10
D         IF( X3-X1 .EQ. 1 ) THEN
D           IF( Y3-Y1 .EQ. 0 ) THEN
D             IF( Z3-Z1 .NE. -2) THEN
D               WRITE(4,*) 'PROBLEM WITH MOVE 10'
D             END IF
D           ELSE IF( Y3-Y1 .EQ. -2) THEN
D             IF( Z3-Z1 .NE. 0) THEN
D               WRITE(4,*) 'PROBLEM WITH MOVE 10'
D             END IF
D           ELSE
D             WRITE(4,*) 'PROBLEM WITH MOVE 10'
D           END IF
D         ELSE
D           IF( X3-X1 .NE. -1 ) THEN
D             WRITE(4,*) 'PROBLEM WITH MOVE 10'
D           END IF
D           IF( Y3-Y1 .EQ. 0 ) THEN
D             IF( Z3-Z1 .NE. 2) THEN
D               WRITE(4,*) 'PROBLEM WITH MOVE 10'
D             END IF
D           ELSE IF( Y3-Y1 .EQ. 2) THEN
D             IF( Z3-Z1 .NE. 0) THEN
D               WRITE(4,*) 'PROBLEM WITH MOVE 10'
D             END IF
D           ELSE
D             WRITE(4,*) 'PROBLEM WITH MOVE 10'
D           END IF
D         END IF

        ELSE IF( ABS(Y3-Y1) .EQ. 1 ) THEN
          XALTERNATIVE1 = X2
          YALTERNATIVE1 = Y2 + 1
          ZALTERNATIVE1 = Z2
          XALTERNATIVE2 = X2
          YALTERNATIVE2 = Y2 - 1
          ZALTERNATIVE2 = Z2
          THISMOVE = 11
D         IF( Y3-Y1 .EQ. 1 ) THEN
D           IF( X3-X1 .EQ. 0 ) THEN
D             IF( Z3-Z1 .NE. -2) THEN
D               WRITE(4,*) 'PROBLEM WITH MOVE 11'
D             END IF
D           ELSE IF( X3-X1 .EQ. -2) THEN
D             IF( Z3-Z1 .NE. 0) THEN
D               WRITE(4,*) 'PROBLEM WITH MOVE 11'
D             END IF
D           ELSE
D             WRITE(4,*) 'PROBLEM WITH MOVE 11'
D           END IF
D         ELSE
D           IF( Y3-Y1 .NE. -1 ) THEN
D             WRITE(4,*) 'PROBLEM WITH MOVE 11'
D           END IF
D           IF( X3-X1 .EQ. 0 ) THEN
D             IF( Z3-Z1 .NE. 2) THEN
D               WRITE(4,*) 'PROBLEM WITH MOVE 11'
D             END IF
D           ELSE IF( X3-X1 .EQ. 2) THEN
D             IF( Z3-Z1 .NE. 0) THEN
D               WRITE(4,*) 'PROBLEM WITH MOVE 11'
D             END IF
D           ELSE
D             WRITE(4,*) 'PROBLEM WITH MOVE 11'
D           END IF
D         END IF

        ELSE
          XALTERNATIVE1 = X2
          YALTERNATIVE1 = Y2
          ZALTERNATIVE1 = Z2 + 1
          XALTERNATIVE2 = X2
          YALTERNATIVE2 = Y2
          ZALTERNATIVE2 = Z2 - 1
          THISMOVE = 12
D         IF( Z3-Z1 .EQ. 1 ) THEN
D           IF( X3-X1 .EQ. 0 ) THEN
D             IF( Y3-Y1 .NE. -2) THEN
D               WRITE(4,*) 'PROBLEM WITH MOVE 12'
D             END IF
D           ELSE IF( X3-X1 .EQ. -2) THEN
D             IF( Y3-Y1 .NE. 0) THEN
D               WRITE(4,*) 'PROBLEM WITH MOVE 12'
D             END IF
D           ELSE
D             WRITE(4,*) 'PROBLEM WITH MOVE 12'
D           END IF
D         ELSE
D           IF( Z3-Z1 .NE. -1 ) THEN
D             WRITE(4,*) 'PROBLEM WITH MOVE 12'
D           END IF
D           IF( X3-X1 .EQ. 0 ) THEN
D             IF( Y3-Y1 .NE. 2) THEN
D               WRITE(4,*) 'PROBLEM WITH MOVE 12'
D             END IF
D           ELSE IF( X3-X1 .EQ. 2) THEN
D             IF( Y3-Y1 .NE. 0) THEN
D               WRITE(4,*) 'PROBLEM WITH MOVE 12'
D             END IF
D           ELSE
D             WRITE(4,*) 'PROBLEM WITH MOVE 12'
D           END IF
D         END IF
        END IF

      ELSE
D       IF( QSUMSQUARE .NE. 6 ) THEN
D         WRITE(4,*) 'PROBLEM 1 WITH QSUMSQUARE = 6'
D         WRITE(4,*) QSUMSQUARE
D       ENDIF
        IF( ABS( X3-X1 ) .EQ. 2 ) THEN
          XALTERNATIVE1 = X2
          YALTERNATIVE1 = Y2 + 1
          ZALTERNATIVE1 = Z2 - 1
          XALTERNATIVE2 = X2
          YALTERNATIVE2 = Y2 - 1
          ZALTERNATIVE2 = Z2 + 1
          THISMOVE = 13
D         IF( X3-X1 .EQ. -2) THEN
D           IF( ( Y3-Y1 .NE. 1 ) .OR.
D    >          ( Z3-Z1 .NE. 1 ) ) THEN
D               WRITE(4,*) 'PROBLEM WITH MOVE 13'
D           END IF
D         ELSE
D           IF( ( X3-X1 .NE. 2 ) .OR.
D    >          ( Y3-Y1 .NE. -1 ) .OR.
D    >          ( Z3-Z1 .NE. -1 ) ) THEN
D             WRITE(4,*) 'PROBLEM WITH MOVE 13'
D           END IF
D         END IF

        ELSE IF( ABS ( Y3-Y1 ) .EQ. 2 ) THEN
          XALTERNATIVE1 = X2 + 1
          YALTERNATIVE1 = Y2
          ZALTERNATIVE1 = Z2 - 1
          XALTERNATIVE2 = X2 - 1
          YALTERNATIVE2 = Y2
          ZALTERNATIVE2 = Z2 + 1
          THISMOVE = 14
D         IF( Y3-Y1 .EQ. -2) THEN
D           IF( ( X3-X1 .NE. 1 ) .OR.
D    >          ( Z3-Z1 .NE. 1 ) ) THEN
D               WRITE(4,*) 'PROBLEM WITH MOVE 14'
D           END IF
D         ELSE
D           IF( ( Y3-Y1 .NE. 2 ) .OR.
D    >          ( X3-X1 .NE. -1 ) .OR.
D    >          ( Z3-Z1 .NE. -1 ) ) THEN
D             WRITE(4,*) 'PROBLEM WITH MOVE 14'
D           END IF
D         END IF

        ELSE
          XALTERNATIVE1 = X2 + 1
          YALTERNATIVE1 = Y2 - 1
          ZALTERNATIVE1 = Z2
          XALTERNATIVE2 = X2 - 1
          YALTERNATIVE2 = Y2 + 1
          ZALTERNATIVE2 = Z2
          THISMOVE = 15
D         IF( ABS (Z3-Z1) .NE. 2 ) THEN
D           WRITE(4,*) 'PROBLEM 2 WITH QSUMSQUARE = 6'
D           WRITE(4,*) QSUMSQUARE
D         END IF
D         IF( Z3-Z1 .EQ. -2) THEN
D           IF( ( X3-X1 .NE. 1 ) .OR.
D    >          ( Y3-Y1 .NE. 1 ) ) THEN
D               WRITE(4,*) 'PROBLEM WITH MOVE 15'
D           END IF
D         ELSE
D           IF( ( Z3-Z1 .NE. 2 ) .OR.
D    >          ( X3-X1 .NE. -1 ) .OR.
D    >          ( Y3-Y1 .NE. -1 ) ) THEN
D             WRITE(4,*) 'PROBLEM WITH MOVE 15'
D           END IF
D         END IF
        END IF
      END IF

C  *****************************************************

D     IF( THISMOVE .EQ. 0 ) THEN
D       WRITE(4,*) 'PROBLEM WITH ASSIGNING ANY MOVE'
D     END IF

C  *****************************************************

C  CHECK WHICH OF THE TWO ALTERNATIVES IS THE CORRECT ONE:
      VECTORX = XALTERNATIVE1 - X1
      VECTORY = YALTERNATIVE1 - Y1
      VECTORZ = ZALTERNATIVE1 - Z1
      DO I = 1, 12
        IF( ( VECTORX.EQ.BONDVECTOR(I,1) ) .AND.
     >      ( VECTORY.EQ.BONDVECTOR(I,2) ) .AND.
     >      ( VECTORZ.EQ.BONDVECTOR(I,3) ) ) THEN
          XALTERNATIVE = XALTERNATIVE1
          YALTERNATIVE = YALTERNATIVE1
          ZALTERNATIVE = ZALTERNATIVE1
D         GO TO 60
          GO TO 70
        ENDIF
      ENDDO

      XALTERNATIVE = XALTERNATIVE2
      YALTERNATIVE = YALTERNATIVE2
      ZALTERNATIVE = ZALTERNATIVE2
C  CHECK WHICH OF THE TWO ALTERNATIVES IS THE CORRECT ONE:

C  *****************************************************

C  TEST - IF ALTERNATIVE = ALTERNATIVE2 BOND-CHECK HAS NOT BEEN PERFORMED
D     VECTORX = XALTERNATIVE - X1
D     VECTORY = YALTERNATIVE - Y1
D     VECTORZ = ZALTERNATIVE - Z1
D     DO I = 1, 12
D       IF( ( VECTORX.EQ.BONDVECTOR(I,1) ) .AND.
D    >      ( VECTORY.EQ.BONDVECTOR(I,2) ) .AND.
D    >      ( VECTORZ.EQ.BONDVECTOR(I,3) ) ) THEN
D         GO TO 60
D       ENDIF
D     ENDDO
D     WRITE(4,*) 'PROBLEM 1 WITH CHOOSING ONE OF THE TWO ALTERNATIVES'
C  TEST - IF ALTERNATIVE = ALTERNATIVE2 BOND-CHECK HAS NOT BEEN PERFORMED

C  *****************************************************

C  TEST - BOND-CHECK TO SECOND BEAD
D  60 VECTORX = XALTERNATIVE - X3
D     VECTORY = YALTERNATIVE - Y3
D     VECTORZ = ZALTERNATIVE - Z3
D     DO I = 1, 12
D       IF( ( VECTORX.EQ.BONDVECTOR(I,1) ) .AND.
D    >      ( VECTORY.EQ.BONDVECTOR(I,2) ) .AND.
D    >      ( VECTORZ.EQ.BONDVECTOR(I,3) ) ) THEN
D         GO TO 70
D       ENDIF
D     ENDDO
D     WRITE(4,*) 'PROBLEM 2 WITH CHOOSING ONE OF THE TWO ALTERNATIVES'
C  TEST - BOND-CHECK TO SECOND BEAD

C  *****************************************************

   70 RETURN
      END
