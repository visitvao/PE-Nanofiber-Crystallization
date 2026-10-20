
C
C --------------------------------------------------------------------
      SUBROUTINE SUBSNNDRISSTATE(X1,Y1,Z1,X2,Y2,Z2,THISVECTORTYPE)
C --------------------------------------------------------------------
C
C  SUBROUTINE TO ASSIGN THE TYPE OF VECTOR THAT IS DESCRIBED BY THE
C  TWO POINTS ( X1,Y1,Z1 ) AND ( X2,Y2,Z2 ).
C
C  THE VECTOR TYPES ARE THE FOLLOWING:
C          INDICES:  1: THE BEADS 1-3 ARE 5.00 ANGSTROM APART (tt)
C                    2: THE BEADS 1-3 ARE 4.33 ANGSTROM APART (t+, t-, +t, -t)
C                    3: THE BEADS 1-3 ARE 3.53 ANGSTROM APART (++, --)
C                    4: THE BEADS 1-3 ARE 2.50 ANGSTROM APART (+-, -+)


      INTEGER  X1, Y1, Z1, X2, Y2, Z2
      INTEGER  QSUM, QSUMSQUARE
      INTEGER  THISVECTORTYPE


      QSUMSQUARE = ( X2 - X1 )**2 + ( Y2 - Y1 )**2 + ( Z2 - Z1 )**2

      IF( QSUMSQUARE .EQ. 1 ) THEN
        THISVECTORTYPE = 4

      ELSE IF( QSUMSQUARE .EQ. 2 ) THEN
        QSUM = ( X2 - X1 ) + ( Y2 - Y1 ) + ( Z2 - Z1 )
        IF( QSUM .EQ. 0 ) THEN
          THISVECTORTYPE = 4
        ELSE
          IF( ABS(QSUM) .NE. 2 ) THEN
            WRITE(15,*) 'PROBLEM WITH QSUMSQUARE=2 IN ENERGY CALC.'
            WRITE(15,*) QSUMSQUARE
            STOP 5
          ENDIF
          THISVECTORTYPE = 2
        ENDIF

      ELSE IF( QSUMSQUARE .EQ. 5 ) THEN
        THISVECTORTYPE = 2

      ELSE IF( QSUMSQUARE .EQ. 6 ) THEN
        THISVECTORTYPE = 2

      ELSE IF( QSUMSQUARE .EQ. 3 ) THEN
        THISVECTORTYPE = 3

      ELSE
        IF( (QSUMSQUARE .NE. 4) .AND. (QSUMSQUARE .NE. 8) ) THEN
          WRITE(15,*) 'PROBELM WITH QSUMSQUARE IN SUBSNNDRISSTATE'
          WRITE(15,*) QSUMSQUARE
          STOP 6
        ENDIF
        THISVECTORTYPE = 1

      END IF

      RETURN
      END
