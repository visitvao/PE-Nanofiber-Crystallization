
C  *******************************************************************
      SUBROUTINE SUBDELTAENERGY(IPRESENTCHAIN, IPRESENTBEAD,
     >  XALTERNATIVE, YALTERNATIVE, ZALTERNATIVE,
     >  BONDVECTOR, BONDVECTORSECONDSHELL, bondvectorthirdshell,
     >  A, B, E, F, ae, bf)

      INTEGER  CHAINL, LATTICES, CHAINN

      PARAMETER (CHAINL=100,LATTICES=100,CHAINN=200)

      INTEGER  I
      INTEGER  A, B, E, F, ae, bf
      INTEGER  IPRESENTCHAIN, IPRESENTBEAD
      INTEGER  CHAINLENGTH, CHAINNUMBER
      INTEGER  LATTICESIZEX, LATTICESIZEY, LATTICESIZEZ
      INTEGER  XTEST2, YTEST2, ZTEST2
      INTEGER  XALTERNATIVE, YALTERNATIVE, ZALTERNATIVE
      INTEGER  BONDVECTOR(12,3)
      INTEGER  BONDVECTORSECONDSHELL(42,3)
      integer  bondvectorthirdshell(92,3)
      INTEGER  X(CHAINN,-2:CHAINL+3), Y(CHAINN,-2:CHAINL+3),
     >         Z(CHAINN,-2:CHAINL+3)
      INTEGER  XP(CHAINN,-2:CHAINL+3), YP(CHAINN,-2:CHAINL+3),
     >         ZP(CHAINN,-2:CHAINL+3)
      INTEGER  BEADOCCUPANCY(0:LATTICES-1,0:LATTICES-1,0:LATTICES-1)
      integer  xtemp,ytemp,ztemp,xtempm2,ytempm2,ztempm2
      integer  xtempp2,ytempp2,ztempp2

      COMMON  /SYSTEMSIZE/      CHAINLENGTH, CHAINNUMBER, 
     >			LATTICESIZEX, LATTICESIZEY, LATTICESIZEZ
      COMMON  /BEADOCCUPANCY/   BEADOCCUPANCY
      COMMON  /COORDINATES/     X, Y, Z, XP, YP, ZP

c  --------------------------------
c  List the iterated variables here
c  --------------------------------
      xtemp = XP(IPRESENTCHAIN,IPRESENTBEAD)
      ytemp = YP(IPRESENTCHAIN,IPRESENTBEAD)
      ztemp = ZP(IPRESENTCHAIN,IPRESENTBEAD)

      xtempm2 = XP(IPRESENTCHAIN,IPRESENTBEAD-2)
      ytempm2 = YP(IPRESENTCHAIN,IPRESENTBEAD-2)
      ztempm2 = ZP(IPRESENTCHAIN,IPRESENTBEAD-2)

      xtempp2 = XP(IPRESENTCHAIN,IPRESENTBEAD+2)
      ytempp2 = YP(IPRESENTCHAIN,IPRESENTBEAD+2)
      ztempp2 = ZP(IPRESENTCHAIN,IPRESENTBEAD+2)
c  --------------------------------
c  List the iterated variables here
c  --------------------------------
     
      A = 0
      DO I = 1, 12
        XTEST2 = xtemp + BONDVECTOR(I,1)
        YTEST2 = ytemp + BONDVECTOR(I,2)
        ZTEST2 = ztemp + BONDVECTOR(I,3)
        IF( IPRESENTBEAD .GE. 3 ) THEN
          IF( (XTEST2 .EQ. xtempm2) .AND.
     >        (YTEST2 .EQ. ytempm2) .AND.
     >        (ZTEST2 .EQ. ztempm2) ) THEN
            A = A - 1
          ENDIF
        ENDIF
        IF( IPRESENTBEAD .LE. CHAINLENGTH-2 ) THEN
          IF( (XTEST2 .EQ. xtempp2) .AND.
     >        (YTEST2 .EQ. ytempp2) .AND.
     >        (ZTEST2 .EQ. ztempp2) ) THEN
            A = A - 1
          ENDIF
        ENDIF
c        CALL SUBFOLDEDPOINT(XTEST2,YTEST2,ZTEST2)

c-- in-line adoption of the folded point calculation
  103 IF( XTEST2 .GE. LATTICESIZEX ) THEN
        XTEST2 = XTEST2 - LATTICESIZEX
        GO TO 103
      END IF
  104 IF( XTEST2 .LT. 0 ) THEN
        XTEST2 = XTEST2 + LATTICESIZEX
        GO TO 104
      END IF

  105 IF( YTEST2 .GE. LATTICESIZEY ) THEN
        YTEST2 = YTEST2 - LATTICESIZEY
        GO TO 105
      END IF
  106 IF( YTEST2 .LT. 0 ) THEN
        YTEST2 = YTEST2 + LATTICESIZEY
        GO TO 106
      END IF

  107 IF( ZTEST2 .GE. LATTICESIZEZ ) THEN
        ZTEST2 = ZTEST2 - LATTICESIZEZ
        GO TO 107
      END IF
  108 IF( ZTEST2 .LT. 0 ) THEN
        ZTEST2 = ZTEST2 + LATTICESIZEZ
        GO TO 108
      END IF
c-- in-line adoption of the folded point calculation

        IF(BEADOCCUPANCY(XTEST2,YTEST2,ZTEST2) .EQ. 1) A = A + 1
      ENDDO

      B = 0
      DO I = 1, 12
        XTEST2 = XALTERNATIVE + BONDVECTOR(I,1)
        YTEST2 = YALTERNATIVE + BONDVECTOR(I,2)
        ZTEST2 = ZALTERNATIVE + BONDVECTOR(I,3)
        IF( IPRESENTBEAD .GE. 3 ) THEN
          IF( (XTEST2 .EQ. xtempm2) .AND.
     >        (YTEST2 .EQ. ytempm2) .AND.
     >        (ZTEST2 .EQ. ztempm2) ) THEN
            B = B - 1
          ENDIF
        ENDIF
        IF( IPRESENTBEAD .LE. CHAINLENGTH-2 ) THEN
          IF( (XTEST2 .EQ. xtempp2) .AND.
     >        (YTEST2 .EQ. ytempp2) .AND.
     >        (ZTEST2 .EQ. ztempp2) ) THEN
            B = B - 1
          ENDIF
        ENDIF
c        CALL SUBFOLDEDPOINT(XTEST2,YTEST2,ZTEST2)

c-- in-line adoption of the folded point calculation
  203 IF( XTEST2 .GE. LATTICESIZEX ) THEN
        XTEST2 = XTEST2 - LATTICESIZEX
        GO TO 203
      END IF
  204 IF( XTEST2 .LT. 0 ) THEN
        XTEST2 = XTEST2 + LATTICESIZEX
        GO TO 204
      END IF

  205 IF( YTEST2 .GE. LATTICESIZEY ) THEN
        YTEST2 = YTEST2 - LATTICESIZEY
        GO TO 205
      END IF
  206 IF( YTEST2 .LT. 0 ) THEN
        YTEST2 = YTEST2 + LATTICESIZEY
        GO TO 206
      END IF

  207 IF( ZTEST2 .GE. LATTICESIZEZ ) THEN
        ZTEST2 = ZTEST2 - LATTICESIZEZ
        GO TO 207
      END IF
  208 IF( ZTEST2 .LT. 0 ) THEN
        ZTEST2 = ZTEST2 + LATTICESIZEZ
        GO TO 208
      END IF
c-- in-line adoption of the folded point calculation

        IF(BEADOCCUPANCY(XTEST2,YTEST2,ZTEST2) .EQ. 1) B = B + 1
      ENDDO


      E = 0
      DO I = 1, 42
        XTEST2 = xtemp +
     >           BONDVECTORSECONDSHELL(I,1)
        YTEST2 = ytemp +
     >           BONDVECTORSECONDSHELL(I,2)
        ZTEST2 = ztemp +
     >           BONDVECTORSECONDSHELL(I,3)
        IF( IPRESENTBEAD .GE. 3 ) THEN
          IF( (XTEST2 .EQ. xtempm2) .AND.
     >        (YTEST2 .EQ. ytempm2) .AND.
     >        (ZTEST2 .EQ. ztempm2) ) THEN
            E = E - 1
          ENDIF
        ENDIF
        IF( IPRESENTBEAD .LE. CHAINLENGTH-2 ) THEN
          IF( (XTEST2 .EQ. xtempp2) .AND.
     >        (YTEST2 .EQ. ytempp2) .AND.
     >        (ZTEST2 .EQ. ztempp2) ) THEN
            E = E - 1
          ENDIF
        ENDIF
c        CALL SUBFOLDEDPOINT(XTEST2,YTEST2,ZTEST2)

c-- in-line adoption of the folded point calculation
  303 IF( XTEST2 .GE. LATTICESIZEX ) THEN
        XTEST2 = XTEST2 - LATTICESIZEX
        GO TO 303
      END IF
  304 IF( XTEST2 .LT. 0 ) THEN
        XTEST2 = XTEST2 + LATTICESIZEX
        GO TO 304
      END IF

  305 IF( YTEST2 .GE. LATTICESIZEY ) THEN
        YTEST2 = YTEST2 - LATTICESIZEY
        GO TO 305
      END IF
  306 IF( YTEST2 .LT. 0 ) THEN
        YTEST2 = YTEST2 + LATTICESIZEY
        GO TO 306
      END IF

  307 IF( ZTEST2 .GE. LATTICESIZEZ ) THEN
        ZTEST2 = ZTEST2 - LATTICESIZEZ
        GO TO 307
      END IF
  308 IF( ZTEST2 .LT. 0 ) THEN
        ZTEST2 = ZTEST2 + LATTICESIZEZ
        GO TO 308
      END IF
c-- in-line adoption of the folded point calculation

        IF(BEADOCCUPANCY(XTEST2,YTEST2,ZTEST2) .EQ. 1) E = E + 1
      ENDDO

      F = 0
      DO I = 1, 42
        XTEST2 = XALTERNATIVE + BONDVECTORSECONDSHELL(I,1)
        YTEST2 = YALTERNATIVE + BONDVECTORSECONDSHELL(I,2)
        ZTEST2 = ZALTERNATIVE + BONDVECTORSECONDSHELL(I,3)
        IF( IPRESENTBEAD .GE. 3 ) THEN
          IF( (XTEST2 .EQ. xtempm2) .AND.
     >        (YTEST2 .EQ. ytempm2) .AND.
     >        (ZTEST2 .EQ. ztempm2) ) THEN
            F = F - 1
          ENDIF
        ENDIF
        IF( IPRESENTBEAD .LE. CHAINLENGTH-2 ) THEN
          IF( (XTEST2 .EQ. xtempp2) .AND.
     >        (YTEST2 .EQ. ytempp2) .AND.
     >        (ZTEST2 .EQ. ztempp2) ) THEN
            F = F - 1
          ENDIF
        ENDIF
c        CALL SUBFOLDEDPOINT(XTEST2,YTEST2,ZTEST2)

c-- in-line adoption of the folded point calculation
  403 IF( XTEST2 .GE. LATTICESIZEX ) THEN
        XTEST2 = XTEST2 - LATTICESIZEX
        GO TO 403
      END IF
  404 IF( XTEST2 .LT. 0 ) THEN
        XTEST2 = XTEST2 + LATTICESIZEX
        GO TO 404
      END IF

  405 IF( YTEST2 .GE. LATTICESIZEY ) THEN
        YTEST2 = YTEST2 - LATTICESIZEY
        GO TO 405
      END IF
  406 IF( YTEST2 .LT. 0 ) THEN
        YTEST2 = YTEST2 + LATTICESIZEY
        GO TO 406
      END IF

  407 IF( ZTEST2 .GE. LATTICESIZEZ ) THEN
        ZTEST2 = ZTEST2 - LATTICESIZEZ
        GO TO 407
      END IF
  408 IF( ZTEST2 .LT. 0 ) THEN
        ZTEST2 = ZTEST2 + LATTICESIZEZ
        GO TO 408
      END IF
c-- in-line adoption of the folded point calculation

        IF(BEADOCCUPANCY(XTEST2,YTEST2,ZTEST2) .EQ. 1) F = F + 1
      ENDDO


      ae = 0
      DO I = 1, 92
        XTEST2 = xtemp +
     >           bondvectorthirdshell(I,1)
        YTEST2 = ytemp +
     >           bondvectorthirdshell(I,2)
        ZTEST2 = ztemp +
     >           bondvectorthirdshell(I,3)
c        IF( IPRESENTBEAD .GE. 3 ) THEN
c          IF( (XTEST2 .EQ. XP(IPRESENTCHAIN,IPRESENTBEAD-2)) .AND.
c     >        (YTEST2 .EQ. YP(IPRESENTCHAIN,IPRESENTBEAD-2)) .AND.
c     >        (ZTEST2 .EQ. ZP(IPRESENTCHAIN,IPRESENTBEAD-2)) ) THEN
c            ae = ae - 1
c          ENDIF
c        ENDIF
c        IF( IPRESENTBEAD .LE. CHAINLENGTH-2 ) THEN
c          IF( (XTEST2 .EQ. XP(IPRESENTCHAIN,IPRESENTBEAD+2)) .AND.
c     >        (YTEST2 .EQ. YP(IPRESENTCHAIN,IPRESENTBEAD+2)) .AND.
c     >        (ZTEST2 .EQ. ZP(IPRESENTCHAIN,IPRESENTBEAD+2)) ) THEN
c            ae = ae - 1
c          ENDIF
c        ENDIF
c        CALL SUBFOLDEDPOINT(XTEST2,YTEST2,ZTEST2)

c-- in-line adoption of the folded point calculation
  503 IF( XTEST2 .GE. LATTICESIZEX ) THEN
        XTEST2 = XTEST2 - LATTICESIZEX
        GO TO 503
      END IF
  504 IF( XTEST2 .LT. 0 ) THEN
        XTEST2 = XTEST2 + LATTICESIZEX
        GO TO 504
      END IF

  505 IF( YTEST2 .GE. LATTICESIZEY ) THEN
        YTEST2 = YTEST2 - LATTICESIZEY
        GO TO 505
      END IF
  506 IF( YTEST2 .LT. 0 ) THEN
        YTEST2 = YTEST2 + LATTICESIZEY
        GO TO 506
      END IF

  507 IF( ZTEST2 .GE. LATTICESIZEZ ) THEN
        ZTEST2 = ZTEST2 - LATTICESIZEZ
        GO TO 507
      END IF
  508 IF( ZTEST2 .LT. 0 ) THEN
        ZTEST2 = ZTEST2 + LATTICESIZEZ
        GO TO 508
      END IF
c-- in-line adoption of the folded point calculation

        IF(BEADOCCUPANCY(XTEST2,YTEST2,ZTEST2) .EQ. 1) ae = ae + 1
      ENDDO

      bf = 0
      DO I = 1, 92
        XTEST2 = XALTERNATIVE + bondvectorthirdshell(I,1)
        YTEST2 = YALTERNATIVE + bondvectorthirdshell(I,2)
        ZTEST2 = ZALTERNATIVE + bondvectorthirdshell(I,3)
c        IF( IPRESENTBEAD .GE. 3 ) THEN
c          IF( (XTEST2 .EQ. XP(IPRESENTCHAIN,IPRESENTBEAD-2)) .AND.
c     >        (YTEST2 .EQ. YP(IPRESENTCHAIN,IPRESENTBEAD-2)) .AND.
c     >        (ZTEST2 .EQ. ZP(IPRESENTCHAIN,IPRESENTBEAD-2)) ) THEN
c            bf = bf - 1
c          ENDIF
c        ENDIF
c        IF( IPRESENTBEAD .LE. CHAINLENGTH-2 ) THEN
c          IF( (XTEST2 .EQ. XP(IPRESENTCHAIN,IPRESENTBEAD+2)) .AND.
c     >        (YTEST2 .EQ. YP(IPRESENTCHAIN,IPRESENTBEAD+2)) .AND.
c     >        (ZTEST2 .EQ. ZP(IPRESENTCHAIN,IPRESENTBEAD+2)) ) THEN
c            bf = bf - 1
c          ENDIF
c        ENDIF
c        CALL SUBFOLDEDPOINT(XTEST2,YTEST2,ZTEST2)

c-- in-line adoption of the folded point calculation
  603 IF( XTEST2 .GE. LATTICESIZEX ) THEN
        XTEST2 = XTEST2 - LATTICESIZEX
        GO TO 603
      END IF
  604 IF( XTEST2 .LT. 0 ) THEN
        XTEST2 = XTEST2 + LATTICESIZEX
        GO TO 604
      END IF

  605 IF( YTEST2 .GE. LATTICESIZEY ) THEN
        YTEST2 = YTEST2 - LATTICESIZEY
        GO TO 605
      END IF
  606 IF( YTEST2 .LT. 0 ) THEN
        YTEST2 = YTEST2 + LATTICESIZEY
        GO TO 606
      END IF

  607 IF( ZTEST2 .GE. LATTICESIZEZ ) THEN
        ZTEST2 = ZTEST2 - LATTICESIZEZ
        GO TO 607
      END IF
  608 IF( ZTEST2 .LT. 0 ) THEN
        ZTEST2 = ZTEST2 + LATTICESIZEZ
        GO TO 608
      END IF
c-- in-line adoption of the folded point calculation

        IF(BEADOCCUPANCY(XTEST2,YTEST2,ZTEST2) .EQ. 1) bf = bf + 1
      ENDDO

      RETURN
      END
