      SUBROUTINE SUBMATMAT(A,B,AB,MFORM,NFORM,LFORM)
C
      INTEGER     MFORM, NFORM, LFORM
      INTEGER     I, J, K

      REAL        A(MFORM,NFORM), B(NFORM,LFORM), AB(MFORM,LFORM)


      DO 5 I = 1, MFORM
        DO 6 J = 1, LFORM
          AB(I,J) = 0.0
    6   CONTINUE
    5 CONTINUE

      DO 10 I = 1, MFORM
	DO 20 J = 1, LFORM
          DO 30 K = 1, NFORM
            AB(I,J) = AB(I,J) + ( A(I,K) * B(K,J) )
   30     CONTINUE
   20   CONTINUE
   10 CONTINUE

      RETURN
      END


      SUBROUTINE SUBMATVEC(A, V, AV, MFORM, NFORM)

      INTEGER      MFORM, NFORM
      INTEGER      I, J

      REAL         A(MFORM,NFORM), V(NFORM), AV(MFORM)

      DO 5 I = 1, MFORM
	AV(I) = 0.0
    5 CONTINUE

      DO 10 I = 1, MFORM
	DO 20 J = 1, NFORM
          AV(I) = AV(I) + ( A(I,J) * V(J) )
   20   CONTINUE
   10 CONTINUE

      RETURN
      END

      SUBROUTINE SUBMATVEC2(A, V, AV, MFORM, NFORM, INDEX)

      INTEGER      MFORM, NFORM, INDEX
      INTEGER      I, J

      REAL         A(MFORM,NFORM,4), V(NFORM), AV(MFORM)

      DO 5 I = 1, MFORM
	AV(I) = 0.0
    5 CONTINUE

      DO 10 I = 1, MFORM
	DO 20 J = 1, NFORM
          AV(I) = AV(I) + ( A(I,J,INDEX) * V(J) )
   20   CONTINUE
   10 CONTINUE

      RETURN
      END


      SUBROUTINE SUBVECMAT(V, A, VA, MFORM, NFORM)
C
      INTEGER      MFORM, NFORM
      INTEGER      I, J

      REAL         A(MFORM,NFORM), V(MFORM), VA(NFORM)

      DO 5 I = 1, NFORM
        VA(I) = 0.0
    5 CONTINUE

      DO 10 I = 1, NFORM
        DO 20 J = 1, MFORM
          VA(I) = VA(I) + ( V(J) * A(J,I) )
   20   CONTINUE
   10 CONTINUE

      RETURN
      END

      SUBROUTINE SUBVECMAT2(V, A, VA, MFORM, NFORM, INDEX)
C
      INTEGER      MFORM, NFORM, INDEX
      INTEGER      I, J

      REAL         A(MFORM,NFORM,4), V(MFORM), VA(NFORM)

      DO 5 I = 1, NFORM
        VA(I) = 0.0
    5 CONTINUE

      DO 10 I = 1, NFORM
        DO 20 J = 1, MFORM
          VA(I) = VA(I) + ( V(J) * A(J,I,INDEX) )
   20   CONTINUE
   10 CONTINUE

      RETURN
      END

      SUBROUTINE SUBMATMATR8(A,B,AB,MFORM,NFORM,LFORM)
C
      INTEGER     MFORM, NFORM, LFORM
      INTEGER     I, J, K

      REAL*8      A(MFORM,NFORM), B(NFORM,LFORM), AB(MFORM,LFORM)


      DO 5 I = 1, MFORM
        DO 6 J = 1, LFORM
          AB(I,J) = 0.0
    6   CONTINUE
    5 CONTINUE

      DO 10 I = 1, MFORM
	DO 20 J = 1, LFORM
          DO 30 K = 1, NFORM
            AB(I,J) = AB(I,J) + ( A(I,K) * B(K,J) )
   30     CONTINUE
   20   CONTINUE
   10 CONTINUE

      RETURN
      END


      SUBROUTINE SUBMATVECR8(A, V, AV, MFORM, NFORM)

      INTEGER      MFORM, NFORM
      INTEGER      I, J

      REAL*8       A(MFORM,NFORM), V(NFORM), AV(MFORM)

      DO 5 I = 1, MFORM
	AV(I) = 0.0
    5 CONTINUE

      DO 10 I = 1, MFORM
	DO 20 J = 1, NFORM
          AV(I) = AV(I) + ( A(I,J) * V(J) )
   20   CONTINUE
   10 CONTINUE

      RETURN
      END


      SUBROUTINE SUBVECMATR8(V, A, VA, MFORM, NFORM)
C
      INTEGER      MFORM, NFORM
      INTEGER      I, J

      REAL*8       A(MFORM,NFORM), V(MFORM), VA(NFORM)

      DO 5 I = 1, NFORM
        VA(I) = 0.0
    5 CONTINUE

      DO 10 I = 1, NFORM
        DO 20 J = 1, MFORM
          VA(I) = VA(I) + ( V(J) * A(J,I) )
   20   CONTINUE
   10 CONTINUE

      RETURN
      END
