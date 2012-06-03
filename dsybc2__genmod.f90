        !COMPILER-GENERATED INTERFACE MODULE: Fri May 11 21:04:10 2012
        MODULE DSYBC2__genmod
          INTERFACE 
            SUBROUTINE DSYBC2(UPLO,A,DESCA,AR,DESCAR,AC,DESCAC,CHKVEC,NB&
     &,INFO)
              INTEGER(KIND=8) :: NB
              INTEGER(KIND=8) :: DESCA(*)
              CHARACTER(LEN=1) :: UPLO
              REAL(KIND=8) :: A(DESCA((9)),*)
              REAL(KIND=8) ,POINTER :: AR(:,:)
              INTEGER(KIND=8) :: DESCAR(*)
              REAL(KIND=8) ,POINTER :: AC(:,:)
              INTEGER(KIND=8) :: DESCAC(*)
              REAL(KIND=8) :: CHKVEC(NB)
              INTEGER(KIND=8) :: INFO
            END SUBROUTINE DSYBC2
          END INTERFACE 
        END MODULE DSYBC2__genmod
