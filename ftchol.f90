MODULE FTCHOL
IMPLICIT NONE
INTEGER            BLOCK_CYCLIC_2D, CSRC_, CTXT_, DLEN_, DTYPE_,  &
                      LLD_, MB_, M_, NB_, N_, RSRC_                  
PARAMETER          ( BLOCK_CYCLIC_2D = 1, DLEN_ = 9, DTYPE_ = 1,  &
                        CTXT_ = 2, M_ = 3, N_ = 4, MB_ = 5, NB_ = 6, &
                        RSRC_ = 7, CSRC_ = 8, LLD_ = 9 )             
DOUBLE PRECISION   ONE, ZERO 
PARAMETER          ( ONE = 1.0D+0, ZERO = 0.0D+0 ) 

DOUBLE PRECISION, ALLOCATABLE :: AR(:, :), AC(:, :), CHKVEC(:)
INTEGER :: DESCAC( DLEN_ ), DESCAR( DLEN_ )

CONTAINS 
subroutine pmat2(name, A, DA)
implicit none
character(len=*) name
character(len=20) editdesc

integer DA( * ), m, n, mb, nb, nprow, npcol,myrow, mycol, &
         li, lj, prow, pcol, i, j, k
double precision A( DA( lld_ ), * )


m = DA( M_ ); n = DA( N_ ); mb = DA( MB_ ); nb = DA( NB_ )

write (editdesc, '(A, I1, A)') '(', nb, 'F8.3)'

if ( mod(m, mb) .NE. 0 ) then 
   print *, 'm%mb!=0'
   return
end if
if ( mod(n, nb) .NE. 0 ) then 
   print *, 'n%nb!=0'
   return
end if

call blacs_gridinfo( DA( CTXT_ ), nprow, npcol, myrow, mycol )

do i = 1, m, mb
   do j =1, (i/mb)*nb+1, nb
      call infog2l( i, j, DA, nprow, npcol, myrow, mycol, li, lj, prow, pcol )
      if ( myrow.eq.prow .and. mycol.eq.pcol ) then
         print '(A, A, I3,I3, A, I3, I3)', name,' proc ', myrow, mycol, ' bi, bj ', i/mb+1, j/nb+1
         do k=li,li+mb-1
            print editdesc, A(k, lj:lj+nb-1) 
         end do
         call flush
      end if
      call blacs_barrier(DA(CTXT_), 'A')
   end do
end do


end subroutine pmat2

subroutine pmat(name, A, dA) 
   !use ifcore
   implicit none
   character(len=*) name
   !integer m_, n_, mb_, nb_, lld_, dlen_, ctxt_, rsrc_, csrc_
   !parameter ( ctxt_ = 2, rsrc_ = 7, csrc_ = 8 )
   !parameter ( m_ = 3, n_ = 4, mb_ = 5, nb_ = 6, lld_ = 9 , dlen_ = 9 )

   integer           m, n, i, dA( dlen_ ), ii,jj
   double precision  A( dA( lld_ ), * )
   character(20)     editdesc
   integer           myrow, mycol, nprow, npcol
   integer           numroc
   logical           res
   external          numroc

   call blacs_gridinfo( dA(ctxt_), nprow, npcol, myrow, mycol )
   m = numroc( dA(m_), dA(mb_), myrow, dA(rsrc_), nprow )
   n = numroc( dA(n_), dA(nb_), mycol, dA(csrc_), npcol )
   write (editdesc, '(A, I1, A)') '(', n, 'F8.3)'

   do ii = 0, nprow-1
      do jj = 0, npcol-1
         if ( myrow .eq. ii .and. mycol .eq. jj ) then
            print '(A, A, I3, I3, A, I3, I3)', name, ' from process ',myrow, mycol, ' dim ',m, n
            !call flush(6)
            do i = 1,m
               print editdesc, A(i,1:n)
            end do
            !call flush(6)
            !res = commitqq(6)
         end if
         call flush
         call blacs_barrier( dA( ctxt_ ), 'A')
      end do
   end do
end subroutine pmat

SUBROUTINE DSYBC2(UPLO, A, DESCA,  NB, INFO)
   !USE CHKVAR
   IMPLICIT NONE
   ! -- FT-ScaLAPACK routine --
   ! Panruo Wu( pwu@mines.edu )
   ! Colorado School of Mines
   ! 2012-03-15
   !
   !
   INTEGER           DESCA( * )
   DOUBLE PRECISION  A( DESCA(9), * )
   CHARACTER         UPLO

   INTEGER           NB, INFO
   !DOUBLE PRECISION  CHKVEC(NB)

   !DOUBLE PRECISION, ALLOCATABLE :: AR(:, :), AC(:, :)
   !INTEGER  :: DESCAC(*), DESCAR(*)
   
   ! Purpose
   ! =======
   ! 
   ! DSYBC2 builds blockwise double checksum of 2d cyclic distributed matrix
   ! A into AR(for row checksums) and AC(for column checksums).
   !
   ! AR and AC are distributed very similar to the way A is distributed, say
   ! 2d-cyclic distribution. 
   !
   !  Notes                                                                
   !  =====                                                                
   !                                                                       
   !  Each global data object is described by an associated description    
   !  vector.  This vector stores the information required to establish    
   !  the mapping between an object element and its corresponding process  
   !  and memory location.                                                 
   !                                                                       
   !  Let A be a generic term for any 2D block cyclicly distributed array. 
   !  Such a global array has an associated description vector DESCA.      
   !  In the following comments, the character _ should be read as         
   !  "of the global array".                                               
   !                                                                       
   !  NOTATION        STORED IN      EXPLANATION                           
   !  --------------- -------------- --------------------------------------
   !  DTYPE_A(global) DESCA( DTYPE_ )The descriptor type.  In this case,   
   !                                 DTYPE_A = 1.                          
   !  CTXT_A (global) DESCA( CTXT_ ) The BLACS context handle, indicating  
   !                                 the BLACS process grid A is distribu- 
   !                                 ted over. The context itself is glo-  
   !                                 bal, but the handle (the integer      
   !                                 value) may vary.                      
   !  M_A    (global) DESCA( M_ )    The number of rows in the global      
   !                                 array A.                              
   !  N_A    (global) DESCA( N_ )    The number of columns in the global   
   !                                 array A.                              
   !  MB_A   (global) DESCA( MB_ )   The blocking factor used to distribute
   !                                 the rows of the array.                
   !  NB_A   (global) DESCA( NB_ )   The blocking factor used to distribute
   !                                 the columns of the array.             
   !  RSRC_A (global) DESCA( RSRC_ ) The process row over which the first  
   !                                 row of the array A is distributed.    
   !  CSRC_A (global) DESCA( CSRC_ ) The process column over which the     
   !                                 first column of the array A is        
   !                                 distributed.                          
   !  LLD_A  (local)  DESCA( LLD_ )  The leading dimension of the local    
   !                                 array.  LLD_A >= MAX(1,LOCr(M_A)).    

   ! Arguments 
   ! =========
   ! UPLO      (global input) 'L': reference lower part of A; 'U' upper part.
   !
   ! A         (local input) local pointer to an array of dimension (LLD_A, *).
   !
   ! DESCA     (global and local input) descriptor for A
   !
   ! AR        (local output) local pointer to array that stores the row
   !           checksum matrix. On exit, the required memory is allocated
   !           and pointed by this pointer. On entry, this is an allocatable
   !           array pointer.
   !
   ! DESCAR    (global and local output) descriptor for AR
   !
   ! AC        (local output) local pointer to array that stores the column
   !           checksum matrix. On exit, the required memory is allocated
   !           and pointed by this pointer. On entry, this is an allocatable
   !           array pointer.
   !
   ! DESCAC    (global and local output) descriptor for AC
   !
   ! INFO      (global output) INTEGER
   !           = 0: successful exit
   !
   ! Assumptions & Limitations:
   ! ==========================
   !  1. Only support the whole A, not its submatrix sub( A )
   !  2. RSRC_A = CSRC_A = 0
   !  3. process grid is square (NPROW = NPCOL)
   !  
   !  This limitations are basically to ease implementation;
   !  they can be lifted once we have more time.
   !
   ! =============================================================================

   !INTEGER            BLOCK_CYCLIC_2D, CSRC_, CTXT_, DLEN_, DTYPE_,  &
      !&                   LLD_, MB_, M_, NB_, N_, RSRC_                  
   !PARAMETER          ( BLOCK_CYCLIC_2D = 1, DLEN_ = 9, DTYPE_ = 1,  &
      !&                     CTXT_ = 2, M_ = 3, N_ = 4, MB_ = 5, NB_ = 6, &
      !&                     RSRC_ = 7, CSRC_ = 8, LLD_ = 9 )             
   !DOUBLE PRECISION   ONE, ZERO 
   !PARAMETER          ( ONE = 1.0D+0, ZERO = 0.0D+0 ) 

   INTEGER           I, J, IB, JB, MYROW, MYCOL, NPROW, NPCOL, ICTXT, GI, GJ
   INTEGER           LMA, LNA, LMAR, LNAC, STAT, II, JJ, ALLOSTAT, LLDA
   INTEGER           IPROC, JPROC, IIR, IIC, M, N

   LOGICAL           LSAME
   INTEGER           ICEIL, NUMROC
   EXTERNAL          NUMROC, ICEIL, DGEMV, LSAME

   DOUBLE PRECISION  WORK( NB ), ONES( NB )
   LOGICAL           UPPER

   UPPER = LSAME( UPLO, 'U' )


   ICTXT = DESCA( CTXT_ )
   CALL BLACS_GRIDINFO( ICTXT, NPROW, NPCOL, MYROW, MYCOL )

   !print *, 'DSYBC2: GRIDINFO', ICTXT, NPROW, NPCOL, MYROW, MYCOL

   ! calculate memory size for local AC and AR
   ! TODO: add allocation error checking here
   
   M = DESCA( M_ ); N = DESCA( N_ )
   LMA = NUMROC( DESCA( M_ ), DESCA( MB_ ), MYROW, DESCA( RSRC_ ), NPROW )
   LNA = NUMROC( DESCA( N_ ), DESCA( NB_ ), MYCOL, DESCA( CSRC_ ), NPCOL )
   LLDA = DESCA( LLD_ )

   !print *, 'DSYBC2: LMA, LNA:', LMA, LNA

   LMAR = ICEIL( LMA, NB ) * 2
   LNAC = ICEIL( LNA, NB ) * 2

   !print *, 'DSYBC2: LMAR, LNAC', LMAR, LNAC
   IF (.NOT.ALLOCATED(CHKVEC)) THEN
      ALLOCATE(CHKVEC(NB))
      DO I=1,NB
         CHKVEC(I) = I
      END DO
   END IF
   !PRINT *,'CHKVEC:',CHKVEC

   IF (ALLOCATED(AR)) DEALLOCATE(AR)
   IF (ALLOCATED(AC)) DEALLOCATE(AC)
   ALLOCATE(AR(LMAR, LNA), STAT=ALLOSTAT)
   ALLOCATE(AC(LMA, LNAC), STAT=ALLOSTAT)
   IF ( ALLOSTAT .NE. 0 ) PRINT *, 'ALLOSTAT', ALLOSTAT
   
   DESCAR( 1:DLEN_ ) = DESCA( 1:DLEN_ )
   DESCAR( M_ ) = ICEIL( DESCA( M_ ), DESCA( MB_ ) ) * 2
   DESCAR( MB_ ) = 2
   DESCAR( LLD_ ) = LMAR

   DESCAC( 1:DLEN_ ) = DESCA( 1:DLEN_ )
   DESCAC( N_ ) = ICEIL( DESCA( N_ ), DESCA( NB_ ) ) * 2
   DESCAC( NB_ ) = 2
   DESCAC( LLD_ ) = LMA

   ONES = ONE


   IF ( UPPER ) THEN
      ! TODO: upper triangle case
   ELSE
      DO I = 1, M, NB
         IB = MIN( M-I+1, NB )
         DO J = 1, I-IB, NB
            JB = NB
            CALL INFOG2L(I, J, DESCA, NPROW, NPCOL, MYROW, MYCOL, &
                  II, JJ, IPROC, JPROC)
            !PRINT *,'PROC', MYROW, MYCOL,'I,J', I, J
            IF ( MYROW .EQ. IPROC .AND. MYCOL .EQ. JPROC ) THEN
               !PRINT  '(A, I3, I3, A, I3, I3, A, I3, I3, A, I3, I3)', &
                        !'PROC ', MYROW, MYCOL, 'I, J', I, J, 'II, JJ', II, JJ,  &
                        !'IB, JB', IB, JB
               CALL DGEMV( 'Trans', IB, JB, ONE, A(II,JJ), LLDA, ONES, 1, ZERO, WORK, 1 )
               AR( II/NB*2 + 1, JJ:JJ+JB-1 ) = WORK(1:JB)
               CALL DGEMV( 'Trans', IB, JB, ONE, A(II,JJ), LLDA, CHKVEC, 1, ZERO, WORK, 1 )
               AR( II/NB*2 + 2, JJ:JJ+JB-1 ) = WORK( 1:JB )
               CALL DGEMV( 'Notrans', IB, JB, ONE, A(II,JJ), LLDA, ONES, 1, ZERO, WORK, 1 )
               AC( II:II+IB-1, JJ/NB*2 + 1 ) = WORK( 1:IB )
               CALL DGEMV( 'Notrans', IB, JB, ONE, A(II,JJ), LLDA, CHKVEC, 1, ZERO, WORK, 1 )
               AC( II:II+IB-1, JJ/NB*2 + 2 ) = WORK( 1:IB )
            END IF
         END DO
         J = I
         IB = MIN( M-I+1, NB)
         JB = IB
         CALL INFOG2L(I, J, DESCA, NPROW, NPCOL, MYROW, MYCOL, &
            II, JJ, IPROC, JPROC)
         IF ( MYROW .EQ. IPROC .AND. MYCOL .EQ. JPROC ) THEN
            CALL DSYMV( UPLO, IB, ONE, A(II,JJ), LLDA, ONES, 1, ZERO, WORK, 1 )
            AR( II/NB*2 + 1, JJ:JJ+JB-1 ) = WORK(1:JB)
            AC( II:II+IB-1, JJ/NB*2 + 1 ) = WORK(1:JB)
            CALL DSYMV( UPLO, IB, ONE, A(II,JJ), LLDA, CHKVEC, 1, ZERO, WORK, 1 )
            AR( II/NB*2 + 2, JJ:JJ+JB-1 ) = WORK( 1:JB )
            AC( II:II+JB-1, JJ/NB*2 + 2 ) = WORK( 1:JB )
         END IF
      END DO
               

   END IF


   ! TODO: provide proper error code
   INFO = 0

END SUBROUTINE DSYBC2

SUBROUTINE CHK1( UPLO, A, IA, JA, DA,  INFO)
! Check after PDPOTF2 procedure
IMPLICIT NONE
CHARACTER      UPLO
INTEGER        N, IA, JA, INFO
INTEGER        DA( DLEN_ ), DAR( DLEN_ ), DAC( DLEN_ )
DOUBLE PRECISION A( DA( LLD_ ), * )

INTEGER        MYROW, MYCOL, NPROW, NPCOL, IPROC, JPROC
INTEGER        I, J, MB, NB, II, JJ
!DOUBLE PRECISION  ZERO
!PARAMETER         ( ZERO = 0.0D+0 )
LOGICAL        UPPER, LSAME
INTEGER        ICEIL
EXTERNAL       ICEIL, LSAME
DOUBLE PRECISION  WORK( DA( MB_ ) )

MB = DA( MB_ )
NB = DA( NB_ )
DAR = DESCAR(1:DLEN_)
DAC = DESCAC(1:DLEN_)

CALL BLACS_GRIDINFO( DA( CTXT_ ), NPROW, NPCOL, MYROW, MYCOL ) 
CALL INFOG2L( IA, JA, DA, NPROW, NPCOL, MYROW, MYCOL, &
           & I, J, IPROC, JPROC )

! check if the requested submatrix is within one block
IF ( IA+NB-1.GT.DA( M_ ) .OR. JA+NB-1.GT.DA( N_ ) ) THEN
   IF ( MYROW.EQ.0 .AND. MYCOL.EQ.0 ) PRINT *, 'CHK1: exceeds boundary'
   RETURN
END IF
IF ( MOD( IA-1, MB ).NE.0 .OR. MOD( JA-1, NB ).NE.0 ) THEN
   IF ( MYROW.EQ.0 .AND. MYCOL.EQ.0 ) PRINT *, 'CHK1: IA, JA not aligned'
   RETURN
END IF
IF ( MB.NE.NB ) THEN
   IF ( MYROW.EQ.0 .AND. MYCOL.EQ.0 ) PRINT*, 'CHK1: MB, NB are not equal'
   RETURN
END IF


UPPER = LSAME( UPLO, 'U' )
IF ( UPPER ) THEN
   IF ( MYROW.EQ.0 .AND. MYCOL.EQ.0 ) PRINT *, 'CHK1: UPPER not supported yet.'
ELSE
   IF ( MYROW.EQ.IPROC .AND. MYCOL.EQ.JPROC ) THEN
      DO JJ = J, J+NB-1
         AR( I/NB*2 + 1, JJ ) = SUM( A(I+JJ-J:I+NB-1, JJ) )
         AR( I/NB*2 + 2, JJ ) = DOT_PRODUCT( A(I+JJ-J:I+NB-1, JJ), &
            CHKVEC(JJ-J+1:NB) )
      END DO
      AC( I:I+NB-1, J/NB*2+1:J/NB*2+2 ) = ZERO
   END IF
END IF


END SUBROUTINE CHK1

SUBROUTINE CHK2(UPLO, M, N, A, IA, JA, DESCA, INFO)
IMPLICIT NONE
CHARACTER   UPLO
INTEGER  M, N, IA, JA, DESCA(*), INFO
DOUBLE PRECISION A( DESCA( LLD_ ), * )
INTEGER  IAR, MAR, JAC

INTEGER NB, I, J, K
INTEGER NPROW, NPCOL, MYROW, MYCOL, LRIND, LCIND, ROW, COL

NB = DESCA(NB_)


IF ( MOD( IA-1, DESCA(MB_) ).NE.0 ) THEN 
   INFO = -5
   RETURN
END IF
IF ( MOD( JA-1, DESCA(NB_) ).NE.0 ) THEN
   INFO = -6
   RETURN
END IF
IF ( MOD( M, DESCA(MB_) ).NE.0 ) THEN
   INFO = -2
   RETURN
END IF
IF ( N.NE.DESCA(NB_) ) THEN
   INFO = -3
   RETURN
END IF

IAR = (IA+NB)/NB*2 + 1
MAR = (M/NB)*2
JAC = JA/NB*2 + 1
!PRINT *, 'IAR, MAR, JAC', IAR, MAR, JAC

CALL PDTRSM('Right', UPLO, 'Transpose', 'Non-Unit', &
      MAR, N, ONE, A, IA, JA, DESCA, AR, IAR, JA, DESCAR)

! checksum checking phase
!CALL BLACS_GRIDINFO( DESCA(CTXT_), NPROW, NPCOL, MYROW, MYCOL )
!DO I=IA+NB:IA+NB+M-1,NB
   !CALL INFOG2L( I, JA, 


END SUBROUTINE CHK2

SUBROUTINE CHK3(UPLO,  N, A, IA, JA, DESCA, INFO)
IMPLICIT NONE
CHARACTER   UPLO
INTEGER  N, IA, JA, DESCA(*), INFO
DOUBLE PRECISION A( DESCA( LLD_ ), * )
INTEGER  IAR, MAR, JAC

INTEGER NB

NB = DESCA(NB_)


IF ( MOD( IA-1, NB ).NE.0 ) THEN 
   INFO = -5
   RETURN
END IF
IF ( MOD( JA-1, NB ).NE.0 ) THEN
   INFO = -6
   RETURN
END IF
IF ( MOD( N, NB ).NE.0 ) THEN
   INFO = -3
   RETURN
END IF


IAR = (IA+NB)/NB*2 + 1
JAC = (JA+NB)/NB*2 + 1
MAR = (N/NB)*2

!PRINT '(9I3)', MAR, NB, N, IAR, JA, IA+NB, JA, IAR, JA+NB
CALL PDGEMM('N', 'T', MAR, N, NB,  &
   -ONE, AR, IAR, JA, DESCAR, &
   A, IA+NB, JA, DESCA, &
   ONE, AR, IAR, JA+NB, DESCAR)
CALL PDGEMM('N', 'T', N, MAR,NB, &
   -ONE, A, IA+NB, JA, DESCA, &
   AR, IAR, JA, DESCAR, &
   ONE, AC, IA+NB, JAC, DESCAC)

   

END SUBROUTINE CHK3

END MODULE FTCHOL