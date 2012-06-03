SUBROUTINE DSYBC2(UPLO, A, DESCA, AR, DESCAR, AC, DESCAC, CHKVEC, NB, INFO)
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
   DOUBLE PRECISION  CHKVEC(NB)

   DOUBLE PRECISION, POINTER:: AR(:, :), AC(:, :)
   INTEGER  :: DESCAC(*), DESCAR(*)
   
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

   INTEGER            BLOCK_CYCLIC_2D, CSRC_, CTXT_, DLEN_, DTYPE_,  &
      &                   LLD_, MB_, M_, NB_, N_, RSRC_                  
   PARAMETER          ( BLOCK_CYCLIC_2D = 1, DLEN_ = 9, DTYPE_ = 1,  &
      &                     CTXT_ = 2, M_ = 3, N_ = 4, MB_ = 5, NB_ = 6, &
      &                     RSRC_ = 7, CSRC_ = 8, LLD_ = 9 )             
   DOUBLE PRECISION   ONE, ZERO 
   PARAMETER          ( ONE = 1.0D+0, ZERO = 0.0D+0 ) 

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
         DO J = 1, I-IB, NB
            IB = MIN( M-I+1, NB )
            JB = NB
            CALL INFOG2L(I, J, DESCA, NPROW, NPCOL, MYROW, MYCOL, &
                  II, JJ, IPROC, JPROC)
            IF ( MYROW .EQ. IPROC .AND. MYCOL .EQ. JPROC ) THEN
               PRINT  '(A, I3, I3, A, I3, I3, A, I3, I3, A, I3, I3)', &
                        'PROC ', MYROW, MYCOL, 'I, J', I, J, 'II, JJ', II, JJ,  &
                        'IB, JB', IB, JB
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


