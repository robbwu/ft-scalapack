program pdpotrf_t1
use ftchol
implicit none

double precision, allocatable :: A(:, :), B(:, :)
integer N, NB, maxlld, nprow, npcol, ictxt, myrow, mycol, info
integer DA(DLEN_), DB(DLEN_)
integer IA, JA, JB
CHARACTER          COLBTOP, ROWBTOP , UPLO

integer ICEIL
external ICEIL


nprow = 7; npcol = 7
N = 1200; NB = 20 
maxlld = ICEIL(N/NB, nprow) * NB


call sl_init( ictxt, nprow, npcol )
call blacs_gridinfo( ictxt, nprow, npcol, myrow, mycol )

allocate(A(maxlld, maxlld), B(maxlld, maxlld))
call descinit(DA, N, N, NB, NB, 0, 0, ictxt, maxlld, info)
call descinit(DB, N, N, NB, NB, 0, 0, ictxt, maxlld, info)
call set_lower_one(B, DB)

!call pmat2('B', B, DB)
call PDSYRK('L', 'N', N, N, ONE, B, 1, 1, DB, ZERO, A, 1, 1, DA)



!CALL DSYBC2('L', A, DA, DA( MB_ ), INFO)
!!call pmat2('A', A, DA)
!!call pmat2('AR', AR, DESCAR)

!CALL PB_TOPGET( ICTXT, 'Broadcast', 'Rowwise', ROWBTOP ) 
!CALL PB_TOPGET( ICTXT, 'Broadcast', 'Columnwise', COLBTOP ) 


!CALL PB_TOPSET( ICTXT, 'Broadcast', 'Rowwise', 'S-ring' ) 
!CALL PB_TOPSET( ICTXT, 'Broadcast', 'Columnwise', ' ' ) 

!IA = 1; JA = 1
!UPLO = 'L'
!JB = NB

!!JN = MIN( ICEIL( JA, DESCA( NB_ ) )*DESCA( NB_ ), JA+N-1 ) 
!!JB = JN - JA + 1 
!!                                                                       
!!        Perform unblocked Cholesky factorization on JB block           
!!                                                                       
!DO IA = 1, N, NB
   !JA = IA
   !!PRINT *, IA
   !CALL PDPOTF2( UPLO, NB, A, IA, JA, DA, INFO ) 
   !!CALL pmat2('A', A, DA)
   !IF( INFO.NE.0 )  PRINT *,'PDPOTF2 INFO=', INFO                                              
   !!SUBROUTINE CHK1( UPLO, A, IA, JA, DA,  INFO)
   !CALL CHK1( UPLO,  A, IA, JA, DA, INFO )
   !!CALL pmat2('AR', AR, DESCAR)

   !IF ( IA+NB.LE.N ) THEN
      !CALL PDTRSM( 'Right', UPLO, 'Transpose', 'Non-Unit',        &
      !&                   N-NB*(IA/NB+1), JB, ONE, A, IA, JA, DA, A, IA+JB, JA, &
      !&                   DA )                                        
      !CALL CHK2(UPLO, N-JB*(IA/NB+1), JB,  A, IA, JA, DA, INFO)
      !!CALL pmat2('A', A, DA)
      !!CALL pmat2('AR', AR, DESCAR)

      !CALL PDSYRK( UPLO, 'No Transpose', N-JB*(IA/NB+1), JB, -ONE, A, IA+JB,&
      !&                   JA, DA, ONE, A, IA+JB, JA+JB, DA )       
      !CALL CHK3(UPLO, N-JB*(IA/NB+1),  A, IA, JA, DA, INFO)
      !IF ( IA.EQ.4 ) THEN
         !!CALL pmat2('A', A, DA)
         !CALL pmat2('AR', AR, DESCAR)
         !CALL pmat2('AC', AC, DESCAC)
      !END IF
   !END IF
!END DO
!!CALL pmat2('A', A, DA)
!!CALL pmat2('AR', AR, DESCAR)
!!CALL pmat2('AC', AC, DESCAC)

!CALL DSYDC2



! FT_DPOTRF testing
! =================

!SUBROUTINE FT_PDPOTRF( UPLO, N, A, IA, JA, DESCA, INFO ) 
CALL FT_PDPOTRF('L', N, A, 1, 1, DA, INFO )
!PRINT '(9I3)', DESCAR
!PRINT *, SHAPE(AR)
!CALL pmat2('A', A, DA)
!CALL pmat2('AR', AR, DESCAR)
!CALL pmat2('AC', AC, DESCAC)
CALL DSYDC2
deallocate(A, B)

call blacs_gridexit( ictxt )
call blacs_exit( 0 )
contains

subroutine set_lower_one(A, DA)
implicit none

integer DA(*),  i, j, lld, nprow, npcol, myrow, mycol,  n, nb, &
               li, lj, prow, pcol, k, l
double precision A( DA(LLD_), *)

call blacs_gridinfo( DA( CTXT_ ), nprow, npcol, myrow, mycol )

n = DA(N_); nb = DA(NB_)
A( 1:DA(LLD_), 1:DA(LLD_) ) = ZERO
do i = 1, n, nb
   do j =1, n, nb
      call infog2l( i, j, DA, nprow, npcol, myrow, mycol, li, lj, prow, pcol )
      if ( myrow.eq.prow .and. mycol.eq.pcol ) then
         if ( i.eq.j ) then
            do k=1, nb
               do l=1,k
                  A(li+k-1, lj+l-1) = ONE
               end do
            end do
         else if ( i.gt.j ) then
            A(li:li+nb-1, lj:lj+nb-1) = ONE
         end if
      end if
   end do
end do
end subroutine set_lower_one

end program
