
program chk2_t
use ftchol
implicit none


integer            m , n, mxllda, mxlocc, mb, nb
parameter         (  m = 6, n = 6, mb = 3, nb = 3)
parameter         ( mxllda = 3, mxlocc = 3 )

!double precision  one
!parameter         ( one = 1.0d+0 )


integer           ictxt, info, mycol, myrow, npcol, nprow
integer           desca( dlen_ )

double precision  A( mxllda, mxlocc )
!double precision, allocatable :: ac(:,:), ar(:,:)

integer           i, j, dbg


nprow = 2; npcol = 2
call sl_init( ictxt, nprow, npcol )
call blacs_gridinfo( ictxt, nprow, npcol, myrow, mycol )
if ( myrow.eq.-1 ) print *, 'grid error'

call descinit( desca, m, n, mb, nb, 0, 0, ictxt, mxllda, info)
if ( info .ne. 0 ) print *, 'info', info
!call random_number( A )
do i = 1, mxllda
   do j = 1, mxlocc
      A( i, j )  = mxlocc * (i-1) + j
   end do
end do
A = A + 10 * myrow + mycol
if ( myrow == 0 .and. mycol == 0 ) print *, '===== Matrix A =====' 
call pmat2('aa',A, desca)

call dsybc2('L', A, desca,  nb, info)
!SUBROUTINE CHK1( UPLO, A, IA, JA, DA,  INFO)
call chk1('L', A, 1, 1, DESCA, INFO)
!SUBROUTINE CHK2(UPLO, M, N, A, IA, JA, DESCA, INFO)
call chk2('L', 3, 3, A, 1, 1, DESCA, INFO)



!if ( myrow == 0 .and. mycol == 0 ) print *, '==== Matrix AR ===='
call pmat2('ar',ar, descar)

!if ( myrow == 0 .and. mycol == 0 ) print *, '==== Matrix AC ===='
call pmat2('ac',ac, descac)


deallocate(ac, ar)
call blacs_gridexit( ictxt )
call blacs_exit( 0 )

contains 
   ! print a small matrix on screen
   ! n < 10

end program chk2_t

