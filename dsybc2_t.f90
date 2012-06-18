program dsybc2_t
! unit test program for dsybc2
!
!use ifport
!use ifcore
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
call pmat(A, desca)

!do i=1, mb
   !chkvec(i) = real(i, kind(one))
!end do
call dsybc2('L', A, desca,  nb, info)


if ( myrow == 0 .and. mycol == 0 ) print *, '==== Matrix AR ===='
call pmat(ar, descar)

if ( myrow == 0 .and. mycol == 0 ) print *, '==== Matrix AC ===='
call pmat(ac, descac)


deallocate(ac, ar)
call blacs_gridexit( ictxt )
call blacs_exit( 0 )

contains 
   ! print a small matrix on screen
   ! n < 10
   subroutine pmat(A, dA) 
      !use ifcore
      implicit none
      integer m_, n_, mb_, nb_, lld_, dlen_, ctxt_, rsrc_, csrc_
      parameter ( ctxt_ = 2, rsrc_ = 7, csrc_ = 8 )
      parameter ( m_ = 3, n_ = 4, mb_ = 5, nb_ = 6, lld_ = 9 , dlen_ = 9 )

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
               print '(A, I3, I3, A, I3, I3)', 'from process ',myrow, mycol, ' dim ',m, n
               !call flush(6)
               do i = 1,m
                  print editdesc, A(i,1:n)
               end do
               !call flush(6)
               !res = commitqq(6)
            end if
            call blacs_barrier( dA( ctxt_ ), 'A')
         end do
      end do
   end subroutine pmat
end program dsybc2_t



