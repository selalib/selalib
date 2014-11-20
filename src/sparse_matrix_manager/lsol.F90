!----------------------------------------------------------------------c
! 2)     T R I A N G U L A R    S Y S T E M    S O L U T I O N S       c
!----------------------------------------------------------------------c
      subroutine lsol (n,x,y,al,jal,ial)
      integer n, jal(*),ial(n+1)
      real(8)  x(n), y(n), al(*)
!-----------------------------------------------------------------------
!   solves    L x = y ; L = lower unit triang. /  CSR format
!-----------------------------------------------------------------------
! solves a unit lower triangular system by standard (sequential )
! forward elimination - matrix stored in CSR format.
!-----------------------------------------------------------------------
!
! On entry:
!----------
! n      = integer. dimension of problem.
! y      = real array containg the right side.
!
! al,
! jal,
! ial,    = Lower triangular matrix stored in compressed sparse row
!          format.
!
! On return:
!-----------
!      x  = The solution of  L x  = y.
!--------------------------------------------------------------------
! local variables
!
      integer k, j
      real(8)  t
!-----------------------------------------------------------------------
      x(1) = y(1)
      do 150 k = 2, n
         t = y(k)
         do 100 j = ial(k), ial(k+1)-1
            t = t-al(j)*x(jal(j))
 100     continue
         x(k) = t
 150  continue
!
      return
!----------end-of-lsol--------------------------------------------------
!-----------------------------------------------------------------------
      end subroutine lsol
