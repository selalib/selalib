!------------------------------------------------------------------------
      subroutine getl (n,a,ja,ia,ao,jao,iao)
      integer n, ia(*), ja(*), iao(*), jao(*)
      real*8 a(*), ao(*)
!------------------------------------------------------------------------
! this subroutine extracts the lower triangular part of a matrix
! and writes the result ao, jao, iao. The routine is in place in
! that ao, jao, iao can be the same as a, ja, ia if desired.
!-----------
! on input:
!
! n     = dimension of the matrix a.
! a, ja,
!    ia = matrix stored in compressed sparse row format.
! On return:
! ao, jao,
!    iao = lower triangular matrix (lower part of a)
!      stored in a, ja, ia, format
! note: the diagonal element is the last element in each row.
! i.e. in  a(ia(i+1)-1 )
! ao, jao, iao may be the same as a, ja, ia on entry -- in which case
! getl will overwrite the result on a, ja, ia.
!
!------------------------------------------------------------------------
! local variables
      real*8 t
      integer ko, kold, kdiag, k, i
!
! inititialize ko (pointer for output matrix)
!
      ko = 0
      do  7 i=1, n
         kold = ko
         kdiag = 0
         do 71 k = ia(i), ia(i+1) -1
            if (ja(k)  .gt. i) goto 71
            ko = ko+1
            ao(ko) = a(k)
            jao(ko) = ja(k)
            if (ja(k)  .eq. i) kdiag = ko
 71      continue
         if (kdiag .eq. 0 .or. kdiag .eq. ko) goto 72
!
!     exchange
!
         t = ao(kdiag)
         ao(kdiag) = ao(ko)
         ao(ko) = t
!
         k = jao(kdiag)
         jao(kdiag) = jao(ko)
         jao(ko) = k
 72      iao(i) = kold+1
 7    continue
!     redefine iao(n+1)
      iao(n+1) = ko+1
      return
!----------end-of-getl -------------------------------------------------
      end subroutine getl 
