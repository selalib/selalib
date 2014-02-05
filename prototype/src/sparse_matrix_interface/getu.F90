!-----------------------------------------------------------------------
      subroutine getu (n,a,ja,ia,ao,jao,iao)
      integer n, ia(*), ja(*), iao(*), jao(*)
      real(8) a(*), ao(*)
!------------------------------------------------------------------------
! this subroutine extracts the upper triangular part of a matrix
! and writes the result ao, jao, iao. The routine is in place in
! that ao, jao, iao can be the same as a, ja, ia if desired.
!-----------
! on input:
!
! n     = dimension of the matrix a.
! a, ja,
!    ia = matrix stored in a, ja, ia, format
! On return:
! ao, jao,
!    iao = upper triangular matrix (upper part of a)
!      stored in compressed sparse row format
! note: the diagonal element is the last element in each row.
! i.e. in  a(ia(i+1)-1 )
! ao, jao, iao may be the same as a, ja, ia on entry -- in which case
! getu will overwrite the result on a, ja, ia.
!
!------------------------------------------------------------------------
! local variables
      real(8) t
      integer ko, k, i, kdiag, kfirst
      ko = 0
      do  7 i=1, n
         kfirst = ko+1
         kdiag = 0
         do 71 k = ia(i), ia(i+1) -1
            if (ja(k)  .lt. i) goto 71
            ko = ko+1
            ao(ko) = a(k)
            jao(ko) = ja(k)
            if (ja(k)  .eq. i) kdiag = ko
 71      continue
         if (kdiag .eq. 0 .or. kdiag .eq. kfirst) goto 72
!     exchange
         t = ao(kdiag)
         ao(kdiag) = ao(kfirst)
         ao(kfirst) = t
!
         k = jao(kdiag)
         jao(kdiag) = jao(kfirst)
         jao(kfirst) = k
 72      iao(i) = kfirst
 7    continue
!     redefine iao(n+1)
      iao(n+1) = ko+1
      return
!----------end-of-getu -------------------------------------------------
!-----------------------------------------------------------------------
      end subroutine getu 
