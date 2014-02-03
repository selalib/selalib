!-----------------------------------------------------------------------
      subroutine getdia (nrow,ncol,job,a,ja,ia,len,diag,idiag,ioff)
      real*8 diag(*),a(*)
      integer nrow, ncol, job, len, ioff, ia(*), ja(*), idiag(*)
!-----------------------------------------------------------------------
! this subroutine extracts a given diagonal from a matrix stored in csr
! format. the output matrix may be transformed with the diagonal removed
! from it if desired (as indicated by job.)
!-----------------------------------------------------------------------
! our definition of a diagonal of matrix is a vector of length nrow
! (always) which contains the elements in rows 1 to nrow of
! the matrix that are contained in the diagonal offset by ioff
! with respect to the main diagonal. if the diagonal element
! falls outside the matrix then it is defined as a zero entry.
! thus the proper definition of diag(*) with offset ioff is
!
!     diag(i) = a(i,ioff+i) i=1,2,...,nrow
!     with elements falling outside the matrix being defined as zero.
!
!-----------------------------------------------------------------------
!
! on entry:
!----------
!
! nrow      = integer. the row dimension of the matrix a.
! ncol      = integer. the column dimension of the matrix a.
! job   = integer. job indicator.  if job = 0 then
!         the matrix a, ja, ia, is not altered on return.
!         if job.ne.0  then getdia will remove the entries
!         collected in diag from the original matrix.
!         this is done in place.
!
! a,ja,
!    ia = matrix stored in compressed sparse row a,ja,ia,format
! ioff  = integer,containing the offset of the wanted diagonal
!        the diagonal extracted is the one corresponding to the
!        entries a(i,j) with j-i = ioff.
!        thus ioff = 0 means the main diagonal
!
! on return:
!-----------
! len   = number of nonzero elements found in diag.
!         (len .le. min(nrow,ncol-ioff)-max(1,1-ioff) + 1 )
!
! diag  = real*8 array of length nrow containing the wanted diagonal.
!        diag contains the diagonal (a(i,j),j-i = ioff ) as defined
!         above.
!
! idiag = integer array of  length len, containing the poisitions
!         in the original arrays a and ja of the diagonal elements
!         collected in diag. a zero entry in idiag(i) means that
!         there was no entry found in row i belonging to the diagonal.
!
! a, ja,
!    ia = if job .ne. 0 the matrix is unchanged. otherwise the nonzero
!         diagonal entries collected in diag are removed from the
!         matrix and therefore the arrays a, ja, ia will change.
!        (the matrix a, ja, ia will contain len fewer elements)
!
!----------------------------------------------------------------------c
!     Y. Saad, sep. 21 1989 - modified and retested Feb 17, 1996.      c
!----------------------------------------------------------------------c
!     local variables
      integer istart, max, iend, i, kold, k, kdiag, ko
!
      istart = max(0,-ioff)
      iend = min(nrow,ncol-ioff)
      len = 0
      do 1 i=1,nrow
         idiag(i) = 0
       diag(i) = 0.0d0
 1    continue
!
!     extract  diagonal elements
!
      do 6 i=istart+1, iend
         do 51 k= ia(i),ia(i+1) -1
            if (ja(k)-i .eq. ioff) then
               diag(i)= a(k)
               idiag(i) = k
               len = len+1
               goto 6
            endif
 51      continue
 6    continue
      if (job .eq. 0 .or. len .eq.0) return
!
!     remove diagonal elements and rewind structure
!
      ko = 0
      do  7 i=1, nrow
         kold = ko
         kdiag = idiag(i)
         do 71 k= ia(i), ia(i+1)-1
            if (k .ne. kdiag) then
               ko = ko+1
               a(ko) = a(k)
               ja(ko) = ja(k)
            endif
 71      continue
         ia(i) = kold+1
 7    continue
!
!     redefine ia(nrow+1)
!
      ia(nrow+1) = ko+1
      return
!------------end-of-getdia----------------------------------------------
!-----------------------------------------------------------------------
      end subroutine getdia 
