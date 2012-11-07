!-------------------------------------
! subroutines for MPI matrix transpose
! Remark: in order to do the transpose a local reordering needs
! to be done before the send and a local transpose after the
! receive. This can be handled using MPI data types and the
! MPI-2 routine MPI_ALLTOALLW when this becomes available
!------------------------------------
!=========================================
!    
!    File:          transpose.f90
!    Project:       vlasov
!    Author(s):     Eric Sonnendrucker
!    Creation:      12.03.1999
!    Last modified: 18.03.1999
!    
!=========================================

subroutine transpose(a,at,nrow,ncol,nproc)
  !-----------------------------------------------------------------------
  !   performs a parallel transpose of a matrix A of total size (nrow,ncol)
  !   Both A and AT are split on the processors according to the columns
  !------------------------------------------------------------------------
  use used_precision
  use sll_collective
  implicit none
  
  real(wp), dimension(nrow,ncol/nproc) :: a   ! matrix
  real(wp), dimension(ncol,nrow/nproc) :: at  ! transpose matrix
  integer :: nrow, ncol  ! number of rows and columns of whole matrix
  integer :: nproc       ! number of processors


  !local variables
  integer :: ierr, icount, ipiecesize,jpiecesize

  integer :: comm

  comm   = sll_world_collective%comm

  ipiecesize = nrow/nproc
  jpiecesize = ncol/nproc
  icount=ipiecesize*jpiecesize

  ! test for wrong dimension
  if (nproc*ipiecesize.ne.nrow) then
     Write(*,*) 'Error in TRANSPOSE'
     write(*,*) 'number of rows on each processor must be the same'
     call MPI_FINALIZE(ierr)
  end if

  if (nproc*jpiecesize.ne.ncol) then
     write(*,*) 'Error in TRANSPOSE'
     write(*,*) 'number of columns on each processor must be the same'
     call MPI_FINALIZE(ierr)
  end if

  call reorder(a,at,ipiecesize,jpiecesize,nproc)

  call mpi_alltoall(at,icount,MPI_REAL8,a,icount,MPI_REAL8, &
       comm,ierr)

  if (mod(ipiecesize,64).eq.0) then
     call bloctransp(a,at,ipiecesize,ncol)
  else
     call loctransp(a,at,ipiecesize,ncol)
  end if
end subroutine transpose

subroutine reorder(f,ft,ipiecesize,jpiecesize,nproc)
  use used_precision
  implicit none
  integer :: ipiecesize,jpiecesize,nproc
  real(wp),dimension(ipiecesize*nproc,jpiecesize) ::f
  real(wp),dimension(ipiecesize,jpiecesize*nproc) ::ft
  ! local variables
  integer :: n,i,j

  do n=0,nproc-1
     do j=1,jpiecesize
        do i=1,ipiecesize
           ft(i,n*jpiecesize+j)=f(n*ipiecesize+i,j)
        end do
     end do
  end do
end subroutine reorder

subroutine loctransp(f,ft,sizex,sizey)
  use used_precision
  implicit none
  integer  :: sizex, sizey   ! local size of arrays
  real(wp),dimension(sizex,sizey) ::f
  real(wp),dimension(sizey,sizex) ::ft
  ! local variables	
  !integer :: i,j, i0,j0 , nblocs

  ft=transpose(f)
!!$  nblocs=2
!!$  i0=sizex/nblocs
!!$  j0=sizey/nblocs
!!$  do i=0,nblocs-1
!!$     do j=0,nblocs-1
!!$        ft(j*j0+1:(j+1)*j0,i*i0+1:(i+1)*i0) = &
!!$             transpose(f(i*i0+1:(i+1)*i0,j*j0+1:(j+1)*j0))
!!$     end do
!!$  end do
end subroutine loctransp

subroutine bloctransp(f,ft,sizex,sizey)
  use used_precision
  implicit none
  integer  :: sizex, sizey   ! local size of arrays
  real(wp),dimension(sizex*sizey) ::f
  real(wp),dimension(sizey*sizex) ::ft
  ! local variables	
  integer :: nblocx,nblocy

  integer :: i,j, ib,jb ,  indx
  integer :: ibloc,jbloc

  ibloc=min(64,sizex)
  jbloc=min(64,sizey)

  nblocx=sizex/ibloc
  nblocy=sizey/jbloc

  if ((ibloc*nblocx.ne.sizex).or.(ibloc*nblocx.ne.sizex)) then
     print*,'bloctranspose: verify size of blocs'
     stop
  end if

  indx = 0
  do jb=0,nblocy-1
     do ib=0,nblocx-1
        do j=0,jbloc-1
           do i=1,ibloc
              indx=indx+1
              ft(indx)=f(ib*ibloc+i+(jb*jbloc+j)*sizex)
           end do
        end do
     end do
  end do

  do jb=0,nblocy-1
     do ib=0,nblocx-1
        call loctransp(ft(ibloc*jbloc*ib+nblocx*ibloc*jbloc*jb+1:), &
             f(ibloc*jbloc*ib+nblocx*ibloc*jbloc*jb+1:),ibloc,jbloc)
     end do
  end do
  do jb=0,nblocy-1
     do ib=0,nblocx-1
        do i= 1,ibloc
           do j=1,jbloc
              ft(jb*jbloc+j+(ib*ibloc+i-1)*sizey)= &
                   f(ibloc*jbloc*ib+nblocx*ibloc*jbloc*jb +j +(i-1)*jbloc)
           end do
        end do
     end do
  end do

end subroutine bloctransp
 
recursive subroutine rectrans(f,ft,sizex,sizey)
  use used_precision
  implicit none
  integer  :: sizex, sizey   ! local size of arrays
  real(wp), dimension(sizex,sizey), intent(in)  :: f
  real(wp), dimension(sizey,sizex), intent(out) :: ft

  integer :: i0,j0

  i0=sizex/2
  j0=sizey/2
  if ((i0.le.512).and.(j0.le.512)) then
     ft(1:j0,1:i0)=transpose(f(1:i0,1:j0))
     ft(j0+1:sizey,1:i0)=transpose(f(1:i0,j0+1:sizey))
     ft(1:j0,i0+1:sizex)=transpose(f(i0+1:sizex,1:j0))
     ft(j0+1:sizey,i0+1:sizex)=transpose(f(i0+1:sizex,j0+1:sizey))
  else
     call rectrans(f(1:i0,1:j0),ft(1:j0,1:i0),i0,j0) 
     call rectrans(f(1:i0,j0+1:sizey),ft(j0+1:sizey,1:i0),i0,j0)
     call rectrans(f(i0+1:sizex,1:j0),ft(1:j0,i0+1:sizex),i0,j0)
     call rectrans(f(i0+1:sizex,j0+1:sizey),ft(j0+1:sizey,i0+1:sizex),i0,j0)
  end if
  
end subroutine rectrans
