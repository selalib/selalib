subroutine mgdsetf(sxf,exf,syf,eyf,szf,ezf,phi,rhs,phif,rhsf,IOUT)

use mpi
implicit none 
integer :: sxf,exf,syf,eyf,szf,ezf,IOUT
real(8) :: phi(sxf-1:exf+1,syf-1:eyf+1,szf-1:ezf+1)
real(8) :: rhs(sxf-1:exf+1,syf-1:eyf+1,szf-1:ezf+1)
real(8) :: phif(sxf-1:exf+1,syf-1:eyf+1,szf-1:ezf+1)
real(8) :: rhsf(sxf-1:exf+1,syf-1:eyf+1,szf-1:ezf+1)
!------------------------------------------------------------------------
! Set the fine grid values in the work vector
!
! Code      : mgd3, 3-D parallel multigrid solver
! Author    : Bernard Bunner (bunner@engin.umich.edu), January 1998
! Called in : mgdsolver
! Calls     : --
!------------------------------------------------------------------------
integer :: i,j,k

do k=szf-1,ezf+1
  do j=syf-1,eyf+1
    do i=sxf-1,exf+1
      phi(i,j,k)=phif(i,j,k)
      rhs(i,j,k)=rhsf(i,j,k)
    end do
  end do
end do

return
end
      subroutine mgdsetf(sxf,exf,syf,eyf,phi,rhs,phif,rhsf,IOUT)
# include "compdir.inc"
      include "mpif.h"
      integer sxf,exf,syf,eyf,IOUT
      REALN phi(sxf-1:exf+1,syf-1:eyf+1),rhs(sxf-1:exf+1,syf-1:eyf+1)
      REALN phif(sxf-1:exf+1,syf-1:eyf+1),rhsf(sxf-1:exf+1,syf-1:eyf+1)
c------------------------------------------------------------------------
c Set the fine grid values in the work vector
c
c Code      : mgd2, 2-D parallel multigrid solver
c Author    : Bernard Bunner (bunner@engin.umich.edu), January 1998
c Called in : mgdsolver
c Calls     : --
c------------------------------------------------------------------------
      integer i,j
# if cdebug
      double precision tinitial
      tinitial=MPI_WTIME()
# endif
c
      do j=syf-1,eyf+1
        do i=sxf-1,exf+1
          phi(i,j)=phif(i,j)
          rhs(i,j)=rhsf(i,j)
        end do
      end do
c
# if cdebug
      timing(88)=timing(88)+MPI_WTIME()-tinitial
# endif
      return
      end
