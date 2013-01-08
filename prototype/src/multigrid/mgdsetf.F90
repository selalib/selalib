module mgdsetf
#include "sll_working_precision.h"
use mpi
implicit none 

contains

subroutine mgdsetf_3d(sxf,exf,syf,eyf,szf,ezf,phi,rhs,phif,rhsf)

integer :: sxf,exf,syf,eyf,szf,ezf
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

end subroutine


subroutine mgdsetf_2d(sxf,exf,syf,eyf,phi,rhs,phif,rhsf)
#include "mgd2.h"
integer sxf,exf,syf,eyf
REALN phi(sxf-1:exf+1,syf-1:eyf+1),rhs(sxf-1:exf+1,syf-1:eyf+1)
REALN phif(sxf-1:exf+1,syf-1:eyf+1),rhsf(sxf-1:exf+1,syf-1:eyf+1)
!------------------------------------------------------------------------
! Set the fine grid values in the work vector
!
! Code      : mgd2, 2-D parallel multigrid solver
! Author    : Bernard Bunner (bunner@engin.umich.edu), January 1998
! Called in : mgdsolver
! Calls     : --
!------------------------------------------------------------------------
integer i,j
# if cdebug
double precision tinitial
tinitial=MPI_WTIME()
# endif

do j=syf-1,eyf+1
  do i=sxf-1,exf+1
    phi(i,j)=phif(i,j)
    rhs(i,j)=rhsf(i,j)
  end do
end do

# if cdebug
timing(88)=timing(88)+MPI_WTIME()-tinitial
# endif

end subroutine
end module
