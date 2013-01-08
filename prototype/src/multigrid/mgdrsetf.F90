module mgdrsetf
#include "sll_working_precision.h"
use mpi
implicit none 

contains

subroutine mgdrsetf_3d(sxf,exf,syf,eyf,szf,ezf,rf,r)

integer :: sxf,exf,syf,eyf,szf,ezf
real(8) :: rf(sxf-1:exf+1,syf-1:eyf+1,szf-1:ezf+1)
real(8) :: r(sxf-1:exf+1,syf-1:eyf+1,szf-1:ezf+1)
!------------------------------------------------------------------------
! Set the fine grid values of the density in the work vector
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
      rf(i,j,k)=r(i,j,k)
    end do
  end do
end do

end subroutine

subroutine mgdrsetf_2d(sxf,exf,syf,eyf,rf,r)
#include "mgd2.h"
integer sxf,exf,syf,eyf
REALN rf(sxf-1:exf+1,syf-1:eyf+1),r(sxf-1:exf+1,syf-1:eyf+1)
!------------------------------------------------------------------------
! For the old version of the multigrid code, set the fine grid values 
! of the density in the work vector
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
    rf(i,j)=r(i,j)
  end do
end do

# if cdebug
timing(85)=timing(85)+MPI_WTIME()-tinitial
# endif

end subroutine
end module
