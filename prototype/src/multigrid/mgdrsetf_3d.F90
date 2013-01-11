subroutine mgdrsetf(sxf,exf,syf,eyf,szf,ezf,rf,r,IOUT)

use mpi
implicit none 
integer :: sxf,exf,syf,eyf,szf,ezf,IOUT
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

return
end
