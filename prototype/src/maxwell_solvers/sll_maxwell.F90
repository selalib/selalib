!>
!>@namespace sll_maxwell_2d
!>
!> @author
!> Pierre Navaro Philippe Helluy
!>
!
! DESCRIPTION: 
!
!> @brief
!> Implements the Maxwell solver in 2D
!>
!>@details
!>This module depends on:
!> - memory
!> - precision
!> - assert 
!> - numerical_utilities
!> - constants
!> - sll_utilities
!>
! REVISION HISTORY:
! 03 02 2012 - Initial Version  (fevrier)
! TODO_dd_mmm_yyyy - TODO_describe_appropriate_changes - TODO_name
!------------------------------------------------------------------------------

module sll_maxwell

#include "sll_working_precision.h"
#include "sll_memory.h"
#include "sll_assert.h"

use numeric_constants

implicit none

enum, bind(C)
   enumerator :: TE_POLARIZATION = 0, TM_POLARIZATION = 1
end enum

enum, bind(C)
   enumerator :: NORTH = 0, EAST = 1, SOUTH = 2, WEST = 3
end enum

!private
!
!interface new
! module procedure new_maxwell_2d
!end interface
!interface solve
! module procedure solve_maxwell_2d
!end interface
!interface delete
! module procedure delete_maxwell_2d
!end interface
!
!public :: new, solve, delete

!> Object with data to solve Maxwell equation on 2d domain
!> Maxwell in TE mode: (Ex,Ey,Hz)
!type, public :: maxwell_2d
!  sll_real64 :: c_light
!  sll_real64 :: epsilon_0
!  sll_int32  :: ix, jx, iy, jy
!  sll_real64 :: dx, dy
!end type maxwell_2d
!
!enum, bind(C)
!   enumerator :: NORTH = 0, EAST = 1, SOUTH = 2, WEST = 3
!end enum
!
!enum, bind(C)
!   enumerator :: FDTD = 0, PSTD = 1
!end enum

contains

!subroutine new_maxwell_2d(this, ix, jx, iy, jy, dx, dy, METHOD )
!
!   type(maxwell_2d) :: this
!   sll_int32        :: ix, jx, iy, jy
!   sll_real64       :: dx, dy
!   sll_int32        :: error
!   sll_int32        :: METHOD
!
!   select case(METHOD)
!   case(FDTD)
!      call new_maxwell_2d_fdtd
!   end select
!
!end subroutine new_maxwell_2d
!
!subroutine solve_maxwell_2d(this, ex, ey, bz, dt)
!
!   type(maxwell_2d)          :: this
!   sll_real64 , intent(inout), dimension(:,:)   :: ex, ey, bz
!   sll_real64 , intent(in)   :: dt
!
!   !B(n-1/2)--> B(n+1/2) sur les pts interieurs   
!   call faraday(this, ex, ey, bz, dt)   
!
!   call cl_periodiques(this, ex, ey, bz, dt)
!
!   !E(n)-->E(n+1) sur les pts interieurs
!   call ampere_maxwell(this, ex, ey, bz, dt) 
!
!end subroutine solve_maxwell_2d
!
!subroutine delete_maxwell_2d(this)
!   type(maxwell_2d), pointer :: this
!  
!end subroutine delete_maxwell_2d

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine plot_fields(fname, f1, f2, iplot, time )

sll_real64, dimension(:,:), intent(in) :: f1, f2
integer :: iplot, i, j
sll_real64, intent(in) :: time
character(len=*) :: fname
character(len=4) :: cplot
sll_int32 :: nx, ny

nx = size(f1,1)
ny = size(f1,2)

SLL_ASSERT(nx == size(f2,1))
SLL_ASSERT(ny == size(f2,2))

call int2string(iplot, cplot)

!write domains
open( 80, file = fname//cplot//".dat" )
   do i=1,nx
      do j=1,ny
         write(80,*) i, j, sngl(f1(i,j)), sngl(f2(i,j))
      end do
      write(80,*) 
   end do
close(80)
   
open( 90, file = fname//'plots.gnu', position="append" )
  if ( iplot == 1 ) then
     rewind(90)
     !write(90,*)"set xr[-0.1:1.1]"
     !write(90,*)"set yr[-0.1:1.1]"
     !write(90,*)"set zr[-1.1:1.1]"
     !write(90,*)"set cbrange[-1:1]"
     !write(90,*)"set pm3d"
     write(90,*)"set surf"
     write(90,*)"set term x11"
  end if
  write(90,*)"set title 'Time = ",time,"'"
  write(90,"(a)",advance='no')"splot '"//fname//cplot//".dat' w lines"
  write(90,"(a)",advance='no')",'"//fname//cplot//".dat' u 1:2:4 w lines"
  close(90)

end subroutine plot_fields

end module sll_maxwell
