!**************************************************************
!  Copyright INRIA
!  Authors : 
!     Pierre Navaro 
!  
!  This code SeLaLib (for Semi-Lagrangian-Library) 
!  is a parallel library for simulating the plasma turbulence 
!  in a tokamak.
!  
!  This software is governed by the CeCILL-B license 
!  under French law and abiding by the rules of distribution 
!  of free software.  You can  use, modify and redistribute 
!  the software under the terms of the CeCILL-B license as 
!  circulated by CEA, CNRS and INRIA at the following URL
!  "http://www.cecill.info". 
!**************************************************************

!> \author
!> Pierre Navaro 
!> Common data for Maxwell solvers
module sll_maxwell_solvers

#include "sll_working_precision.h"
#include "sll_memory.h"
#include "sll_assert.h"
#include "sll_constants.h"
#include "sll_maxwell_solvers_macros.h"

implicit none

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

!Object with data to solve Maxwell equation on 2d domain
!Maxwell in TE mode: (Ex,Ey,Hz)
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

!> write files to visualize 2d fields with gnuplot
subroutine plot_fields(fname, f1, f2, iplot, time )

sll_real64, dimension(:,:), intent(in) :: f1 !< first field 2d
sll_real64, dimension(:,:), intent(in) :: f2 !< second field 2d
integer :: iplot !< plot counter
integer :: i, j
sll_real64, intent(in) :: time !< time
character(len=*) :: fname !< output file name
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

end module sll_maxwell_solvers
