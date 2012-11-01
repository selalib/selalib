!>
!>@namespace sll_maxwell_fdtd_2d
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

module sll_maxwell_2d

#include "sll_working_precision.h"
#include "sll_memory.h"
#include "sll_assert.h"

use numeric_constants
use sll_maxwell_fdtd_2d

implicit none
private

interface new
 module procedure new_maxwell_2d
end interface
interface solve
 module procedure solve_maxwell_2d
end interface
interface delete
 module procedure delete_maxwell_2d
end interface

public :: new, solve, delete

!> Object with data to solve Maxwell equation on 2d domain
!> Maxwell in TE mode: (Ex,Ey,Hz)
type, public :: maxwell_2d
  sll_real64 :: c_light
  sll_real64 :: epsilon_0
  sll_int32  :: ix, jx, iy, jy
  sll_real64 :: dx, dy
end type maxwell_2d

enum, bind(C)
   enumerator :: NORTH = 0, EAST = 1, SOUTH = 2, WEST = 3
end enum

enum, bind(C)
   enumerator :: FDTD = 0, PSTD = 1
end enum

contains

subroutine new_maxwell_2d(this, ix, jx, iy, jy, dx, dy, METHOD )

   type(maxwell_2d) :: this
   sll_int32        :: ix, jx, iy, jy
   sll_real64       :: dx, dy
   sll_int32        :: error
   sll_int32        :: METHOD

   select case(METHOD)
   case(FDTD)
      call new_maxwell_2d_fdtd

end subroutine new_maxwell_2d

subroutine solve_maxwell_2d(this, ex, ey, bz, dt)

   type(maxwell_2d)          :: this
   sll_real64 , intent(inout), dimension(:,:)   :: ex, ey, bz
   sll_real64 , intent(in)   :: dt

   !B(n-1/2)--> B(n+1/2) sur les pts interieurs   
   call faraday(this, ex, ey, bz, dt)   

   call cl_periodiques(this, ex, ey, bz, dt)

   !E(n)-->E(n+1) sur les pts interieurs
   call ampere_maxwell(this, ex, ey, bz, dt) 

end subroutine solve_maxwell_2d

subroutine delete_maxwell_2d(this)
   type(maxwell_2d), pointer :: this
  
end subroutine delete_maxwell_2d

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine faraday( this, ex, ey, bz, dt )

type(maxwell_2d) :: this
sll_real64, dimension(:,:), intent(in)    :: ex, ey
sll_real64, dimension(:,:), intent(inout) :: bz
sll_real64, intent(in) :: dt
sll_int32 :: i, j
sll_real64 :: dex_dy, dey_dx
sll_real64 :: dx, dy
sll_int32  :: ix, jx, iy, jy

csq = this%c_light * this%c_light
dx  = this%dx
dy  = this%dy

!*** On utilise l'equation de Faraday sur un demi pas
!*** de temps pour le calcul du champ magnetique  Bz 
!*** a l'instant n puis n+1/2 

ix = this%ix
jx = this%jx
iy = this%iy
jy = this%jy

do i=ix,jx-1
do j=iy,jy-1
   dex_dy  = (ex(i,j+1)-ex(i,j)) / dy
   dey_dx  = (ey(i+1,j)-ey(i,j)) / dx
   bz(i,j) = bz(i,j) + dt * (dex_dy - dey_dx)
end do
end do

end subroutine faraday

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine ampere_maxwell( this, ex, ey, bz, dt )

type(maxwell_2d) :: this
sll_int32 :: ix, jx, iy, jy
sll_real64, dimension(:,:), intent(inout) :: ex, ey
sll_real64, dimension(:,:), intent(in)    :: bz
sll_int32 :: i, j
sll_real64 :: dex_dx, dey_dy
sll_real64 :: dex_dy, dey_dx
sll_real64 :: dbz_dx, dbz_dy
sll_real64, intent(in) :: dt
sll_real64 :: dx, dy

ix = this%ix
jx = this%jx
iy = this%iy
jy = this%jy

csq = this%c_light * this%c_light
dx  = this%dx
dy  = this%dy

!*** Calcul du champ electrique E au temps n+1
!*** sur les points internes du maillage
!*** Ex aux points (i+1/2,j)
!*** Ey aux points (i,j+1/2)

do i=ix,jx
do j=iy+1,jy
   dbz_dy  = (bz(i,j)-bz(i,j-1)) / dy
   ex(i,j) = ex(i,j) + dt*csq*dbz_dy 
end do
end do

do i=ix+1,jx
do j=iy,jy
   dbz_dx  = (bz(i,j)-bz(i-1,j)) / dx
   ey(i,j) = ey(i,j) - dt*csq*dbz_dx 
end do
end do

end subroutine ampere_maxwell

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine cl_periodiques(this, ex, ey, bz, dt)

type(maxwell_2d) :: this
sll_int32 :: ix, jx, iy, jy
sll_real64, dimension(:,:), intent(inout) :: ex, ey, bz
sll_real64 :: dex_dx, dey_dy
sll_real64 :: dex_dy, dey_dx
sll_real64 :: dbz_dx, dbz_dy
sll_real64, intent(in) :: dt
sll_real64 :: dx, dy
sll_int32 :: i, j

ix = this%ix
jx = this%jx
iy = this%iy
jy = this%jy

csq = this%c_light * this%c_light
dx  = this%dx
dy  = this%dy

do i = ix, jx-1
   bz(i,jy) = bz(i,iy)
end do1G1G
do j = iy, jy-1
   bz(jx,j) = bz(ix,j)
end do

bz(jx,jy) = bz(ix,iy)

do i = ix, jx
   dbz_dy = (bz(i,iy)-bz(i,jy-1)) / dy
   ex(i,iy) = ex(i,jy) + dt*csq*dbz_dy 
end do
     
do j = iy, jy
   dbz_dx = (bz(ix,j)-bz(jx-1,j)) / dx
   ey(ix,j) = ey(jx,j) - dt*csq*dbz_dx 
end do

end subroutine cl_periodiques

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine cl_condparfait(this, ex, ey, bz, t log side)

type(maxwell_2d) :: this
sll_int32, intent(in) :: side
sll_real64, dimension(:,:), intent(inout) :: ex, ey, bz
sll_int32 :: ix, jx, iy, jy
sll_int32 :: i, j

ix = this%ix
jx = this%jx
iy = this%iy
jy = this%jy

select case(side)
case(SOUTH)
   do i = ix, jx
      ex(i,iy) = 0.d0
   end do
case(NORTH)
   do i = ix, jx
      ex(i,jy) = 0.d0
      bz(i,jy) = bz(i,jy-1)
   end do
case(WEST)
   do j = iy, jy
      ey(ix,j) = 0.d0
   end do
case(EAST)
   do j = iy, jy
      ey(jx,j) = 0.d0
      bz(jx,j) = bz(jx-1,j)
   end do
end select

end subroutine cl_condparfait


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine silver_muller( this, ex, ey, bz, ccall, dt )

type(maxwell_2d) :: this
sll_int32, intent(in) :: ccall
sll_int32 :: ix, jx, iy, jy
sll_real64 :: a11,a12,a21,a22,b1,b2,dis
sll_int32 :: i, j
sll_real64, intent(in) :: dt
sll_real64 :: dx, dy, c, csq
sll_real64, dimension(:,:), pointer :: ex, ey, bz

c   = this%c_light 
csq = c * c
dx  = this%dx
dy  = this%dy

ix = this%ix
jx = this%jx
iy = this%iy
jy = this%jy

!Conditions de Silver-Muller
!------------------------------------
!Ey = -c Bz sur la frontiere ouest
!Ey =  c Bz sur la frontiere est
!Ex = -c Bz sur la frontiere nord
!Ex =  c Bz sur la frontiere sud
   
!On effectue le calcul de B sur les points fictifs du maillage
!simultanement avec la prise en compte des conditions limites sur
!E. Pour cela, on resout sur les points frontieres, l'equation de la
!condition limite en moyennant en temps pour E et en espace pour B puis
!l'equation d'Ampere

select case (ccall)

case (NORTH)
   !Frontiere Nord : Ex = -c Bz 
   do i = ix, jx
         
      a11 = 1.; a12 = + c
      a21 = 1./dt; a22 = - csq / dy
      b1  = - ex(i,jy) - c * bz(i,jy-1)
      b2  =   ex(i,jy)/dt - csq/dy*bz(i,jy-1)
         
      dis = a11*a22-a21*a12 
         
      !ex(i,jy) = (b1*a22-b2*a12)/dis
      bz(i,jy) = (a11*b2-a21*b1)/dis
         
   end do
      
case (SOUTH)

   !Frontiere Sud : Ex =  c Bz
   do i = ix, jx
         
      a11 = 1.; a12 = - c
      a21 = 1./dt; a22 = csq / dy
      b1  = - ex(i,iy) + c * bz(i,iy+1)
      b2  = ex(i,iy)/dt + csq / dy * bz(i,iy+1) 
         
      dis = a11*a22-a21*a12 
         
      ex(i,iy) = (b1*a22-b2*a12)/dis
      !bz(i,iy) = (a11*b2-a21*b1)/dis
         
   end do
      
case (EAST)

   !Frontiere Est : Ey =  c Bz
   do j = iy, jy
         
      a11 = 1.; a12 = - c
      a21 = 1./dt; a22 = + csq / dx
      b1  = - ey(jx,j) + c * bz(jx-1,j)
      b2  = ey(jx,j)/dt + csq/dx*bz(jx-1,j) 
         
      dis = a11*a22-a21*a12 
         
      !ey(jx,j) = (b1*a22-b2*a12)/dis
      bz(jx,j) = (a11*b2-a21*b1)/dis
      
   end do
      
case (WEST)

   !Frontiere Ouest : Ey = -c Bz
   do j = iy, jy
      
      a11 = 1.; a12 = + c
      a21 = 1./dt; a22 = - csq / dx
      b1  = - ey(ix,j) - c * bz(ix+1,j)
      b2  =   ey(ix,j)/dt - csq/dx*bz(ix+1,j) 
      
      dis = a11*a22-a21*a12 
   
      ey(ix,j) = (b1*a22-b2*a12)/dis
      !bz(ix,j) = (a11*b2-a21*b1)/dis
      
   end do

end select

end subroutine silver_muller

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module sll_maxwell_2d
