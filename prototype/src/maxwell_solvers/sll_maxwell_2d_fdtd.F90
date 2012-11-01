!>
!>@namespace sll_maxwell_2d_fdtd
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

module sll_maxwell_2d_fdtd

#include "sll_working_precision.h"
#include "sll_memory.h"
#include "sll_assert.h"

use sll_maxwell
use numeric_constants

implicit none
private

interface initialize
 module procedure new_maxwell_2d_fdtd
end interface
interface solve
 module procedure solve_maxwell_2d_fdtd
end interface
interface free
 module procedure delete_maxwell_2d_fdtd
end interface

public :: initialize, solve, free

!> Object with data to solve Maxwell equation on 2d domain
!> Maxwell in TE mode: (Ex,Ey,Hz)
type, public :: maxwell_fdtd
  sll_real64 :: c
  sll_real64 :: e_0
  sll_int32  :: i1, j1, i2, j2
  sll_real64 :: dx, dy
end type maxwell_fdtd

contains

subroutine new_maxwell_2d_fdtd(this, i1, j1, i2, j2, dx, dy  )

   type(maxwell_fdtd) :: this
   sll_int32        :: i1, j1, i2, j2
   sll_real64       :: dx, dy
   !sll_int32        :: error

   this%c   = 1.0_f64
   this%e_0 = 1.0_f64
   
   this%dx = dx
   this%dy = dy
   this%i1 = i1
   this%j1 = j1
   this%i2 = i2
   this%j2 = j2

end subroutine new_maxwell_2d_fdtd

subroutine solve_maxwell_2d_fdtd(this, ex, ey, bz, dt)

   type(maxwell_fdtd)          :: this
   sll_real64 , intent(inout), dimension(:,:)   :: ex, ey, bz
   sll_real64 , intent(in)   :: dt

   !B(n-1/2)--> B(n+1/2) sur les pts interieurs   
   call faraday(this, ex, ey, bz, dt)   

   call cl_periodiques(this, ex, ey, bz, dt)

   !E(n)-->E(n+1) sur les pts interieurs
   call ampere_maxwell(this, ex, ey, bz, dt) 

end subroutine solve_maxwell_2d_fdtd

subroutine delete_maxwell_2d_fdtd( this)
   type(maxwell_fdtd) :: this
  
end subroutine delete_maxwell_2d_fdtd

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine faraday( this, ex, ey, bz, dt )

type(maxwell_fdtd) :: this
sll_real64, dimension(:,:), intent(in)    :: ex, ey
sll_real64, dimension(:,:), intent(inout) :: bz
sll_real64, intent(in) :: dt
sll_int32 :: i, j
sll_real64 :: dex_dy, dey_dx
sll_real64 :: dx, dy
sll_int32  :: i1, j1, i2, j2

dx  = this%dx
dy  = this%dy

!*** On utilise l'equation de Faraday sur un demi pas
!*** de temps pour le calcul du champ magnetique  Bz 
!*** a l'instant n puis n+1/2 

i1 = this%i1
j1 = this%j1
i2 = this%i2
j2 = this%j2

do i=i1,j1-1
do j=i2,j2-1
   dex_dy  = (ex(i,j+1)-ex(i,j)) / dy
   dey_dx  = (ey(i+1,j)-ey(i,j)) / dx
   bz(i,j) = bz(i,j) + dt * (dex_dy - dey_dx)
end do
end do

end subroutine faraday

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine ampere_maxwell( this, ex, ey, bz, dt, jx, jy )

type(maxwell_fdtd) :: this
sll_int32 :: i1, j1, i2, j2
sll_real64, dimension(:,:), intent(inout) :: ex, ey
sll_real64, dimension(:,:), intent(in)    :: bz
sll_real64, dimension(:,:), intent(in), optional :: jx, jy
sll_int32 :: i, j
sll_real64 :: dbz_dx, dbz_dy
sll_real64, intent(in) :: dt
sll_real64 :: dx, dy
sll_real64 :: csq

i1 = this%i1
j1 = this%j1
i2 = this%i2
j2 = this%j2

csq = this%c * this%c
dx  = this%dx
dy  = this%dy

!*** Calcul du champ electrique E au temps n+1
!*** sur les points internes du maillage
!*** Ex aux points (i+1/2,j)
!*** Ey aux points (i,j+1/2)

do i=i1,j1
do j=i2+1,j2
   dbz_dy  = (bz(i,j)-bz(i,j-1)) / dy
   ex(i,j) = ex(i,j) + dt*csq*dbz_dy 
end do
end do

do i=i1+1,j1
do j=i2,j2
   dbz_dx  = (bz(i,j)-bz(i-1,j)) / dx
   ey(i,j) = ey(i,j) - dt*csq*dbz_dx 
end do
end do

if (present(jx) .and. present(jy)) then

   ex(i1:j1,i2+1:j2) = ex(i1:j1,i2+1:j2) - dt * jx(i1:j1,i2+1:j2) / this%e_0
   ey(i1+1:j1,i2:j2) = ey(i1+1:j1,i2:j2) - dt * jy(i1+1:j1,i2:j2) / this%e_0

endif

end subroutine ampere_maxwell

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine cl_periodiques(this, ex, ey, bz, dt)

type(maxwell_fdtd) :: this
sll_int32 :: i1, j1, i2, j2
sll_real64, dimension(:,:), intent(inout) :: ex, ey, bz
sll_real64 :: dbz_dx, dbz_dy
sll_real64, intent(in) :: dt
sll_real64 :: dx, dy
sll_int32 :: i, j
sll_real64 :: csq

i1 = this%i1
j1 = this%j1
i2 = this%i2
j2 = this%j2

csq = this%c * this%c
dx  = this%dx
dy  = this%dy

do i = i1, j1-1
   bz(i,j2) = bz(i,i2)
end do
do j = i2, j2-1
   bz(j1,j) = bz(i1,j)
end do

bz(j1,j2) = bz(i1,i2)

do i = i1, j1
   dbz_dy = (bz(i,i2)-bz(i,j2-1)) / dy
   ex(i,i2) = ex(i,j2) + dt*csq*dbz_dy 
end do
     
do j = i2, j2
   dbz_dx = (bz(i1,j)-bz(j1-1,j)) / dx
   ey(i1,j) = ey(j1,j) - dt*csq*dbz_dx 
end do

end subroutine cl_periodiques

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine cl_condparfait(this, ex, ey, bz, side)

type(maxwell_fdtd) :: this
sll_int32, intent(in) :: side
sll_real64, dimension(:,:), intent(inout) :: ex, ey, bz
sll_int32 :: i1, j1, i2, j2
sll_int32 :: i, j

i1 = this%i1
j1 = this%j1
i2 = this%i2
j2 = this%j2

select case(side)
case(SOUTH)
   do i = i1, j1
      ex(i,i2) = 0.d0
   end do
case(NORTH)
   do i = i1, j1
      ex(i,j2) = 0.d0
      bz(i,j2) = bz(i,j2-1)
   end do
case(WEST)
   do j = i2, j2
      ey(i1,j) = 0.d0
   end do
case(EAST)
   do j = i2, j2
      ey(j1,j) = 0.d0
      bz(j1,j) = bz(j1-1,j)
   end do
end select

end subroutine cl_condparfait


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine silver_muller( this, ex, ey, bz, ccall, dt )

type(maxwell_fdtd) :: this
sll_int32, intent(in) :: ccall
sll_int32 :: i1, j1, i2, j2
sll_real64 :: a11,a12,a21,a22,b1,b2,dis
sll_int32 :: i, j
sll_real64, intent(in) :: dt
sll_real64 :: dx, dy, c, csq
sll_real64, dimension(:,:), pointer :: ex, ey, bz

c   = this%c
csq = c * c
dx  = this%dx
dy  = this%dy

i1 = this%i1
j1 = this%j1
i2 = this%i2
j2 = this%j2

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
   do i = i1, j1
         
      a11 = 1.; a12 = + c
      a21 = 1./dt; a22 = - csq / dy
      b1  = - ex(i,j2) - c * bz(i,j2-1)
      b2  =   ex(i,j2)/dt - csq/dy*bz(i,j2-1)
         
      dis = a11*a22-a21*a12 
         
      !ex(i,j2) = (b1*a22-b2*a12)/dis
      bz(i,j2) = (a11*b2-a21*b1)/dis
         
   end do
      
case (SOUTH)

   !Frontiere Sud : Ex =  c Bz
   do i = i1, j1
         
      a11 = 1.; a12 = - c
      a21 = 1./dt; a22 = csq / dy
      b1  = - ex(i,i2) + c * bz(i,i2+1)
      b2  = ex(i,i2)/dt + csq / dy * bz(i,i2+1) 
         
      dis = a11*a22-a21*a12 
         
      ex(i,i2) = (b1*a22-b2*a12)/dis
      !bz(i,i2) = (a11*b2-a21*b1)/dis
         
   end do
      
case (EAST)

   !Frontiere Est : Ey =  c Bz
   do j = i2, j2
         
      a11 = 1.; a12 = - c
      a21 = 1./dt; a22 = + csq / dx
      b1  = - ey(j1,j) + c * bz(j1-1,j)
      b2  = ey(j1,j)/dt + csq/dx*bz(j1-1,j) 
         
      dis = a11*a22-a21*a12 
         
      !ey(j1,j) = (b1*a22-b2*a12)/dis
      bz(j1,j) = (a11*b2-a21*b1)/dis
      
   end do
      
case (WEST)

   !Frontiere Ouest : Ey = -c Bz
   do j = i2, j2
      
      a11 = 1.; a12 = + c
      a21 = 1./dt; a22 = - csq / dx
      b1  = - ey(i1,j) - c * bz(i1+1,j)
      b2  =   ey(i1,j)/dt - csq/dx*bz(i1+1,j) 
      
      dis = a11*a22-a21*a12 
   
      ey(i1,j) = (b1*a22-b2*a12)/dis
      !bz(i1,j) = (a11*b2-a21*b1)/dis
      
   end do

end select

end subroutine silver_muller

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module sll_maxwell_2d_fdtd
