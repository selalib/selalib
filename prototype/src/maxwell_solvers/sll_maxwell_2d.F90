!------------------------------------------------------------------------------
! SELALIB
!------------------------------------------------------------------------------
!
! MODULE: sll_maxwell_2d
!
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
!> - mesh_types
!> - diagnostics
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
#include "sll_mesh_types.h"

use numeric_constants

implicit none
private
sll_real64 :: csq

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
  type(field_2d_vec3),      pointer :: fields
  type(field_2d_vec3),      pointer :: dt_fields  
  type(mesh_descriptor_2d), pointer :: descriptor
  sll_real64 :: c_light
  sll_real64 :: epsilon_0
end type maxwell_2d

enum, bind(C)
   enumerator :: NORTH = 0, EAST = 1, SOUTH = 2, WEST = 3
end enum


contains

function new_maxwell_2d(fields)

   type(maxwell_2d),         pointer :: new_maxwell_2d
   type(field_2D_vec3),      pointer :: fields
   type(mesh_descriptor_2d), pointer :: mesh
   sll_int32                         :: nc_eta1
   sll_int32                         :: nc_eta2
   sll_int32                         :: error

   SLL_ASSERT(associated(fields))

   mesh => fields%descriptor

   nc_eta1 = GET_MESH_NC_ETA1(mesh)
   nc_eta2 = GET_MESH_NC_ETA2(mesh)

   SLL_ALLOCATE(new_maxwell_2d, error)
   SLL_ALLOCATE(new_maxwell_2d%descriptor, error)

   new_maxwell_2d%c_light   = 1.0_f64
   new_maxwell_2d%epsilon_0 = 1.0_f64

   SLL_ALLOCATE(new_maxwell_2d%dt_fields,error)
   !new_maxwell_2d%dt_fields => new_field_2D_vec3(mesh)
   new_maxwell_2d%fields => fields

end function new_maxwell_2d

subroutine solve_maxwell_2d(this,dt)

   type(maxwell_2d), pointer :: this
   sll_real64 , intent(in)   :: dt
   sll_int32                 :: nc_eta1
   sll_int32                 :: nc_eta2

   SLL_ASSERT(associated(this))
   SLL_ASSERT(associated(this%fields))
   nc_eta1 = this%fields%descriptor%nc_eta1
   nc_eta2 = this%fields%descriptor%nc_eta2

   !B(n-1/2)--> B(n+1/2) sur les pts interieurs   
   call faraday(this, 1, nc_eta1+1, 1, nc_eta2+1, dt)   

   call cl_periodiques(this, 1, nc_eta1+1, 1, nc_eta2+1, dt)

   !E(n)-->E(n+1) sur les pts interieurs
   call ampere_maxwell(this, 1, nc_eta1+1, 1, nc_eta2+1, dt) 

 end subroutine solve_maxwell_2d

!> Delete the Maxwell object
subroutine delete_maxwell_2d(this)
  
  type(maxwell_2d), pointer :: this

  this => null()

end subroutine delete_maxwell_2d

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine faraday( this, ix, jx, iy, jy, dt )

type(maxwell_2d), pointer :: this
sll_int32, intent(in) :: ix, jx, iy, jy
sll_real64, dimension(:,:), pointer :: ex, ey, bz
sll_int32 :: i, j
sll_real64 :: dex_dy, dey_dx
sll_real64, intent(in) :: dt
sll_real64 :: dx, dy

SLL_ASSERT(associated(this))
SLL_ASSERT(associated(this%fields))

ex => this%fields%data%v1 
ey => this%fields%data%v2 
bz => this%fields%data%v3 

csq = this%c_light * this%c_light
dx  = this%fields%descriptor%delta_eta1
dy  = this%fields%descriptor%delta_eta2

!*** On utilise l'equation de Faraday sur un demi pas
!*** de temps pour le calcul du champ magnetique  Bz 
!*** a l'instant n puis n+1/2 

do i=ix,jx-1
do j=iy,jy-1
   dex_dy  = (ex(i,j+1)-ex(i,j)) / dy
   dey_dx  = (ey(i+1,j)-ey(i,j)) / dx
   bz(i,j) = bz(i,j) + dt * (dex_dy - dey_dx)
end do
end do

end subroutine faraday

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine ampere_maxwell( this, ix, jx, iy, jy, dt )

type(maxwell_2d), pointer :: this
sll_int32, intent(in) :: ix, jx, iy, jy
sll_real64, dimension(:,:), pointer :: ex, ey, bz
sll_int32 :: i, j
sll_real64 :: dex_dx, dey_dy
sll_real64 :: dex_dy, dey_dx
sll_real64 :: dbz_dx, dbz_dy
sll_real64, intent(in) :: dt
sll_real64 :: dx, dy

ex => this%fields%data%v1 
ey => this%fields%data%v2 
bz => this%fields%data%v3 

csq = this%c_light * this%c_light
dx  = this%fields%descriptor%delta_eta1
dy  = this%fields%descriptor%delta_eta2

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

subroutine cl_periodiques(this, ix, jx, iy, jy, dt)

type(maxwell_2d) :: this
sll_int32, intent(in) :: ix, jx, iy, jy
sll_real64, dimension(:,:), pointer :: ex, ey, bz
sll_real64 :: dex_dx, dey_dy
sll_real64 :: dex_dy, dey_dx
sll_real64 :: dbz_dx, dbz_dy
sll_real64, intent(in) :: dt
sll_real64 :: dx, dy
sll_int32 :: i, j

ex => this%fields%data%v1 
ey => this%fields%data%v2 
bz => this%fields%data%v3 

csq = this%c_light * this%c_light
dx  = this%fields%descriptor%delta_eta1
dy  = this%fields%descriptor%delta_eta2

do i = ix, jx-1
   bz(i,jy) = bz(i,iy)
end do
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

subroutine cl_condparfait(this, ix, jx, iy, jy, side)

type(maxwell_2d) :: this
sll_int32, intent(in) :: side
sll_int32, intent(in) :: ix, jx, iy, jy
sll_int32 :: i, j
sll_real64, dimension(:,:), pointer :: ex, ey, bz

ex => this%fields%data%v1 
ey => this%fields%data%v2 
bz => this%fields%data%v3 

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

subroutine silver_muller( this, ix, jx, iy, jy, ccall, dt )

type(maxwell_2d) :: this
sll_int32, intent(in) :: ccall
sll_int32, intent(in) :: ix, jx, iy, jy
sll_real64 :: a11,a12,a21,a22,b1,b2,dis
sll_int32 :: i, j
sll_real64, intent(in) :: dt
sll_real64 :: dx, dy, c, csq
sll_real64, dimension(:,:), pointer :: ex, ey, bz

ex => this%fields%data%v1 
ey => this%fields%data%v2 
bz => this%fields%data%v3 
c   = this%c_light 
csq = c * c
dx  = this%fields%descriptor%delta_eta1
dy  = this%fields%descriptor%delta_eta2

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
         
      !fields%dataex(i,jy) = (b1*a22-b2*a12)/dis
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
