!**************************************************************
!  Copyright INRIA
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
!> @brief
!> Implements the Maxwell solver in 2D with FDTD method
!>
!> @details
!> Equation
!>
!>\f$ \displaystyle \frac{\partial E_x}{\partial t} = {c^2} \frac{\partial B_z}{\partial y} \f$,
!>
!>\f$\displaystyle \frac{\partial E_y}{\partial t} = -{c^2} \frac{\partial B_z}{\partial x} \f$
!>
!>\f$\displaystyle \frac{\partial B_z}{\partial t} =  \frac{\partial E_x}{\partial y} - \frac{\partial E_y}{\partial x}  \f$.
!>
!>FDTD scheme
!>
!>\f$\displaystyle B_{z}^{n+1/2} = B_z^{n-1/2} - \Delta t \big( \frac{\partial E_y^n}{\partial x}
!>- \frac{\partial E_x^n}{\partial y} \big)\f$
!>
!>\f$\displaystyle E_x^{n+1} = E_x^{n} + c^2\Delta t  \frac{\partial B_z^{n+1/2}}{\partial y} \f$
!>
!>\f$\displaystyle E_y^{n+1} = E_y^{n} - c^2\Delta t  \frac{\partial B_z^{n+1/2}}{\partial x} \f$
!>
module sll_maxwell_2d_fdtd

#include "sll_working_precision.h"
#include "sll_memory.h"
#include "sll_assert.h"
#include "sll_maxwell_solvers_macros.h"
#include "sll_constants.h"

implicit none
!private

!> Initialize maxwell solver 2d with FDTD scheme
interface initialize
 module procedure new_maxwell_2d_fdtd
end interface
!> Solve maxwell solver 2d with FDTD scheme
interface solve
 module procedure solve_maxwell_2d_fdtd
end interface
!> Solve Ampere-Maxwell equation
interface ampere
 module procedure ampere_2d_fdtd
end interface
!> Solve Faraday equation
interface faraday
 module procedure ampere_2d_fdtd
end interface

public :: initialize, solve

!> @brief Object with data to solve Maxwell equation 
!> Maxwell in TE mode: (Ex,Ey,Bz)
type, public :: maxwell_2d_fdtd
  sll_real64 :: c            !< light speed
  sll_real64 :: e_0          !< electric conductivity
  sll_int32  :: i1           !< first indice of the block dimension 1
  sll_int32  :: j1           !< last indice of the block dimension 1
  sll_int32  :: i2           !< first indice of the block dimesnion 2
  sll_int32  :: j2           !< last indice of the block dimension 2
  sll_real64 :: dx           !< step size along dimension 1
  sll_real64 :: dy           !< step size along dimension 2
  sll_int32  :: polarization !< polarization type (TE or TM)
end type maxwell_2d_fdtd

contains

!>Initilialize the maxwell solver
subroutine new_maxwell_2d_fdtd(this, i1, j1, i2, j2, dx, dy, polarization )

   type(maxwell_2d_fdtd) :: this           !< maxwell solver object
   sll_int32          :: i1             !< first incidice along x
   sll_int32          :: j1             !< last indice along x
   sll_int32          :: i2             !< first indice along y
   sll_int32          :: j2             !< last indice along y
   sll_real64         :: dx             !< size step along x
   sll_real64         :: dy             !< size step along y
   sll_int32          :: polarization   !< TE or TM
   !sll_int32        :: error           !< error code

   this%c   = 1.0_f64
   this%e_0 = 1.0_f64
   
   this%dx = dx
   this%dy = dy
   this%i1 = i1
   this%j1 = j1
   this%i2 = i2
   this%j2 = j2
   this%polarization = polarization

end subroutine new_maxwell_2d_fdtd

!> this routine exists only for testing purpose. Use ampere and faraday
!> in your appication.
subroutine solve_maxwell_2d_fdtd(this, fx, fy, fz, dt)

   type(maxwell_2d_fdtd)         :: this !< maxwell object
   sll_real64, dimension(:,:) :: fx   !< Ex or Bx
   sll_real64, dimension(:,:) :: fy   !< Ey or By
   sll_real64, dimension(:,:) :: fz   !< Bz or Ez
   sll_real64, intent(in)     :: dt   !< time step

   call faraday_2d_fdtd(this, fx, fy, fz, 0.5*dt)   
   call cl_periodiques_2d_fdtd(this, fx, fy, fz, dt)
   call ampere_2d_fdtd(this, fx, fy, fz, dt) 
   call cl_periodiques_2d_fdtd(this, fx, fy, fz, dt)
   call faraday_2d_fdtd(this, fx, fy, fz, 0.5*dt)   
   call cl_periodiques_2d_fdtd(this, fx, fy, fz, dt)

end subroutine solve_maxwell_2d_fdtd

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!> Solve Faraday equation
subroutine faraday_2d_fdtd( this, fx, fy, fz, dt )

type(maxwell_2d_fdtd)                 :: this !< Maxwell object
sll_real64, dimension(:,:), target :: fx   !< Ex or Bx
sll_real64, dimension(:,:), target :: fy   !< Ey or By
sll_real64, dimension(:,:), target :: fz   !< Bz or Ez
sll_real64, intent(in) :: dt               !< time step
sll_int32 :: i, j
sll_real64 :: dex_dy, dey_dx, dez_dx, dez_dy
sll_real64 :: dx, dy
sll_int32  :: i1, j1, i2, j2
sll_real64, dimension(:,:), pointer :: ex
sll_real64, dimension(:,:), pointer :: ey
sll_real64, dimension(:,:), pointer :: ez
sll_real64, dimension(:,:), pointer :: bx
sll_real64, dimension(:,:), pointer :: by
sll_real64, dimension(:,:), pointer :: bz

dx  = this%dx
dy  = this%dy

!*** On utilise l'equation de Faraday sur un demi pas
!*** de temps pour le calcul du champ magnetique  Bz 
!*** a l'instant n puis n+1/2 

i1 = this%i1
j1 = this%j1
i2 = this%i2
j2 = this%j2

if (this%polarization == TE_POLARIZATION) then

   ex => fx; ey => fy; bz => fz
   do i=i1,j1-1
   do j=i2,j2-1
      dex_dy  = (ex(i,j+1)-ex(i,j)) / dy
      dey_dx  = (ey(i+1,j)-ey(i,j)) / dx
      bz(i,j) = bz(i,j) + dt * (dex_dy - dey_dx)
   end do
   end do

end if

if (this%polarization == TM_POLARIZATION) then

   bx => fx; by => fy; ez => fz
   do i=i1,j1
   do j=i2+1,j2
      dez_dy  = (ez(i,j)-ez(i,j-1)) / dy
      bx(i,j) = bx(i,j) - dt * dez_dy 
   end do
   end do
   
   do i=i1+1,j1
   do j=i2,j2
      dez_dx  = (ez(i,j)-ez(i-1,j)) / dx
      by(i,j) = by(i,j) + dt * dez_dx 
   end do
   end do

end if

end subroutine faraday_2d_fdtd

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!> Solve ampere-maxwell equation with FDTD scheme
subroutine ampere_2d_fdtd( this, fx, fy, fz, dt, jx, jy )

type(maxwell_2d_fdtd) :: this !< Maxwell object
sll_int32 :: i1, j1, i2, j2
sll_real64, dimension(:,:), intent(inout), target :: fx !< Ex or Bx
sll_real64, dimension(:,:), intent(inout), target :: fy !< Ey or By
sll_real64, dimension(:,:), intent(inout), target :: fz !< Bz or Ez
sll_real64, dimension(:,:), intent(in), optional :: jx  !< Jx current
sll_real64, dimension(:,:), intent(in), optional :: jy  !< Jy current
sll_int32 :: i, j
sll_real64 :: dbz_dx, dbz_dy, dbx_dy, dby_dx
sll_real64, intent(in) :: dt !< time step
sll_real64 :: dx, dy
sll_real64 :: csq
sll_real64, dimension(:,:), pointer :: ex
sll_real64, dimension(:,:), pointer :: ey
sll_real64, dimension(:,:), pointer :: ez
sll_real64, dimension(:,:), pointer :: bx
sll_real64, dimension(:,:), pointer :: by
sll_real64, dimension(:,:), pointer :: bz

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

if (this%polarization == TE_POLARIZATION) then

   ex => fx; ey => fy; bz => fz

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

end if

if (this%polarization == TM_POLARIZATION) then

   bx => fx; by => fy; ez => fz

   do i=i1,j1-1
   do j=i2,j2-1
      dbx_dy  = (bx(i,j+1)-bx(i,j)) / dy
      dby_dx  = (by(i+1,j)-by(i,j)) / dx
      ez(i,j) = ez(i,j) - dt * (dbx_dy - dby_dx)
   end do
   end do

   if (present(jx) .and. .not. present(jy)) then

      ez(i1:j1-1,i2:j2-1) = ez(i1:j1-1,i2:j2-1) - dt * jx(i1:j1-1,i2:j2-1) / this%e_0

   endif

end if

end subroutine ampere_2d_fdtd

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!> Set boundary conditions 
subroutine cl_periodiques_2d_fdtd(this, fx, fy, fz, dt)

type(maxwell_2d_fdtd) :: this !< maxwell solver object
sll_int32 :: i1, j1, i2, j2
sll_real64, dimension(:,:), intent(inout) :: fx !< Ex or Bx
sll_real64, dimension(:,:), intent(inout) :: fy !< Ey or By
sll_real64, dimension(:,:), intent(inout) :: fz !< Bz or Ez
sll_real64 :: dbz_dx, dbz_dy, dez_dx, dez_dy
sll_real64, intent(in) :: dt !< time step
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

fz(i1:j1-1,j2) = fz(i1:j1-1,i2)
fz(j1,i2:j2-1) = fz(i1,i2:j2-1)

fz(j1,j2) = fz(i1,i2)

if ( this%polarization == TE_POLARIZATION) then
   do i = i1, j1
      dbz_dy = (fz(i,i2)-fz(i,j2-1)) / dy
      fx(i,i2) = fx(i,j2) + dt*csq*dbz_dy 
   end do
     
   do j = i2, j2
      dbz_dx = (fz(i1,j)-fz(j1-1,j)) / dx
      fy(i1,j) = fy(j1,j) - dt*csq*dbz_dx 
   end do
end if

if ( this%polarization == TM_POLARIZATION) then
   do i = i1, j1
      dez_dy = (fz(i,i2)-fz(i,j2-1)) / dy
      fx(i,i2) = fx(i,j2) - dt*csq*dez_dy 
   end do
     
   do j = i2, j2
      dez_dx = (fz(i1,j)-fz(j1-1,j)) / dx
      fy(i1,j) = fy(j1,j) + dt*csq*dez_dx 
   end do
end if

end subroutine cl_periodiques_2d_fdtd

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!> Set periodic bounday conditions
subroutine cl_condparfait_2d_fdtd(this, fx, fy, fz, side)

type(maxwell_2d_fdtd) :: this !< maxwell object
sll_int32, intent(in) :: side !< which domain edge
sll_real64, dimension(:,:), target  :: fx !< Ex or Bx
sll_real64, dimension(:,:), target  :: fy !< Ey or By
sll_real64, dimension(:,:), target  :: fz !< Bz or Ez
sll_real64, dimension(:,:), pointer :: bx, by, ez
sll_real64, dimension(:,:), pointer :: ex, ey, bz
sll_int32 :: i1, j1, i2, j2

i1 = this%i1
j1 = this%j1
i2 = this%i2
j2 = this%j2

if (this%polarization == TE_POLARIZATION) then
   ex => fx; ey => fy; bz => fz
   select case(side)
   case(SOUTH)
      ex(i1:j1,i2) = 0.d0
   case(NORTH)
      bz(i1:j1,j2) = bz(i1:j1,j2-1)
   case(WEST)
      ey(i1,i2:j2) = 0.d0
   case(EAST)
      bz(j1,i2:j2) = bz(j1-1,i2:j2)
   end select
end if

if (this%polarization == TM_POLARIZATION) then
   bx => fx; by => fy; ez => fz
   select case(side)
   case(SOUTH)
      bx(i1:j1,i2) = 0.0_f64
   case(NORTH)
      ez(i1:j1,j2) = ez(i1:j1,j2-1)
   case(WEST)
      by(i1,i2:j2) = 0.d0
   case(EAST)
      ez(j1,i2:j2) = ez(j1-1,i2:j2)
   end select
end if

end subroutine cl_condparfait_2d_fdtd

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!> Bundary conditions
subroutine silver_muller_2d_fdtd( this, ex, ey, bz, ccall, dt )

type(maxwell_2d_fdtd) :: this !< maxwell object
sll_int32, intent(in) :: ccall !< domain edge (N,S,E,W)
sll_int32 :: i1, j1, i2, j2
sll_real64 :: a11,a12,a21,a22,b1,b2,dis
sll_int32 :: i, j
sll_real64, intent(in) :: dt !< time step
sll_real64 :: dx, dy, c, csq
sll_real64, dimension(:,:), pointer :: ex !< x electric field
sll_real64, dimension(:,:), pointer :: ey !< y electric field
sll_real64, dimension(:,:), pointer :: bz !< z magnetic field

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

end subroutine silver_muller_2d_fdtd

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module sll_maxwell_2d_fdtd
