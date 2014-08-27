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
!
!
!
!  Contact : Pierre Navaro http://wwww-irma.u-strasbg.fr/~navaro
!
!
!**************************************************************

module sll_maxwell_fdtd_hex_mesh

#include "sll_working_precision.h"
#include "sll_memory.h"
#include "sll_assert.h"
#include "sll_maxwell_solvers_macros.h"
#include "sll_constants.h"
use hex_mesh

implicit none
!private

!> Initialize maxwell solver 2d with FDTD scheme
interface initialize
 module procedure initialize_maxwell_hex_mesh_fdtd
end interface
!> Solve Ampere-Maxwell equation
interface ampere
 module procedure ampere_hex_mesh_fdtd
end interface
!> Solve Faraday equation
interface faraday
 module procedure ampere_hex_mesh_fdtd
end interface

public :: initialize, ampere, faraday

!> @brief Object with data to solve Maxwell equation 
!> Maxwell in TE mode: (Ex,Ey,Bz)
type, public :: maxwell_hex_mesh_fdtd

  type(hex_mesh_2d), pointer :: mesh          !< hexagonal mesh
  sll_real64                 :: c             !< light speed
  sll_real64                 :: e_0           !< electric conductivity
  sll_int32                  :: polarization  !< polarization type (TE or TM)

end type maxwell_hex_mesh_fdtd

contains

!>Initilialize the maxwell solver
subroutine initialize_maxwell_hex_mesh_fdtd(this, mesh, polarization )

   type(maxwell_hex_mesh_fdtd) :: this         !< maxwell solver object
   type(hex_mesh_2d)           :: mesh         !< hexagonal mesh
   sll_int32                   :: polarization !< TE or TM

   this%mesh = mesh
   this%c    = 1.0_f64
   this%e_0  = 1.0_f64
   
   this%polarization = polarization

end subroutine initialize_maxwell_hex_mesh_fdtd

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!> Solve Faraday equation
subroutine faraday_hex_mesh_fdtd( this, fx, fy, fz, dt )

type(maxwell_hex_mesh_fdtd)                 :: this !< Maxwell object
sll_real64, dimension(:,:), target :: fx   !< Ex or Bx
sll_real64, dimension(:,:), target :: fy   !< Ey or By
sll_real64, dimension(:,:), target :: fz   !< Bz or Ez
sll_real64, intent(in) :: dt               !< time step
sll_int32 :: i, j
sll_real64 :: dex_dy, dey_dx, dez_dx, dez_dy
sll_real64 :: delta
sll_real64, dimension(:,:), pointer :: ex
sll_real64, dimension(:,:), pointer :: ey
sll_real64, dimension(:,:), pointer :: ez
sll_real64, dimension(:,:), pointer :: bx
sll_real64, dimension(:,:), pointer :: by
sll_real64, dimension(:,:), pointer :: bz

delta = this%mesh%delta

!*** On utilise l'equation de Faraday sur un demi pas
!*** de temps pour le calcul du champ magnetique  Bz 
!*** a l'instant n puis n+1/2 

!i1 = this%i1
!j1 = this%j1
!i2 = this%i2
!j2 = this%j2
!
!if (this%polarization == TE_POLARIZATION) then
!
!   ex => fx; ey => fy; bz => fz
!   do i=i1,j1-1
!   do j=i2,j2-1
!      dex_dy  = (ex(i,j+1)-ex(i,j)) / dy
!      dey_dx  = (ey(i+1,j)-ey(i,j)) / dx
!      bz(i,j) = bz(i,j) + dt * (dex_dy - dey_dx)
!   end do
!   end do
!
!end if
!
!if (this%polarization == TM_POLARIZATION) then
!
!   bx => fx; by => fy; ez => fz
!   do i=i1,j1
!   do j=i2+1,j2
!      dez_dy  = (ez(i,j)-ez(i,j-1)) / dy
!      bx(i,j) = bx(i,j) - dt * dez_dy 
!   end do
!   end do
!   
!   do i=i1+1,j1
!   do j=i2,j2
!      dez_dx  = (ez(i,j)-ez(i-1,j)) / dx
!      by(i,j) = by(i,j) + dt * dez_dx 
!   end do
!   end do
!
!end if

end subroutine faraday_hex_mesh_fdtd

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!> Solve ampere-maxwell equation with FDTD scheme
subroutine ampere_hex_mesh_fdtd( this, fx, fy, fz, dt, jx, jy )

type(maxwell_hex_mesh_fdtd) :: this !< Maxwell object
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

!i1 = this%i1
!j1 = this%j1
!i2 = this%i2
!j2 = this%j2
!
!csq = this%c * this%c
!dx  = this%dx
!dy  = this%dy
!
!!*** Calcul du champ electrique E au temps n+1
!!*** sur les points internes du maillage
!!*** Ex aux points (i+1/2,j)
!!*** Ey aux points (i,j+1/2)
!
!if (this%polarization == TE_POLARIZATION) then
!
!   ex => fx; ey => fy; bz => fz
!
!   do i=i1,j1
!   do j=i2+1,j2
!      dbz_dy  = (bz(i,j)-bz(i,j-1)) / dy
!      ex(i,j) = ex(i,j) + dt*csq*dbz_dy 
!   end do
!   end do
!   
!   do i=i1+1,j1
!   do j=i2,j2
!      dbz_dx  = (bz(i,j)-bz(i-1,j)) / dx
!      ey(i,j) = ey(i,j) - dt*csq*dbz_dx 
!   end do
!   end do
!
!   if (present(jx) .and. present(jy)) then
!
!      ex(i1:j1,i2+1:j2) = ex(i1:j1,i2+1:j2) - dt * jx(i1:j1,i2+1:j2) / this%e_0
!      ey(i1+1:j1,i2:j2) = ey(i1+1:j1,i2:j2) - dt * jy(i1+1:j1,i2:j2) / this%e_0
!
!   endif
!
!end if
!
!if (this%polarization == TM_POLARIZATION) then
!
!   bx => fx; by => fy; ez => fz
!
!   do i=i1,j1-1
!   do j=i2,j2-1
!      dbx_dy  = (bx(i,j+1)-bx(i,j)) / dy
!      dby_dx  = (by(i+1,j)-by(i,j)) / dx
!      ez(i,j) = ez(i,j) - dt * (dbx_dy - dby_dx)
!   end do
!   end do
!
!   if (present(jx) .and. .not. present(jy)) then
!
!      ez(i1:j1-1,i2:j2-1) = ez(i1:j1-1,i2:j2-1) - dt * jx(i1:j1-1,i2:j2-1) / this%e_0
!
!   endif
!
!end if

end subroutine ampere_hex_mesh_fdtd

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module sll_maxwell_fdtd_hex_mesh
