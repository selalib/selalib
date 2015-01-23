!******************************************************************
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
!  Contact : Pierre Navaro http://wwww-irma.u-strasbg.fr/~navaro
!
!******************************************************************

!> @ingroup ampere_solvers
!> @brief
!> Implements the Ampere solver in 1D 
!> @details
!> \f[
!> \nabla \times H = j + \frac{\partial D}{\partial t} 
!> \f]
module sll_module_ampere_1d_pstd
#include "sll_working_precision.h"
#include "sll_memory.h"
#include "sll_assert.h"
#include "sll_constants.h"

implicit none

!> @brief 
!> Ampere solver object
!> @details
!> We solve Maxwell system with PSTD numerical method. The type contains
!> information about FFT, mesh and physical properties.
type, public :: sll_ampere_1d_pstd

   private
   sll_int32           :: nc_eta1      !< x cells number
   sll_real64          :: e_0          !< electric conductivity
   sll_real64          :: mu_0         !< magnetic permeability
   sll_real64          :: c            !< speed of light
   sll_real64          :: eta1_min     !< left side 
   sll_real64          :: eta1_max     !< right side
   sll_real64          :: delta_eta1   !< step size

end type sll_ampere_1d_pstd

!> Initialize ampere solver 2d cartesian periodic with PSTD scheme
interface sll_create
 module procedure initialize_ampere_1d_pstd
end interface sll_create

!> Solve ampere solver 2d cartesian periodic with PSTD scheme
interface sll_solve
 module procedure solve_ampere_1d_pstd
end interface sll_solve

!> Delete ampere solver 2d cartesian periodic with PSTD scheme
interface sll_delete
 module procedure free_ampere_1d_pstd
end interface sll_delete

public sll_create
public sll_delete
public sll_solve
public new_ampere_1d_pstd

private

contains

!> Initialize 2d ampere solver on cartesian mesh with PSTD scheme
function new_ampere_1d_pstd(xmin,xmax,nc_x) result(self)

   type(sll_ampere_1d_pstd), pointer :: self    !< Solver object
   sll_real64, intent(in)            :: xmin    !< x min
   sll_real64, intent(in)            :: xmax    !< x max
   sll_int32,  intent(in)            :: nc_x    !< x cells number

   sll_int32                         :: error   !< error code
   sll_real64                        :: dx      !< x space step
   sll_real64                        :: kx0
   sll_int32                         :: i

   self%nc_eta1 = nc_x

   self%e_0  = 1._f64
   self%mu_0 = 1._f64

end function new_ampere_1d_pstd

!> Initialize 2d ampere solver on cartesian mesh with PSTD scheme
subroutine initialize_ampere_1d_pstd(self,xmin,xmax,nc_x) 

   type(sll_ampere_1d_pstd) :: self    !< Solver object
   sll_real64, intent(in)   :: xmin    !< x min
   sll_real64, intent(in)   :: xmax    !< x max
   sll_int32,  intent(in)   :: nc_x    !< x cells number

   sll_int32                :: error   !< error code
   sll_real64               :: dx      !< x space step
   sll_real64               :: kx0
   sll_int32                :: i

   self%nc_eta1 = nc_x

   self%e_0  = 1._f64
   self%mu_0 = 1._f64

end subroutine initialize_ampere_1d_pstd


!> Solve \f$ \frac{\partial E_x}{\partial_t} + J_x = 0 \f$
subroutine solve_ampere_1d_pstd(self, ex, dt, jx)

   type(sll_ampere_1d_pstd),intent(inout)  :: self   !< Solver object
   sll_real64, dimension(:), intent(inout) :: ex     !< E field x
   sll_real64 , intent(in)                 :: dt     !< time step
   sll_real64, dimension(:), optional      :: jx     !< J current x

   sll_real64  :: dt_e

   dt_e = dt / self%e_0

   if (present(jx)) then
      ex = ex - dt_e * jx 
   end if

end subroutine solve_ampere_1d_pstd

!> delete ampere solver object
subroutine free_ampere_1d_pstd(self)
type(sll_ampere_1d_pstd) :: self

end subroutine free_ampere_1d_pstd

end module sll_module_ampere_1d_pstd
