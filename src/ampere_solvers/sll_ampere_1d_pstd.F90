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
#include "sll_fftw.h"

use fftw3

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
   fftw_plan           :: fwx          !< forward fft plan along x
   fftw_plan           :: bwx          !< backward fft plan along x
   fftw_plan           :: p_ek         !< pointer for fft memory allocation
   fftw_plan           :: p_jk         !< pointer for fft memory allocation
   fftw_comp, pointer  :: ek(:)        !< e fft transform
   fftw_comp, pointer  :: jk(:)        !< j fft transform
   fftw_int            :: sz_ek        !< size for memory allocation
   fftw_int            :: sz_jk        !< size for memory allocation

contains

   procedure :: compute_e_from_j => solve_ampere_1d_pstd


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
sll_int32 :: error

contains

!> Initialize 2d ampere solver on cartesian mesh with PSTD scheme
function new_ampere_1d_pstd(xmin,xmax,nc_x) result(self)

   type(sll_ampere_1d_pstd), pointer :: self    !< Solver object
   sll_real64, intent(in)            :: xmin    !< x min
   sll_real64, intent(in)            :: xmax    !< x max
   sll_int32,  intent(in)            :: nc_x    !< x cells number


   SLL_ALLOCATE(self, error)

   call initialize_ampere_1d_pstd(self,xmin,xmax,nc_x) 

end function new_ampere_1d_pstd

!> Initialize 2d ampere solver on cartesian mesh with PSTD scheme
subroutine initialize_ampere_1d_pstd(self,xmin,xmax,nc_x) 

   type(sll_ampere_1d_pstd) :: self    !< Solver object
   sll_real64, intent(in)   :: xmin    !< x min
   sll_real64, intent(in)   :: xmax    !< x max
   sll_int32,  intent(in)   :: nc_x    !< x cells number

   sll_real64, allocatable  :: tmp(:)

   self%nc_eta1 = nc_x
   self%eta1_min = xmin
   self%eta1_max = xmax

   self%e_0  = 1._f64
   self%mu_0 = 1._f64

   allocate(tmp(nc_x))

   FFTW_ALLOCATE(self%ek,(nc_x/2+1),self%sz_ek,self%p_ek)
   FFTW_ALLOCATE(self%jk,(nc_x/2+1),self%sz_jk,self%p_jk)

   NEW_FFTW_PLAN_R2C_1D(self%fwx, nc_x, tmp,  self%ek)
   NEW_FFTW_PLAN_C2R_1D(self%bwx, nc_x, self%jk,  tmp)

   deallocate(tmp)

end subroutine initialize_ampere_1d_pstd

!> Solve \f$ \frac{\partial E_x}{\partial_t} + J_x = 0 \f$
!> using spectral method
subroutine solve_ampere_1d_pstd(self, dt, jx, ex)

   class(sll_ampere_1d_pstd), intent(inout) :: self   !< Solver object
   sll_real64,               intent(in)     :: dt     !< time step
   sll_real64, dimension(:), intent(inout)  :: jx     !< J current x
   sll_real64, dimension(:), intent(inout)  :: ex     !< E field x

   sll_int32   :: nc_x
   sll_real64  :: dt_e

   dt_e = dt / self%e_0
   nc_x = self%nc_eta1 

   call fftw_execute_dft_r2c(self%fwx, ex, self%ek)
   call fftw_execute_dft_r2c(self%fwx, jx, self%jk)
   self%ek = self%ek - dt_e * self%jk
   call fftw_execute_dft_c2r(self%bwx, self%ek, ex)
   call fftw_execute_dft_c2r(self%bwx, self%jk, jx)
   ex = ex / nc_x
   jx = jx / nc_x

end subroutine solve_ampere_1d_pstd

!> delete ampere solver object
subroutine free_ampere_1d_pstd(self)

   type(sll_ampere_1d_pstd) :: self
   
#ifdef FFTW_F2003
   if (c_associated(self%p_ek)) call fftw_free(self%p_ek)
   if (c_associated(self%p_jk)) call fftw_free(self%p_jk)
#endif
   
   call fftw_destroy_plan(self%fwx)
   call fftw_destroy_plan(self%bwx)

end subroutine free_ampere_1d_pstd

end module sll_module_ampere_1d_pstd
