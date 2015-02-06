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
module sll_module_ampere_vlasov_1d
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
type, public :: sll_ampere_1d

   private
   sll_int32           :: nc_eta1      !< x cells number
   sll_int32           :: nc_eta2      !< v cells number
   sll_real64          :: e_0          !< electric conductivity
   sll_real64          :: mu_0         !< magnetic permeability
   sll_real64          :: c            !< speed of light
   sll_real64          :: eta1_min     !< left side 
   sll_real64          :: eta1_max     !< right side
   sll_real64          :: delta_eta1   !< step size
contains
   procedure :: compute_e_from_j => solve_ampere_1d
end type sll_ampere_1d

type, extends(sll_ampere_1d), public :: sll_ampere_vlasov_1d

   sll_real64          :: eta2_min     !< left side 
   sll_real64          :: eta2_max     !< right side
   sll_real64          :: delta_eta2   !< step size
   fftw_plan           :: fwx          !< forward fft plan along x
   fftw_plan           :: bwx          !< backward fft plan along x
   fftw_plan           :: p_ek         !< pointer for fft memory allocation
   fftw_plan           :: p_rk         !< pointer for fft memory allocation
   fftw_plan           :: p_fk         !< pointer for fft memory allocation
   fftw_comp, pointer  :: ek(:)        !< e fft transform
   fftw_comp, pointer  :: rk(:)        !< rho fft transform
   fftw_comp, pointer  :: fk(:)        !< f fft transform
   fftw_int            :: sz_ek        !< size for memory allocation
   fftw_int            :: sz_fk        !< size for memory allocation
   fftw_int            :: sz_rk        !< size for memory allocation
   sll_real64, pointer :: kx(:)        !< wave number
   sll_real64, pointer :: df_dx(:)     !< f derivative x

contains

   !> Compute electric field from current
   procedure :: compute_e_from_f => solve_ampere_vlasov_1d


end type sll_ampere_vlasov_1d

!> Initialize ampere solver 2d cartesian periodic with PSTD scheme
interface sll_create
 module procedure initialize_ampere_1d
 module procedure initialize_ampere_vlasov_1d
end interface sll_create

!> Solve ampere solver 2d cartesian periodic with PSTD scheme
interface sll_solve
 module procedure solve_ampere_1d
 module procedure solve_ampere_vlasov_1d
end interface sll_solve

!> Delete ampere solver 2d cartesian periodic with PSTD scheme
interface sll_delete
 module procedure free_ampere_vlasov_1d
end interface sll_delete

public sll_create
public sll_delete
public sll_solve
public new_ampere_1d
public new_ampere_vlasov_1d

private
sll_int32 :: error !< error code

contains

!> Initialize 1d ampere solver on cartesian mesh 
function new_ampere_1d(xmin,xmax,nc_x) result(self)

   type(sll_ampere_1d), pointer :: self    !< Solver object
   sll_real64, intent(in)       :: xmin    !< x min
   sll_real64, intent(in)       :: xmax    !< x max
   sll_int32,  intent(in)       :: nc_x    !< x cells number

   SLL_ALLOCATE(self, error)

   call initialize_ampere_1d(self,xmin,xmax,nc_x)

end function new_ampere_1d

!> Initialize 1d ampere vlasov solver on cartesian mesh 
function new_ampere_vlasov_1d(xmin,xmax,nc_x, &
                              vmin,vmax,nc_v) result(self)

   type(sll_ampere_vlasov_1d), pointer :: self    !< Solver object
   sll_real64, intent(in)              :: xmin    !< x min
   sll_real64, intent(in)              :: xmax    !< x max
   sll_int32,  intent(in)              :: nc_x    !< x cells number
   sll_real64, intent(in)              :: vmin    !< v min
   sll_real64, intent(in)              :: vmax    !< v max
   sll_int32,  intent(in)              :: nc_v    !< v cells number

   SLL_ALLOCATE(self, error)

   call initialize_ampere_vlasov_1d(self,             &
                                    xmin,xmax,nc_x,   &
                                    vmin,vmax,nc_v) 

end function new_ampere_vlasov_1d

!> Initialize 1d ampere solver on cartesian mesh 
subroutine initialize_ampere_1d(self, xmin, xmax, nc_x)

   type(sll_ampere_1d)     :: self    !< Solver object
   sll_real64, intent(in)  :: xmin    !< x min
   sll_real64, intent(in)  :: xmax    !< x max
   sll_int32,  intent(in)  :: nc_x    !< x cells number

   self%nc_eta1  = nc_x
   self%eta1_min = xmin
   self%eta1_max = xmax

   self%e_0  = 1._f64
   self%mu_0 = 1._f64

end subroutine initialize_ampere_1d

!> Initialize 1d ampere vlasov solver on cartesian mesh
subroutine initialize_ampere_vlasov_1d(self,             &
                                       xmin, xmax, nc_x, & 
                                       vmin, vmax, nc_v  )

   type(sll_ampere_vlasov_1d) :: self    !< Solver object
   sll_real64, intent(in)     :: xmin    !< x min
   sll_real64, intent(in)     :: xmax    !< x max
   sll_int32,  intent(in)     :: nc_x    !< x cells number
   sll_real64, intent(in)     :: vmin    !< v min
   sll_real64, intent(in)     :: vmax    !< v max
   sll_int32,  intent(in)     :: nc_v    !< v cells number

   sll_real64, allocatable    :: tmp(:)
   sll_int32                  :: i
   sll_real64                 :: kx0

   self%nc_eta1  = nc_x
   self%eta1_min = xmin
   self%eta1_max = xmax

   self%nc_eta2  = nc_v
   self%eta2_min = vmin
   self%eta2_max = vmax

   self%e_0  = 1._f64
   self%mu_0 = 1._f64

   FFTW_ALLOCATE(self%ek,(nc_x/2+1),self%sz_ek,self%p_ek)
   FFTW_ALLOCATE(self%rk,(nc_x/2+1),self%sz_rk,self%p_rk)
   FFTW_ALLOCATE(self%fk,(nc_x/2+1),self%sz_fk,self%p_fk)

   SLL_CLEAR_ALLOCATE(self%df_dx(1:nc_x), error)

   NEW_FFTW_PLAN_R2C_1D(self%fwx, nc_x, self%df_dx,  self%ek)
   NEW_FFTW_PLAN_C2R_1D(self%bwx, nc_x, self%rk,  self%df_dx)

   SLL_CLEAR_ALLOCATE(self%kx(1:nc_x/2+1), error)
    
   kx0 = 2.0_f64*sll_pi/(xmax-xmin)

   self%kx(1) = 1.0_f64
   do i=2,nc_x/2+1
      self%kx(i) = (i-1)*kx0
   end do
        

end subroutine initialize_ampere_vlasov_1d

!> @brief
!> Solve \f$ \frac{\partial E_x}{\partial_t} + J_x = 0 \f$
!> @details
!> Simplest time scheme
subroutine solve_ampere_1d(self, dt, jx, ex)

   class(sll_ampere_1d),     intent(inout) :: self   !< Solver object
   sll_real64,               intent(in)    :: dt     !< time step
   sll_real64, dimension(:), intent(inout) :: jx     !< J current x
   sll_real64, dimension(:), intent(inout) :: ex     !< E field x

   sll_real64  :: dt_e

   dt_e = dt / self%e_0

   ex = ex - dt_e * jx

end subroutine solve_ampere_1d

!> @brief
!> Solve Vlasov-Ampere equation in 1D using discrete Fourier transforms
!> @details
!> Solve 
!> \f[ f_k^{n+1}(v_j) = \exp{-2ikv\Delta t}/L f_k^n(v_j) \f]
!> \f[ \rho_k^{n+1}   = \Delta v \sum_j f_k^{n+1}(v_j) \f]
!> \f[ E^{n+1}_k      = \rho_k^{n+1} (eL) / (2ik \epsilon_0 \f]
subroutine solve_ampere_vlasov_1d(self, dt, f, ex)

   class(sll_ampere_vlasov_1d), intent(inout) :: self !< Solver object
   sll_real64,                  intent(in)    :: dt   !< time step
   sll_real64, dimension(:,:),  intent(inout) :: f    !< f distribution function
   sll_real64, dimension(:),    intent(inout) :: ex   !< E field x

   sll_int32   :: j
   sll_int32   :: nc_x, nc_y
   sll_real64  :: dt_e, a

   dt_e = dt / self%e_0
   nc_x = self%nc_eta1 
   nc_y = self%nc_eta2

   !self%rk = cmplx(0.0,0.0,kind=f64)
   do j = 1, nc_y+1
     a = (self%eta2_min + (j-1) * self%delta_eta2)*dt
     self%df_dx = f(1:nc_x,j)
     call fftw_execute_dft_r2c(self%fwx, self%df_dx, self%fk)
     !self%rk = self%rk + self%delta_eta2 * self%fk
     self%fk = self%fk*cmplx(cos(self%kx*a),-sin(self%kx*a),kind=f64)
     call fftw_execute_dft_c2r(self%bwx, self%fk, self%df_dx)
     f(1:nc_x,j) = self%df_dx / nc_x
     f(nc_x+1,j) = f(1,j)
   end do


   !self%ek = self%rk / cmplx(cos(self%kx*self%e_0),sin(self%kx*self%e_0),kind=f64)
   !call fftw_execute_dft_c2r(self%bwx, self%ek, ex)
   !ex = ex / nc_x

end subroutine solve_ampere_vlasov_1d

!> delete ampere solver object
subroutine free_ampere_vlasov_1d(self)

   type(sll_ampere_vlasov_1d) :: self
   
#ifdef FFTW_F2003
   if (c_associated(self%p_ek)) call fftw_free(self%p_ek)
   if (c_associated(self%p_rk)) call fftw_free(self%p_rk)
   if (c_associated(self%p_fk)) call fftw_free(self%p_fk)
#endif
   
   call fftw_destroy_plan(self%fwx)
   call fftw_destroy_plan(self%bwx)

end subroutine free_ampere_vlasov_1d

end module sll_module_ampere_vlasov_1d
