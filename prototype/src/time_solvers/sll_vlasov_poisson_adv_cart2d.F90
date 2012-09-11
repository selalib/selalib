module sll_vp_cartesian_2d
#include "sll_working_precision.h"
#include "sll_memory.h"
#include "sll_assert.h"
  use sll_time_splitting
  use sll_module_interpolators_1d_base
  use distribution_function
  use sll_poisson_1d_periodic
  implicit none

  type, extends(time_splitting) :: vp_cartesian_2d
     
     type(sll_distribution_function_2d), pointer   :: f
     type(poisson_1d_periodic), pointer            :: poisson_1d
     sll_int32 :: n1, n2
     class(sll_interpolator_1d_base), pointer    :: interp1, interp2
   contains
     procedure, pass(this) :: operator1 => advection_x
     procedure, pass(this) :: operator2 => advection_v
  end type const_coef_advection_2d

contains
  subroutine initialize(this, data, n1, n2, a1, a2, interp1, interp2)
    type(const_coef_advection_2d) :: this 
    sll_real64, dimension(:,:), target :: data
    sll_int32 :: n1, n2
    sll_real64 :: a1, a2
    class(sll_interpolator_1d_base), pointer    :: interp1, interp2
    this%data => data
    this%n1 = n1
    this%n2 = n2
    this%a1 = a1
    this%a2 = a2
    this%interp1 => interp1
    this%interp2 => interp2
  end subroutine initialize

  subroutine advection_x(this, dt)
    class(vp_cartesian_2d) :: this 
    sll_real64 :: dt
    ! local variables
    sll_real64, dimension(:), pointer :: f1d
    sll_real64 :: displacement
    sll_int32 :: j
    
    do j = 1, this%n2
       displacement = this%a1 * dt
       f1d => this%data(:,j)
       f1d = this%interp1%interpolate_array_disp(this%n1, f1d, displacement)
    end do
  end subroutine advection_x

  subroutine advection_v(this, dt)
    class(vp_cartesian_2d) :: this 
    sll_real64 :: dt
    ! local variables
    sll_real64, dimension(:), pointer :: f1d
    sll_real64 :: displacement
    sll_int32 :: i,j
    
    ! compute electric field
    rho = 1.0_f64 - this%dist_func%compute_rho()
    call solve(poisson_1d, this%efield, this%rho)
     if (driven) then
        call PFenvelope(adr, istep*dt, tflat, tL, tR, twL, twR, &
             t0, turn_drive_off)
        do i = 1, Ncx + 1
           E_app(i) = Edrmax * adr * kmode * sin(kmode * (i-1) * delta_x - omegadr*istep*dt)
        enddo
     endif
    ! do advection for given electric field
    do i = 1, Ncx+1
        displacement = -(efield(i)+e_app(i)) * 0.5_f64 * dt
        f1d => FIELD_DATA(f) (i,:) 
        f1d = interp_v%interpolate_array_disp(Ncv+1, f1d, displacement)
        if (is_delta_f) then
           ! add equilibrium contribution
           do j=1, Ncv + 1
              v = vmin + (j-1) * delta_v
              f1d(j) = f1d(j) + f_equilibrium(v-alpha) - f_equilibrium(v)
           end do
        end if
     end do
  end subroutine advection_v
end module vp_cartesian_2d
