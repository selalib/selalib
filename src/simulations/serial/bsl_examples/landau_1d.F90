program landau_1d
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#include "sll_memory.h"
#include "sll_working_precision.h"

  use sll_m_boundary_condition_descriptors, only: &
    sll_p_periodic

  use sll_m_constants, only: &
    sll_p_pi

  use sll_m_cubic_spline_interpolator_1d, only: &
    sll_t_cubic_spline_interpolator_1d

  use sll_m_interpolators_1d_base, only: &
    sll_c_interpolator_1d

  use sll_m_poisson_1d_periodic, only: &
    sll_o_initialize, &
    sll_t_poisson_1d_periodic, &
    sll_o_solve

  implicit none
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  
sll_int32  :: i_step, n_step, i, j, j_step
sll_int32  :: nc_eta1, nc_eta2

sll_real64 :: eta1_min, eta1_max, delta_eta1
sll_real64 :: eta2_min, eta2_max, delta_eta2
sll_real64 :: delta_t, time

sll_real64, dimension(:), allocatable :: nrj
sll_real64, dimension(:,:), allocatable :: df

sll_int32  :: error

sll_real64, dimension(:)   , allocatable :: eta1, eta2

class(sll_c_interpolator_1d), pointer    :: interp_x
class(sll_c_interpolator_1d), pointer    :: interp_v

type(sll_t_cubic_spline_interpolator_1d), target  :: spline_x
type(sll_t_cubic_spline_interpolator_1d), target  :: spline_v

sll_real64, dimension(:), allocatable :: ex
sll_real64, dimension(:), allocatable :: rho
type (sll_t_poisson_1d_periodic)            :: poisson
sll_real64 :: eps, kx

interp_x => spline_x
interp_v => spline_v

eta1_min = 0.0_f64
eta1_max = 4.0_f64*sll_p_pi
nc_eta1  = 256

eta2_min = -6.0_f64
eta2_max = +6.0_f64
nc_eta2  = 256

SLL_ALLOCATE(rho(nc_eta1+1),error)
SLL_ALLOCATE(ex(nc_eta1+1),error)
SLL_ALLOCATE(df(nc_eta1+1,nc_eta2+1),error)

delta_eta1 = (eta1_max-eta1_min)/nc_eta1
delta_eta2 = (eta2_max-eta2_min)/nc_eta2

call sll_o_initialize(poisson, eta1_min, eta2_max, nc_eta1, error) 
call spline_x%initialize(nc_eta1+1, eta1_min, eta1_max, sll_p_periodic )
call spline_v%initialize(nc_eta2+1, eta2_min, eta2_max, sll_p_periodic )

SLL_ALLOCATE(eta1(nc_eta1+1),error)
SLL_ALLOCATE(eta2(nc_eta2+1),error)

!sll_o_initialize distribution function
eps = 0.05_f64
kx  = 0.5_f64 
do j=1, nc_eta2+1
   eta2(j) = eta2_min + (j-1)*delta_eta2
   do i=1, nc_eta1+1
      eta1(i) = eta1_min + (i-1)*delta_eta1
      df(i,j)=(1.0_f64+eps*cos(kx*eta1(i)))/(2.0_f64*sll_p_pi) &
                * exp(-0.5_f64*eta2(j)*eta2(j))
   end do
end do

time = 0.0_f64
n_step = 1000
delta_t = 0.05_f64
SLL_ALLOCATE(nrj(n_step), error)
print*, "set title'", delta_eta1, delta_eta2,"'"
print*, 'set term x11'

call advection_x(0.5*delta_t)

do i_step = 1, n_step

   do i = 1, nc_eta1+1
      rho(i) = sum(df(i,:))*delta_eta2
   end do

   call sll_o_solve(poisson, ex , rho)

   call advection_v(delta_t)
   call advection_x(delta_t)
   
   nrj(i_step) = 0.5_f64*log(sum(ex*ex)*delta_eta1)
   
   open(11, file='thf_2d.dat', position='append')
   if (i_step == 1) rewind(11)
   write(11,*) time, nrj(i_step)
   close(11)

   time = time + delta_t

   write(*,100) .0,n_step*delta_t,-29.5,0.5
   do j_step = 1, i_step
      print*, (j_step-1)*delta_t, nrj(j_step)
   end do
   print*, 'e'

end do

print*,'PASSED'
100 format('p [',f5.1,':',f5.1,'][',f6.1,':',f6.1,'] ''-'' w l')

contains

   subroutine advection_x(dt)
    sll_real64, intent(in) :: dt
    do j = 1, nc_eta2+1
      call interp_x%interpolate_array_disp_inplace(nc_eta1+1,df(:,j),-dt*eta2(j))
    end do
   end subroutine advection_x

   subroutine advection_v(dt)
    sll_real64, intent(in) :: dt
    do i = 1, nc_eta1+1
      call interp_v%interpolate_array_disp_inplace(nc_eta2+1,df(i,:),-dt*ex(i))
    end do
   end subroutine advection_v

end program landau_1d
