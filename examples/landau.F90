!> @brief This is an example of selalib 
!> @details We use selalib to solve the Vlasov-Poisson system
!> for the Landau test case
!>
!> \f[ f(x,v,t=0) = (\epsilon + \sin(t)) \exp^{k_x} \f]
!>
program landau_1d
#include "selalib.h"

implicit none
  
integer :: i_step, n_step, i, j, j_step
integer :: nc_eta1, nc_eta2

real(8) :: eta1_min, eta1_max, delta_eta1
real(8) :: eta2_min, eta2_max, delta_eta2
real(8) :: delta_t, time

real(8), dimension(:),   allocatable :: nrj
real(8), dimension(:,:), allocatable :: df

integer  :: error

real(8), dimension(:), allocatable :: eta1
real(8), dimension(:), allocatable :: eta2

class(sll_interpolator_1d_base), pointer :: interp_x
class(sll_interpolator_1d_base), pointer :: interp_v

type(cubic_spline_1d_interpolator), target :: spline_x
type(cubic_spline_1d_interpolator), target :: spline_v

real(8), dimension(:), allocatable :: ex
real(8), dimension(:), allocatable :: rho
real(8)                            :: eps
real(8)                            :: kx
type (poisson_1d_periodic)         :: poisson

interp_x => spline_x
interp_v => spline_v

eta1_min = 0.0_f64
eta1_max = 4.0_f64*sll_pi
nc_eta1  = 256

eta2_min = -6.0_f64
eta2_max = +6.0_f64
nc_eta2  = 256

allocate(rho(nc_eta1+1))
allocate(ex(nc_eta1+1))
allocate(df(nc_eta1+1,nc_eta2+1))

delta_eta1 = (eta1_max-eta1_min)/nc_eta1
delta_eta2 = (eta2_max-eta2_min)/nc_eta2

call initialize(poisson, eta1_min, eta2_max, nc_eta1, error) 
call spline_x%initialize(nc_eta1+1, eta1_min, eta1_max, SLL_PERIODIC )
call spline_v%initialize(nc_eta2+1, eta2_min, eta2_max, SLL_HERMITE )

allocate(eta1(nc_eta1+1))
allocate(eta2(nc_eta2+1))

!Initialize distribution function
eps = 0.05_f64
kx  = 0.5_f64 
do j=1, nc_eta2+1
   eta2(j) = eta2_min + (j-1)*delta_eta2
   do i=1, nc_eta1+1
      eta1(i) = eta1_min + (i-1)*delta_eta1
      df(i,j)=(1.0_f64+eps*cos(kx*eta1(i)))/(2.0_f64*sll_pi) &
                * exp(-0.5_f64*eta2(j)*eta2(j))
   end do
end do

time = 0.0_f64
n_step = 1000
delta_t = 0.05_f64
allocate(nrj(n_step))
print'(a)', 'set term x11'

call advection_x(0.5*delta_t)

do i_step = 1, n_step

   do i = 1, nc_eta1+1
      rho(i) = sum(df(i,:))*delta_eta2
   end do

   call solve(poisson, ex , rho)

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
      print'(2e15.3)', (j_step-1)*delta_t, nrj(j_step)
   end do
   print'(a)','e'

end do

100 format('p [',f5.1,':',f5.1,'][',f6.1,':',f6.1,'] ''-'' w l')

contains

   subroutine advection_x(dt)
    real(8), intent(in) :: dt
    do j = 1, nc_eta2+1
      df(:,j) = interp_x%interpolate_array_disp(nc_eta1+1,df(:,j),dt*eta2(j))
    end do
   end subroutine advection_x

   subroutine advection_v(dt)
    real(8), intent(in) :: dt
    do i = 1, nc_eta1+1
      df(i,:) = interp_v%interpolate_array_disp(nc_eta2+1,df(i,:),dt*ex(i))
    end do
   end subroutine advection_v

end program landau_1d
