!
!  Contact : Pierre Navaro http://wwww-irma.u-strasbg.fr/~navaro
!
!-------------------------------------------------------------------
!  test 1D Ampere Vlasov solver based on FFT
!-------------------------------------------------------------------
program test_ampere_vlasov_1d
#include "sll_working_precision.h"
#include "sll_memory.h"
#include "sll_assert.h"
#include "sll_constants.h"
#include "sll_file_io.h"
#include "sll_interpolators.h"
#include "sll_poisson_solvers.h"

use sll_module_ampere_vlasov_1d
use sll_poisson_1d_periodic  
use sll_module_poisson_1d_periodic_solver

implicit none

sll_real64 :: x_max, x_min, delta_x
sll_int32  :: nc_x
sll_int32  :: error

type(sll_ampere_1d)        :: ampere
type(sll_ampere_vlasov_1d) :: ampere_vlasov
sll_int32                  :: i
sll_int32                  :: j
sll_real64                 :: w
sll_real64                 :: time
sll_int32                  :: istep, nstep
sll_real64                 :: dt
sll_real64                 :: cfl = 0.5

sll_int32,  parameter      :: mode = 2

sll_real64, dimension(:), allocatable :: x
sll_real64, dimension(:), allocatable :: ex
sll_real64, dimension(:), allocatable :: ex_exact
sll_real64, dimension(:), allocatable :: jx
  
sll_int32  :: i_step, n_step, j_step
sll_int32  :: nc_eta1, nc_eta2

sll_real64 :: eta1_min, eta1_max, delta_eta1
sll_real64 :: eta2_min, eta2_max, delta_eta2
sll_real64 :: delta_t

sll_real64, dimension(:), allocatable :: nrj
sll_real64, dimension(:,:), allocatable :: df

sll_real64, dimension(:)   , allocatable :: eta1, eta2

class(sll_interpolator_1d_base), pointer    :: interp_x
class(sll_interpolator_1d_base), pointer    :: interp_v

type(sll_cubic_spline_interpolator_1d), target  :: spline_x
type(sll_cubic_spline_interpolator_1d), target  :: spline_v

sll_real64, dimension(:), allocatable :: rho
class(sll_poisson_1d_base), pointer   :: poisson 
sll_real64 :: eps, kx, t

sll_real64 :: tstart, tend

call cpu_time(tstart)

x_min =  0.0_f64; x_max = 1.0_f64

nc_x = 64

delta_x = (x_max-x_min)/nc_x

SLL_ALLOCATE(x(nc_x+1),error)

do i = 1, nc_x+1
  x(i) = x_min + (i-1)*delta_x
end do

dt = cfl  / sqrt(1./(delta_x*delta_x))
nstep = 1000

time  = 0.

w = 1.0_f64

SLL_CLEAR_ALLOCATE(ex(1:nc_x+1),       error)
SLL_CLEAR_ALLOCATE(jx(1:nc_x+1),       error)
SLL_CLEAR_ALLOCATE(ex_exact(1:nc_x+1), error)

call sll_create(ampere, x_min, x_max, nc_x)

do istep = 1, nstep !*** Loop over time

  ex_exact = w*sin(time*w)*sin(sll_pi*x)
  time     = time + 0.5_f64*dt
  jx       = -w*w*cos(time*w)*sin(sll_pi*x)

  !write(*,"(10x,' istep = ',I6)",advance="no") istep
  !write(*,"(' time = ',g12.3,' sec')",advance="no") time
  !write(*,"(' error max = ',g15.5)") maxval(abs(ex - ex_exact))

  call sll_solve(ampere, dt, jx, ex)
  ex(nc_x+1) = ex(1)
  
  call sll_gnuplot_1d( 'ex', ex, ex_exact, istep)

  time = time + 0.5_f64*dt

end do ! next time step

deallocate(ex)

interp_x => spline_x
interp_v => spline_v

eta1_min = 0.0_f64
eta1_max = 4.0_f64*sll_pi
nc_eta1  = 64

eta2_min = -6.0_f64
eta2_max = +6.0_f64
nc_eta2  = 64

SLL_ALLOCATE(rho(nc_eta1+1),error)
SLL_ALLOCATE(ex(nc_eta1+1),error)
SLL_ALLOCATE(df(nc_eta1+1,nc_eta2+1),error)

delta_eta1 = (eta1_max-eta1_min)/nc_eta1
delta_eta2 = (eta2_max-eta2_min)/nc_eta2

poisson => new_poisson_1d_periodic_solver(eta1_min,eta1_max,nc_eta1)

call spline_x%initialize(nc_eta1+1, eta1_min, eta1_max, SLL_PERIODIC )
call spline_v%initialize(nc_eta2+1, eta2_min, eta2_max, SLL_PERIODIC )

SLL_ALLOCATE(eta1(nc_eta1+1),error)
SLL_ALLOCATE(eta2(nc_eta2+1),error)

!Initialize distribution function
eps = 0.001_f64
kx  = 0.5_f64 
do j=1, nc_eta2+1
   eta2(j) = eta2_min + (j-1)*delta_eta2
   do i=1, nc_eta1+1
      eta1(i) = eta1_min + (i-1)*delta_eta1
      df(i,j)=(1.0_f64+eps*cos(kx*eta1(i)))/sqrt(2.0_f64*sll_pi) &
                * exp(-0.5_f64*eta2(j)*eta2(j))
   end do
end do

time = 0.0_f64
n_step = 600
delta_t = 0.1_f64
SLL_ALLOCATE(nrj(n_step), error)
print*, "set title'", delta_eta1, delta_eta2,"'"
print*, 'set term x11'
print*, 'set log y'

call sll_create(ampere_vlasov,               &
                eta1_min, eta1_max, nc_eta1, &
                eta2_min, eta2_max, nc_eta2)
              
call advection_poisson(0.5*delta_t)

do i_step = 1, n_step

   call advection_v(delta_t)

   call advection_ampere(delta_t) 
   !call advection_poisson(delta_t)
   
   nrj(i_step) = 0.5_f64*sum(ex(1:nc_eta1)*ex(1:nc_eta1))*delta_eta1
   
   open(11, file='thf_2d.dat', position='append')
   if (i_step == 1) rewind(11)
   write(11,*) time, sqrt(nrj(i_step))
   close(11)

   time = time + delta_t

   !write(*,100) .0,n_step*delta_t,-35.,-10.
   !do j_step = 1, i_step
   !   t =  (j_step-1)*delta_t
   !   print*, t, log(nrj(j_step))
   !end do
   !print*, 'e'

end do

print*,"plot (2.6e-3*exp(-0.1533*x)),'thf_2d.ref' w l, 'thf_2d.dat' w l"
!100 format('p [t=',f5.1,':',f5.1,'][',f6.1,':',f6.1 &
!           ,'] 2.6e-3*exp(-0.1533*t), ''thf_2d.dat'' w l')


call cpu_time(tend)
print"('#CPU time : ',g15.3)", tend-tstart
print*,'#PASSED'

deallocate(ex,jx,ex_exact)

contains

   subroutine advection_poisson(dt)
    sll_real64, intent(in) :: dt
    do j = 1, nc_eta2+1
      df(:,j) = interp_x%interpolate_array_disp(nc_eta1+1,df(:,j),dt*eta2(j))
    end do
    do i = 1, nc_eta1+1
      rho(i) = sum(df(i,:))*delta_eta2
    end do
    !rho = 1.0_f64 - rho
    call poisson%compute_e_from_rho( ex, rho )
   end subroutine advection_poisson

   subroutine advection_ampere(dt)
   sll_real64, intent(in) :: dt

   call sll_solve(ampere_vlasov, dt, df, ex) 

   end subroutine advection_ampere

   subroutine advection_v(dt)
    sll_real64, intent(in) :: dt
    do i = 1, nc_eta1+1
      df(i,:) = interp_v%interpolate_array_disp(nc_eta2+1,df(i,:),dt*ex(i))
    end do
   end subroutine advection_v

end program test_ampere_vlasov_1d
