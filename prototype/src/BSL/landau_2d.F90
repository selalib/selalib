program landau_4d

#include "sll_assert.h"
#include "sll_working_precision.h"
#include "sll_memory.h"
#include "sll_poisson_solvers.h"

use sll_constants
use sll_module_interpolators_1d_base
use sll_cubic_spline_interpolator_1d
use sll_utilities, only: int2string

implicit none
  
!Geometry
sll_real64 :: eta1, eta2, eta3, eta4
sll_real64 :: eta1_min, eta1_max, eta2_min, eta2_max
sll_real64 :: eta3_min, eta3_max, eta4_min, eta4_max

!Grids
sll_int32  :: nc_eta1, nc_eta2, nc_eta3, nc_eta4
sll_real64 :: delta_eta1, delta_eta2, delta_eta3, delta_eta4

!Time domain
sll_int32  :: i_step, n_step
sll_real64 :: delta_t
sll_real64 :: time

!Distribution function 4D
sll_real64 :: eps, kx, ky, v2
sll_real64, dimension(:,:,:,:), allocatable :: f

!Electric fields and charge density
sll_real64, dimension(:,:), allocatable :: ex
sll_real64, dimension(:,:), allocatable :: ey
sll_real64, dimension(:,:), allocatable :: phi
sll_real64, dimension(:,:), allocatable :: rho

type(poisson_2d_periodic)               :: poisson

class(sll_interpolator_1d_base), pointer    :: interp_1
class(sll_interpolator_1d_base), pointer    :: interp_2
class(sll_interpolator_1d_base), pointer    :: interp_3
class(sll_interpolator_1d_base), pointer    :: interp_4

type(cubic_spline_1d_interpolator), target  :: spl_eta1
type(cubic_spline_1d_interpolator), target  :: spl_eta2
type(cubic_spline_1d_interpolator), target  :: spl_eta3
type(cubic_spline_1d_interpolator), target  :: spl_eta4

!Diagnostics and errors
sll_int32                             :: error
sll_real64, dimension(:), allocatable :: nrj

!Local indices
sll_int32  :: i1, i2, i3, i4

!x domain
eta1_min =  0.0_f64; eta1_max =  4.0_f64 * sll_pi
eta2_min =  0.0_f64; eta2_max =  4.0_f64 * sll_pi

nc_eta1 = 31; nc_eta2 = 31

delta_eta1 = (eta1_max-eta1_min)/nc_eta1
delta_eta2 = (eta2_max-eta2_min)/nc_eta2

!v domain
eta3_min = -6.0_f64; eta3_max = 6.0_f64 
eta4_min = -6.0_f64; eta4_max = 6.0_f64 

nc_eta3 = 31; nc_eta4 = 31

delta_eta3 = (eta3_max-eta3_min)/nc_eta3
delta_eta4 = (eta4_max-eta4_min)/nc_eta4

SLL_ALLOCATE(f(nc_eta1+1,nc_eta2+1,nc_eta3+1,nc_eta4+1),error)
SLL_ALLOCATE(phi(nc_eta1+1,nc_eta2+1),error)
SLL_ALLOCATE(rho(nc_eta1+1,nc_eta2+1),error)
SLL_ALLOCATE(ex(nc_eta1+1,nc_eta2+1),error)
SLL_ALLOCATE(ey(nc_eta1+1,nc_eta2+1),error)


call initialize(poisson, eta1_min, eta1_max, nc_eta1, &
                         eta2_min, eta2_max, nc_eta2, error)

call spl_eta1%initialize(nc_eta1+1, eta1_min, eta1_max, SLL_PERIODIC )
call spl_eta2%initialize(nc_eta2+1, eta2_min, eta2_max, SLL_PERIODIC )
call spl_eta3%initialize(nc_eta3+1, eta3_min, eta3_max, SLL_PERIODIC )
call spl_eta4%initialize(nc_eta4+1, eta4_min, eta4_max, SLL_PERIODIC )

interp_1 => spl_eta1
interp_2 => spl_eta2
interp_3 => spl_eta3
interp_4 => spl_eta4

eps = 0.05
kx  = 2*sll_pi/(eta1_max-eta1_min)
ky  = 2*sll_pi/(eta2_max-eta2_min)

eta4 = eta4_min
do i4=1,nc_eta4+1
   eta3 = eta3_min
   do i3=1,nc_eta3+1
      eta2 = eta2_min
      v2 = eta3*eta3+eta4*eta4
      do i2=1,nc_eta2+1
         eta1 = eta1_min
         do i1=1,nc_eta1+1
            f(i1,i2,i3,i4)=(1+eps*cos(kx*eta1)*cos(ky*eta2))/(2*sll_pi)*exp(-.5*v2)
            eta1 = eta1 + delta_eta1
         end do
         eta2 = eta2 + delta_eta2
      end do
      eta3 = eta3 + delta_eta3
   end do
   eta4 = eta4 + delta_eta4
end do

n_step = 1000
SLL_CLEAR_ALLOCATE(nrj(1:n_step), error)
delta_t = .01_f64
time = 0.0_f32

if ( delta_t > 0.5/sqrt(1./(delta_eta1*delta_eta1)+1./(delta_eta2*delta_eta2))) &
  stop 'Warning CFL'


call advection_x1(0.5_f64*delta_t)
call advection_x2(0.5_f64*delta_t)

do i_step = 1, n_step !Loop over time

   time  = time + 0.5 * delta_t

   call compute_rho()
   call solve(poisson,ex,ey,rho,nrj(i_step))
   call online_plot() 


   call advection_v1(delta_t)
   call advection_v2(delta_t)

   if (i_step == 1 .or. mod(i_step, 10) == 0) then
      call plot_field(ex,"ex",i_step/10)
      call plot_field(rho,"rho",i_step/10)
      !call plot_field(f(:,:,1,1),"f",i_step/10)
   end if

   time  = time + 0.5 * delta_t

   call advection_x1(delta_t)
   call advection_x2(delta_t)


end do !next time step


contains

subroutine compute_rho()

   do i2 = 1, nc_eta2+1
      do i1 = 1, nc_eta1+1
         rho(i1,i2) = sum(f(i1,i2,:,:))
      end do
   end do
   rho = rho * delta_eta3 * delta_eta4

end subroutine compute_rho

subroutine online_plot()

   sll_real64 :: cell_volume  
   open(11, file='thf_4d.dat', position='append')
   if (i_step == 1) rewind(11)
   write(11,*) time, nrj(i_step)
   close(11)

   cell_volume = delta_eta1 * delta_eta2 * delta_eta3 * delta_eta4
   print"(5f10.3)", time, nrj(i_step), &
                   cell_volume * sum(f), &
                   cell_volume * sum(abs(f)), &
                   cell_volume * sum(f*f) 

   !write(*,100) .0,10.,-15.,1.
   !do j_step = 1, i_step
   !   print*, (j_step-1)*delta_t, nrj(j_step)
   !end do
   !print*, 'e'

   !100 format('p [',f5.1,':',f5.1,'][',f6.1,':',f6.1,'] ''-'' w l')

end subroutine online_plot

subroutine advection_x1(dt)

   sll_real64, intent(in) :: dt
   sll_real64             :: eta3

   do i4 = 1, nc_eta4+1
   eta3 = eta3_min
   do i3 = 1, nc_eta3+1
   do i2 = 1, nc_eta2+1
      f(:,i2,i3,i4) = interp_1%interpolate_array_disp(nc_eta1+1, &
                      f(:,i2,i3,i4),dt*eta3)
   end do
   eta3 = eta3 + delta_eta3
   end do
   end do

end subroutine advection_x1

subroutine advection_x2(dt)

   sll_real64, intent(in) :: dt
   sll_real64             :: eta4

   eta4 = eta4_min
   do i4 = 1, nc_eta4+1
   do i3 = 1, nc_eta3+1
   do i1 = 1, nc_eta1+1
            f(i1,:,i3,i4) = interp_2%interpolate_array_disp(nc_eta2+1, &
                            f(i1,:,i3,i4),dt*eta4)
   end do
   end do
   eta4 = eta4 + delta_eta4
   end do


end subroutine advection_x2

subroutine advection_v1(dt)
   sll_real64, intent(in) :: dt

   do i4 = 1, nc_eta4+1
   do i2 = 1, nc_eta2+1
   do i1 = 1, nc_eta1+1
      f(i1,i2,:,i4) = interp_3%interpolate_array_disp(nc_eta3+1, &
                      f(i1,i2,:,i4),ex(i1,i2)*dt)
   end do
   end do
   end do

end subroutine advection_v1

subroutine advection_v2(dt)

   sll_real64, intent(in) :: dt

   do i3 = 1, nc_eta3+1
   do i2 = 1, nc_eta2+1
   do i1 = 1, nc_eta1+1
      f(i1,i2,i3,:) = interp_4%interpolate_array_disp(nc_eta4+1, &
                      f(i1,i2,i3,:),ey(i1,i2)*dt)
   end do
   end do
   end do

end subroutine advection_v2

subroutine plot_field(f, fname, iplot)

   sll_int32 :: iplot, i, j
   sll_real64 :: x, y
   sll_real64, dimension(:,:) :: f
   character(len=*) :: fname
   character(len=4) :: cplot
 
   call int2string(iplot,cplot)

   open(11, file=fname//cplot//".dat")
   do i = 1, size(f,1)
      do j = 1, size(f,2)
         x = eta1_min + (i-1)*delta_eta1
         y = eta2_min + (j-1)*delta_eta2
         write(11,*) x,y,f(i,j)
      end do
      write(11,*)
   end do
   close(11)
   
   open( 90, file = fname//'.gnu', position="append" )
   if ( iplot == 1 ) then
      rewind(90)
      !write(90,*)"set cbrange[-1:1]"
      !write(90,*)"set pm3d"
      write(90,*)"set surf"
      write(90,*)"set term x11"
   end if

   write(90,*)"set title 'step = ",iplot,"'"
   write(90,"(a)")"splot '"//fname//cplot//".dat' u 1:2:3 w lines"
   close(90)

end subroutine plot_field

end program landau_4d
