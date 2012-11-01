program unit_test_1d
#include "sll_working_precision.h"
#include "sll_memory.h"
#include "sll_assert.h"

  use numeric_constants
  use sll_module_interpolators_1d_base
  use sll_cubic_spline_interpolator_1d

  implicit none

  sll_int32  :: nc_x, nc_v
  sll_int32  :: i, j, it, n_steps
  sll_real64 :: x_min, x_max, v_min, v_max
  sll_real64 :: delta_t, error
  sll_int32  :: info

  sll_real64, dimension(:)   , allocatable :: x, dx
  sll_real64, dimension(:)   , allocatable :: v, dv
  sll_real64, dimension(:,:) , allocatable :: df

  sll_real64 :: advfield_x, advfield_v

  class(sll_interpolator_1d_base), pointer    :: interp_x
  class(sll_interpolator_1d_base), pointer    :: interp_v

  type(cubic_spline_1d_interpolator), target  :: spline_x
  type(cubic_spline_1d_interpolator), target  :: spline_v

  print*,'*********************'
  print*,' 1D case             '
  print*,' 1D in x and 1D in v '
  print*,'*********************'

  print*, 'set domain size'
  x_min =  -5.0_f64; x_max =  5.0_f64
  v_min =  -5.0_f64; v_max =  5.0_f64 

  nc_x = 100; nc_v = 100

  SLL_ALLOCATE(x(nc_x), info)
  SLL_ALLOCATE(v(nc_v), info)

  do i = 1, nc_x
     x(i) = x_min + (i-1)*(x_max-x_min)/(nc_x-1)
  end do

  do j = 1, nc_v
     v(j) = v_min + (j-1)*(v_max-v_min)/(nc_v-1)
  end do

  SLL_ALLOCATE(df(nc_x,nc_v), info)

  do j = 1, nc_v
     do i = 1, nc_x
        df(i,j) =  exp(-(x(i)**2+v(j)**2))
     end do
  end do

  advfield_x = 1_f64 
  advfield_v = 0.0 

  print*, 'initialize 2d distribution function f(x,v) gaussian'

  Print*, 'checking advection of a Gaussian in a uniform field'
  
  call spline_x%initialize(nc_x, x_min, x_max, PERIODIC_SPLINE )
  call spline_v%initialize(nc_v, v_min, v_max, PERIODIC_SPLINE )

  interp_x => spline_x
  interp_v => spline_v

  SLL_ALLOCATE(dx(nc_x),info)
  SLL_ALLOCATE(dv(nc_v),info)

  ! run BSL method using 10 time steps
  n_steps = 100
  delta_t = 10.0_f64/n_steps
  do it = 1, n_steps

     !call plot_df( df, it )

     call advection_x(df, interp_x, 0.5*delta_t)
     call advection_v(df, interp_v, delta_t)
     call advection_x(df, interp_x, 0.5*delta_t)

     ! compute error when Gaussian arrives at center (t=1)

  end do

  do j = 1, nc_v
     do i = 1, nc_x
        error = max(error,abs(df(i,j)-exp(-(x(i)**2+v(j)**2))))
     end do
  end do

  print*, ' 100 nodes, ', it, ' time steps. Error= ', error

  print *, 'Successful, exiting program.'
  print *, 'PASSED'

contains

   subroutine advection_x(df, interp_x, dt)
   class(sll_interpolator_1d_base), pointer  :: interp_x
   sll_real64, intent(inout), dimension(:,:) :: df
   sll_real64, intent(in) :: dt

     do i = 1, nc_x
        dx(i) = x(1) + modulo(x(i)-x(1)-dt*advfield_x,x(nc_x)-x(1))
     end do

     do j = 1, nc_v
        df(:,j) = interp_x%interpolate_array( nc_x, df(:,j), dx )
     end do

   end subroutine advection_x

   subroutine advection_v(df, interp_v, dt)
   class(sll_interpolator_1d_base), pointer    :: interp_v
   sll_real64, intent(inout), dimension(:,:) :: df
   sll_real64, intent(in) :: dt

     do j = 1, nc_v
        dv(j) = v(1) + modulo(v(j)-v(1)-dt*advfield_v,v(nc_v)-v(1))
     end do

     do i = 1, nc_x
        df(i,:) = interp_v%interpolate_array( nc_v, df(i,:), dv )
     end do

   end subroutine advection_v


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine plot_df( df, iplot )

sll_real64, dimension(:,:), intent(in) :: df
integer :: iplot, i, j
integer :: kk0, kk1, kk2, kk3, kk4
character(len=4) :: fin
character(len=1) :: aa,bb,cc,dd

kk0 = iplot
kk1 = kk0/1000
aa  = char(kk1 + 48)
kk2 = (kk0 - kk1*1000)/100
bb  = char(kk2 + 48)
kk3 = (kk0 - (kk1*1000) - (kk2*100))/10
cc  = char(kk3 + 48)
kk4 = (kk0 - (kk1*1000) - (kk2*100) - (kk3*10))/1
dd  = char(kk4 + 48)
fin = aa//bb//cc//dd

open(11, file="df-"//fin//".dat")
do i = 1, size(x)
   do j = 1, size(v)
      write(11,*) x(i),v(j),df(i,j)
   end do
   write(11,*)
end do
close(11)
   
open( 90, file = 'df.gnu', position="append" )
if ( iplot == 1 ) then
   rewind(90)
   !write(90,*)"set cbrange[-1:1]"
   !write(90,*)"set pm3d"
   write(90,*)"set surf"
   write(90,*)"set term x11"
end if

write(90,*)"set title 'step = ",iplot,"'"
write(90,"(a)")"splot 'df-"//fin//".dat' u 1:2:3 w lines"
close(90)

end subroutine plot_df

end program unit_test_1d
