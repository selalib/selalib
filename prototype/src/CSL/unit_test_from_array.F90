program unit_test_from_array
#include "sll_working_precision.h"
#include "sll_mesh_types.h"
#include "sll_memory.h"
  use sll_constants
  use distribution_function
  use sll_diagnostics
  use sll_csl
  implicit none

  sll_int32 :: nc_eta1, nc_eta2
  sll_int32 :: i1, i2, it, n_steps, ierr
  sll_real64 :: eta1_min, eta1_max,  eta2_min, eta2_max, eta1, eta2, eta1c, eta2c
  sll_real64 :: deltat, val, error, error1 
  sll_real64 :: alpha1, alpha2
  sll_real64 :: r_min, r_max, delta_r, delta_theta, delta_eta1, delta_eta2
  sll_real64, dimension(:,:), pointer :: x1c_array, x2c_array, jac_array
  sll_real64, dimension(:,:), pointer :: x1n_array, x2n_array
  sll_int32 :: k1, k2,visu_step
  sll_real64 :: xx1, xx2
  type(geometry_2D), pointer :: geom
  type(mesh_descriptor_2D), pointer :: mesh
  type(sll_distribution_function_2D_t), pointer :: dist_func
  type(field_2D_vec1), pointer :: rotating_field
  type(field_2D_vec1), pointer :: uniform_field
  type(csl_workspace), pointer :: csl_work
  character(32), parameter  :: name = 'distribution_function'
  logical, parameter :: read_from_file = .false.

  r_min = 0._f64
  r_max = 8.0_f64
  eta1_min =  0.0_f64
  eta1_max =  1.0_f64
  eta2_min =  0.0_f64
  eta2_max =  1.0_f64
  nc_eta1 = 100
  nc_eta2 = 100
  delta_r = (r_max-r_min) / nc_eta1
  delta_theta = 2*sll_pi / nc_eta2
  delta_eta1 = 1.0_f64 / nc_eta1
  delta_eta2 = 1.0_f64 / nc_eta2
  SLL_ALLOCATE(x1n_array(nc_eta1+1, nc_eta2+1), ierr)
  SLL_ALLOCATE(x2n_array(nc_eta1+1, nc_eta2+1), ierr)
  SLL_ALLOCATE(x1c_array(nc_eta1+1, nc_eta2+1), ierr)
  SLL_ALLOCATE(x2c_array(nc_eta1+1, nc_eta2+1), ierr)
  SLL_ALLOCATE(jac_array(nc_eta1+1, nc_eta2+1), ierr)
  if (read_from_file) then
     ! read array from file
     open(1, FILE='../mesh_types/xxvv_node.dat')
     open(2, FILE='../mesh_types/xxvv_cell.dat') 
     open(3, FILE='../mesh_types/jac_cell.dat') 
     do i1=1, nc_eta1+1
        do i2=1, nc_eta2+1
           read(1,*) k1, k2, xx1, xx2 
           x1n_array(k1,k2) = xx1
           x2n_array(k1,k2) = xx2
!           if (k2 == 1) then
!              x1n_array(nc_eta1+1,k2) = xx1
!              x2n_array(nc_eta1+1,k2) = xx2
!           endif
        end do
        
     end do
     do i1=1, nc_eta1
        do i2=1, nc_eta2
           read(2,*) k1, k2, xx1, xx2 
           x1c_array(k1,k2) = xx1
           x2c_array(k1,k2) = xx2
           read(3,*) k1, k2, xx1
           jac_array(k1,k2) = xx1
        end do
     end do
  else 
     ! test polar coordinates
     eta2 = 0.0_f64 
     !eta2c = eta2 + 0.5_f64*delta_theta
     eta2c = eta2 + 0.5_f64*delta_eta2
     do i2= 1, nc_eta2 + 1
        !eta1 = r_min
        !eta1c = eta1 + 0.5_f64*delta_r
        eta1 = 0.0
        eta1c = 0.5_f64*delta_eta1
        do i1 = 1, nc_eta1 + 1
           !x1n_array(i1,i2) = eta1 * cos(eta2)
           !x2n_array(i1,i2) = eta1 * sin(eta2)
           !x1c_array(i1,i2) = eta1c * cos(eta2c)
           !x2c_array(i1,i2) = eta1c * sin(eta2c)
           !jac_array(i1,i2) = eta1c
           x1n_array(i1,i2) = eta1 + 0.1_f64 * sin(2*sll_pi*eta1) * sin(2*sll_pi*eta2)
           x2n_array(i1,i2) = eta2 + 0.1_f64 * sin(2*sll_pi*eta1) * sin(2*sll_pi*eta2)
           x1c_array(i1,i2) = eta1c + 0.1_f64 * sin(2*sll_pi*eta1c) * sin(2*sll_pi*eta2c)
           x2c_array(i1,i2) = eta2c + 0.1_f64 * sin(2*sll_pi*eta1c) * sin(2*sll_pi*eta2c)
           jac_array(i1,i2) = (1.0_f64 + 0.2_f64 *sll_pi * cos (2*sll_pi*eta1c) * sin (2*sll_pi*eta2c)) * &
         (1.0_f64 + 0.2_f64 * sll_pi * sin (2*sll_pi*eta1c) * cos (2*sll_pi*eta2c)) - &
         0.2_f64 *sll_pi * sin (2*sll_pi*eta1c) * cos (2*sll_pi*eta2c) * &
         0.2_f64 * sll_pi * cos (2*sll_pi*eta1c) * sin (2*sll_pi*eta2c)
           !eta1 = eta1 + delta_r
           !eta1c = eta1c + delta_r
           eta1 = eta1 + delta_eta1
           eta1c = eta1c + delta_eta1
        end do
        !eta2 = eta2 + delta_theta
        !eta2c = eta2c + delta_theta
        eta2 = eta2 + delta_eta2
        eta2c = eta2c + delta_eta2
     end do
  endif

! geom => new_geometry_2D ('from_array',nc_eta1+1,nc_eta2+1, &
!       x1n_array, x2n_array, x1c_array, x2c_array, jac_array,PERIODIC,COMPACT)
!  geom => new_geometry_2D ('from_array',nc_eta1+1,nc_eta2+1, &
!       x1n_array, x2n_array, x1c_array, x2c_array, jac_array,COMPACT, PERIODIC)
  geom => new_geometry_2D ('from_array',nc_eta1+1,nc_eta2+1, &
       x1n_array, x2n_array, x1c_array, x2c_array, jac_array,PERIODIC, PERIODIC)


!  mesh => new_mesh_descriptor_2D(eta1_min, eta1_max, nc_eta1, &
!       COMPACT, eta2_min, eta2_max, nc_eta2, PERIODIC, geom)
  mesh => new_mesh_descriptor_2D(eta1_min, eta1_max, nc_eta1, &
       PERIODIC, eta2_min, eta2_max, nc_eta2, PERIODIC, geom)
 ! mesh => new_mesh_descriptor_2D(r_min, r_max, nc_eta1, &
 !      COMPACT, 0.0_f64, 2*sll_pi, nc_eta2, PERIODIC, geom)
  dist_func => sll_new_distribution_function_2D(mesh,CELL_CENTERED_DF, name)


  call sll_init_distribution_function_2D( dist_func, GAUSSIAN )
  call write_mesh_2D(mesh)

  call write_distribution_function ( dist_func )

  Print*, 'checking advection of a Gaussian in a uniform field'
  !    no splitting error. First and second order splitting should be same.
  !    only interpolation error

  ! define uniform field on coarse_mesh (using stream function)
  uniform_field => new_field_2D_vec1(mesh)
  ! components of field
  alpha1 = 0.0_f64
  alpha2 = 10.0_f64
  do i1 = 1, nc_eta1+1 
     do i2 = 1, nc_eta2+1
        FIELD_2D_AT_I( uniform_field, i1, i2 ) = alpha1 * x2n_array(i1,i2) &
             - alpha2 * x1n_array(i1,i2)
     end do
  end do
  ! initialize CSL  
  csl_work => new_csl_workspace( dist_func )
  ! run CSL method for 10 time steps
  n_steps = 100
  visu_step = 10
  deltat = 0.1_f64/real(n_steps,f64)!0._f64!10.0_f64/n_steps
  !deltat = 0.4_f64
  do it = 1, n_steps
     !print*, 'iteration=',it
     print*, 't=',real(it,f64)*deltat
     !call csl_first_order(csl_work, dist_func, uniform_field, deltat)
     call csl_second_order(csl_work, dist_func, uniform_field, uniform_field, deltat)
     if(mod(it,visu_step)==0)then
       call write_distribution_function ( dist_func )
     endif  
  end do
  ! compute error when Gaussian arrives at center (t=1)
  error = 0.0_f64
  do i1 = 1, nc_eta1
     do i2 = 1, nc_eta2
        val = sll_get_df_val(dist_func, i1, i2) / jac_array(i1,i2)
        !if (val > 1) then
        !   print*, 'val ',val, i1,i2
        !end if
        !error = max (error, abs(val - exp(-0.5_f64*40*((x1c_array(i1,i2)+.5)**2+(x2c_array(i1,i2)-.5)**2))))
        error = max (error, abs(val - exp(-0.5_f64*40*((x1c_array(i1,i2)-.5)**2+(x2c_array(i1,i2)-.5)**2))))
     end do
  end do
  print*, ' 1st order splitting, 100 cells, 10 time steps. Error= ', error
  ! reinitialize distribution function

  print*, 'checking advection in rotating field' 
  ! define rotating field
  rotating_field => new_field_2D_vec1(mesh)
  eta1 = eta1_min 
  do i1 = 1, nc_eta1+1
     eta2 = eta2_min 
     do i2 = 1, nc_eta2+1
        FIELD_2D_AT_I( rotating_field, i1, i2 ) = 0.5_f64*(x1n_array(i1,i2)**2 &
             + x2n_array(i1,i2)**2)
     end do
  end do

  ! reinitialize distribution function
  call sll_init_distribution_function_2D( dist_func, GAUSSIAN )
  ! run CSL method
  n_steps = 120
  deltat = 12*0.5_f64*sll_pi/n_steps  ! do one quarter turn
  do it = 1, n_steps
     call csl_first_order(csl_work, dist_func, rotating_field, deltat)
     !call write_distribution_function ( dist_func )
  end do
  ! compute error after one quarter turn
  error = 0.0_f64
  do i1 = 1, nc_eta1
     do i2 = 1, nc_eta2
        val = sll_get_df_val(dist_func, i1, i2) / jac_array(i1,i2)
        error = max (error, abs(val - exp(-0.5_f64*((x1c_array(i1,i2)-0.0_f64)**2 &
             + (x2c_array(i1,i2)+0.0_f64)**2))))
     end do
  end do
  print*, '    fine mesh, 1st order splitting, 200 cells, 10 time steps,  Error=', error
  error1 = error

  print *, 'Successful, exiting program.'

end program unit_test_from_array
