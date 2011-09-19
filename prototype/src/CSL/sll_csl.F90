module sll_csl
#include "sll_working_precision.h"
#include "sll_memory.h"
#include "sll_assert.h"
#include "sll_mesh_types.h"

  use numeric_constants
  use sll_splines
  use ode_solvers
  use distribution_function
  implicit none

  type csl_workspace
     type (sll_spline_1D), pointer :: spl_eta1
     type (sll_spline_1D), pointer :: spl_eta2
  end type csl_workspace

  sll_real64 :: global_eta1, global_eta2
  procedure(scalar_function_2D), pointer        :: x1
  procedure(scalar_function_2D), pointer        :: jac12
  procedure(scalar_function_2D), pointer        :: jac21
  procedure(scalar_function_2D), pointer        :: jac22
  procedure(scalar_function_2D), pointer        :: jac
contains
  
  function jac_first_direction( eta1 )
    use sll_working_precision
    sll_real64 :: jac_first_direction
    sll_real64, intent(in)  :: eta1

    jac_first_direction = jac (eta1, global_eta2)
  end function jac_first_direction

  function jac_second_direction( eta2 )
    use sll_working_precision
    sll_real64 :: jac_second_direction
    sll_real64, intent(in)  :: eta2

    jac_second_direction = jac (global_eta1, eta2)
  end function jac_second_direction
  
  function new_csl_workspace(dist_func_2D)
    type (csl_workspace), pointer :: new_csl_workspace
    type (sll_distribution_function_2D_t), pointer  :: dist_func_2D 
    sll_int32  :: ierr
    sll_int32  :: nc_eta1
    sll_int32  :: nc_eta2
    sll_real64 :: eta1_min
    sll_real64 :: eta1_max
    sll_real64 :: eta2_min
    sll_real64 :: eta2_max
    sll_int32  :: boundary1_type
    sll_int32  :: boundary2_type
    sll_real64 :: delta_eta2

    ! allocate pointer
    SLL_ALLOCATE(new_csl_workspace,ierr)

    ! get dimensions
    nc_eta1    = get_df_nc_eta1( dist_func_2D ) 
    eta1_min   = get_df_eta1_min( dist_func_2D )
    eta1_max   = get_df_eta1_max( dist_func_2D )
    nc_eta2    = get_df_nc_eta2( dist_func_2D ) 
    eta2_min   = get_df_eta2_min( dist_func_2D )
    eta2_max   = get_df_eta2_max( dist_func_2D )
    boundary1_type = get_df_boundary1_type( dist_func_2D )
    boundary2_type = get_df_boundary2_type( dist_func_2D )

    ! initialize splines
    if (boundary1_type == PERIODIC) then
       new_csl_workspace%spl_eta1 => new_spline_1D( nc_eta1+1,          &
                                                 eta1_min,         &
                                                 eta1_max,         &
                                                 PERIODIC_SPLINE )
    else if (boundary1_type == COMPACT) then
       new_csl_workspace%spl_eta1 => new_spline_1D( nc_eta1+1,          &
                                                 eta1_min,         &
                                                 eta1_max,         &
                                                 HERMITE_SPLINE )
    else
       print*, 'sll_csl.F90: new_csl_workspace. boundary1_type ', boundary1_type, ' not implemented'
       stop
    end if
    if (boundary2_type == PERIODIC) then
       new_csl_workspace%spl_eta2 => new_spline_1D( nc_eta2+1,        &
                                                 eta2_min,         &
                                                 eta2_max,         &
                                                 PERIODIC )  
    else if (boundary2_type == COMPACT) then
       new_csl_workspace%spl_eta2 => new_spline_1D( nc_eta2+1,        &
                                                 eta2_min,         &
                                                 eta2_max,         &
                                                 HERMITE_SPLINE )   
    else
       print*, 'sll_csl.F90: new_csl_workspace. boundary2_type ', boundary2_type, ' not implemented'
       stop
    end if

  end function new_csl_workspace

  subroutine delete_csl_workspace(csl_worksp)
    type (csl_workspace), pointer :: csl_worksp
    sll_int32   :: ierr

    if( .not. (associated(csl_worksp))) then
       write (*,'(a)') 'ERROR: delete_csl_workspace(), not associated argument.'
       STOP
    end if
    nullify(csl_worksp%spl_eta1)
    nullify(csl_worksp%spl_eta2)
    SLL_DEALLOCATE(csl_worksp, ierr)
  end subroutine delete_csl_workspace

  ! the code between first and second order is very close. 
  ! direction can be decoupled and should maybe be put in different 
  ! subroutines.

  ! subroutine csl_first_order
  ! Advances the distribution function on a time step deltat using a first 
  ! order time split conservative semi-Lagrangian scheme
  subroutine csl_first_order(csl_work, dist_func_2D, advfield, deltat)
    type (csl_workspace), pointer                   :: csl_work
    type (sll_distribution_function_2D_t), pointer  :: dist_func_2D  
    type (field_2D_vec1), pointer                   :: advfield ! advection field defined by its stream function
    sll_real64  ::  deltat  ! time step

    sll_int32  :: order 
    sll_real64, dimension(:), pointer  ::  advfield_1D_1
    sll_real64, dimension(:), pointer  ::  advfield_1D_2
    sll_real64, dimension(:), pointer  ::  primitive1
    sll_real64, dimension(:), pointer  ::  primitive2
    sll_real64, dimension(:), pointer  ::  eta1_out 
    sll_real64, dimension(:), pointer  ::  eta2_out
    sll_real64, dimension(:), pointer  ::  zeros1
    sll_real64, dimension(:), pointer  ::  zeros2
    procedure(scalar_function_1D), pointer  :: jacobian
    sll_int32  :: i1
    sll_int32  :: i2
    sll_int32  :: ierr
    sll_int32  :: nc_eta1
    sll_int32  :: nc_eta2
    sll_real64 :: delta_eta1
    sll_real64 :: delta_eta2
    sll_real64 :: eta1_min
    sll_real64 :: eta2_min
    sll_real64 :: eta1_max
    sll_real64 :: eta2_max
    sll_int32  :: boundary1_type
    sll_int32  :: boundary2_type

    sll_real64 :: val
    sll_real64 :: avg
    sll_real64 :: eta1
    sll_real64 :: eta2

    ! order of scheme
    order = 1

    ! parameter checking. We need to consider putting these checks on a more
    ! permanent basis, with their corresponding if-statements, instead of
    ! the assertions, which can be turned off.
    SLL_ASSERT(associated(csl_work))
    SLL_ASSERT(associated(dist_func_2D))
    SLL_ASSERT(associated(advfield))

    ! get dimensions
    nc_eta1    = get_df_nc_eta1( dist_func_2D ) 
    delta_eta1 = get_df_delta_eta1( dist_func_2D )
    eta1_min   = get_df_eta1_min( dist_func_2D )
    eta1_max   = get_df_eta1_max( dist_func_2D )
    nc_eta2    = get_df_nc_eta2( dist_func_2D ) 
    delta_eta2 = get_df_delta_eta2( dist_func_2D )
    eta2_min   = get_df_eta2_min( dist_func_2D )
    eta2_max   = get_df_eta2_max( dist_func_2D )
    boundary1_type = get_df_boundary1_type( dist_func_2D )
    boundary2_type = get_df_boundary2_type( dist_func_2D )
    x1 = get_df_x1( dist_func_2D )
    jac12 = get_df_jac12( dist_func_2D )
    jac21 = get_df_jac21( dist_func_2D )
    jac22 = get_df_jac22( dist_func_2D )
    jac = get_df_jac( dist_func_2D )
    
    ! allocation
    SLL_ALLOCATE(zeros1(nc_eta1+1),ierr)
    SLL_ALLOCATE(zeros2(nc_eta2+1),ierr)
    SLL_ALLOCATE(advfield_1D_1(nc_eta1+1),ierr)
    SLL_ALLOCATE(advfield_1D_2(nc_eta2+1),ierr)
    SLL_ALLOCATE(primitive1(nc_eta1+1),ierr)
    SLL_ALLOCATE(primitive2(nc_eta2+1),ierr)
    SLL_ALLOCATE(eta1_out(nc_eta1+1),ierr)
    SLL_ALLOCATE(eta2_out(nc_eta2+1),ierr)

   
    zeros1(:) = 0.0_f64
    zeros2(:) = 0.0_f64
    
    ! advection along the first direction 
    eta2 = eta2_min + 0.5_f64*delta_eta2  ! at cell center in this direction
    do i2=1, nc_eta2
       eta1 = eta1_min  ! at nodes
       primitive1 (1) = 0.0_f64  ! set primitive to 0 on left boundary 
       !advfield_1D_1 ( 1 ) = jac11(eta1,eta2)*FIELD_2D_AT_I_V1( advfield, 1, i2 ) &
       !                    + jac12(eta1,eta2)*FIELD_2D_AT_I_V2( advfield, 1, i2 )
       advfield_1D_1 ( 1 ) = (FIELD_2D_AT_I( advfield, 1, i2+1 ) - FIELD_2D_AT_I( advfield, 1, i2 )) / &
            ( delta_eta2 )!* jac(eta1,eta2) )
       do i1 = 2, nc_eta1+1
          ! extract subarray from advection field
          !advfield_1D_1 ( i1 ) = jac11(eta1,eta2)*FIELD_2D_AT_I_V1( advfield, i1, i2 ) &
          !                     + jac12(eta1,eta2)*FIELD_2D_AT_I_V2( advfield, i1, i2 )
           advfield_1D_1 ( i1 ) = (FIELD_2D_AT_I( advfield, i1, i2+1 ) - FIELD_2D_AT_I( advfield, i1, i2 )) / &
            ( delta_eta2)! * jac(eta1,eta2) )
          ! compute primitive of distribution function along this line
          primitive1 ( i1 ) = primitive1 ( i1-1 ) &
               + delta_eta1 * sll_get_df_val( dist_func_2D, i1-1, i2 ) 
          eta1 = eta1 + delta_eta1
       end do
       global_eta2 = eta2
       jacobian => jac_first_direction
       call advance_1D(primitive1, advfield_1D_1, zeros1, order, deltat, &
            eta1_min, nc_eta1, delta_eta1, boundary1_type, csl_work%spl_eta1, &
            eta1_out, jacobian) 

       ! update average value of distribution function in cell using 
       ! difference of primitive
       do i1 = 1, nc_eta1 
          val = (primitive1 ( i1+1 ) - primitive1 ( i1 )) / delta_eta1  
          call sll_set_df_val( dist_func_2D, i1, i2, val )
       end do
       eta2 = eta2 + delta_eta2
    end do
    ! advection along the second direction
    eta1 = eta1_min + 0.5_f64*delta_eta1
    do i1=1, nc_eta1
       eta2 = eta2_min 
       primitive2 (1) = 0.0_f64  ! set primitive to 0 on left boundary 
       !advfield_1D_2(1) =  jac21(eta1,eta2)*FIELD_2D_AT_I_V1( advfield, i1, 1 ) &
       !                        + jac22(eta1,eta2)*FIELD_2D_AT_I_V2( advfield, i1, 1 )
       advfield_1D_2 ( 1 ) = (FIELD_2D_AT_I( advfield, i1, 1 ) - FIELD_2D_AT_I( advfield, i1+1, 1 )) / &
            ( delta_eta1)! * jac(eta1,eta2) )
       do i2 = 2, nc_eta2+1
          eta2 = eta2 + delta_eta2
          ! extract subarray from advection field
          !advfield_1D_2(i2) =  jac21(eta1,eta2)*FIELD_2D_AT_I_V1( advfield, i1, i2 ) &
          !                     + jac22(eta1,eta2)*FIELD_2D_AT_I_V2( advfield, i1, i2 )
          advfield_1D_2 ( i2 ) = (FIELD_2D_AT_I( advfield, i1, i2 ) - FIELD_2D_AT_I( advfield, i1+1, i2 )) / &
            ( delta_eta1)! * jac(eta1,eta2) )
          ! compute primitive of distribution function along this line
          primitive2 (i2) = primitive2 (i2-1) &
               + delta_eta2 * sll_get_df_val( dist_func_2D, i1, i2-1 ) 
       end do
       global_eta1 = eta1
       jacobian => jac_second_direction
       call advance_1D( primitive2, advfield_1D_2, zeros2, order, deltat, &
                        eta2_min, nc_eta2, delta_eta2,          &
                        boundary2_type, csl_work%spl_eta2, eta2_out, jacobian ) 
       ! update average value of distribution function in cell using 
       ! difference of primitives
       do i2 = 1, nc_eta2 
          val = ( primitive2 ( i2+1 ) - primitive2 ( i2 )  ) / delta_eta2
          call sll_set_df_val( dist_func_2D, i1, i2, val )
       end do
       eta1 = eta1 + delta_eta1
    end do

  end subroutine csl_first_order

  ! subroutine csl_second_order
  ! Advances the distribution function on a time step deltat using a second 
  ! order time split (Strang splitting)
  ! conservative semi-Lagrangian scheme
  subroutine csl_second_order( csl_work,       &
                               dist_func_2D,   &
                               advfield_old,   &
                               advfield_new,   &
                               deltat )
    type (csl_workspace), pointer :: csl_work
    type (sll_distribution_function_2D_t), pointer  :: dist_func_2D
    type (field_2D_vec1), pointer  :: advfield_old   ! adv. field at (t)
    type (field_2D_vec1), pointer  :: advfield_new   ! adv. field at (t+dt)
    sll_real64  ::  deltat                           ! dt

    sll_int32  :: order 
    procedure(scalar_function_1D), pointer  :: jacobian
    sll_real64, dimension(:), pointer  ::  advfield_1D_1_old
    sll_real64, dimension(:), pointer  ::  advfield_1D_2_old
    sll_real64, dimension(:), pointer  ::  advfield_1D_1_new
    sll_real64, dimension(:), pointer  ::  advfield_1D_2_new
    sll_real64, dimension(:), pointer  ::  primitive1
    sll_real64, dimension(:), pointer  ::  primitive2
    sll_real64, dimension(:), pointer  ::  eta1_out 
    sll_real64, dimension(:), pointer  ::  eta2_out
    sll_int32  :: i1
    sll_int32  :: i2
    sll_int32  :: ierr
    sll_int32  :: nc_eta1
    sll_int32  :: nc_eta2
    sll_real64 :: delta_eta1
    sll_real64 :: delta_eta2
    sll_real64 :: eta1_min
    sll_real64 :: eta2_min
    sll_real64 :: eta1_max
    sll_real64 :: eta2_max
    sll_int32  :: boundary1_type
    sll_int32  :: boundary2_type
    sll_real64 :: val
    sll_real64 :: avg
    sll_real64 :: eta1
    sll_real64 :: eta2

    ! order of scheme
    order = 2
    ! parameter checking
    SLL_ASSERT(associated(csl_work))
    SLL_ASSERT(associated(dist_func_2D))
    SLL_ASSERT(associated(advfield_old))
    SLL_ASSERT(associated(advfield_new))

    ! get dimensions
    nc_eta1    = get_df_nc_eta1( dist_func_2D ) 
    delta_eta1 = get_df_delta_eta1( dist_func_2D )
    eta1_min   = get_df_eta1_min( dist_func_2D )
    eta1_max   = get_df_eta1_max( dist_func_2D )
    nc_eta2    = get_df_nc_eta2( dist_func_2D ) 
    delta_eta2 = get_df_delta_eta2( dist_func_2D )
    eta2_min   = get_df_eta2_min( dist_func_2D )
    eta2_max   = get_df_eta2_max( dist_func_2D )
    boundary1_type = get_df_boundary1_type( dist_func_2D )
    boundary2_type = get_df_boundary2_type( dist_func_2D )
    
    ! allocation
    SLL_ALLOCATE(advfield_1D_1_old(nc_eta1+1),ierr)
    SLL_ALLOCATE(advfield_1D_2_old(nc_eta2+1),ierr)
    SLL_ALLOCATE(advfield_1D_1_new(nc_eta1+1),ierr)
    SLL_ALLOCATE(advfield_1D_2_new(nc_eta2+1),ierr)
    SLL_ALLOCATE(primitive1(nc_eta1+1),ierr)
    SLL_ALLOCATE(primitive2(nc_eta2+1),ierr)
    SLL_ALLOCATE(eta1_out(nc_eta1+1),ierr)
    SLL_ALLOCATE(eta2_out(nc_eta2+1),ierr)
    
    ! half advection along the first direction 
    eta2 = eta2_min + 0.5_f64*delta_eta2  ! at cell center in this direction
    do i2 = 1, nc_eta2
       eta1 = eta1_min  ! at nodes
       primitive1 (1) = 0.0_f64  ! set primitive to 0 on left boundary 
       !advfield_1D_1_old ( 1 ) = FIELD_2D_AT_I_V1( advfield_old, 1, i2 )
       !advfield_1D_1_new ( 1 ) = FIELD_2D_AT_I_V1( advfield_new, 1, i2 )
       advfield_1D_1_old ( 1 ) = (FIELD_2D_AT_I( advfield_old, 1, i2+1 ) - FIELD_2D_AT_I( advfield_old, 1, i2 )) / &
            ( delta_eta2)! * jac(eta1,eta2) )
       advfield_1D_1_new ( 1 ) = (FIELD_2D_AT_I( advfield_new, 1, i2+1 ) - FIELD_2D_AT_I( advfield_new, 1, i2 )) / &
            ( delta_eta2)! * jac(eta1,eta2) )
       do i1 = 2, nc_eta1+1
          ! extract subarray from advection field
          !advfield_1D_1_old ( i1 ) = FIELD_2D_AT_I_V1( advfield_old, i1, i2 )
          !advfield_1D_1_new ( i1 ) = FIELD_2D_AT_I_V1( advfield_new, i1, i2 )
          advfield_1D_1_old ( i1 ) = (FIELD_2D_AT_I( advfield_old, i1, i2+1 ) &
               - FIELD_2D_AT_I( advfield_old, i1, i2 )) / ( delta_eta2)! * jac(eta1,eta2) )
          advfield_1D_1_new ( i1 ) = (FIELD_2D_AT_I( advfield_new, i1, i2+1 ) &
               - FIELD_2D_AT_I( advfield_new, i1, i2 )) / ( delta_eta2)! * jac(eta1,eta2) )
          ! compute primitive of distribution function along this line
          primitive1 ( i1 ) = primitive1 ( i1-1 ) &
               + delta_eta1 * sll_get_df_val( dist_func_2D, i1-1, i2 )
          eta1 = eta1 + delta_eta1
       end do
       global_eta2 = eta2
       jacobian => jac_first_direction
       call advance_1D( primitive1,        &
                        advfield_1D_1_old, &
                        advfield_1D_1_new, &
                        order,             &
                        0.5_f64*deltat,    &
                        eta1_min,          &
                        nc_eta1,           &
                        delta_eta1,        &
                        boundary1_type,    &
                        csl_work%spl_eta1, &
                        eta1_out,          & 
                        jacobian) 
       ! update average value of distribution function in cell using 
       ! difference of primitives
       do i1 = 1, nc_eta1 
          val = (primitive1 ( i1+1 ) - primitive1 ( i1 )) / delta_eta1
          call sll_set_df_val( dist_func_2D, i1, i2, val )
       end do
    eta2 = eta2 + delta_eta2
    end do
    ! advection along the second direction
    eta1 = eta1_min + 0.5_f64*delta_eta1 ! cell centered
    do i1=1, nc_eta1
       eta2 = eta2_min  ! node centered
       primitive2 (1) = 0.0_f64  ! set primitive to 0 on left boundary 
       !advfield_1D_2_old(1) = FIELD_2D_AT_I_V2( advfield_old, i1, 1 )
       !advfield_1D_2_new(1) = FIELD_2D_AT_I_V2( advfield_new, i1, 1 )
       advfield_1D_2_old ( 1 ) = (FIELD_2D_AT_I( advfield_old, i1, 1 ) &
               - FIELD_2D_AT_I( advfield_old, i1+1, 1 )) / ( delta_eta1)! * jac(eta1,eta2) )
       advfield_1D_2_new ( 1 ) = (FIELD_2D_AT_I( advfield_new, i1, 1 ) &
            - FIELD_2D_AT_I( advfield_new, i1+1, 1 )) / ( delta_eta1)! * jac(eta1,eta2) )
       do i2 = 2, nc_eta2+1
          eta2 = eta2 + delta_eta2
          ! extract subarray from advection field
          !advfield_1D_2_old(i2) = FIELD_2D_AT_I_V2( advfield_old, i1, i2 )
          !advfield_1D_2_new(i2) = FIELD_2D_AT_I_V2( advfield_new, i1, i2 )
          advfield_1D_2_old ( i2 ) = (FIELD_2D_AT_I( advfield_old, i1, i2 ) &
               - FIELD_2D_AT_I( advfield_old, i1+1, i2 )) / ( delta_eta1)! * jac(eta1,eta2) )
          advfield_1D_2_new ( i2 ) = (FIELD_2D_AT_I( advfield_new, i1, i2 ) &
               - FIELD_2D_AT_I( advfield_new, i1+1, i2 )) / ( delta_eta1)! * jac(eta1,eta2) )
          ! compute primitive of distribution function along this line
          primitive2 (i2) = primitive2 (i2-1) &
               + delta_eta2 * sll_get_df_val( dist_func_2D, i1, i2-1 )
       end do
       global_eta1 = eta1
       jacobian => jac_second_direction
       call advance_1D( primitive2,        &
                        advfield_1D_2_old, &
                        advfield_1D_2_new, &
                        order,             &
                        deltat,            &
                        eta2_min,          &
                        nc_eta2,           &
                        delta_eta2,        &
                        boundary2_type,    &
                        csl_work%spl_eta2, &
                        eta2_out,          &
                        jacobian) 
       ! update average value of distribution function in cell using 
       ! difference of primitives
       do i2 = 1, nc_eta2 
          val = (primitive2(i2+1) - primitive2(i2))/delta_eta2 
          call sll_set_df_val( dist_func_2D, i1, i2, val )
       end do
       eta1 = eta1 + delta_eta1
    end do
    ! half advection along the first direction 
    ! this is exactly the same code as the first half advection.
    ! might be coded better
    eta2 = eta2_min + 0.5_f64*delta_eta2  ! at cell center in this direction
    do i2 = 1, nc_eta2
       eta1 = eta1_min  ! at nodes
       primitive1 (1) = 0.0_f64  ! set primitive to 0 on left boundary 
       !advfield_1D_1_old ( 1 ) = FIELD_2D_AT_I_V1( advfield_old, 1, i2 )
       !advfield_1D_1_new ( 1 ) = FIELD_2D_AT_I_V1( advfield_new, 1, i2 )
       advfield_1D_1_old ( 1 ) = (FIELD_2D_AT_I( advfield_old, 1, i2+1 ) - FIELD_2D_AT_I( advfield_old, 1, i2 )) / &
            ( delta_eta2)! * jac(eta1,eta2) )
       advfield_1D_1_new ( 1 ) = (FIELD_2D_AT_I( advfield_new, 1, i2+1 ) - FIELD_2D_AT_I( advfield_new, 1, i2 )) / &
            ( delta_eta2)! * jac(eta1,eta2) )
       do i1 = 2, nc_eta1+1
          ! extract subarray from advection field
          !advfield_1D_1_old ( i1 ) = FIELD_2D_AT_I_V1( advfield_old, i1, i2 )
          !advfield_1D_1_new ( i1 ) = FIELD_2D_AT_I_V1( advfield_new, i1, i2 )
          advfield_1D_1_old ( i1 ) = (FIELD_2D_AT_I( advfield_old, i1, i2+1 ) &
               - FIELD_2D_AT_I( advfield_old, i1, i2 )) / ( delta_eta2)! * jac(eta1,eta2) )
          advfield_1D_1_new ( i1 ) = (FIELD_2D_AT_I( advfield_new, i1, i2+1 ) &
               - FIELD_2D_AT_I( advfield_new, i1, i2 )) / ( delta_eta2)! * jac(eta1,eta2) )
          ! compute primitive of distribution function along this line
          primitive1 ( i1 ) = primitive1 ( i1-1 ) &
               + delta_eta1 * sll_get_df_val( dist_func_2D, i1-1, i2 )
          eta1 = eta1 + delta_eta1
       end do
       global_eta2 = eta2
       jacobian => jac_first_direction
       call advance_1D( primitive1,        &
                        advfield_1D_1_old, &
                        advfield_1D_1_new, &
                        order,             &
                        0.5_f64*deltat,    &
                        eta1_min,          &
                        nc_eta1,           &
                        delta_eta1,        &
                        boundary1_type,    &
                        csl_work%spl_eta1, &
                        eta1_out,          &
                        jacobian) 
       ! update average value of distribution function in cell using 
       ! difference of primitives
       do i1 = 1, nc_eta1 
          val = (primitive1 ( i1+1 ) - primitive1 ( i1 )) / delta_eta1
          call sll_set_df_val( dist_func_2D, i1, i2, val )
       end do
       eta2 = eta2 + delta_eta2
    end do
  end subroutine csl_second_order


  subroutine advance_1D( primitive,     &
                         fieldn,        &
                         fieldnp1,      &
                         order,         &
                         deltat,        &
                         eta_min,       &
                         nc_eta,        & 
                         delta_eta,     &
                         boundary_type, &
                         spline,        &
                         eta_out,   &
                         jacobian)  
    sll_real64, dimension(:), pointer, intent(inout) :: primitive
    sll_real64, dimension(:), pointer, intent(in)    :: fieldn
    sll_real64, dimension(:), pointer, intent(in)    :: fieldnp1
    sll_int32                               :: order
    sll_real64                              :: deltat
    sll_real64                              :: eta_min
    sll_int32                               :: nc_eta
    sll_real64                              :: delta_eta
    sll_int32                               :: boundary_type
    type (sll_spline_1D), pointer           :: spline
    sll_real64, dimension(:), intent(out)   :: eta_out
    procedure(scalar_function_1D), pointer  :: jacobian
    
    ! local variables
    sll_real64  :: avg
    sll_real64  :: eta_max
    sll_real64  :: eta
    sll_real64  :: eta_new
    sll_real64  :: lperiod
    sll_int32   :: i

    ! check array dimensions
    SLL_ASSERT(size(primitive) >= nc_eta+1)
    SLL_ASSERT(size(fieldn)    >= nc_eta+1)
    SLL_ASSERT(size(fieldnp1)  >= nc_eta+1)
    SLL_ASSERT(size(eta_out)   >= nc_eta+1)

    select case (boundary_type)
    case (PERIODIC)
       eta_max = eta_min + nc_eta * delta_eta 
       ! average of dist func along the line
       avg = primitive ( nc_eta+1 ) / (eta_max - eta_min) 
       ! modify primitive so that it becomes periodic
       eta = eta_min
       do i = 2, nc_eta+1
          eta = eta + delta_eta
          primitive ( i ) = primitive ( i ) - avg * (eta-eta_min)
       end do
!!$       call implicit_ode( order,        &
!!$                          deltat,       &
!!$                          eta_min,      &
!!$                          nc_eta,       &
!!$                          delta_eta,    &
!!$                          PERIODIC_ODE, &
!!$                          eta_out,      &
!!$                          fieldn,       &
!!$                          fieldnp1 )
       call rk4( 5,        &
                -deltat,       &
                 eta_min,      &
                 nc_eta,       &
                 delta_eta,    &
                 PERIODIC_ODE, &
                 eta_out,      &
                 fieldn,       &
                 jacobian)
       call compute_spline_1D_periodic( primitive, spline )
       ! interpolate primitive at origin of characteritics
       call interpolate_array_values( eta_out, primitive, nc_eta+1, spline )
       ! come back to real primitive by adding average
       if (eta_out(1) > eta_min + 0.5_f64*(eta_max-eta_min)) then
          lperiod = -1.0_f64
       else
          lperiod = 0.0_f64
       end if
       eta_new = eta_out(1) +  lperiod*(eta_max-eta_min)
       primitive ( 1 ) = primitive ( 1 ) + avg * (eta_new-eta_min)
       !if ((eta_new > eta_max) .or. (eta_new <eta_min)) then
       !   print*, 1, eta_new, eta_out(1), primitive(1)
       !end if
       eta = eta_min
       do i = 2, nc_eta+1
          eta = eta_min + delta_eta
          ! eta_out(i) is an increasing sequence because grid points cannot catch up with each other
          ! We need here to find the points where it has been modified by periodicity
          if (eta_out(i) < eta_out(i-1)) then
             lperiod = lperiod + 1.0_f64
          end if
          eta_new = eta_out(i) +  lperiod*(eta_max-eta_min)
          primitive ( i ) = primitive ( i ) + avg * (eta_new-eta_min)
          !if ((eta_new > eta_max) .or. (eta_new <eta_min)) then
          !   print*, i, eta_new, eta_out(i), primitive(i)
          !end if
       end do

    case (COMPACT)
!!$       call implicit_ode( order,       &
!!$                          deltat,      &
!!$                          eta_min,     &
!!$                          nc_eta,      &
!!$                          delta_eta,   &
!!$                          COMPACT_ODE, &
!!$                          eta_out,     &
!!$                          fieldn,      &
!!$                          fieldnp1 ) 
       call rk4( 5,        &
                -deltat,       &
                 eta_min,      &
                 nc_eta,       &
                 delta_eta,    &
                 COMPACT_ODE, &
                 eta_out,      &
                 fieldn,       &
                 jacobian)
       call compute_spline_1D_hermite( primitive, spline )
       ! interpolate primitive at origin of characteritics
       call interpolate_array_values( eta_out, primitive, nc_eta+1, spline )
    end select
  end subroutine advance_1D



end module sll_csl
