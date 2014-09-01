!> \brief Implements the conservative semi-Lagrangian algorithm
!>
!> The algorithm is described in Crouseilles, Mehrenberger, Sonnendrucker, J. Comput. Phys. 229, pp 1927-1953, 2010
!> http://dx.doi.org/10.1016/j.jcp.2009.11.007, preprint is available on HAL http://hal.inria.fr/hal-00363643/en

module sll_csl
#include "sll_working_precision.h"
#include "sll_memory.h"
#include "sll_assert.h"
!#include "sll_mesh_types.h"
#include "sll_field_2d.h"

  use sll_constants
  !use sll_splines
  use cubic_non_uniform_splines
  use ode_solvers
  use distribution_function
  implicit none
  integer, parameter :: PERIODIC = 0,COMPACT = 1

  type csl_workspace
     type (cubic_nonunif_spline_1D), pointer :: spl_eta1
     type (cubic_nonunif_spline_1D), pointer :: spl_eta2
  end type csl_workspace

contains
    
!> initialize opaque pointer that contains information used by the algorithm
!> \param[in] dist_func_2D 2D distribution function object
!> \return pointer to opaque data type 
  function new_csl_workspace(dist_func_2D)
    type (csl_workspace), pointer :: new_csl_workspace
    !type (sll_distribution_function_2d), pointer  :: dist_func_2D 
    type (sll_distribution_function_2d) :: dist_func_2D 
    sll_int32  :: ierr
    sll_int32  :: nc_eta1
    sll_int32  :: nc_eta2
    sll_real64 :: eta1_min
    sll_real64 :: eta1_max
    sll_real64 :: eta2_min
    sll_real64 :: eta2_max
    sll_int32  :: boundary1_type
    sll_int32  :: boundary2_type

    ! allocate pointer
    SLL_ALLOCATE(new_csl_workspace,ierr)

    ! get dimensions
    nc_eta1    = GET_FIELD_NC_ETA1( dist_func_2D ) 
    eta1_min   = 0._f64!get_df_eta1_min( dist_func_2D )
    eta1_max   = 1._f64!get_df_eta1_max( dist_func_2D )
    nc_eta2    = GET_FIELD_NC_ETA2( dist_func_2D ) 
    eta2_min   = 0._f64!get_df_eta2_min( dist_func_2D )
    eta2_max   = 1._f64!get_df_eta2_max( dist_func_2D )
    boundary1_type = COMPACT!get_df_boundary1_type( dist_func_2D )
    boundary2_type = PERIODIC!get_df_boundary2_type( dist_func_2D )

    ! initialize splines
    if (boundary1_type == PERIODIC) then
       new_csl_workspace%spl_eta1 => new_cubic_nonunif_spline_1D( nc_eta1, SLL_PERIODIC)
    else if (boundary1_type == COMPACT) then
       new_csl_workspace%spl_eta1 => new_cubic_nonunif_spline_1D( nc_eta1, SLL_HERMITE)
    else
       print*, 'sll_csl.F90: new_csl_workspace. boundary1_type ', boundary1_type, ' not implemented'
       stop
    end if
    if (boundary2_type == PERIODIC) then
       new_csl_workspace%spl_eta2 => new_cubic_nonunif_spline_1D( nc_eta2, SLL_PERIODIC)
    else if (boundary2_type == COMPACT) then
       new_csl_workspace%spl_eta2 => new_cubic_nonunif_spline_1D( nc_eta2, SLL_HERMITE)  
    else
       print*, 'sll_csl.F90: new_csl_workspace. boundary2_type ', boundary2_type, ' not implemented'
       stop
    end if

  end function new_csl_workspace

!> delete pointer on opaque data type
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

  !> Advances the distribution function on a time step deltat using a first 
  !> order time split conservative semi-Lagrangian scheme
  !> \param[in] csl_work pointer on CSL opaque object
  !> \param[in,out] dist_func_2D distribution function which is advanced
  !> \param[in] advfield advection field used for advancing distribution function
  !> \param[in] deltat time step on which distribution function is advanced
  subroutine csl_first_order(csl_work, dist_func_2D, advfield, deltat)
    type (csl_workspace), pointer                   :: csl_work
    type (sll_distribution_function_2D), pointer  :: dist_func_2D  
    type (scalar_field_2d), pointer                   :: advfield ! advection field defined by its stream function
    sll_real64, intent(in)  ::  deltat  ! time step
    ! local variables
    sll_int32, parameter   :: order = 1 ,field_order=4   ! order of scheme

    call csl_advance_1(csl_work, dist_func_2D, advfield, advfield, deltat, order,field_order)
    call csl_advance_2(csl_work, dist_func_2D, advfield, advfield, deltat, order,field_order)
  end subroutine csl_first_order

  !> Advances the distribution function on a time step deltat using a second
  !> order time split conservative semi-Lagrangian scheme (Strang splitting)
  !> \param[in] csl_work pointer on CSL opaque object
  !> \param[in,out] dist_func_2D distribution function which is advanced
  !> \param[in] advfield_old advection field at t used for advancing distribution function
  !> \param[in] advfield_new advection field at t+dt used for advancing distribution function
  !> \param[in] deltat time step on which distribution function is advanced
  subroutine csl_second_order(csl_work, dist_func_2D, advfield_old, advfield_new, deltat)
    type (csl_workspace), pointer                   :: csl_work
    type (sll_distribution_function_2D), pointer  :: dist_func_2D  
    type (scalar_field_2d), pointer                   :: advfield_old ! advection field at t
    type (scalar_field_2d), pointer                   :: advfield_new ! advection field at t+dt
    sll_real64, intent(in)  ::  deltat  ! time step
    ! local variables
    sll_int32, parameter   :: order = 2 ,field_order=4   ! order of scheme
 
    call csl_advance_1(csl_work, dist_func_2D, advfield_old, advfield_new, 0.5_f64*deltat, order,field_order)
    call csl_advance_2(csl_work, dist_func_2D, advfield_old, advfield_new, deltat, order,field_order)
    call csl_advance_1(csl_work, dist_func_2D, advfield_old, advfield_new, 0.5_f64*deltat, order,field_order)
  end subroutine csl_second_order


  !> Part of a split semi-Lagrangian algorithm that advances the distribution function in the first direction
  !> on one time step
  !> \param[in] csl_work pointer on CSL opaque object
  !> \param[in,out] dist_func_2D distribution function which is advanced
  !> \param[in] advfield_old advection field at t used for advancing distribution function
  !> \param[in] advfield_new advection field at t+dt used for advancing distribution function
  !> \param[in] deltat time step on which distribution function is advanced
  !> \param[in] order order of time scheme 1 or 2 needed by advance_1D
  subroutine csl_advance_1( csl_work,       &
                            dist_func_2D,   &
                            advfield_old,   &
                            advfield_new,   &
                            deltat,         &
                            order,field_order)
    type (csl_workspace), pointer :: csl_work
    type (sll_distribution_function_2D), pointer  :: dist_func_2D
    type (scalar_field_2d), pointer  :: advfield_old   ! adv. field at (t)
    type (scalar_field_2d), pointer  :: advfield_new   ! adv. field at (t+dt)
    sll_real64, intent(in)  ::  deltat                           ! dt
    sll_int32, intent(in)  :: order ,field_order

    sll_real64, dimension(:), pointer  ::  advfield_1D_1_old
    sll_real64, dimension(:), pointer  ::  advfield_1D_1_new
    sll_real64, dimension(:), pointer  ::  primitive1
    sll_real64, dimension(:), pointer  ::  eta1_out 
    sll_real64, dimension(:), pointer  ::  jacobian
    sll_real64, dimension(:,:), pointer  ::  df_jac_at_i
    
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
    sll_real64 :: val
    sll_real64 :: eta1
    sll_real64 :: eta2

    ! parameter checking
    SLL_ASSERT(associated(csl_work))
    SLL_ASSERT(associated(dist_func_2D))
    SLL_ASSERT(associated(advfield_old))
    SLL_ASSERT(associated(advfield_new))

    ! get dimensions
    nc_eta1    = GET_FIELD_NC_ETA1( dist_func_2D ) 
    delta_eta1 = GET_FIELD_DELTA_ETA1( dist_func_2D )
    eta1_min   = 0._f64!get_df_eta1_min( dist_func_2D )
    eta1_max   = 1._f64!get_df_eta1_max( dist_func_2D )
    nc_eta2    = GET_FIELD_NC_ETA2( dist_func_2D ) 
    delta_eta2 = GET_FIELD_DELTA_ETA2( dist_func_2D )
    eta2_min   = 0._f64!get_df_eta2_min( dist_func_2D )
    eta2_max   = 1._f64!get_df_eta2_max( dist_func_2D )
    boundary1_type = COMPACT!get_df_boundary1_type( dist_func_2D )
    
    ! allocation
    SLL_ALLOCATE(advfield_1D_1_old(nc_eta1+1),ierr)
    SLL_ALLOCATE(advfield_1D_1_new(nc_eta1+1),ierr)
    SLL_ALLOCATE(primitive1(nc_eta1+1),ierr)
    SLL_ALLOCATE(eta1_out(nc_eta1+1),ierr)
    SLL_ALLOCATE(jacobian(nc_eta1+1),ierr)
    SLL_ALLOCATE(df_jac_at_i(nc_eta1, nc_eta2),ierr)

    !df_jac_at_i => FIELD_JACOBIAN_CELL_DATA( dist_func_2D )
    !df_jac_at_i => FIELD_JACOBIAN_CELL_DATA( dist_func_2D )
    do i2 = 1, nc_eta2
       do i1 = 1, nc_eta1
          df_jac_at_i(i2,i1) = dist_func_2D%mesh%jacobian_at_cell(i2,i1)
       end do
    end do
    !do i2=1,nc_eta2
    !  do i1=1,nc_eta1
    !    df_jac_at_i(i1,i2) = FIELD_2D_JACOBIAN_AT_I( dist_func_2D ,i1,i2)
    !  enddo
    !enddo

    ! advection along the first direction 
    eta2 = eta2_min + 0.5_f64*delta_eta2  ! at cell center in this direction
    do i2 = 1, nc_eta2
       eta1 = eta1_min  ! at nodes
       primitive1 (1) = 0.0_f64  ! set primitive to 0 on left boundary 
       !advfield_1D_1_old ( 1 ) = FIELD_2D_AT_I_V1( advfield_old, 1, i2 )
       !advfield_1D_1_new ( 1 ) = FIELD_2D_AT_I_V1( advfield_new, 1, i2 )
       advfield_1D_1_old ( 1 ) = (FIELD_2D_AT_I( advfield_old, 1, i2+1 ) - FIELD_2D_AT_I( advfield_old, 1, i2 )) / &
            ( delta_eta2)
       if(field_order==4)then 
       
              advfield_1D_1_old ( 1 ) = 9._f64/8._f64*advfield_1D_1_old ( 1 )&
-1._f64/24._f64*(FIELD_2D_AT_I( advfield_old, 1,modulo(i2+2-1,nc_eta2)+1) &
               - FIELD_2D_AT_I( advfield_old, 1,modulo(i2-1-1+nc_eta2,nc_eta2)+1 )) / ( delta_eta2 )
       endif
            
            
            
       if (order == 1) then
          advfield_1D_1_new ( 1 ) = 0.0_f64
       else
          advfield_1D_1_new ( 1 ) = (FIELD_2D_AT_I( advfield_new, 1, i2+1 )   &
               - FIELD_2D_AT_I( advfield_new, 1, i2 )) / ( delta_eta2)
       if(field_order==4)then 

              advfield_1D_1_new ( 1 ) = 9._f64/8._f64*advfield_1D_1_new ( 1 )&
-1._f64/24._f64*(FIELD_2D_AT_I( advfield_new, 1,modulo(i2+2-1,nc_eta2)+1) &
               - FIELD_2D_AT_I( advfield_new, 1,modulo(i2-1-1+nc_eta2,nc_eta2)+1 )) / ( delta_eta2 )

    endif
       end if
       do i1 = 2, nc_eta1+1
          ! extract subarray from advection field
          !advfield_1D_1_old ( i1 ) = FIELD_2D_AT_I_V1( advfield_old, i1, i2 )
          !advfield_1D_1_new ( i1 ) = FIELD_2D_AT_I_V1( advfield_new, i1, i2 )
          advfield_1D_1_old ( i1 ) = (FIELD_2D_AT_I( advfield_old, i1, i2+1 ) &
               - FIELD_2D_AT_I( advfield_old, i1, i2 )) / delta_eta2
        if(field_order==4)then 
              
                        advfield_1D_1_old ( i1 ) = 9._f64/8._f64*advfield_1D_1_old ( i1 )&
-1._f64/24._f64*(FIELD_2D_AT_I( advfield_old, i1,modulo(i2+2-1,nc_eta2)+1) &
               - FIELD_2D_AT_I( advfield_old, i1,modulo(i2-1-1+nc_eta2,nc_eta2)+1 )) / ( delta_eta2 )
     
          endif     
          if (order == 1) then
             advfield_1D_1_new ( i1 ) = 0.0_f64
          else
             advfield_1D_1_new ( i1 ) = (FIELD_2D_AT_I( advfield_new, i1, i2+1 ) &
                  - FIELD_2D_AT_I( advfield_new, i1, i2 )) / delta_eta2
       if(field_order==4)then 

                        advfield_1D_1_new ( i1 ) = 9._f64/8._f64*advfield_1D_1_new ( i1 )&
-1._f64/24._f64*(FIELD_2D_AT_I( advfield_new, i1,modulo(i2+2-1,nc_eta2)+1) &
               - FIELD_2D_AT_I( advfield_new, i1,modulo(i2-1-1+nc_eta2,nc_eta2)+1 )) / ( delta_eta2 )

endif
          end if
          ! compute primitive of distribution function along this line
          primitive1 ( i1 ) = primitive1 ( i1-1 ) &
               + delta_eta1 * FIELD_2D_AT_I( dist_func_2D, i1-1, i2 )
          eta1 = eta1 + delta_eta1
          jacobian(i1) = df_jac_at_i( i1-1, i2 )
       end do
       call advance_1D_nonuniform( primitive1,        &
                        advfield_1D_1_old, &
                        advfield_1D_1_new, &
                        jacobian,          &
                        order,             &
                        deltat,            &
                        nc_eta1,           &
                        delta_eta1,        &
                        boundary1_type,    &
                        csl_work%spl_eta1, &
                        eta1_out) 
       ! update average value of distribution function in cell using 
       ! difference of primitives
       do i1 = 1, nc_eta1 
          val = (primitive1 ( i1+1 ) - primitive1 ( i1 )) / delta_eta1
          call sll_set_df_val( dist_func_2D, i1, i2, val )
          !if (val/df_jac_at_i(i1,i2)>1.) then
          !   print*, 'val', i1,i2, val, primitive1(i1) , primitive1(i1+1), df_jac_at_i(i1,i2), delta_eta1
          !end if
       end do
    eta2 = eta2 + delta_eta2
    end do
  end subroutine csl_advance_1



  !> Part of a split semi-Lagrangian algorithm that advances the distribution function in the second direction
  !> on one time step
  !> \param[in] csl_work pointer on CSL opaque object
  !> \param[in,out] dist_func_2D distribution function which is advanced
  !> \param[in] advfield_old advection field at t used for advancing distribution function
  !> \param[in] advfield_new advection field at t+dt used for advancing distribution function
  !> \param[in] deltat time step on which distribution function is advanced
  !> \param[in] order order of time scheme 1 or 2 needed by advance_1D
  subroutine csl_advance_2( csl_work,       &
                            dist_func_2D,   &
                            advfield_old,   &
                            advfield_new,   &
                            deltat,         &
                            order ,field_order)
    type (csl_workspace), pointer :: csl_work
    type (sll_distribution_function_2D), pointer  :: dist_func_2D
    type (scalar_field_2d), pointer  :: advfield_old   ! adv. field at (t)
    type (scalar_field_2d), pointer  :: advfield_new   ! adv. field at (t+dt)
    sll_real64, intent(in)  :: deltat                            ! dt
    sll_int32, intent(in)   :: order ,field_order                            ! order

    ! local variables
    sll_real64, dimension(:), pointer  ::  advfield_1D_2_old
    sll_real64, dimension(:), pointer  ::  advfield_1D_2_new
    sll_real64, dimension(:), pointer  ::  primitive2
    sll_real64, dimension(:), pointer  ::  eta2_out
    sll_real64, dimension(:,:), pointer  :: df_jac_at_i, x1c_at_i, x2c_at_i
    sll_real64, dimension(:), pointer  ::  jacobian
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
    sll_int32  :: boundary2_type
    sll_real64 :: val
    sll_real64 :: eta1
    sll_real64 :: eta2

    ! parameter checking
    SLL_ASSERT(associated(csl_work))
    SLL_ASSERT(associated(dist_func_2D))
    SLL_ASSERT(associated(advfield_old))
    SLL_ASSERT(associated(advfield_new))

    ! get dimensions
    nc_eta1    = GET_FIELD_NC_ETA1( dist_func_2D ) 
    delta_eta1 = GET_FIELD_DELTA_ETA1( dist_func_2D )
    eta1_min   = 0._f64!get_df_eta1_min( dist_func_2D )
    eta1_max   = 1._f64!get_df_eta1_max( dist_func_2D )
    nc_eta2    = GET_FIELD_NC_ETA2( dist_func_2D ) 
    delta_eta2 = GET_FIELD_DELTA_ETA2( dist_func_2D )
    eta2_min   = 0._f64!get_df_eta2_min( dist_func_2D )
    eta2_max   = 1._f64!get_df_eta2_max( dist_func_2D )
    boundary2_type = PERIODIC!get_df_boundary2_type( dist_func_2D )
    
    ! allocation
    SLL_ALLOCATE(advfield_1D_2_old(nc_eta2+1),ierr)
    SLL_ALLOCATE(advfield_1D_2_new(nc_eta2+1),ierr)
    SLL_ALLOCATE(primitive2(nc_eta2+1),ierr)
    SLL_ALLOCATE(eta2_out(nc_eta2+1),ierr)
    SLL_ALLOCATE(jacobian(nc_eta2+1),ierr)
    SLL_ALLOCATE(df_jac_at_i(nc_eta1, nc_eta2),ierr)
    SLL_ALLOCATE(x1c_at_i(nc_eta1, nc_eta2),ierr)
    SLL_ALLOCATE(x2c_at_i(nc_eta1, nc_eta2),ierr)

    !df_jac_at_i => FIELD_JACOBIAN_CELL_DATA( dist_func_2D )
    do i2 = 1, nc_eta2
       do i1 = 1, nc_eta1
          df_jac_at_i(i2,i1) = dist_func_2D%mesh%jacobian_at_cell(i2,i1)
          x1c_at_i(i2,i1) = dist_func_2D%mesh%x1_at_cell(i2,i1)
          x2c_at_i(i2,i1) = dist_func_2D%mesh%x2_at_cell(i2,i1)
       end do
    end do
    !x1c_at_i => get_df_x1c_at_i( dist_func_2D )        
    !x2c_at_i => get_df_x2c_at_i( dist_func_2D )
    !x1c_at_i => FIELD_X1_CELL( dist_func_2D )
    !x2c_at_i => FIELD_X2_CELL( dist_func_2D )
    
    ! advection along the second direction
    eta1 = eta1_min + 0.5_f64*delta_eta1 ! cell centered
    do i1=1, nc_eta1
       eta2 = eta2_min  ! node centered
       primitive2 (1) = 0.0_f64  ! set primitive to 0 on left boundary 
       !advfield_1D_2_old(1) = FIELD_2D_AT_I_V2( advfield_old, i1, 1 )
       !advfield_1D_2_new(1) = FIELD_2D_AT_I_V2( advfield_new, i1, 1 )
       advfield_1D_2_old ( 1 ) = (FIELD_2D_AT_I( advfield_old, i1, 1 ) &
               - FIELD_2D_AT_I( advfield_old, i1+1, 1 )) / ( delta_eta1 )
       if(field_order==4)then 
       advfield_1D_2_old ( 1 ) = 9._f64/8._f64*advfield_1D_2_old ( 1 )&
+1._f64/24._f64*(FIELD_2D_AT_I( advfield_old, modulo(i1+2-1,nc_eta1)+1, 1 ) &
               - FIELD_2D_AT_I( advfield_old, modulo(i1-1-1+nc_eta1,nc_eta1)+1, 1 )) / ( delta_eta1 )
       endif
       if (order == 1) then
          advfield_1D_2_new ( 1 ) = 0.0_f64
       else
          advfield_1D_2_new ( 1 ) = (FIELD_2D_AT_I( advfield_new, i1, 1 ) &
               - FIELD_2D_AT_I( advfield_new, i1+1, 1 )) / delta_eta1
       if(field_order==4)then 
       advfield_1D_2_new ( 1 ) = 9._f64/8._f64*advfield_1D_2_new ( 1 )&
+1._f64/24._f64*(FIELD_2D_AT_I( advfield_new, modulo(i1+2-1,nc_eta1)+1, 1 ) &
               - FIELD_2D_AT_I( advfield_new, modulo(i1-1-1+nc_eta1,nc_eta1)+1, 1 )) / ( delta_eta1 )
       endif
       end if
       do i2 = 2, nc_eta2+1
          eta2 = eta2 + delta_eta2
          ! extract subarray from advection field
          !advfield_1D_2_old(i2) = FIELD_2D_AT_I_V2( advfield_old, i1, i2 )
          !advfield_1D_2_new(i2) = FIELD_2D_AT_I_V2( advfield_new, i1, i2 )
          advfield_1D_2_old ( i2 ) = (FIELD_2D_AT_I( advfield_old, i1, i2 ) &
               - FIELD_2D_AT_I( advfield_old, i1+1, i2 )) / ( delta_eta1)
       if(field_order==4)then 

                 advfield_1D_2_old ( i2 ) = 9._f64/8._f64*advfield_1D_2_old ( i2 )&
+1._f64/24._f64*(FIELD_2D_AT_I( advfield_old, modulo(i1+2-1,nc_eta1)+1, i2 ) &
               - FIELD_2D_AT_I( advfield_old, modulo(i1-1-1+nc_eta1,nc_eta1)+1, i2 )) / ( delta_eta1 )
     
       endif               
               
               
          if (order == 1) then
             advfield_1D_2_new ( i2 ) = 0.0_f64
          else
             advfield_1D_2_new ( i2 ) = (FIELD_2D_AT_I( advfield_new, i1, i2 ) &
                  - FIELD_2D_AT_I( advfield_new, i1+1, i2 )) / delta_eta1
       if(field_order==4)then 

            advfield_1D_2_new ( i2 ) = 9._f64/8._f64*advfield_1D_2_new ( i2 )&
+1._f64/24._f64*(FIELD_2D_AT_I( advfield_new, modulo(i1+2-1,nc_eta1)+1, i2 ) &
               - FIELD_2D_AT_I( advfield_new, modulo(i1-1-1+nc_eta1,nc_eta1)+1, i2 )) / ( delta_eta1 )

        endif          
          end if
          ! compute primitive of distribution function along this line
          primitive2 (i2) = primitive2 (i2-1) &
               + delta_eta2 * FIELD_2D_AT_I( dist_func_2D, i1, i2-1 )
          jacobian(i2) = df_jac_at_i( i1, i2-1 )
       end do
       !i2 = nc_eta2+1
       !print*, 'new', i1,i2, sqrt(x1c_at_i(i1,10)**2+x2c_at_i(i1,10)**2), primitive2 (i2), jacobian(i2)
       call advance_1D_nonuniform( primitive2,        &
                        advfield_1D_2_old, &
                        advfield_1D_2_new, &
                        jacobian,          &
                        order,             &
                        deltat,            &
                        nc_eta2,           &
                        delta_eta2,        &
                        boundary2_type,    &
                        csl_work%spl_eta2, &
                        eta2_out) 
       ! update average value of distribution function in cell using 
       ! difference of primitives
       do i2 = 1, nc_eta2 
          val = (primitive2(i2+1) - primitive2(i2))/delta_eta2 
          call sll_set_df_val( dist_func_2D, i1, i2, val )
          !if (val/df_jac_at_i(i1,i2)>1.) then
          !   print*, 'val', i1,i2, val, primitive2(i2) , primitive2(i2+1), df_jac_at_i(i1,i2), delta_eta2
          !end if
       end do
       eta1 = eta1 + delta_eta1
    end do
  end subroutine csl_advance_2

  !> Performs a 1D conservative semi-Lagragian algorithm
  !> Might be a better option to remove geometry information that should be known by spline object
  !> \param[in] primitive  primitive of functions to be advanced
  !> \param[in] fieldn advection field at time t_n
  !> \param[in] fieldnp1 advection field at time t_n+1
  !> \param[in] order order of time scheme only order 1 and 2 implemented
  !> \param[in] deltat time step
  !> \param[in] eta_min first point in eta grid
  !> \param[in] nc_eta number of cells in eta grid
  !> \param[in] delta_eta cell size of eta grid
  !> \param[in] boundary_type type of the boundary for spline interpolation (periodic or hermite implemented)
  !> \param[in] spline spline object (previously initialized)
  !> \param[out] eta_out origin of cells ending at grid points
  subroutine advance_1D_nonuniform( primitive,     &
                                    fieldn,        &
                                    fieldnp1,      &
                                    jacobian,      &
                                    order,         &
                                    deltat,        &
                                    nc_eta,        & 
                                    delta_eta,     &
                                    boundary_type, &
                                    spline,        &
                                    xi_out)  
!    sll_real64, dimension(:), pointer, intent(inout) :: primitive
!    sll_real64, dimension(:), pointer, intent(in)    :: fieldn
!    sll_real64, dimension(:), pointer, intent(in)    :: fieldnp1
!    sll_real64, dimension(:), pointer, intent(in)    :: jacobian
!    sll_real64, dimension(:), allocatable            :: xi
    sll_real64, dimension(:), pointer                :: primitive
    sll_real64, dimension(:), pointer                :: fieldn
    sll_real64, dimension(:), pointer                :: fieldnp1
    sll_real64, dimension(:), pointer                :: jacobian
    sll_real64, dimension(:), pointer                :: xi
    sll_int32, intent(in)                            :: order
    sll_real64, intent(in)                           :: deltat
    sll_int32, intent(in)                            :: nc_eta
    sll_real64, intent(in)                           :: delta_eta
    sll_int32, intent(in)                            :: boundary_type
    type (cubic_nonunif_spline_1D), pointer           :: spline
    sll_real64, dimension(:), intent(out)   :: xi_out
    
    ! local variables
    sll_real64  :: avg
    sll_real64  :: xi_max
    sll_real64  :: xi_new
    sll_real64  :: lperiod
    sll_int32   :: i
    sll_int32   :: ierr

    ! check array dimensions
    SLL_ASSERT(size(primitive) >= nc_eta+1)
    SLL_ASSERT(size(fieldn)    >= nc_eta+1)
    SLL_ASSERT(size(fieldnp1)  >= nc_eta+1)
    SLL_ASSERT(size(xi_out)   >= nc_eta+1)

    ! array allocation
    SLL_ALLOCATE(xi(nc_eta+1),ierr)

    ! compute xi (primitive of jacobian)
    xi(1) = 0_f64  
    do i = 2, nc_eta+1
       xi(i) = xi(i-1) + jacobian(i)* delta_eta
       !print*, 'xi ', xi(i), xi(i) - xi(i-1)
    end do
    
    select case (boundary_type)
    case (PERIODIC)
       ! average of dist func along the line
       avg = primitive ( nc_eta+1 ) / xi(nc_eta+1)
       ! modify primitive so that it becomes periodic
       do i = 2, nc_eta+1
          primitive ( i ) = primitive ( i ) - avg * xi(i)
       end do
       xi_max = xi(nc_eta + 1)

       !print*, 'avg', avg, xi_max
       !print*, 'a',fieldn

       call implicit_ode_nonuniform( order,   &
                          deltat,             &
                          xi,                 &
                          nc_eta,             &
                          PERIODIC_ODE,       &
                          xi_out,             &
                          fieldn,             &
                          fieldnp1 ) 
       !call compute_spline_1D_periodic( primitive, spline )
       call compute_spline_nonunif( primitive, spline, xi)
       ! interpolate primitive at origin of characteritics
       !call interpolate_array_values( eta_out, primitive, nc_eta+1, spline )
       call interpolate_array_value_nonunif( xi_out, primitive, nc_eta+1, spline)
       ! come back to real primitive by adding average
       if (xi_out(1) > 0.5_f64*(xi_max)) then
          lperiod = -1.0_f64
       else
          lperiod = 0.0_f64
       end if
       xi_new = xi_out(1) +  lperiod*xi_max
       primitive ( 1 ) = primitive ( 1 ) + avg * xi_new
       !if ((xi_new > xi_max) .or. (xi_new <xi(1))) then
       !   print*, 1, xi_new, xi_out(1), primitive(1)
       !end if
       do i = 2, nc_eta+1
          ! We need here to find the points where it has been modified by periodicity
          if (xi_out(i) < xi_out(i-1)) then
             lperiod = lperiod + 1.0_f64
          end if
          xi_new = xi_out(i) +  lperiod*xi_max
          primitive ( i ) = primitive ( i ) + avg * xi_new
          !if (i>98) then
          !   print*, 'iii', i, xi_new, xi_out(i),xi_max, primitive(i), primitive(i-1)
          !endif
       end do

    case (COMPACT)

       call implicit_ode_nonuniform( order,       &
                          deltat,      &
                          xi,     &
                          nc_eta,      &
                          COMPACT_ODE, &
                          xi_out,     &
                          fieldn,      &
                          fieldnp1 ) 
       !call compute_spline_1D_hermite( primitive, spline )
       call compute_spline_nonunif( primitive, spline, xi)
       ! interpolate primitive at origin of characteritics
       call interpolate_array_value_nonunif( xi_out, primitive, nc_eta+1, spline)
       !do i = 2, nc_eta+1
       !   print*, 'iii', i, xi(1), xi_out(i),xi(nc_eta+1), primitive(i), primitive(i-1)
       !end do
    end select
  end subroutine advance_1D_nonuniform

end module sll_csl
