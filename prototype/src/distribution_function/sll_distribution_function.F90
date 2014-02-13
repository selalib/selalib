!------------------------------------------------------------------------------
! SELALIB
!------------------------------------------------------------------------------
!
! MODULE: sll_distribution_function
!
!> @author
!> - Eric
!> - Michel
!> - Pierre
!> - Edwin
!>
!
! DESCRIPTION: 
!
!> @brief
!> Implements the distribution function types
!>
!>@details
!>
!> This module depends on:
!>    - memory
!>    - precision
!>    - assert
!>    - utilities
!>    - constants
!>    - diagnostics
!>    - splines
!>    - mesh_types
!
! REVISION HISTORY:
! DD Mmm YYYY - Initial Version
! TODO_dd_mmm_yyyy - TODO_describe_appropriate_changes - TODO_name
!------------------------------------------------------------------------------
module distribution_function
#include "sll_working_precision.h"
#include "sll_memory.h"
#include "sll_assert.h"
#include "sll_field_2d.h"
  use sll_constants
  use sll_utilities   ! for int2string
  use sll_scalar_field_initializers_base
  implicit none

!#ifdef STDF95
!  type  :: sll_distribution_function_2D
   !  type(sll_mapped_mesh_2d_discrete), pointer  :: mesh
   !  type(cubic_spline_1d_interpolator), pointer :: eta1_interpolator
   !  type(cubic_spline_1d_interpolator), pointer :: eta2_interpolator
   !  sll_real64, dimension(:,:), pointer      :: data
   !  sll_int32                                :: data_position
   !  character(len=64)                        :: name
   !  sll_int32                                :: plot_counter
!     type(scalar_field_2d) :: extend_type
!     sll_real64        :: pmass
!     sll_real64        :: pcharge           
!     sll_real64        :: average 
!  end type  sll_distribution_function_2D
!#else  
#define NEW_TYPE_FOR_DF( new_df_type, extended_type)                 \
  type, extends(extended_type) :: new_df_type;                       \
     sll_real64      :: pmass;                                       \
     sll_real64      :: pcharge;                                     \
     sll_real64      :: average;                                     \
  end type new_df_type

!NEW_TYPE_FOR_DF(sll_distribution_function_2D_t, scalar_field_2d)
!NEW_TYPE_FOR_DF(sll_distribution_function_4D_t, scalar_field_4d)

NEW_TYPE_FOR_DF( sll_distribution_function_2d, scalar_field_2d )


!!$  interface write_distribution_function
!!$     module procedure write_distribution_function_2D, &
!!$                      write_distribution_function_4D
!!$  end interface

#undef NEW_TYPE_FOR_DF
!#endif

contains

#if 1
  subroutine sll_new_distribution_function_2d( &
    this, &
    transf, &
    data_position, &
    name, &
    data_func ) 
    
    class(sll_distribution_function_2D)   :: this
    class(sll_coordinate_transformation_2d_base), pointer :: transf
    procedure(scalar_function_2D)           :: data_func
    sll_int32, intent(in)                   :: data_position
    character(len=*), intent(in)            :: name
    sll_int32                         :: ierr
    sll_int32  :: i1, i2
    sll_real64 :: eta1, eta2
    sll_real64 :: delta1, delta2
    type(sll_logical_mesh_2d), pointer :: mesh

    this%transf => transf
    this%plot_counter = 0
    this%name = name
    this%data_position = data_position
    this%pcharge = 1.0_f64
    this%pmass = 1.0_f64
    mesh => transf%get_logical_mesh()

    if (data_position == NODE_CENTERED_FIELD) then
       SLL_ALLOCATE(this%data(mesh%num_cells1+1,mesh%num_cells2+1), ierr)
       do i2 = 1, mesh%num_cells2+1
          do i1 = 1, mesh%num_cells1+1
             this%data(i1,i2) = data_func(transf%x1_at_node(i1,i2), &
                  transf%x2_at_node(i1,i2))
          end do
       end do
    else if (data_position == CELL_CENTERED_FIELD) then
       SLL_ALLOCATE(this%data(mesh%num_cells1+1,mesh%num_cells2+1), ierr)
       delta1 = 1.0_f64/mesh%num_cells1
       delta2 = 1.0_f64/mesh%num_cells2
       eta2 = 0.5_f64 * delta2
       do i2 = 1, mesh%num_cells2
          eta1 = 0.5_f64 * delta1
          do i1 = 1, mesh%num_cells1
             this%data(i1,i2) = data_func(transf%x1(eta1,eta2), &
                  transf%x2(eta1,eta2)) * transf%jacobian(eta1,eta2)
             eta1 = eta1 + delta1
          end do
          eta2 = eta2 + delta2
       end do
    endif
  end subroutine sll_new_distribution_function_2d

  subroutine initialize_distribution_function_2d( &
    this, &
    mass, &
    charge, &
    field_name, &
    transf, &
    data_position, &
    eta1_interpolator, &
    eta2_interpolator, &
    initializer )

    class(sll_coordinate_transformation_2d_base), pointer :: transf
    class(sll_interpolator_1d_base), pointer            :: eta1_interpolator
    class(sll_interpolator_1d_base), pointer            :: eta2_interpolator
    class(scalar_field_2d_initializer_base), pointer, optional :: initializer
    type(sll_distribution_function_2d), intent(inout)   :: this
    sll_real64, intent(in)                              :: mass
    sll_real64, intent(in)                              :: charge
    character(len=*), intent(in)                        :: field_name
    sll_int32, intent(in)                               :: data_position

    this%pmass = mass
    this%pcharge = charge

    call initialize_scalar_field_2d( &
         this, &
         field_name, &
         transf, &
         data_position, &
         eta1_interpolator, &
         eta2_interpolator, &
         initializer )
  end subroutine initialize_distribution_function_2d
#endif



!!$  function sll_new_distribution_function_4D( mesh_descriptor_x,  &
!!$                                             mesh_descriptor_v,  &
!!$                                             center, name ) 
!!$
!!$    type(sll_distribution_function_4D_t), pointer :: &
!!$         sll_new_distribution_function_4D
!!$    type(mesh_descriptor_2D), pointer :: mesh_descriptor_x
!!$    type(mesh_descriptor_2D), pointer :: mesh_descriptor_v
!!$    sll_int32, intent(in)             :: center
!!$    character(len=*), intent(in)      :: name
!!$    sll_int32                         :: ierr
!!$    !
!!$    SLL_ASSERT(associated(mesh_descriptor_x))
!!$    SLL_ASSERT(associated(mesh_descriptor_v))
!!$    SLL_ALLOCATE(sll_new_distribution_function_4D, ierr)
!!$    sll_new_distribution_function_4D%field => &
!!$         new_field_4D_vec1( mesh_descriptor_x, mesh_descriptor_v )
!!$    sll_new_distribution_function_4D%pcharge = 1.0_f64
!!$    sll_new_distribution_function_4D%pmass = 1.0_f64
!!$    sll_new_distribution_function_4D%plot_counter = 0
!!$    sll_new_distribution_function_4D%center = center
!!$    sll_new_distribution_function_4D%name = name
!!$
!!$  end function sll_new_distribution_function_4D


!!$  subroutine sll_delete_distribution_function( f )
!!$    type(sll_distribution_function_2D_t), pointer      :: f
!!$    sll_int32 :: ierr
!!$    call delete_field_2D_vec1(f%field)
!!$    SLL_DEALLOCATE(f, ierr)
!!$  end subroutine sll_delete_distribution_function

!!$#define NEW_ACCESS_FUNCTION( func_name, typ, slot ) \
!!$  function func_name( f ) ; \
!!$    typ :: func_name ;\
!!$    type(sll_distribution_function_2D_t), pointer :: f ;\
!!$    func_name = f%field%descriptor%slot ;\
!!$  end function func_name
!!$
!!$  NEW_ACCESS_FUNCTION(get_df_nc_eta1,     sll_int32,  nc_eta1)
!!$  NEW_ACCESS_FUNCTION(get_df_eta1_min,   sll_real64, eta1_min)
!!$  NEW_ACCESS_FUNCTION(get_df_eta1_max,   sll_real64, eta1_max)
!!$  NEW_ACCESS_FUNCTION(get_df_delta_eta1, sll_real64, delta_eta1)   
!!$  NEW_ACCESS_FUNCTION(get_df_nc_eta2,     sll_int32,  nc_eta2)
!!$  NEW_ACCESS_FUNCTION(get_df_eta2_min,   sll_real64, eta2_min)
!!$  NEW_ACCESS_FUNCTION(get_df_eta2_max,   sll_real64, eta2_max)
!!$  NEW_ACCESS_FUNCTION(get_df_delta_eta2, sll_real64, delta_eta2)   
!!$  NEW_ACCESS_FUNCTION(get_df_boundary1_type, sll_int32, boundary1_type)
!!$  NEW_ACCESS_FUNCTION(get_df_boundary2_type, sll_int32, boundary2_type)
!!$#undef NEW_ACCESS_FUNCTION
!!$
!!$  function get_df_x1 ( f )
!!$    procedure(scalar_function_2D), pointer        :: get_df_x1
!!$    type(sll_distribution_function_2D_t), pointer :: f
!!$    get_df_x1 => f%mesh%x1
!!$  end function get_df_x1
!!$
!!$  function get_df_x2 ( f )
!!$    procedure(scalar_function_2D), pointer        :: get_df_x2
!!$    type(sll_distribution_function_2D_t), pointer :: f
!!$    get_df_x2 => f%field%descriptor%geom%x2
!!$  end function get_df_x2
!!$
!!$  function get_df_eta1 ( f )
!!$    procedure(scalar_function_2D), pointer        :: get_df_eta1
!!$    type(sll_distribution_function_2D_t), pointer :: f
!!$    get_df_eta1 => f%field%descriptor%geom%eta1
!!$  end function get_df_eta1
!!$
!!$  function get_df_eta2 ( f )
!!$    procedure(scalar_function_2D), pointer        :: get_df_eta2
!!$    type(sll_distribution_function_2D_t), pointer :: f
!!$    get_df_eta2 => f%field%descriptor%geom%eta2
!!$  end function get_df_eta2
!!$
!!$  function get_df_jac11 ( f )
!!$    procedure(scalar_function_2D), pointer        :: get_df_jac11
!!$    type(sll_distribution_function_2D_t), pointer :: f
!!$    get_df_jac11 => f%field%descriptor%geom%Jacobian11
!!$  end function get_df_jac11
!!$
!!$  function get_df_jac12 ( f )
!!$    procedure(scalar_function_2D), pointer        :: get_df_jac12
!!$    type(sll_distribution_function_2D_t), pointer :: f
!!$    get_df_jac12 => f%field%descriptor%geom%Jacobian12
!!$  end function get_df_jac12
!!$
!!$  function get_df_jac21 ( f )
!!$    procedure(scalar_function_2D), pointer        :: get_df_jac21
!!$    type(sll_distribution_function_2D_t), pointer :: f
!!$    get_df_jac21 => f%field%descriptor%geom%Jacobian21
!!$  end function get_df_jac21
!!$
!!$  function get_df_jac22 ( f )
!!$    procedure(scalar_function_2D), pointer        :: get_df_jac22
!!$    type(sll_distribution_function_2D_t), pointer :: f
!!$    get_df_jac22 => f%field%descriptor%geom%Jacobian22
!!$  end function get_df_jac22  
!!$
!!$  function get_df_jac ( f )
!!$    procedure(scalar_function_2D), pointer        :: get_df_jac
!!$    type(sll_distribution_function_2D_t), pointer :: f
!!$    get_df_jac => f%field%descriptor%geom%Jacobian
!!$  end function get_df_jac
!!$
!!$  function get_df_x1_at_i ( f )
!!$    sll_real64, dimension(:,:), pointer        :: get_df_x1_at_i
!!$    type(sll_distribution_function_2D_t), pointer :: f
!!$    get_df_x1_at_i => f%field%descriptor%geom%x1_at_i
!!$  end function get_df_x1_at_i
!!$
!!$  function get_df_x2_at_i ( f )
!!$    sll_real64, dimension(:,:), pointer        :: get_df_x2_at_i
!!$    type(sll_distribution_function_2D_t), pointer :: f
!!$    get_df_x2_at_i => f%field%descriptor%geom%x2_at_i
!!$  end function get_df_x2_at_i
!!$
!!$  function get_df_x1c_at_i ( f )
!!$    sll_real64, dimension(:,:), pointer        :: get_df_x1c_at_i
!!$    type(sll_distribution_function_2D_t), pointer :: f
!!$    get_df_x1c_at_i => f%field%descriptor%geom%x1c_at_i
!!$  end function get_df_x1c_at_i
!!$
!!$  function get_df_x2c_at_i ( f )
!!$    sll_real64, dimension(:,:), pointer        :: get_df_x2c_at_i
!!$    type(sll_distribution_function_2D_t), pointer :: f
!!$    get_df_x2c_at_i => f%field%descriptor%geom%x2c_at_i
!!$  end function get_df_x2c_at_i
!!$
!!$  function get_df_jac_at_i ( f )
!!$    sll_real64, dimension(:,:), pointer        :: get_df_jac_at_i
!!$    type(sll_distribution_function_2D_t), pointer :: f
!!$    get_df_jac_at_i => f%field%descriptor%geom%Jacobian_at_i
!!$  end function get_df_jac_at_i
!!$
!!$
!!$  function sll_get_df_val( f, i, j )
!!$    sll_real64 :: sll_get_df_val
!!$    type(sll_distribution_function_2D_t), pointer      :: f
!!$    sll_int32 :: i, j
!!$    sll_get_df_val = f%field%data(i,j)
!!$  end function sll_get_df_val
!!$
!!$  subroutine sll_set_df_val( f, i, j, val )
!!$    type(sll_distribution_function_2D_t), pointer      :: f
!!$    sll_int32 :: i, j
!!$    sll_real64 :: val
!!$    f%field%data(i,j) = val
!!$  end subroutine sll_set_df_val
!!$
!!$  subroutine sll_init_distribution_function_2D( dist_func_2D, test_case)
!!$    type(sll_distribution_function_2D_t), pointer      :: dist_func_2D
!!$    sll_int32  :: test_case
!!$    ! local variables
!!$    procedure(scalar_function_2D), pointer :: x1, x2, jac
!!$    sll_real64, dimension(:,:), pointer :: x1c_at_i, x2c_at_i, jac_at_i
!!$    sll_int32 :: nc_eta1, nc_eta2, i1, i2
!!$    sll_real64 :: delta_eta1, delta_eta2,  eta1_min, eta2_min
!!$    sll_real64 :: x, vx, xx, vv, eps, kx, xi, v0, fval, xoffset, voffset
!!$    sll_real64 :: eta1
!!$    sll_real64 :: eta2
!!$    sll_real64 :: avg, avg_jac
!!$    
!!$    nc_eta1 = get_df_nc_eta1( dist_func_2D ) 
!!$    delta_eta1 = get_df_delta_eta1( dist_func_2D )
!!$    eta1_min = get_df_eta1_min( dist_func_2D )
!!$    nc_eta2 = get_df_nc_eta2( dist_func_2D ) 
!!$    delta_eta2 = get_df_delta_eta2( dist_func_2D )
!!$    eta2_min = get_df_eta2_min( dist_func_2D )
!!$    x1 => get_df_x1 ( dist_func_2D )
!!$    x2 => get_df_x2 ( dist_func_2D )
!!$    jac => get_df_jac ( dist_func_2D )
!!$    x1c_at_i => get_df_x1c_at_i ( dist_func_2D )
!!$    x2c_at_i => get_df_x2c_at_i ( dist_func_2D )
!!$    jac_at_i => get_df_jac_at_i ( dist_func_2D )
!!$
!!$    if (dist_func_2D%center==CELL_CENTERED_DF) then ! half cell offset
!!$       eta1_min = eta1_min + 0.5_f64 * delta_eta1
!!$       eta2_min = eta2_min + 0.5_f64 * delta_eta2
!!$    end if
!!$
!!$   
!!$  end subroutine sll_init_distribution_function_2D
!!$  
!!$  subroutine sll_init_distribution_function_4D( dist_func_4D, test_case)
!!$
!!$    type(sll_distribution_function_4D_t), pointer :: dist_func_4D
!!$    sll_int32  :: test_case
!!$
!!$    sll_int32 :: nnode_x1, nnode_x2, nnode_v1, nnode_v2
!!$    sll_real64 :: delta_x1, delta_x2,  x1_min, x2_min
!!$    sll_real64 :: delta_v1, delta_v2,  v1_min, v2_min
!!$    sll_real64 :: x1, v1, x2, v2, eps, kx, ky, xi, vsq
!!$    sll_int32  :: ix, jx, iv, jv
!!$    
!!$    x1_min =   dist_func_4D%field%descriptor_x%eta1_min
!!$    nnode_x1 = dist_func_4D%field%descriptor_x%nc_eta1+1
!!$    delta_x1 = dist_func_4D%field%descriptor_x%delta_eta1
!!$
!!$    x2_min =   dist_func_4D%field%descriptor_x%eta2_min
!!$    nnode_x2 = dist_func_4D%field%descriptor_x%nc_eta2+1
!!$    delta_x2 = dist_func_4D%field%descriptor_x%delta_eta2
!!$
!!$    v1_min =   dist_func_4D%field%descriptor_v%eta1_min
!!$    nnode_v1 = dist_func_4D%field%descriptor_v%nc_eta1+1
!!$    delta_v1 = dist_func_4D%field%descriptor_v%delta_eta1
!!$
!!$    v2_min =   dist_func_4D%field%descriptor_v%eta2_min
!!$    nnode_v2 = dist_func_4D%field%descriptor_v%nc_eta2+1
!!$    delta_v2 = dist_func_4D%field%descriptor_v%delta_eta2
!!$
!!$    select case (test_case)
!!$    case (LANDAU)
!!$       xi  = 0.9
!!$       eps = 0.05
!!$       kx  = 2.0_f64*sll_pi/(nnode_x1*nnode_x2)
!!$       ky  = 2.0_f64*sll_pi/(nnode_v1*nnode_v2)
!!$       do jv=1, nnode_v2
!!$          v2 = v2_min+(jv-1)*delta_v2
!!$          do iv=1, nnode_v1
!!$             v1 = v1_min+(iv-1)*delta_v1
!!$             vsq = v1*v1+v2*v2
!!$             do jx=1, nnode_x2
!!$                x2=x2_min+(jx-1)*delta_x2
!!$                do ix=1, nnode_x1
!!$                   x1=x1_min+(ix-1)*delta_x1
!!$                   dist_func_4D%field%data(ix,jx,iv,jv)= &
!!$                      (1+eps*cos(kx*x1))*1/(2.*sll_pi)*exp(-.5*vsq)
!!$                end do
!!$             end do
!!$          end do
!!$       end do
!!$    end select
!!$  end subroutine sll_init_distribution_function_4D
!!$
!!$    ! compute integral of f with respect to x2 (-> rho)
!!$    ! using a trapezoidal rule on a uniform grid of physical space
!!$  subroutine compute_rho(dist_func_2D,rho,npoints)
!!$    type(sll_distribution_function_2D_t), pointer      :: dist_func_2D
!!$    type(field_1D_vec1), pointer                       :: rho 
!!$    sll_int32                                          :: npoints ! number of integration points
!!$
!!$    ! local variables
!!$    procedure(scalar_function_2D), pointer :: x1f
!!$    procedure(scalar_function_2D), pointer :: x2f
!!$    procedure(scalar_function_2D), pointer :: eta1f
!!$    procedure(scalar_function_2D), pointer :: eta2f
!!$    sll_int32 :: nc_eta1
!!$    sll_int32 :: nc_eta2
!!$    sll_int32 :: nc_rho
!!$    sll_real64 :: eta1_min
!!$    sll_real64 :: eta2_min
!!$    sll_real64 :: delta_eta1
!!$    sll_real64 :: delta_eta2
!!$    sll_real64 :: x1min_rho
!!$    sll_real64 :: delta_rho
!!$    sll_int32 :: i
!!$    sll_int32 :: i1
!!$    sll_int32 :: i2
!!$    sll_real64 :: x2min
!!$    sll_real64 :: x2max
!!$    sll_real64 :: x1
!!$    sll_real64 :: x2
!!$    sll_real64 :: delta_int
!!$    sll_real64 :: sum
!!$    sll_real64 :: eta1
!!$    sll_real64 :: eta2
!!$
!!$    ! get mesh data attached to f
!!$    nc_eta1 = get_df_nc_eta1( dist_func_2D ) 
!!$    delta_eta1 = get_df_delta_eta1( dist_func_2D )
!!$    eta1_min = get_df_eta1_min( dist_func_2D )
!!$    nc_eta2 = get_df_nc_eta2( dist_func_2D ) 
!!$    delta_eta2 = get_df_delta_eta2( dist_func_2D )
!!$    eta2_min = get_df_eta2_min( dist_func_2D )
!!$    x1f => get_df_x1 ( dist_func_2D )
!!$    x2f => get_df_x2 ( dist_func_2D )
!!$    eta1f => get_df_eta1 ( dist_func_2D )
!!$    eta2f => get_df_eta2 ( dist_func_2D )
!!$
!!$    ! get mesh data attached to rho
!!$    nc_rho = GET_FIELD_NC_ETA1( rho )
!!$    x1min_rho = GET_FIELD_ETA1_MIN( rho )
!!$    delta_rho = GET_FIELD_DELTA_ETA1( rho )
!!$
!!$    ! find minimal and maximal values of x2
!!$    x2min = x2f(1.0_f64,1.0_f64)
!!$    x2max = x2min
!!$    eta1 = eta1_min
!!$    do i1 = 1, nc_eta1 + 1
!!$       eta2 = eta2_min
!!$       do i2 = 1, nc_eta2 + 1
!!$          x2 = x2f(eta1,eta2)
!!$          if ( x2 >  x2max ) then
!!$             x2max = x2
!!$          else if ( x2 <  x2min ) then
!!$             x2min = x2
!!$          end if
!!$          eta2 = eta2 + delta_eta2
!!$       end do
!!$       eta1 = eta1 + delta_eta1
!!$    end do
!!$    ! set delta_int the subdivision step
!!$    delta_int = (x2max-x2min)/npoints
!!$    x1 = x1min_rho
!!$    do i1 = 1, nc_rho + 1
!!$       sum = 0.0_f64
!!$       x2 = x2min
!!$       do i2 = 1, npoints
!!$          !sum = sum + FIELD_2D_AT_X( dist_func_2D, eta1f(x1, x2), eta2f(x1,x2) )
!!$          ! FIELD_2D_AT_X needs to be defined and implemented using linear or 2D cubic spline interpolation 
!!$          ! beware of case where eta1f and eta2f fall outside the grid (in this case 0 should be returned)
!!$          x2 = x2 + delta_int
!!$       end do
!!$       SET_FIELD_1D_VALUE_AT_I( rho, i1, delta_int * sum ) 
!!$       x1 = x1 + delta_rho 
!!$    end do
!!$
!!$  end subroutine compute_rho
!!$
!!$  subroutine write_distribution_function_2D ( f )
!!$    type(sll_distribution_function_2D_t), pointer      :: f
!!$    character(len=4) :: counter
!!$    character(64) :: name
!!$    logical, parameter   :: jacobian = .true.
!!$    call int2string(f%plot_counter,counter)
!!$    name = trim(f%name)//counter
!!$    call write_field_2d_vec1 ( f%field, name, jacobian, f%average, f%center )
!!$    f%plot_counter = f%plot_counter + 1
!!$  end subroutine write_distribution_function_2D
!!$
!!$
!!$  subroutine write_distribution_function_4D ( f )
!!$    type(sll_distribution_function_4D_t), pointer      :: f
!!$    character(len=4) :: counter
!!$    character(64) :: name
!!$    type(mesh_descriptor_2D), pointer :: mesh_x
!!$    type(mesh_descriptor_2D), pointer :: mesh_v
!!$    sll_int32  :: nnode_x1, nnode_x2
!!$    sll_int32  :: nnode_v1, nnode_v2
!!$    sll_real64, dimension(:,:), pointer :: val
!!$    sll_int32  :: error
!!$    sll_int32  :: ix
!!$    sll_int32  :: jx
!!$
!!$    call int2string(f%plot_counter,counter)
!!$    name = trim(f%name)//counter
!!$
!!$    mesh_x => f%field%descriptor_x
!!$    mesh_v => f%field%descriptor_v
!!$
!!$    SLL_ASSERT(associated(mesh_x))
!!$    SLL_ASSERT(associated(mesh_v))
!!$
!!$    nnode_x1 = mesh_x%nc_eta1 + 1
!!$    nnode_x2 = mesh_x%nc_eta2 + 1
!!$    nnode_v1 = mesh_v%nc_eta1 + 1
!!$    nnode_v2 = mesh_v%nc_eta2 + 1
!!$
!!$    SLL_ALLOCATE(val(nnode_x1,nnode_x2), error)
!!$    do ix = 1, nnode_x1
!!$       do jx = 1, nnode_x2
!!$          val(ix,jx) = sum(f%field%data(ix,jx,:,:))
!!$       end do
!!$    end do
!!$
!!$    call write_vec1d(val,nnode_x1,nnode_x2,name,"mesh_x",f%center)
!!$    
!!$    f%plot_counter = f%plot_counter + 1
!!$  end subroutine write_distribution_function_4D

end module distribution_function
