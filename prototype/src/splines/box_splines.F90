!> @brief  
!> Box-Splines of degree 2 using formula from : 
!> @Condat2006 :Three-directional box splines
!> 
module box_splines
#include "sll_working_precision.h"
#include "sll_memory.h"
#include "sll_assert.h"
#include "sll_splines.h"
#include "sll_utilities.h"
use hex_logical_meshes!, only:find_neighbour


  implicit none

  !> Spline object
  type linear_box_spline_2D
      sll_int32  :: num_pts    !< number of points in any direction from origin
      sll_real64 :: radius     ! distance between origin and external vertex
      sll_real64 :: center_x1  ! x1 cartesian coordinate of the origin
      sll_real64 :: center_x2  ! x2 cartesian coordinate of the origin
      sll_real64 :: delta      ! cell spacing
      !generator vectors (r1, r2, r3) coordinates
      sll_real64 :: r1_x1 
      sll_real64 :: r1_x2 
      sll_real64 :: r2_x1 
      sll_real64 :: r2_x2 
      sll_real64 :: r3_x1 
      sll_real64 :: r3_x2 
      sll_int32 SLL_PRIV                           :: x1_bc_type
      sll_int32 SLL_PRIV                           :: x2_bc_type
      !< spline coefficients
      sll_real64, dimension(:), pointer :: coeffs
      ! Boundary conditions
      sll_real64, dimension(:), pointer SLL_PRIV   :: x1_min_slopes => null()
      sll_real64, dimension(:), pointer SLL_PRIV   :: x1_max_slopes => null()
      sll_real64, dimension(:), pointer SLL_PRIV   :: x2_min_slopes => null()
      sll_real64, dimension(:), pointer SLL_PRIV   :: x2_max_slopes => null()
  end type linear_box_spline_2D

  ! Some useful macros that should probably be put in a different file to 
  ! make them more widely accessible. Here we use them to compute default
  ! values for the slopes. 
  ! Actually the proper way to factor this is with a function that takes 
  ! an array segment of a certain length and returns the slope. Such
  ! function(s) could be made to work in 1D or 2D by passing the right
  ! array segment. Keep this in mind...

    ! (-25/12, 4, -3, 4/3, -1/4) stencil
#define FORWARD_FD_5PT( f, r_delta ) \
    r_delta*(-(25.0_f64/12.0_f64)*f(1) + 4.0_f64*f(2) -3.0_f64*f(3) + (4.0_f64/3.0_f64)*f(4) - 0.25_f64*f(5))

    ! (1/4, -4/3, 3, -4, 25/12) stencil
#define BACKWARD_FD_5PT( f, r_delta, np ) \
    r_delta*(0.25_f64*f(np-4) - (4.0_f64/3.0_f64)*f(np-3) + 3.0_f64*f(np-2) - 4.0_f64*f(np-1) + (25.0_f64/12.0_f64)*f(np) )

    ! (-3/2, 2, -1/2 ) stencil
#define FORWARD_FD_3PT( f, r_delta ) \
    r_delta*(-1.5_f64*f(1) + 2.0_f64*f(2) - 0.5_f64*f(3))

    ! ( 1/2, -2, 3/2 ) stencil
#define BACKWARD_FD_3PT( f, r_delta, np ) \
    r_delta*(0.5_f64*f(np-2)-2.0_f64*f(np-1) +1.5_f64*f(np))

  
  
contains  ! ****************************************************************


#define MAKE_GET_SLOT_FUNCTION( fname, datatype, slot, ret_type )    \
  function fname( spline_obj ) result(val);                \
    type(datatype), pointer :: spline_obj;                 \
    ret_type :: val;                                       \
    val = spline_obj%slot;                                 \
  end function fname

MAKE_GET_SLOT_FUNCTION(get_radius_lbs2d, linear_box_spline_2d, radius, sll_real64 )
MAKE_GET_SLOT_FUNCTION(get_center_x1_lbs2d, linear_box_spline_2d, center_x1, sll_real64 )
MAKE_GET_SLOT_FUNCTION(get_center_x2_lbs2d, linear_box_spline_2d, center_x2, sll_real64 )
MAKE_GET_SLOT_FUNCTION(get_r1_x1_lbs2d, linear_box_spline_2d, r1_x1, sll_real64 )
MAKE_GET_SLOT_FUNCTION(get_r1_x2_lbs2d, linear_box_spline_2d, r1_x2, sll_real64 )
MAKE_GET_SLOT_FUNCTION(get_r2_x1_lbs2d, linear_box_spline_2d, r2_x1, sll_real64 )
MAKE_GET_SLOT_FUNCTION(get_r2_x2_lbs2d, linear_box_spline_2d, r2_x2, sll_real64 )
MAKE_GET_SLOT_FUNCTION(get_r3_x1_lbs2d, linear_box_spline_2d, r3_x1, sll_real64 )
MAKE_GET_SLOT_FUNCTION(get_r3_x2_lbs2d, linear_box_spline_2d, r3_x2, sll_real64 )



  function new_linear_box_spline_2D( &
    num_pts,    &
    radius,       &
    x1_bc_type,   &
    x2_bc_type,   &
    const_slope_x1_min, &
    const_slope_x1_max, &
    const_slope_x2_min, &
    const_slope_x2_max, &
    x1_min_slopes, &
    x1_max_slopes, &
    x2_min_slopes, &
    x2_max_slopes)

    type(linear_box_spline_2D), pointer  :: new_linear_box_spline_2D
    sll_int32,  intent(in)               :: num_pts
    sll_real64, intent(in)               :: radius
    sll_int32,  intent(in)               :: x1_bc_type
    sll_int32,  intent(in)               :: x2_bc_type
    sll_real64, intent(in), optional     :: const_slope_x1_min
    sll_real64, intent(in), optional     :: const_slope_x1_max
    sll_real64, intent(in), optional     :: const_slope_x2_min
    sll_real64, intent(in), optional     :: const_slope_x2_max
    sll_real64, intent(in), dimension(:), optional :: x1_min_slopes
    sll_real64, intent(in), dimension(:), optional :: x1_max_slopes
    sll_real64, intent(in), dimension(:), optional :: x2_min_slopes
    sll_real64, intent(in), dimension(:), optional :: x2_max_slopes
    sll_int32                            :: bc_selector
    sll_int32                            :: ierr
    sll_int32                            :: num_pts_tot


    SLL_ALLOCATE( new_linear_box_spline_2D, ierr )
    new_linear_box_spline_2D%num_pts = num_pts
    new_linear_box_spline_2D%radius    = radius
    new_linear_box_spline_2D%delta     = radius/real((num_pts-1),f64)
    new_linear_box_spline_2D%x1_bc_type = x1_bc_type
    new_linear_box_spline_2D%x2_bc_type = x2_bc_type


!    if (num_pts .le. NUM_TERMS) then
!      print *, 'ERROR, new_linear_box_spline_2D: Because of the algorithm used, this ', &
!       'function is meant to be used with arrays that are at least of size = 28'
!       STOP 'new_linear_box_spline_2D()'
!    end if
    if (radius .le. 0.)  then
       print *, 'ERROR , new_linear_box_spline_2D : negative radius '
       STOP
    end if

    ! Check that slope arrays are of the right size. Consider making this 
    ! something more permanent than an assertion. Does this compile if the
    ! assertions are turned off???
    if( present(x1_min_slopes) ) then
       SLL_ASSERT(size(x1_min_slopes) .ge. num_pts )
    end if
    if( present(x1_max_slopes) ) then
       SLL_ASSERT(size(x1_max_slopes) .ge. num_pts )
    end if
    if( present(x2_min_slopes) ) then
       SLL_ASSERT(size(x2_min_slopes) .ge. num_pts )
    end if
    if( present(x2_max_slopes) ) then
       SLL_ASSERT(size(x2_max_slopes) .ge. num_pts )
    end if


    ! Treat the bc_selector variable essentially like a bit field, to 
    ! accumulate the information on the different boundary conditions
    ! given. This scheme allows to add more types of boundary conditions
    ! if necessary.
    bc_selector = 0
    if( x1_bc_type .eq. SLL_NEUMANN ) then 
       bc_selector = bc_selector + 1
    end if
    if( x1_bc_type .eq. SLL_PERIODIC ) then
       bc_selector = bc_selector + 2
    end if
    if( x2_bc_type .eq. SLL_NEUMANN ) then 
       bc_selector = bc_selector + 4
    end if

    select case (bc_selector)
    case ( 6 ) 
       ! Periodic on x1 and Neumann on x2
       if( present(x1_min_slopes) .or. present(x1_max_slopes) .or. &
           present(const_slope_x1_min) .or. present(const_slope_x1_max) )then

          print *, 'new_linear_box_spline_2D(): ', &
               'it is not allowed to specify the end slopes at x1', &
               'in the case of periodic-neumann boundary conditions.', &
               'Exiting program...'
          STOP 'new_linear_box_spline_2D'
       else
          ! X1 slopes are not needed
          new_linear_box_spline_2d%x1_min_slopes => null()
          new_linear_box_spline_2d%x1_max_slopes => null()
          ! But X2 slopes are.
          SLL_ALLOCATE(new_linear_box_spline_2d%x2_min_slopes(num_pts),ierr)
          SLL_ALLOCATE(new_linear_box_spline_2d%x2_max_slopes(num_pts),ierr)
       end if
    case ( 5 ) 
       !  Neumann on x1 and x2 
       if( &
          present(x2_min_slopes) .or. present(x2_max_slopes) .or. &
          present(const_slope_x2_min) .or. present(const_slope_x2_max) ) then
          print *, 'new_linear_box_spline_2D(): hermite-periodic case, it is not ', &
               'allowed to specify the end slopes in the case of periodic ', &
               'boundary conditions.', &
               'Exiting program...'
          STOP 'new_linear_box_spline_2D'
       end if

       ! X1 and X2 slopes are needed
       SLL_ALLOCATE(new_linear_box_spline_2d%x1_min_slopes(num_pts),ierr)
       SLL_ALLOCATE(new_linear_box_spline_2d%x1_max_slopes(num_pts),ierr)
       SLL_ALLOCATE(new_linear_box_spline_2d%x2_min_slopes(num_pts),ierr)
       SLL_ALLOCATE(new_linear_box_spline_2d%x2_max_slopes(num_pts),ierr)


    case default
       print *, 'ERROR: new_linear_box_spline_2D(): ', &
            'did not recognize given boundary conditions.'
       STOP
    end select

    ! rk: num_pts-1 = number of nested hexagones
    num_pts_tot = 6*num_pts*(num_pts-1)/2+1 
    SLL_ALLOCATE( new_linear_box_spline_2D%coeffs(0:num_pts_tot), ierr )
  end function new_linear_box_spline_2D


  ! Pre-filter to compute the box splines coefficients
  ! Reference : @Condat and Van De Ville (2007)
  !             "Quasi-interpolating spline models 
  !             for hexagonally-sampled data."
  function pre_filter_piir2(local_index) result(weight)
      sll_int32, intent(in)     :: local_index
      sll_real64                :: weight

      if (local_index .eq. 0) then
         weight = real((0),f64)
      else if (local_index .lt. 7) then
         weight = real((256./6912.), f64)
      else if (local_index .lt. 13) then
         if (modulo(local_index, 2) .eq. 1) then
            weight = real((1./13824.), f64)
         else
            weight = real((11./6912.), f64)
         end if
      else
         weight = 0.
      end if
   end function pre_filter_piir2



  subroutine compute_linear_box_spline_2D( data, spline )
    sll_real64, dimension(:), intent(in), target :: data  ! data to be fit
    type(linear_box_spline_2D), pointer      :: spline
    sll_int32  :: bc1
    sll_int32  :: bc2
    sll_int32  :: num_pts_tot
    sll_int32  :: bc_selector
    sll_int32  :: i,k
    sll_int32  :: nei
    if( .not. associated(spline) ) then
       ! FIXME: THROW ERROR
       print *, 'ERROR: compute_linear_box_spline_2D(): ', &
            'uninitialized spline object passed as argument. '
       print *, "Exiting..."
       STOP
    end if
 
    bc1 = spline%x1_bc_type
    bc2 = spline%x2_bc_type

    ! Treat the bc_selector variable essentially like a bit field, to 
    ! accumulate the information on the different boundary conditions
    ! given. This scheme allows to add more types of boundary conditions
    ! if necessary.
    bc_selector = 0

    ! We make every case explicit to facilitate adding more BC types in
    ! the future.
    if( spline%x1_bc_type .eq. SLL_NEUMANN ) then 
       bc_selector = bc_selector + 1
    end if
    if( spline%x1_bc_type .eq. SLL_PERIODIC ) then
       bc_selector = bc_selector + 2
    end if
    if( spline%x2_bc_type .eq. SLL_NEUMANN ) then 
       bc_selector = bc_selector + 4
    end if

    select case (bc_selector)
       case ( 5 ) 
          ! both boundary condition types are neumann
          call compute_spline_2D_neum_neum( data, spline )
       case ( 6 )
          ! periodic in x1 and neumann in x2
          call compute_spline_2D_prdc_neum( data, spline )
       case default
          print *, 'ERROR: compute_linear_box_spline_2D(): ', &
            'did not recognize given boundary condition combination.'
       STOP
    end select

  
    num_pts_tot = 6*spline%num_pts*(spline%num_pts-1)/2+1 
    do i = 0, num_pts_tot
       spline%coeffs(i) = real(0,f64)
       do k = 0, 18
          nei = find_neighbour(i,k)
          if (nei .le. num_pts_tot) then
             spline%coeffs(i) = spline%coeffs(i) + data(nei) * & 
                                pre_filter_piir2(k)
          !TODO : ELSE ??? BC here !!
          end if
       end do
    end do
   
    !END DO : QUE FAIRE AVEC LES CONDITIONS LIMITES

  end subroutine compute_linear_box_spline_2D


  function chi2(spline, x1, x2)
    ! Reference : @Condat and Van De Ville (2006)
    !             "Three directional Box Splines:
    !             Characterization and Efficient Evaluation."
    sll_real64                          :: chi2
    sll_real64, intent(in)              :: x1
    sll_real64, intent(in)              :: x2
    type(linear_box_spline_2D), pointer :: spline
    sll_real64              :: u
    sll_real64              :: v
    sll_real64              :: delta
    sll_real64              :: g

    delta = spline%r1_x1 * spline%r2_x2 - spline%r2_x1 * spline%r1_x2
    u = (-spline%r2_x2 * abs(x1) - spline%r2_x1 * abs(x2))/delta
    v = ( spline%r1_x1 * abs(x2) + spline%r1_x2 * abs(x1))/delta
    
    if (u.lt.0.0) then
       ! r2 symmetry
       u =   - u
       v = v + u
    end if

    if (2.*u.lt.v) then
       !r2 + r3 symmetry
       u = v - u
    end if

    g = u-v/2.0
      
    if (v.gt.2.0) then
      chi2 = 0.0
    else if (v.lt.1.0) then
       chi2 = 0.5 + ((5/3. - v/8.0)*v-3.0)*v*v/4.0 + &
              ((1.-v/4.)*v+g*g/6.-1.)*g*g
    else if (u.gt.1.0) then
       chi2 = (v-2.)*(v-2.)*(v-2.)*(g-1.)/6.0
    else 
       chi2 = 5./6. + ((1.+ (1./3. - v/8.0)*v)*v/4.0 - 1.0)*v + &
              ((1-v/4.)*v + g*g/6.0 - 1.)*g*g
    end if
    end function chi2

  function hex_interpolate_value(mesh, x1, x2, spline) result(val)
    type(hex_logical_mesh_2d), pointer  :: mesh
    sll_real64                          :: val
    intrinsic                           :: associated, int, real
    sll_real64, intent(in)              :: x1
    sll_real64, intent(in)              :: x2
    type(linear_box_spline_2D), pointer :: spline
    sll_real64              :: xm1
    sll_real64              :: xm2
    sll_real64              :: r1_x1, r1_x2
    sll_real64              :: r2_x1, r2_x2
    sll_int32               :: k1
    sll_int32               :: k2
    sll_int32               :: i

    r1_x1 = spline%r1_x1
    r1_x2 = spline%r1_x2
    r2_x1 = spline%r2_x1
    r2_x2 = spline%r2_x2

    !Lower left corner of encapsulating rhomboid
    k1 = from_cart_index_k1(mesh, x1, x2)
    k2 = from_cart_index_k2(mesh, x1, x2)
    i  = global_index(k1, k2) 

    xm1 = x1-r1_x1*k1-r2_x1*k2
    xm2 = x2-r1_x2*k1-r2_x2*k2
    val = spline%coeffs(i)*chi2(spline,xm1, xm2)

    !Lower right corner of encapsulating rhomboid
    k1 = from_cart_index_k1(mesh, x1, x2) + 1
    k2 = from_cart_index_k2(mesh, x1, x2)
    i  = global_index(k1, k2) 

    xm1 = x1-r1_x1*k1-r2_x1*k2 - r1_x1
    xm2 = x2-r1_x2*k1-r2_x2*k2 - r1_x2
    val = val + spline%coeffs(i)*chi2(spline,xm1, xm2)
    
    !Upper left corner of encapsulating rhomboid
    k1 = from_cart_index_k1(mesh, x1, x2) 
    k2 = from_cart_index_k2(mesh, x1, x2) + 1
    i  = global_index(k1, k2) 

    xm1 = x1-r1_x1*k1-r2_x1*k2 - r2_x1
    xm2 = x2-r1_x2*k1-r2_x2*k2 - r2_x2
    val = val + spline%coeffs(i)*chi2(spline,xm1, xm2)

    !Upper right corner of encapsulating rhomboid
    k1 = from_cart_index_k1(mesh, x1, x2) + 1
    k2 = from_cart_index_k2(mesh, x1, x2) + 1
    i  = global_index(k1, k2) 

    xm1 = x1-r1_x1*k1-r2_x1*k2 - r1_x1 - r2_x1
    xm2 = x2-r1_x2*k1-r2_x2*k2 - r1_x2 - r2_x2
    val = val + spline%coeffs(i)*chi2(spline,xm1, xm2)

  end function hex_interpolate_value

end module box_splines

