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
use hex_mesh!, only:find_neighbour


  implicit none

  !> Spline object
  type linear_box_spline_2d
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

  end type linear_box_spline_2d


    
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



  function new_linear_box_spline_2d( &
    num_pts,    &
    radius,       &
    x1_bc_type,   &
    x2_bc_type)

    type(linear_box_spline_2d), pointer  :: new_linear_box_spline_2d
    sll_int32,  intent(in)               :: num_pts
    sll_real64, intent(in)               :: radius
    sll_int32,  intent(in)               :: x1_bc_type
    sll_int32,  intent(in)               :: x2_bc_type
    sll_int32                            :: bc_selector
    sll_int32                            :: ierr
    sll_int32                            :: num_pts_tot


    SLL_ALLOCATE( new_linear_box_spline_2d, ierr )
    new_linear_box_spline_2d%num_pts = num_pts
    new_linear_box_spline_2d%radius    = radius
    new_linear_box_spline_2d%delta     = radius/real((num_pts-1),f64)
    new_linear_box_spline_2d%x1_bc_type = x1_bc_type
    new_linear_box_spline_2d%x2_bc_type = x2_bc_type


!    if (num_pts .le. NUM_TERMS) then
!      print *, 'ERROR, new_linear_box_spline_2d: Because of the algorithm used, this ', &
!       'function is meant to be used with arrays that are at least of size = 28'
!       STOP 'new_linear_box_spline_2d()'
!    end if
    if (radius .le. 0.)  then
       print *, 'ERROR , new_linear_box_spline_2d : negative radius '
       STOP
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

    ! rk: num_pts-1 = number of nested hexagones
    num_pts_tot = 6*num_pts*(num_pts-1)/2+1 
    SLL_ALLOCATE( new_linear_box_spline_2d%coeffs(0:num_pts_tot), ierr )
  end function new_linear_box_spline_2d


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



  subroutine compute_linear_box_spline_2d( data, spline )
    sll_real64, dimension(:), intent(in), target :: data  ! data to be fit
    type(linear_box_spline_2d), pointer      :: spline
    sll_int32  :: bc1
    sll_int32  :: bc2
    sll_int32  :: bc_selector


    if( .not. associated(spline) ) then
       ! FIXME: THROW ERROR
       print *, 'ERROR: compute_linear_box_spline_2d(): ', &
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
          call compute_box_spline_2d_neum_neum( data, spline )
       case ( 6 )
          ! periodic in x1 and neumann in x2
          call compute_box_spline_2d_prdc_neum( data, spline )
       case default
          print *, 'ERROR: compute_linear_box_spline_2d(): ', &
            'did not recognize given boundary condition combination.'
       STOP
    end select

  
  end subroutine compute_linear_box_spline_2d


  subroutine compute_box_spline_2d_neum_neum( data, spline )
    sll_real64, dimension(:), intent(in), target :: data  ! data to be fit
    type(linear_box_spline_2d), pointer      :: spline
    sll_int32  :: num_pts_tot
    sll_int32  :: i,k
    sll_int32  :: nei
    sll_int32  :: k1
    sll_int32  :: k2
    sll_int32  :: temp
    sll_int32  :: hex_num
    sll_int32  :: num_pts

    num_pts = spline%num_pts
    num_pts_tot = 6*num_pts*(num_pts-1)/2+1 
    do i = 0, num_pts_tot-1
       spline%coeffs(i+1) = real(0,f64)
       do k = 0, 18
          nei = find_neighbour(i,k)
          if (nei .lt. num_pts_tot) then
             spline%coeffs(i+1) = spline%coeffs(i+1) + data(nei+1) * & 
                                pre_filter_piir2(k)
          else
             k1 = from_global_index_k1(nei)
             k2 = from_global_index_k2(nei)
             ! We compute the hexagon-ring number
             if (k1*k2 .ge. 0.) then
                hex_num = max( abs(k1), abs(k2))
             else 
                hex_num = abs(k1) + abs(k2)
             end if

             if (k1 .eq. hex_num)  then
                ! index out on first edge
                k1 = 2*(num_pts - 1) - k1
                k2 = k2 + k1 - num_pts + 1
             elseif (k2 .eq. hex_num) then
                ! index out on second edge
                k2 = 2*(num_pts - 1) - k2
                k1 = k2 + k1 - num_pts + 1
             elseif (( k1 .lt. 0) .and. (k2 .gt. 0)) then
                ! index out on third edge
                temp = k2 - num_pts + 1
                k2   = k1 + num_pts - 1
                k1   = temp
             elseif (k1 .eq. -hex_num) then
                ! index out on fourth edge
                k1 = -2*(num_pts - 1) - k1
                k2 =  k2 + k1 + num_pts - 1
             elseif (k2 .eq. -hex_num) then
                ! index out on fifth edge
                k2 = -2*(num_pts - 1) - k2
                k1 = k2 + k1 + num_pts - 1
             elseif ((k1 .gt. 0).and.(k2 .lt. 0)) then
                ! index out on sixth edge
                temp = k2 + num_pts - 1
                k2   = k1 - num_pts + 1
                k1   = temp
             end if
             nei = global_index(k1,k2)
             spline%coeffs(i+1) = spline%coeffs(i+1) + data(nei+1) * & 
                                pre_filter_piir2(k)             
          end if
       end do
    end do

  end subroutine compute_box_spline_2d_neum_neum


  subroutine compute_box_spline_2d_prdc_neum( data, spline )
    sll_real64, dimension(:), intent(in), target :: data  ! data to be fit
    type(linear_box_spline_2d), pointer      :: spline
    sll_int32  :: num_pts_tot
    sll_int32  :: i
    sll_int32  :: num_pts

    print *, ' WARNING : BOUNDARY CONDITIONS PERIODIC_NEUMANN NOT &
      & YET IMPLEMENTED'
    num_pts = spline%num_pts
    num_pts_tot = 6*num_pts*(num_pts-1)/2+1 
    do i = 0, num_pts_tot
       spline%coeffs(i) = real(0,f64)*data(i)
    end do

  end subroutine compute_box_spline_2d_prdc_neum




  function chi2(spline, x1, x2)
    ! Reference : @Condat and Van De Ville (2006)
    !             "Three directional Box Splines:
    !             Characterization and Efficient Evaluation."
    sll_real64                          :: chi2
    sll_real64, intent(in)              :: x1
    sll_real64, intent(in)              :: x2
    type(linear_box_spline_2d), pointer :: spline
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
    type(hex_mesh_2d), pointer  :: mesh
    sll_real64                          :: val
    intrinsic                           :: associated, int, real
    sll_real64, intent(in)              :: x1
    sll_real64, intent(in)              :: x2
    type(linear_box_spline_2d), pointer :: spline
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

