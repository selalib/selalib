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
    mesh,         &
    x1_bc_type,   &
    x2_bc_type)

    type(linear_box_spline_2d), pointer  :: new_linear_box_spline_2d
    type(hex_mesh_2d), pointer           :: mesh
    sll_int32,  intent(in)               :: x1_bc_type
    sll_int32,  intent(in)               :: x2_bc_type
    sll_int32                            :: bc_selector
    sll_int32                            :: ierr
    sll_int32                            :: num_pts_tot


    SLL_ALLOCATE( new_linear_box_spline_2d, ierr )
    new_linear_box_spline_2d%num_pts   = mesh%num_cells + 1
    new_linear_box_spline_2d%radius    = mesh%radius
    new_linear_box_spline_2d%delta     = mesh%delta
    new_linear_box_spline_2d%x1_bc_type = x1_bc_type
    new_linear_box_spline_2d%x2_bc_type = x2_bc_type

    new_linear_box_spline_2d%r1_x1 = mesh%r1_x1
    new_linear_box_spline_2d%r1_x2 = mesh%r1_x2
    new_linear_box_spline_2d%r2_x1 = mesh%r2_x1
    new_linear_box_spline_2d%r2_x2 = mesh%r2_x2


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
    num_pts_tot = mesh%num_pts_tot
    SLL_ALLOCATE( new_linear_box_spline_2d%coeffs(0:num_pts_tot-1), ierr )
!    new_linear_box_spline_2d%coeffs(:) = 0._f64
  end function new_linear_box_spline_2d


  ! Pre-filter to compute the box splines coefficients
  ! Reference : @Condat and Van De Ville (2007)
  !             "Quasi-interpolating spline models 
  !             for hexagonally-sampled data."
  function pre_filter_piir2(local_index) result(weight)
      sll_int32, intent(in)     :: local_index
      sll_real64                :: weight

      if (local_index .eq. 0) then
         weight = real((1775./2304.),f64)
      else if (local_index .lt. 7) then
         weight = real((256./6912.), f64)
      else if (local_index .lt. 19) then
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

  function twelve_fold_symmetry(out_index, num_pts) result(in_index)
      sll_int32    :: out_index
      sll_int32    :: num_pts
      sll_int32    :: in_index
      sll_int32    :: k1
      sll_int32    :: k2
      sll_int32    :: temp
      sll_int32    :: hex_num


      k1 = from_global_index_k1(out_index)
      k2 = from_global_index_k2(out_index)
      ! We compute the hexagon-ring number
      if (k1*k2 .ge. 0.) then
         hex_num = max( abs(k1), abs(k2))
      else 
         hex_num = abs(k1) + abs(k2)
      end if

      if (k1 .eq. hex_num)  then
         ! index out on first edge : symmetry by r2_ext
         k1 = 2*(num_pts - 1) - k1
         k2 = k2 + k1 - num_pts + 1
      elseif (k2 .eq. hex_num) then
         ! index out on second edge : symmetry by r1_ext
         k2 = 2*(num_pts - 1) - k2
         k1 = k2 + k1 - num_pts + 1
      elseif (( k1 .lt. 0) .and. (k2 .gt. 0)) then
         ! index out on third edge : symmetry by r3_ext
         temp = k2 - num_pts + 1
         k2   = k1 + num_pts - 1
         k1   = temp
      elseif (k1 .eq. -hex_num) then
         ! index out on fourth edge : symmetry by r2_ext
         k1 = -2*(num_pts - 1) - k1
         k2 =  k2 + k1 + num_pts - 1
      elseif (k2 .eq. -hex_num) then
         ! index out on fifth edge : symmetry by r1_ext
         k2 = -2*(num_pts - 1) - k2
         k1 = k2 + k1 + num_pts - 1
      elseif ((k1 .gt. 0).and.(k2 .lt. 0)) then
         ! index out on sixth edge : symmetry by r3_ext
         temp = k2 + num_pts - 1
         k2   = k1 - num_pts + 1
         k1   = temp
      end if
      in_index = global_index(k1,k2)
  end function twelve_fold_symmetry



  subroutine compute_box_spline_2d_neum_neum( data, spline )
    sll_real64, dimension(:), intent(in), target :: data  ! data to be fit
    type(linear_box_spline_2d), pointer      :: spline
    sll_int32  :: num_pts_tot
    sll_int32  :: i,k
    sll_int32  :: nei
    sll_int32  :: num_pts

    num_pts = spline%num_pts
    num_pts_tot = 6*num_pts*(num_pts-1)/2+1 
    

    do i = 0, num_pts_tot-1
       spline%coeffs(i) = real(0,f64)
       do k = 0, 18
          nei = find_neighbour(i,k)
          if (nei .lt. num_pts_tot) then
             spline%coeffs(i) = spline%coeffs(i) + data(nei+1) * & 
                                pre_filter_piir2(k)
          else
             nei = twelve_fold_symmetry(nei, num_pts)
             spline%coeffs(i) = spline%coeffs(i) + data(nei+1) * & 
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
    do i = 0, num_pts_tot-1
       spline%coeffs(i) = real(0,f64)*data(i)
    end do

  end subroutine compute_box_spline_2d_prdc_neum

  function factorial (n) result (res)
 
    implicit none
    integer, intent (in) :: n
    integer :: res
    integer :: i
 
    res = product ((/(i, i = 1, n)/))
 
  end function factorial
 
  function choose (n, k) result (res)
 
    implicit none
    integer, intent (in) :: n
    integer, intent (in) :: k
    integer :: res
 
    res = factorial (n) / (factorial (k) * factorial (n - k))
 
  end function choose


  function chi_gen(spline,x1,x2,deg) result(val)
    ! Reference : @Condat and Van De Ville (2006)
    !             "Three directional Box Splines:
    !             Characterization and Efficient Evaluation."
    sll_real64, dimension(:), allocatable   :: val
    sll_real64, dimension(:)                :: x1
    sll_real64, dimension(:)                :: x2
    sll_int32, intent(in)                   :: deg
    type(linear_box_spline_2d), pointer     :: spline
    sll_real64, dimension(:), allocatable   :: u
    sll_real64, dimension(:), allocatable   :: v
    sll_real64, dimension(:), allocatable   :: aux
    sll_real64, dimension(:), allocatable   :: aux2
    sll_real64              :: det
    sll_real64              :: coeff
    sll_int32               :: num_pts_tot
    sll_int32               :: K
    sll_int32               :: L
    sll_int32               :: i
    sll_int32               :: d
    sll_int32               :: ierr


    num_pts_tot = 6*spline%num_pts*(spline%num_pts-1)/2+1 
    SLL_ALLOCATE( val(num_pts_tot), ierr )
    SLL_ALLOCATE( u(num_pts_tot), ierr )
    SLL_ALLOCATE( v(num_pts_tot), ierr )
    SLL_ALLOCATE( aux(num_pts_tot), ierr )
    SLL_ALLOCATE( aux2(num_pts_tot), ierr )
    val(:) = 0._f64


    det = spline%r1_x1 * spline%r2_x2 - spline%r2_x1 * spline%r1_x2

    x1 = -abs(x1); x2 = abs(x2);
    u(:) = 1./det*( spline%r2_x2*x1(:) - spline%r2_x1*x2(:))
    v(:) = 1./det*(-spline%r1_x2*x1(:) + spline%r1_x1*x2(:))

    WHERE(v.gt.0) v = -v
    WHERE(v.gt.0) u = u + v
    WHERE(v.gt.u/2) v = u-v


    do K = -deg, CEILING(MAXVAL(u))-1
       do L = -deg, CEILING(MAXVAL(v))-1
            do i = 0,min(deg+K, deg+L)
                coeff = (-1.0_f64)**(K+L+i)*choose(deg,i-K)*choose(deg,i-L)*choose(deg,i)
                do d = 0,deg-1
                    aux=abs(2./sqrt(3.)*(spline%r1_x2*u+spline%r2_x2*v)+K-L)
                    aux2=(2.*(spline%r1_x1*u+spline%r2_x1*v)-K-L-aux)/2
                    WHERE(aux2.lt.0) aux2 = 0
                    val(:) = val(:) + coeff*choose(deg-1+d,d)   &
                        /factorial(2*deg-1+d)/factorial(deg-1-d)  &
                        * aux**(deg-1-d)\
                        * aux2**(2*deg-1+d)
                 end do
            end do
         end do
      end do
      
  end function chi_gen


  function chi_gen_val(spline,x1,x2,deg) result(val)
    ! Reference : @Condat and Van De Ville (2006)
    !             "Three directional Box Splines:
    !             Characterization and Efficient Evaluation."
    sll_real64    :: val
    sll_real64    :: x1
    sll_real64    :: x2
    sll_int32, intent(in)     :: deg
    type(linear_box_spline_2d), pointer     :: spline
    sll_real64    :: u
    sll_real64    :: v
    sll_real64    :: aux
    sll_real64    :: aux2
    sll_real64    :: det
    sll_real64    :: coeff
    sll_int32     :: K
    sll_int32     :: L
    sll_int32     :: i
    sll_int32     :: d


    val = 0._f64
    det = spline%r1_x1 * spline%r2_x2 - spline%r2_x1 * spline%r1_x2

    x1 = -abs(x1); x2 = abs(x2);
    u = 1./det*( spline%r2_x2*x1 - spline%r2_x1*x2)
    v = 1./det*(-spline%r1_x2*x1 + spline%r1_x1*x2)

    if (v.gt.0) then
      v = -v
      u = u + v
    end if
    if(v.gt.u/2) then
      v = u-v
    end if

    do K = -deg, CEILING(u)-1
       do L = -deg, CEILING(v)-1
            do i = 0,min(deg+K, deg+L)
                coeff = (-1.0_f64)**(K+L+i)*choose(deg,i-K)*choose(deg,i-L)*choose(deg,i)
                do d = 0,deg-1
                    aux=abs(2./sqrt(3.)*(spline%r1_x2*u+spline%r2_x2*v)+K-L)
                    aux2=(2.*(spline%r1_x1*u+spline%r2_x1*v)-K-L-aux)/2
                    if(aux2.lt.0) then
                       aux2 = 0
                    end if
                    val = val + coeff*choose(deg-1+d,d)   &
                        /factorial(2*deg-1+d)/factorial(deg-1-d)  &
                        * aux**(deg-1-d)\
                        * aux2**(2*deg-1+d)
                 end do
            end do
         end do
      end do
      
  end function chi_gen_val


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
    sll_int32               :: num_pts
    sll_int32               :: num_pts_tot

    r1_x1 = spline%r1_x1
    r1_x2 = spline%r1_x2
    r2_x1 = spline%r2_x1
    r2_x2 = spline%r2_x2

    num_pts = mesh%num_cells + 1
    num_pts_tot = mesh%num_pts_tot

    !Lower left corner of encapsulating rhomboid
    k1 = from_cart_index_k1(mesh, x1, x2)
    k2 = from_cart_index_k2(mesh, x1, x2)
    i  = global_index(k1, k2) 
    if (i .ge. num_pts_tot)  then 
       i = twelve_fold_symmetry(i, num_pts)
    end if
    
    xm1 = x1-r1_x1*k1-r2_x1*k2
    xm2 = x2-r1_x2*k1-r2_x2*k2
    val = spline%coeffs(i)*chi_gen_val(spline,xm1, xm2,1)

    !Lower right corner of encapsulating rhomboid
    k1 = from_cart_index_k1(mesh, x1, x2) + 1
    k2 = from_cart_index_k2(mesh, x1, x2)
    i  = global_index(k1, k2) 
    if (i .ge. num_pts_tot)  then 
       i = twelve_fold_symmetry(i, num_pts)
    end if

    xm1 = x1-r1_x1*k1-r2_x1*k2 - r1_x1
    xm2 = x2-r1_x2*k1-r2_x2*k2 - r1_x2
    val = val + spline%coeffs(i)*chi2(spline,xm1, xm2)
    
    !Upper left corner of encapsulating rhomboid
    k1 = from_cart_index_k1(mesh, x1, x2) 
    k2 = from_cart_index_k2(mesh, x1, x2) + 1
    i  = global_index(k1, k2) 
    if (i .ge. num_pts_tot)  then 
       i = twelve_fold_symmetry(i, num_pts)
    end if

    xm1 = x1-r1_x1*k1-r2_x1*k2 - r2_x1
    xm2 = x2-r1_x2*k1-r2_x2*k2 - r2_x2
    val = val + spline%coeffs(i)*chi2(spline,xm1, xm2)

    !Upper right corner of encapsulating rhomboid
    k1 = from_cart_index_k1(mesh, x1, x2) + 1
    k2 = from_cart_index_k2(mesh, x1, x2) + 1
    i  = global_index(k1, k2) 
    if (i .ge. num_pts_tot)  then 
       i = twelve_fold_symmetry(i, num_pts)
    end if

    xm1 = x1-r1_x1*k1-r2_x1*k2 - r1_x1 - r2_x1
    xm2 = x2-r1_x2*k1-r2_x2*k2 - r1_x2 - r2_x2
    val = val + spline%coeffs(i)*chi2(spline,xm1, xm2)

  end function hex_interpolate_value

end module box_splines

