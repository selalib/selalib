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
    new_linear_box_spline_2d%r3_x1 = mesh%r3_x1
    new_linear_box_spline_2d%r3_x2 = mesh%r3_x2

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
         if (modulo(local_index, 2) .eq. 0) then
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
      sll_int32,intent(in)    :: out_index
      sll_int32,intent(in)    :: num_pts
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
    type(linear_box_spline_2d), pointer          :: spline
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
!             if ((i .eq. 0).and.(k.eq.0)) then
!                print*,"data =",data(nei+1)
!                print*,"res  =",spline%coeffs(i)
!                print*,""
!             end if
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
    sll_int32, intent (in) :: n
    sll_int32 :: res
    sll_int32 :: i
 
    res = product ((/(i, i = 1, n)/))
 
  end function factorial
 
  function choose (n, k) result (res)
 
    implicit none
    sll_int32, intent (in) :: n
    sll_int32, intent (in) :: k
    sll_real64 :: res
    if (k.lt.0) then
       res = 0.
    else if (n .lt. k) then
       res = 0.
    else 
       res = factorial (n) / (factorial (k) * factorial (n - k))
    end if
  end function choose


  function chi_gen_val(spline,x1_in,x2_in,deg) result(val)
    ! Reference : @Condat and Van De Ville (2006)
    !             "Three directional Box Splines:
    !             Characterization and Efficient Evaluation."
    sll_real64, intent(in)               :: x1_in
    sll_real64, intent(in)               :: x2_in
    sll_int32,  intent(in)               :: deg
    type(linear_box_spline_2d), pointer  :: spline
    sll_real64    :: x1
    sll_real64    :: x2
    sll_real64    :: val
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

    x1 = -abs(x1_in)
    x2 =  abs(x2_in)
    u = 1./det*( spline%r2_x2*x1 - spline%r2_x1*x2)/(real(deg+1,f64))
      !*spline%radius/real(spline%num_pts,f64)
    v = 1./det*(-spline%r1_x2*x1 + spline%r1_x1*x2)/(real(deg+1,f64))
      !*spline%radius/real(spline%num_pts,f64)

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
                coeff = (-1.0_f64)**(K+L+i)* &
                        choose(deg,i-K)*     &
                        choose(deg,i-L)*     &
                        choose(deg,i)
                do d = 0,deg-1
                    aux=abs(v-L-u+K)  
                    aux2=(u-K+v-L-aux)/2._f64
                    if(aux2.lt.0.) then
                       aux2 = 0._f64
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



  function hex_interpolate_value(mesh, x1, x2, spline) result(val)
    type(hex_mesh_2d), pointer, intent(in)  :: mesh
    sll_real64, intent(in)                  :: x1
    sll_real64, intent(in)                  :: x2
    type(linear_box_spline_2d), pointer     :: spline
    sll_real64                          :: val
    intrinsic                           :: associated, int, real
    sll_real64              :: xm1
    sll_real64              :: xm2
    sll_real64              :: x1_center
    sll_real64              :: x2_center
    sll_real64              :: r1_x1, r1_x2
    sll_real64              :: r2_x1, r2_x2
    sll_real64              :: r3_x1, r3_x2
    sll_int32               :: k1
    sll_int32               :: k2
    sll_int32               :: i
    sll_int32               :: deg = 1
    sll_int32               :: num_pts
    sll_int32               :: num_pts_tot
!    sll_int32               :: dr1,dr2

    r1_x1 = spline%r1_x1
    r1_x2 = spline%r1_x2
    r2_x1 = spline%r2_x1
    r2_x2 = spline%r2_x2
    r3_x1 = spline%r3_x1
    r3_x2 = spline%r3_x2

    num_pts = mesh%num_cells + 1
    num_pts_tot = mesh%num_pts_tot

    !Lower left corner of encapsulating rhomboid
    k1 = from_cart_index_k1(mesh, x1, x2)
    k2 = from_cart_index_k2(mesh, x1, x2)
    i  = global_index(k1, k2)        

    if (i .ge. num_pts_tot)  then 
       i = twelve_fold_symmetry(i, num_pts)
    end if

    x1_center = x1-r1_x1*k1-r2_x1*k2
    x2_center = x2-r1_x2*k1-r2_x2*k2
    val = spline%coeffs(i)*chi_gen_val(spline, x1_center, x2_center, deg)*spline%delta

    if ((x1.eq.0.).and.(x2.eq.0.)) then
      print*, "delta =", spline%delta
       print *, "i0 = ", i, " coeffs =",spline%coeffs(i), "chi =",chi_gen_val(spline,x1_center, x2_center, deg)
    end if


!    if ((x1.eq.0.).and.(x2.eq.0)) then
!       print *, "val   =", val
!       print *, "coeff =", spline%coeffs(i)
!       print *, "chi   =", chi_gen_val(spline,xm1, xm2, deg)
!    end if
!
    !Lower right corner of encapsulating rhomboid
    k1 = from_cart_index_k1(mesh, x1, x2) + 1
    k2 = from_cart_index_k2(mesh, x1, x2)
    i  = global_index(k1, k2) 
    if (i .ge. num_pts_tot)  then 
       i = twelve_fold_symmetry(i, num_pts)
    end if

    xm1 = x1_center + r1_x1
    xm2 = x2_center + r1_x2
    val = val + spline%coeffs(i)*chi_gen_val(spline,xm1, xm2, deg)
    
    if ((x1.eq.0.).and.(x2.eq.0.)) then
       print *, "i1 = ", i, " coeffs =",spline%coeffs(i), "chi =",chi_gen_val(spline,xm1, xm2, deg)
    end if


    !Upper left corner of encapsulating rhomboid
    k1 = from_cart_index_k1(mesh, x1, x2)+ 1
    k2 = from_cart_index_k2(mesh, x1, x2)+ 1
    i  = global_index(k1, k2) 
    if (i .ge. num_pts_tot)  then 
       i = twelve_fold_symmetry(i, num_pts)
    end if


    xm1 = x1_center + r3_x1
    xm2 = x2_center + r3_x2
    val = val + spline%coeffs(i)*chi_gen_val(spline,xm1, xm2, deg)


    if ((x1.eq.0.).and.(x2.eq.0.)) then
       print *, "i2 = ", i, " coeffs =",spline%coeffs(i), "chi =",chi_gen_val(spline,xm1, xm2, deg)
    end if

    !Upper right corner of encapsulating rhomboid
    k1 = from_cart_index_k1(mesh, x1, x2) + 2
    k2 = from_cart_index_k2(mesh, x1, x2) + 1
    i  = global_index(k1, k2) 
    if (i .ge. num_pts_tot)  then 
       i = twelve_fold_symmetry(i, num_pts)
    end if


    xm1 = x1_center + r1_x1 + r3_x1
    xm2 = x2_center + r1_x2 + r3_x2
    val = val + spline%coeffs(i)*chi_gen_val(spline,xm1, xm2, deg)



    if ((x1.eq.0.).and.(x2.eq.0.)) then
       print *, "i3 = ", i, " coeffs =",spline%coeffs(i), "chi =",chi_gen_val(spline,xm1, xm2, deg)
    end if

!    do dr1 = -deg+1,deg
!       do dr2 = -deg+1,deg
!          if ((x1.eq.0.).and.(x2.eq.0)) then
!             print *, "dr1, dr2 =", dr1, dr2
!
!             ! TODO : a changer pour que cela n'arrive qu'à deg cellules plus loin
!          end if
!       end do
!    end do


  end function hex_interpolate_value

end module box_splines

