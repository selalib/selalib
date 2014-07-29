!> @brief  
!> Box-Splines of degree n using formula from : 
!> Condat2006 : Three-directional box splines
!> @author : Laura Mendoza (mela@ipp.mpg.de)

module box_splines
#include "sll_working_precision.h"
#include "sll_memory.h"
#include "sll_assert.h"
#include "sll_splines.h"
#include "sll_utilities.h"
use hex_pre_filters
use hex_mesh


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
   sll_int32 SLL_PRIV  :: x1_bc_type
   sll_int32 SLL_PRIV  :: x2_bc_type
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

    new_linear_box_spline_2d%r1_x1 = real(0.5,f64)
    new_linear_box_spline_2d%r1_x2 = -sqrt(3.)/2.
    new_linear_box_spline_2d%r2_x1 = real(0.5,f64)
    new_linear_box_spline_2d%r2_x2 = sqrt(3.0)/2.
    new_linear_box_spline_2d%r3_x1 = 1.0_f64
    new_linear_box_spline_2d%r3_x2 = 1.0_f64

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


  subroutine compute_linear_box_spline_2d( data, deg, spline )
    sll_real64, dimension(:), intent(in), target    :: data  ! data to be fit
    sll_int32, intent(in)                           :: deg
    type(linear_box_spline_2d), pointer, intent(in) :: spline
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
          call compute_box_spline_2d_neum_neum( data, deg, spline )
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



  subroutine compute_box_spline_2d_neum_neum( data, deg, spline )
    sll_real64, dimension(:), intent(in), target  :: data  ! data to be fit
    type(linear_box_spline_2d), pointer           :: spline
    sll_int32, intent(in)                         :: deg
    sll_int32  :: num_pts_tot
    sll_int32  :: i,k
    sll_int32  :: nei
    sll_int32  :: num_pts

    num_pts = spline%num_pts
    num_pts_tot = 6*num_pts*(num_pts-1)/2+1 
    

    do i = 0, num_pts_tot-1
       spline%coeffs(i) = real(0,f64)
       do k = 0, 3*(2*deg)*(2*deg+1)
          nei = find_neighbour(i,k)
          if (nei .lt. num_pts_tot) then
             spline%coeffs(i) = spline%coeffs(i) + data(nei+1) * & 
                                pre_filter_pfir(k,deg)
             ! if (i .eq. 35) then
             !    print *,""
             !    print *,"nei  =", nei
             !    print *,"data = ", data(nei+1)
             !    print *,"filt = ", pre_filter_piir1(k,deg)
             !    print *,"coef = ", spline%coeffs(i)
             ! end if
          else
             spline%coeffs(i) = spline%coeffs(i)
             ! nei = twelve_fold_symmetry(nei, num_pts)
             ! spline%coeffs(i) = spline%coeffs(i) + data(nei+1) * & 
             !                    pre_filter_pfir(k,deg)
          end if
       end do
    end do
!   print*, " ++++++++++ WARNING :  coeffs = data +++++++++"
!   spline%coeffs(:) = data(:)
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


  function Factorial (n) result (res)
    sll_int32, intent (in) :: n
    sll_int32 :: res
    sll_int32 :: i
    sll_int32 :: tab(n)

    tab = (/(i, i = 1, n)/)
    res = product (tab)
 
  end function factorial
 
  function choose (n, k) result (res)

    sll_int32, intent (in) :: n
    sll_int32, intent (in) :: k
    sll_real64 :: res
    if (k.lt.0) then
       res = 0.
    else if (n .lt. k) then
       res = 0.
    else 
       res = Factorial(n) / (Factorial (k) * Factorial (n - k))
    end if
  end function choose


  function chi_gen_val(x1_in,x2_in,deg) result(val)
    ! Reference : @Condat and Van De Ville (2006)
    !             "Three directional Box Splines:
    !             Characterization and Efficient Evaluation."
    sll_real64, intent(in)               :: x1_in
    sll_real64, intent(in)               :: x2_in
    sll_int32,  intent(in)               :: deg
    sll_real64    :: x1
    sll_real64    :: x2
    sll_real64    :: val
    sll_real64    :: u
    sll_real64    :: v
    sll_real64    :: aux
    sll_real64    :: aux2
    sll_real64    :: coeff
    sll_real64    :: g
    sll_int32     :: K
    sll_int32     :: L
    sll_int32     :: i
    sll_int32     :: d



    if (deg.ne.2) then
       ! ------------------------------------------------
       ! FOR ARBITRARY ORDER SPLINES :
       ! We take advantage of the fact that the mesh is symetrical
       ! So we constrain the points to the first triangles by
       ! using the symmetry properties
       ! ------------------------------------------------
       ! Transformation to hexagonal coordinates :
       x1 = -abs(x1_in)
       x2 =  abs(x2_in)
       u = x1 - x2/sqrt(3._f64)
       v = x1 + x2/sqrt(3._f64)
       ! First generating vector (r1') symmetry
       if (v.gt.0) then
          v = -v
          u = u + v
       end if
       ! First triangle axe (r2'+r3') symmetry
       if(v.gt.u/2) then
          v = u-v
       end if
       
       val = 0._f64
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
                        * aux**(deg-1-d) &
                        * aux2**(2*deg-1+d)
                end do
             end do
          end do
       end do

    else
       ! ------------------------------------------------
       ! OPTIMIZED ALGORITHM FOR DEGREE TWO SPLINES :
       ! We take advantage of the fact that the mesh is symetrical
       ! So we constrain the points to the first triangles by
       ! using the symmetry properties
       ! ------------------------------------------------
       ! Transformation to hexagonal coordinates :
       x1 = abs(x1_in)
       x2 = abs(x2_in)
       u = x1 - x2/sqrt(3._f64)
       v = x1 + x2/sqrt(3._f64)
       if (u.lt.0) then
          !Symmetry r2
          u = -u
          v = v+u
       end if
       if (2*u.lt.v) then
          !Symmetry r2+r3
          u = v-u
       end if
       g = u - v/2._f64
       if (v.gt.2.0) then
          val = 0.0_f64
       else if (v.lt.1.0) then
          val = 0.5 + ((5._f64/3._f64 - v/8.0_f64)*v-3.0_f64)*v*v/4.0_f64 + &
               ((1._f64 - v/4._f64)*v + g*g/6._f64 - 1._f64)*g*g
       else if (u.gt.1.0) then
          val = (v - 2._f64)*(v - 2._f64)*(v - 2._f64)*(g - 1._f64)/6.0_f64
       else
          val = 5._f64/6._f64 + ((1._f64 + (1._f64/3._f64 - v/8._f64)*v)*v/4._f64-1._f64)*v + &
               ((1._f64 - v/4._f64)*v + g*g/6.0_f64 - 1._f64)*g*g
       end if
    end if

  end function chi_gen_val


  function change_basis_x1(mesh, spline, x1, x2) result(x1_basis)
    type(hex_mesh_2d), pointer, intent(in)  :: mesh
    sll_real64, intent(in)                  :: x1
    sll_real64, intent(in)                  :: x2
    type(linear_box_spline_2d), pointer     :: spline
    sll_real64                              :: delta_q
    sll_real64                              :: k1_basis
    sll_real64                              :: k2_basis
    sll_real64                              :: x1_basis
    sll_real64              :: q11, q12
    sll_real64              :: q21, q22
    sll_real64              :: q31, q32
    sll_real64              :: r11, r12
    sll_real64              :: r21, r22
    sll_real64              :: r31, r32


    r11 = spline%r1_x1
    r12 = spline%r1_x2
    r21 = spline%r2_x1
    r22 = spline%r2_x2
    r31 = spline%r3_x1
    r32 = spline%r3_x2

    q11 = mesh%r1_x1
    q12 = mesh%r1_x2
    q21 = mesh%r2_x1
    q22 = mesh%r2_x2
    q31 = mesh%r3_x1
    q32 = mesh%r3_x2

    !change of basis :
    delta_q  = q11*q22 - q12*q21
    k1_basis = 1./delta_q*(q22*x1 - q21*x2)
    k2_basis = 1./delta_q*(q11*x2 - q12*x1)
    x1_basis = r11*k1_basis+r21*k2_basis
    x1_basis = x1_basis
  end function change_basis_x1


  function change_basis_x2(mesh, spline, x1, x2) result(x2_basis)
    type(hex_mesh_2d), pointer, intent(in)  :: mesh
    sll_real64, intent(in)                  :: x1
    sll_real64, intent(in)                  :: x2
    type(linear_box_spline_2d), pointer     :: spline
    sll_real64                              :: delta_q
    sll_real64                              :: k1_basis
    sll_real64                              :: k2_basis
    sll_real64                              :: x2_basis
    sll_real64              :: q11, q12
    sll_real64              :: q21, q22
    sll_real64              :: q31, q32
    sll_real64              :: r11, r12
    sll_real64              :: r21, r22
    sll_real64              :: r31, r32

    r11 = spline%r1_x1
    r12 = spline%r1_x2
    r21 = spline%r2_x1
    r22 = spline%r2_x2
    r31 = spline%r3_x1
    r32 = spline%r3_x2

    q11 = mesh%r1_x1
    q12 = mesh%r1_x2
    q21 = mesh%r2_x1
    q22 = mesh%r2_x2
    q31 = mesh%r3_x1
    q32 = mesh%r3_x2

    !change of basis :
    delta_q  = q11*q22 - q12*q21
    k1_basis = 1./delta_q*(q22*x1 - q21*x2)
    k2_basis = 1./delta_q*(q11*x2 - q12*x1)
    x2_basis = r12*k1_basis+r22*k2_basis
    x2_basis = x2_basis
  end function change_basis_x2

  function hex_interpolate_value(mesh_geom, x1, x2, spline, deg) result(val)
    type(hex_mesh_2d), pointer, intent(in)  :: mesh_geom
    type(linear_box_spline_2d), pointer     :: spline
    sll_int32,  intent(in)                  :: deg
    sll_real64, intent(in)                  :: x1
    sll_real64, intent(in)                  :: x2
    sll_real64              :: val
    sll_real64              :: xm1
    sll_real64              :: xm2
    sll_real64              :: x1_basis
    sll_real64              :: x2_basis
    sll_real64              :: r11, r12
    sll_real64              :: r21, r22
    sll_real64              :: r31, r32
    sll_int32               :: k1_asso
    sll_int32               :: k2_asso
    sll_int32               :: k1
    sll_int32               :: k2
    sll_int32               :: ind
    sll_int32               :: ki, kj
    sll_int32               :: num_pts
    sll_int32               :: num_pts_tot
    sll_real64              :: delta
    sll_real64              :: k1_temp
    sll_real64              :: k2_temp

    val = 0._f64

    r11 = mesh_geom%r1_x1
    r12 = mesh_geom%r1_x2
    r21 = mesh_geom%r2_x1
    r22 = mesh_geom%r2_x2
    r31 = mesh_geom%r3_x1
    r32 = mesh_geom%r3_x2

    num_pts = mesh_geom%num_cells + 1
    num_pts_tot = mesh_geom%num_pts_tot

    ! First we need to compute the coordinates of 
    ! the closest mesh point associated to (x1,x2)
    k1_asso = from_cart_index_k1(mesh_geom, x1, x2)
    k2_asso = from_cart_index_k2(mesh_geom, x1, x2)

    ! Then we will do a loop for all the points 
    ! on the envelopping rhomboid of radius=deg
    do ki= 1-deg, deg
      do kj= 1-deg, deg
        
        k1  = k1_asso + ki
        k2  = k2_asso + kj
        ind = global_index(k1, k2)

        
        do while (ind .ge. num_pts_tot)
          ! If point got out of the domain
          ind = twelve_fold_symmetry(ind, num_pts)
        end do

        ! We centralize and shift the coordinates
        ! i.e. centralize : xm = x - Rk
        !      shifting   : xm to the associated rhomboid point
        xm1 = x1 - r11*k1_asso - r21*k2_asso - ki*r11 - kj*r21
        xm2 = x2 - r12*k1_asso - r22*k2_asso - ki*r12 - kj*r22

        ! change of basis : geometrical basis => spline basis
        x1_basis = change_basis_x1(mesh_geom, spline, xm1, xm2)
        x2_basis = change_basis_x2(mesh_geom, spline, xm1, xm2)

        val = val + spline%coeffs(ind) * &
                    chi_gen_val(x1_basis, x2_basis, deg)

      ! ! if ((x1.eq.0.0).and.(x2.eq.0.0)) then
      !   if ((k1_asso.eq.-10).and.(k2_asso.eq.-2)) then
      !   ! if (ind.eq.245) then
      !      if(chi_gen_val(x1_basis, x2_basis, deg).gt.0.) then
      !         print*," "
      !         print*," indice =", ind
      !         print*," chi_gv =", chi_gen_val(x1_basis, x2_basis, deg)
      !         print*," coeffs =", spline%coeffs(ind)
      !         print*," k1, k2 =", k1_asso, k2_asso
      !         print*," ki, kj =", ki, kj
      !         ! print*," x , y  =", x1, x2
      !         ! print*," shifts =", ki*r11 + kj*r21, ki*r12 + kj*r22
      !         ! print*," modifi =", xm1, xm2
      !         ! print*," basiss =", x1_basis, x2_basis
      !      else
      !         print*,""
      !         print*," indice =", ind
      !         ! print*," chi_gv =", chi_gen_val(x1_basis, x2_basis, deg)
      !         ! print*," coeffs =", spline%coeffs(ind)
      !         ! print*," shifts =", ki*r11 + kj*r21, ki*r12 + kj*r22
      !         ! print*," modifi =", xm1, xm2
      !         ! print*," basiss =", x1_basis, x2_basis

      !      end if
      !   end if

      end do
   end do

  end function hex_interpolate_value

end module box_splines

