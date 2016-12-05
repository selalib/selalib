!> @ingroup splines
!> @author Laura S. Mendoza
!> @brief Provides capabilities for values and derivatives
!> interpolation with box splines on a hexagonal mesh
!> @details This modules contains the computation of boxsplines
!> of arbitrary degree. There is a special optimized algorithm
!> for degree 2. The boxsplines here are defined only on hexagonal
!> meshes (regular equilateral triangles elements).
!> The only boundary condition supported right now
!> is dirichlet.
!>  Reference :
!>     \@Condat2006 "Three-directional box splines"

module sll_m_box_splines
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#include "sll_memory.h"
#include "sll_working_precision.h"

  use sll_m_boundary_condition_descriptors, only: &
    sll_p_dirichlet, &
    sll_p_neumann, &
    sll_p_periodic

  use sll_m_constants, only: &
    sll_p_epsilon_0, &
    sll_p_sqrt3

  use sll_m_hex_pre_filters, only: &
       sll_s_pre_filter_pfir, &
       sll_f_pre_filter_int, &
       sll_f_pre_filter_piir2, &
       sll_f_pre_filter_piir1

  use sll_m_hexagonal_meshes, only: &
    sll_f_cart_to_hex1, &
    sll_f_cart_to_hex2, &
    sll_f_cells_to_origin, &
    sll_f_change_elements_notation, &
    sll_o_delete, &
    sll_s_get_cell_vertices_index, &
    sll_f_local_to_global, &
    sll_f_new_hex_mesh_2d, &
    sll_t_hex_mesh_2d, &
    sll_s_write_caid_files

  use sll_m_utilities, only: &
    sll_o_factorial

  implicit none

  public :: &
    sll_f_boxspline_x1_derivative, &
    sll_f_boxspline_x2_derivative, &
    sll_f_compute_box_spline, &
    sll_f_hex_interpolate_value, &
    sll_f_new_box_spline_2d, &
    sll_f_boxspline_val_der, &
    sll_t_box_spline_2d, &
    sll_o_delete, &
    sll_s_write_connectivity

  private
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  !> @brief
  !> basic type for 2 dimensional box splines dara
  !> @details
  !> 2D Box spline type, containing the mesh information, the
  !> boundary condition, and the spline coefficients
  type sll_t_box_spline_2d
     type(sll_t_hex_mesh_2d), pointer  :: mesh !< Hexagonal mesh
     sll_int32, private :: bc_type !< Boundary conditions definition
     sll_real64, dimension(:), pointer :: coeffs !< Spline coefficients

   contains
     procedure, pass(spline) :: compute_coeff_box_spline_2d
     procedure, pass(spline) :: compute_coeff_box_spline_2d_diri
     procedure, pass(spline) :: compute_coeff_box_spline_2d_prdc
     procedure, pass(spline) :: compute_coeff_box_spline_2d_neum
  end type sll_t_box_spline_2d

  !> @brief Generic sub-routine defined for 2D box spline types.
  !> Deallocates the memory associated with the given box spline object.
  !> @param[inout] spline_object.
  interface sll_o_delete
     module procedure delete_box_spline_2d
  end interface sll_o_delete

contains  ! ****************************************************************


#define MAKE_GET_SLOT_FUNCTION( fname, datatype, slot, ret_type )    \
  function fname( spline_obj ) result(val);                \
    type(datatype), pointer :: spline_obj;                 \
    ret_type :: val;                                       \
    val = spline_obj%slot;                                 \
  end function fname


  !---------------------------------------------------------------------------
  !> @brief Creates a box spline element
  !> @details Creates a box spline element to the associated mesh and with
  !> the given boundary conditions. Also allocates memory for the coefficients
  !> @param[in] mesh Hexagonal mesh on which the box spline will be defined
  !> @param[in] bc_type Integer defining the boundary condition, usual values
  !> are sll_p_periodic, sll_p_dirichlet, ...
  !> @return box spline element
  function sll_f_new_box_spline_2d( &
       mesh,         &
       bc_type)

    type(sll_t_box_spline_2d), pointer  :: sll_f_new_box_spline_2d
    type(sll_t_hex_mesh_2d),   pointer  :: mesh
    sll_int32,  intent(in)        :: bc_type
    sll_int32                     :: ierr


    SLL_ALLOCATE( sll_f_new_box_spline_2d, ierr )

    sll_f_new_box_spline_2d%mesh => mesh

    sll_f_new_box_spline_2d%bc_type = bc_type

    SLL_ALLOCATE( sll_f_new_box_spline_2d%coeffs(1:mesh%num_pts_tot), ierr )

  end function sll_f_new_box_spline_2d

  !---------------------------------------------------------------------------
  !> @brief Computes box splines coefficients
  !> @details Computes box splines coefficients for a box spline of degree deg
  !> and fitted to the data vector
  !> @param[in] data vector containing the data to be fit
  !> @param[in] deg integer representing the box spline degree
  !> @param[in] spline box spline type element, containting the mesh, bc, ...
  subroutine compute_coeff_box_spline_2d( spline, data, deg )
    class(sll_t_box_spline_2d),         intent(inout) :: spline
    sll_real64, dimension(:), target, intent(in) :: data
    sll_int32,                        intent(in) :: deg

    ! We make every case explicit to facilitate adding more BC types in
    ! the future.
    select case(spline%bc_type)
    case(sll_p_dirichlet)
       ! boundary condition type is dirichlet
       call spline%compute_coeff_box_spline_2d_diri( data, deg)
    case(sll_p_periodic)
       ! boundary condition type is periodic
       call spline%compute_coeff_box_spline_2d_prdc( data, deg)
    case(sll_p_neumann)
       ! boundary condition type is neumann
       call spline%compute_coeff_box_spline_2d_neum( data, deg)
    case default
       print *, 'ERROR: compute_coeff_box_spline_2d(): ', &
            'did not recognize given boundary condition combination.'
       STOP
    end select

  end subroutine compute_coeff_box_spline_2d


  !---------------------------------------------------------------------------
  !> @brief Computes box splines coefficients for dirichlet BC
  !> @details Computes box splines coefficients for a box spline of degree deg
  !> and fitted to the data vector on a hexmesh with dirichlet BC
  !> @param[in] data vector containing the data to be fit
  !> @param[in] deg integer representing the box spline degree
  !> @param[in] spline box spline type element, containting the mesh, bc, ...
  subroutine compute_coeff_box_spline_2d_diri(spline, data, deg)
    class(sll_t_box_spline_2d), intent(inout)          :: spline
    sll_real64, dimension(:), intent(in), target  :: data  ! data to be fit
    sll_int32,                intent(in)          :: deg

    sll_int32  :: num_pts_tot
    sll_int32  :: k1_ref, k2_ref
    sll_int32  :: k
    sll_int32  :: i
    sll_int32  :: ierr
    sll_int32  :: nei
    sll_int32  :: num_pts_radius
    sll_real64 :: filter
    sll_real64, allocatable, dimension(:) :: filter_array

    num_pts_tot = spline%mesh%num_pts_tot
    ! we will work on a radius of 'deg' cells
    ! we compute the number of total points on that radius
    num_pts_radius = 3*deg*(deg+1) + 1

    ! Create a table for the filter values and fill it:

    ! TODO: this should be controlled some other way........ sorry.
    call sll_s_pre_filter_pfir(spline%mesh, deg, filter_array)
    ! If pINT, pIIR1 or  pIIR2:
    ! SLL_ALLOCATE(filter_array(num_pts_radius), ierr)
    ! do k=1, num_pts_radius
       ! filter_array(k) = sll_f_pre_filter_int(spline%mesh, k, deg)
       ! filter_array(k) = sll_f_pre_filter_piir1(spline%mesh, k, deg)
    !    filter_array(k) = sll_f_pre_filter_piir2(spline%mesh, k, deg)
    ! end do

    do i = 1, num_pts_tot

       spline%coeffs(i) = real(0,f64)
       k1_ref = spline%mesh%global_to_hex1(i)
       k2_ref = spline%mesh%global_to_hex2(i)

       ! We don't need to fo through all points, just till a certain radius
       ! which depends on the degree of the spline we are evaluating
       do k = 1, num_pts_radius
          filter = filter_array(k)
          nei = spline%mesh%local_hex_to_global(k1_ref, k2_ref, k)
          if ((nei .le. num_pts_tot).and.(nei .gt. 0)) then
             spline%coeffs(i) = spline%coeffs(i) + data(nei) * filter
          else
             ! Boundary conditions (BC) to be treated here :
             ! With dirichlet boundary conditions data(out_of_domain) = 0
             spline%coeffs(i) = spline%coeffs(i)
          end if
       end do
    end do

  end subroutine compute_coeff_box_spline_2d_diri


  !---------------------------------------------------------------------------
  !> @brief Computes box splines coefficients with periodic BC.
  !> @details NOT YET IMPLEMENTED. Computes periodic box splines coefficients
  !> for a box spline of degree deg and fitted to the data vector
  !> @param[in] data vector containing the data to be fit
  !> @param[in] deg integer representing the box spline degree
  !> @param[in] spline box spline type element, containting the mesh, bc, ...
  subroutine compute_coeff_box_spline_2d_prdc(spline, data, deg)
    class(sll_t_box_spline_2d), intent(inout)          :: spline
    sll_int32,                intent(in)          :: deg
    sll_real64, dimension(:), intent(in), target  :: data  ! data to be fit

    sll_int32  :: num_pts_tot
    sll_int32  :: i

    print *, ' WARNING : BOUNDARY CONDITIONS PERIODIC NOT &
         & YET IMPLEMENTED', deg
    num_pts_tot = spline%mesh%num_pts_tot
    do i = 1, num_pts_tot
       spline%coeffs(i) = real(0,f64)*data(i)
    end do

  end subroutine compute_coeff_box_spline_2d_prdc

  !---------------------------------------------------------------------------
  !> @brief Computes box splines coefficients with neumann BC.
  !> @details NOT YET IMPLEMENTED. Computes neuman box splines coefficients for
  !> a box spline of degree deg and fitted to the data vector
  !> @param[in] data vector containing the data to be fit
  !> @param[in] deg integer representing the box spline degree
  !> @param[in] spline box spline type element, containting the mesh, bc, ...
  subroutine compute_coeff_box_spline_2d_neum(spline, data, deg)
    class(sll_t_box_spline_2d), intent(inout)          :: spline
    sll_real64, dimension(:), intent(in), target  :: data  ! data to be fit
    sll_int32,                intent(in)          :: deg

    sll_int32  :: num_pts_tot
    sll_int32  :: i

    print *, ' WARNING : BOUNDARY CONDITIONS PERIODIC NOT &
         & YET IMPLEMENTED', deg
    num_pts_tot = spline%mesh%num_pts_tot
    do i = 1, num_pts_tot
       spline%coeffs(i) = real(0,f64)*data(i)
    end do

  end subroutine compute_coeff_box_spline_2d_neum

  !---------------------------------------------------------------------------
  !> @brief Computes the binomial coefficient (n, k)
  !> @details Computes the binomial coefficient (n, k)
  !> @param[in] n integer
  !> @param[in] k integer
  !> @result res binomial coefficient \f[ \dfrac{n!}{(k!(n-k)!} \f]
  function choose (n, k) result (res)
    sll_int32, intent (in) :: n
    sll_int32, intent (in) :: k
    sll_real64 :: res
    if (k.lt.0) then
       res = 0._f64
    else if (n .lt. k) then
       res = 0._f64
    else
       res = real(sll_o_factorial(n),f64) / real((sll_o_factorial(k) &
            * sll_o_factorial(n - k)), f64)
    end if
  end function choose


  !---------------------------------------------------------------------------
  !> @brief Computes the box spline of degree deg on (x1, x2)
  !> @details Computes the value of the box spline (chi) of degree deg
  !> on the point of cartesian coordiantes (x1, x2). The algorithm is specific
  !> to the hexagonal mesh and has been optimized for degree 2 splines.
  !> Reference : \@Condat and Van De Ville (2006)
  !>             "Three directional Box Splines:
  !>             Characterization and Efficient Evaluation."
  !> @param[in] x1_in real containing first coordinate of point
  !> where spline is to be evaluated
  !> @param[in] x2_in real containing second coordinate of point
  !> where spline is to be evaluated
  !> @param[in] deg integer containing the degree of the spline
  !> @return chi_gen_val  value of the box spline
  function chi_gen_val(x1_in,x2_in,deg) result(val)
    sll_real64, intent(in) :: x1_in
    sll_real64, intent(in) :: x2_in
    sll_int32,  intent(in) :: deg
    sll_real64 :: x1
    sll_real64 :: x2
    sll_real64 :: val
    sll_real64 :: u
    sll_real64 :: v
    sll_real64 :: aux
    sll_real64 :: aux2
    sll_real64 :: coeff
    sll_real64 :: g
    sll_int32  :: K
    sll_int32  :: L
    sll_int32  :: i
    sll_int32  :: d

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
       u = x1 - x2/sll_p_sqrt3
       v = x1 + x2/sll_p_sqrt3
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
          if ((x1_in.eq.0.).and.(x2_in.eq.0.8)) then
             !             print *, " K = ", K
          end if
          do L = -deg, CEILING(v)-1
             if ((x1_in.eq.0.).and.(x2_in.eq.0.8)) then
                !                print *, "    L = ", L
             end if
             do i = 0,min(deg+K, deg+L)
                if ((x1_in.eq.0.).and.(x2_in.eq.0.8)) then
                   !                   print *, "      i = ", i
                end if
                coeff = (-1.0_f64)**(K+L+i)* &
                     choose(deg,i-K)*     &
                     choose(deg,i-L)*     &
                     choose(deg,i)
                if ((x1_in.eq.0.).and.(x2_in.eq.0.8)) then
                   !                  print *, "      coeff = ", coeff
                end if
                do d = 0,deg-1
                   if ((x1_in.eq.0.).and.(x2_in.eq.0.8)) then
                      !                    print *, "          d = ", d
                   end if
                   aux=abs(v-L-u+K)
                   aux2=(u-K+v-L-aux)/2._f64
                   if(aux2.lt.0.) then
                      aux2 = 0._f64
                   end if
                   val = val + coeff*choose(deg-1+d,d)   &
                        /real(sll_o_factorial(2*deg-1+d), f64) &
                        /real(sll_o_factorial(deg -1 -d), f64) &
                        * aux**(deg-1-d) &
                        * aux2**(2*deg-1+d)
                   if ((x1_in.eq.0.).and.(x2_in.eq.0.8)) then
                      !                   print *, "            aux, aux2, val = ", aux, aux2, val
                   end if
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
       u = x1 - x2/sll_p_sqrt3
       v = x1 + x2/sll_p_sqrt3
       if (u.lt.0) then
          !Symmetry r2
          u = -u
          v = v + u
       end if
       if (u.lt.v/2) then
          !Symmetry r2+r3
          u = v - u
       end if
       g = u - v/2._f64
       if (v.gt.2._f64) then
          val = 0._f64
       else if (v.lt.1._f64) then
          val = 0.5_f64 + &
               ((5._f64/3._f64 - v/8.0_f64)*v - 3.0_f64)*v*v/4.0_f64 + &
               ((1._f64 - v/4._f64)*v + g*g/6._f64 - 1._f64)*g*g
       else if (u.gt.1._f64) then
          val = (v - 2._f64)*(v - 2._f64)*(v - 2._f64)*(g - 1._f64)/6.0_f64
       else
          val = 5._f64/6._f64 + &
               ((1._f64 + (1._f64/3._f64 - v/8._f64)*v)*v/4._f64 - 1._f64)*v + &
               ((1._f64 - v/4._f64)*v + g*g/6.0_f64 - 1._f64)*g*g
       end if
    end if

  end function chi_gen_val


  !---------------------------------------------------------------------------
  !> @brief Computes the value of a box spline
  !> @details This function computes the value of a box spline of degree
  !> deg at the point (x1,x2)
  !> @param[in] spline box spline which contains the reference hexagonal mesh
  !> @param[in] x1 real containing first coordinate of point
  !> @param[in] x2 real containing second coordinate of point
  !> @param[in] deg real containing the degree of the spline to be computed
  !> @return the value of the box spline at (x1,x2)
  function sll_f_compute_box_spline(spline, x1, x2, deg) result(val)
    type(sll_t_box_spline_2d), pointer, intent(in):: spline
    sll_real64, intent(in) :: x1
    sll_real64, intent(in) :: x2
    sll_int32,  intent(in) :: deg
    sll_real64 :: val
    sll_real64 :: x1_basis
    sll_real64 :: x2_basis

    x1_basis = change_basis_x1(spline, x1, x2)
    x2_basis = change_basis_x2(spline, x1, x2)

    val = chi_gen_val(x1_basis, x2_basis, deg)

  end function sll_f_compute_box_spline



  !---------------------------------------------------------------------------
  !> @brief 1st coo. of (x1, x2) in reference hex-mesh coo.
  !> @details This function allows to change a point of coordinates (x1, x2)
  !> on the spline basis to the mesh basis. Gives 1st coordinate.
  !> @param[in] spline box spline who contains the reference hexagonal mesh
  !> @param[in] x1 real containing first coordinate of point
  !> @param[in] x2 real containing second coordinate of point
  !> @return the first coordinate after change of coordinate sys.
  function change_basis_x1(spline, x1, x2) result(x1_basis)
    type(sll_t_box_spline_2d), pointer    :: spline
    sll_real64, intent(in) :: x1
    sll_real64, intent(in) :: x2
    sll_real64             :: delta_q
    sll_real64             :: k1_basis
    sll_real64             :: k2_basis
    sll_real64             :: x1_basis
    sll_real64             :: inv_delta_q
    sll_real64             :: q11, q12
    sll_real64             :: q21, q22
    sll_real64             :: r11, r12
    sll_real64             :: r21, r22

    ! Algorithms basis
    r11 = 0.5_f64
    r12 = -sll_p_sqrt3 * 0.5_f64
    r21 =  r11
    r22 = -r12

    ! Getting mesh generator vectors coordinates
    q11 = spline%mesh%r1_x1
    q12 = spline%mesh%r1_x2
    q21 = spline%mesh%r2_x1
    q22 = spline%mesh%r2_x2

    !change of basis :
    delta_q  = q11*q22 - q12*q21
    inv_delta_q = 1._f64/delta_q
    k1_basis = inv_delta_q*(q22*x1 - q21*x2)
    k2_basis = inv_delta_q*(q11*x2 - q12*x1)
    x1_basis = r11*k1_basis+r21*k2_basis
  end function change_basis_x1

  !---------------------------------------------------------------------------
  !> @brief 2nd coo. of (x1, x2) in reference hex-mesh coo.
  !> @details This function allows to change a point of coordinates (x1, x2)
  !> on the spline basis to the mesh basis. Gives 2nd coordinate.
  !> @param[in] spline box spline who contains the reference hexagonal mesh
  !> @param[in] x1 real containing first coordinate of point
  !> @param[in] x2 real containing second coordinate of point
  !> @return the second coordinate after change of coordinate sys.
  function change_basis_x2(spline, x1, x2) result(x2_basis)
    ! This function allows to change a point of coordinates (x1, x2)
    ! on the spline basis to the mesh basis
    type(sll_t_box_spline_2d), pointer           :: spline
    sll_real64, intent(in) :: x1
    sll_real64, intent(in) :: x2
    sll_real64             :: delta_q
    sll_real64             :: inv_delta_q
    sll_real64             :: k1_basis
    sll_real64             :: k2_basis
    sll_real64             :: x2_basis
    sll_real64             :: q11, q12
    sll_real64             :: q21, q22
    sll_real64             :: r11, r12
    sll_real64             :: r21, r22

    ! Getting spline generator vectors coordinates
    ! Algorithms basis
    r11 =  0.5_f64
    r12 = -0.5_f64 * sll_p_sqrt3
    r21 =  r11
    r22 = -r12

    ! Getting mesh generator vectors coordinates
    q11 = spline%mesh%r1_x1
    q12 = spline%mesh%r1_x2
    q21 = spline%mesh%r2_x1
    q22 = spline%mesh%r2_x2

    !change of basis :
    delta_q  = q11*q22 - q12*q21
    inv_delta_q = 1._f64/delta_q
    k1_basis = inv_delta_q*(q22*x1 - q21*x2)
    k2_basis = inv_delta_q*(q11*x2 - q12*x1)
    x2_basis = r12*k1_basis+r22*k2_basis

  end function change_basis_x2

  !---------------------------------------------------------------------------
  !> @brief Interpolates on point of cartesian coordinates (x1, x2)
  !> @details Interpolation with box splines of degree deg on the point
  !> of cartesian coordinates (x1, x2) on the mesh mesh_geom
  !> @param[in] mesh_geom geometric mesh
  !> @param[in] x1 real containing 1st coordinate of point where to interpolate
  !> @param[in] x2 real containing 2nd coordinate of point where to interpolate
  !> @param[in] spline box spline with coefficients already pre computed
  !> @param[in] deg integer with degree of splines
  !> @return real Interpolated value
  function sll_f_hex_interpolate_value(mesh_geom, x1, x2, spline, deg) result(val)
    type(sll_t_hex_mesh_2d), pointer, intent(in) :: mesh_geom
    type(sll_t_box_spline_2d), pointer           :: spline
    sll_int32,  intent(in)  :: deg
    sll_real64, intent(in)  :: x1,      x2
    sll_real64              :: xm1,     xm2
    sll_real64              :: x1_spl,  x2_spl
    sll_real64              :: r11,     r12
    sll_real64              :: r21,     r22
    sll_real64              :: r31,     r32
    sll_int32               :: k1_asso, k2_asso
    sll_int32               :: k1,      k2
    sll_int32               :: ki,      kj
    sll_int32               :: num_pts_tot
    sll_int32               :: num_cells
    sll_int32               :: distance
    sll_int32               :: ind
    sll_real64              :: val


    val = 0._f64

    r11 = mesh_geom%r1_x1
    r12 = mesh_geom%r1_x2
    r21 = mesh_geom%r2_x1
    r22 = mesh_geom%r2_x2
    r31 = mesh_geom%r3_x1
    r32 = mesh_geom%r3_x2

    num_pts_tot = mesh_geom%num_pts_tot
    num_cells = spline%mesh%num_cells

    ! First we need to compute the coordinates of
    ! the closest mesh point associated to (x1,x2)
    k1_asso = sll_f_cart_to_hex1(mesh_geom, x1, x2)
    k2_asso = sll_f_cart_to_hex2(mesh_geom, x1, x2)

    ! Then we will do a loop for all the points 
    ! on the envelopping rhomboid of radius=deg
    do ki= 1-deg, deg
       do kj= 1-deg, deg

          k1  = k1_asso + ki
          k2  = k2_asso + kj
          distance = sll_f_cells_to_origin(k1, k2)

          ! We test if we are in the domain
          if (distance.le.num_cells) then

             ind = spline%mesh%hex_to_global(k1, k2)

             ! We centralize and shift the coordinates
             ! i.e. centralize : xm = x - Rk
             !      shifting   : xm to the associated rhomboid point
             xm1 = x1 - r11*k1_asso - r21*k2_asso - ki*r11 - kj*r21
             xm2 = x2 - r12*k1_asso - r22*k2_asso - ki*r12 - kj*r22

             ! change of basis : geometrical basis => spline basis
             x1_spl = change_basis_x1(spline, xm1, xm2)
             x2_spl = change_basis_x2(spline, xm1, xm2)

             val = val + spline%coeffs(ind) * &
                  chi_gen_val(x1_spl, x2_spl, deg)
          else
             !! TREAT HERE BC
             ! TODO @LM
             if (spline%bc_type .eq. sll_p_dirichlet) then
                val = val !no update
                ind = spline%mesh%hex_to_global(k1_asso, k2_asso)
             else
                print *, "Error : Boundary condition type not yet implemented"
                STOP
             end if
          end if
       end do
    end do
  end function sll_f_hex_interpolate_value


  !---------------------------------------------------------
  !> @brief Computes indices of non null splines on a given cell
  !> @details The function returns for a given cell and a certain degree
  !> the indices of the splines of that degree that are different to 0
  !> on that cell.
  !> @param[IN] mesh hexagonal mesh, the domain
  !> @param[IN] cell_index index of the cell where we wish to know the
  !> indices of the non vanishing splines
  !> @param[IN] deg integer of the degree of the splines
  !> @param[OUT] index_nZ vector of size 3*deg*deg containing the global
  !> indices of all non-zero splines
  function non_zeros_splines(mesh, cell_index, deg) result(index_nZ)
    type(sll_t_hex_mesh_2d), pointer, intent(in) :: mesh
    sll_int32,  intent(in)  :: deg
    sll_int32,  intent(in)  :: cell_index
    sll_int32               :: index_nZ(3*deg*deg)
    !sll_int32               :: ierr
    sll_int32               :: nei_point
    sll_int32               :: non_Zero
    sll_int32               :: distance
    sll_int32               :: edge1, edge2, edge3
    sll_int32               :: first
    sll_int32               :: last
    sll_int32               :: type
    sll_int32               :: i, j
    sll_int32               :: last_point
    sll_int32               :: current_nZ

    ! Number of non zero splines on a cell:
    non_Zero = 3 * deg * deg
    index_nZ(1:non_Zero) = -1

    !type of cell
    call mesh%cell_type(cell_index, type)

    !Getting the cell vertices which are the 1st indices of the non null splines
    call sll_s_get_cell_vertices_index(mesh%center_cartesian_coord(1,cell_index), &
         mesh%center_cartesian_coord(2,cell_index),&
         mesh, &
         edge1, edge2, edge3)

    ! index_nZ(1) = edge1
    ! index_nZ(2) = edge2
    ! index_nZ(3) = edge3

    if (type.eq.2) then
       index_nZ(1) = edge1
       index_nZ(2) = edge2
       index_nZ(3) = edge3
    else
       index_nZ(1) = edge1
       index_nZ(2) = edge3
       index_nZ(3) = edge2
    end if

    current_nZ = 4
    do distance = 1,deg-1
       first = 1 + (distance-1)*(distance-1)*3
       last = distance * distance * 3
       do i=first,last
          last_point = index_nZ(i)
          if (last_point .eq. -1) then
             print *, "ERROR in non_zero_splines: wrong index, i=", i
          end if
          do j=2,7 ! this represents the direct neighbours for a point
             nei_point = sll_f_local_to_global(mesh, last_point, j)
             if ((nei_point.ne.-1).and.(.not.(ANY(index_nZ==nei_point)))) then
                index_nZ(current_nZ) = nei_point
                current_nZ = current_nZ + 1
             end if
          end do
       end do
    end do
  end function non_zeros_splines


  !---------------------------------------------------------------------------
  !> @brief Computes x-derivative on (x,y)
  !> @details Using the 5 point stencil computes the x-derivative on (x,y)
  !> @param[in] x1 real containing second coordinate of point
  !> @param[in] x2 real containing second coordinate of point
  !> @param[in] deg integer with degree of splines
  !> @return real derivative on x1 of box spline
  function sll_f_boxspline_x1_derivative(x1, x2, deg) result(val)
    sll_int32,  intent(in)  :: deg
    sll_real64, intent(in)  :: x1
    sll_real64, intent(in)  :: x2
    sll_real64 :: h
    sll_real64 :: fm2h
    sll_real64 :: fm1h
    sll_real64 :: fp1h
    sll_real64 :: fp2h
    sll_real64 :: val

    !Computing step on each direction
    h = max(10.*sll_p_epsilon_0*abs(x1), sll_p_epsilon_0)

    ! Finite difference method of order 5
    fm2h = chi_gen_val(x1-2.0_f64*h, x2, deg)
    fm1h = chi_gen_val(x1 - h,       x2, deg)
    fp2h = chi_gen_val(x1+2.0_f64*h, x2, deg)
    fp1h = chi_gen_val(x1 + h,       x2, deg)

    val = ( - fp2h + 8._f64 * fp1h - 8._f64 * fm1h + fm2h) / 12._f64 / h

  end function sll_f_boxspline_x1_derivative


  !---------------------------------------------------------------------------
  !> @brief Computes y-derivative on (x,y)
  !> @details Using the 5 point stencil computes the y-derivative on (x,y)
  !> @param[in] x1 real containing second coordinate of point
  !> @param[in] x2 real containing second coordinate of point
  !> @param[in] deg integer with degree of splines
  !> @return real derivative on x2 of box spline
  function sll_f_boxspline_x2_derivative(x1, x2, deg) result(val)
    sll_int32,  intent(in)  :: deg
    sll_real64, intent(in)  :: x1
    sll_real64, intent(in)  :: x2
    sll_real64 :: h
    sll_real64 :: fm2h
    sll_real64 :: fm1h
    sll_real64 :: fp1h
    sll_real64 :: fp2h
    sll_real64 :: val

    !Computing step on each direction
    h = max(10.*sll_p_epsilon_0*abs(x2), sll_p_epsilon_0)

    ! Finite difference method of order 5
    fm2h = chi_gen_val(x1, x2-2.0*h, deg)
    fm1h = chi_gen_val(x1, x2 - h,   deg)
    fp2h = chi_gen_val(x1, x2+2.0*h, deg)
    fp1h = chi_gen_val(x1, x2 + h,   deg)

    val = 0.25/3._f64/h * ( - fp2h + 8._f64 * fp1h - 8._f64 * fm1h + fm2h)

  end function sll_f_boxspline_x2_derivative


  !---------------------------------------------------------------------------
  !> @brief Computes the values or derivatives of box splines
  !> @details Depending on nderiv1 and nderiv2 will compute box spline or
  !> derivative with respect to x (nderiv1 > 0) or/and to y (nderiv2 > 0)
  !> @param[in] x1 real containing second coordinate of point
  !> @param[in] x2 real containing second coordinate of point
  !> @param[in] deg integer with degree of splines
  !> @param[in] nderiv1 integer number of times to derive on the x direction
  !> @param[in] nderiv2 integer number of times to derive on the y direction
  !> @return real nderiv-derivatives of boxspline
  function sll_f_boxspline_val_der(x1, x2, deg, nderiv1, nderiv2) result(val)
    sll_int32,  intent(in)  :: deg
    sll_int32,  intent(in)  :: nderiv1
    sll_int32,  intent(in)  :: nderiv2
    sll_real64, intent(in)  :: x1
    sll_real64, intent(in)  :: x2
    sll_real64 :: val
    sll_real64 :: x1_basis
    sll_real64 :: x2_basis
    sll_int32  :: ierr
    type(sll_t_hex_mesh_2d),   pointer  :: mesh
    type(sll_t_box_spline_2d), pointer  :: spline

    mesh => sll_f_new_hex_mesh_2d(1)
    spline => sll_f_new_box_spline_2d(mesh, sll_p_dirichlet)
    x1_basis = change_basis_x1(spline, x1, x2)
    x2_basis = change_basis_x2(spline, x1, x2)

    val = 0._f64

    if (nderiv1.eq.0) then
       if (nderiv2.eq.0) then
          !> no derivative to compute
          val = chi_gen_val(x1_basis, x2_basis, deg)
       else if (nderiv2.eq.1) then
          !> derivative with respect to the second coo
          val = sll_f_boxspline_x2_derivative(x1_basis, x2_basis, deg)
       else
          print *, "Error in sll_f_boxspline_val_der : cannot compute this derivative"
       end if
    else if (nderiv1.eq.1) then
       ! derivative with respecto to the first coo
       if (nderiv2.eq.0) then
          val = sll_f_boxspline_x1_derivative(x1_basis, x2_basis, deg)
       else
          print *, "Error in sll_f_boxspline_val_der : cannot compute this derivative"
       end if
    end if

    SLL_DEALLOCATE_ARRAY(spline%coeffs,ierr)
    SLL_DEALLOCATE(spline,ierr)
    call sll_o_delete(mesh)
  end function sll_f_boxspline_val_der

  !---------------------------------------------------------------------------
  !> @brief Writes connectivity for CAID / DJANGO
  !> @details write connectivity info for CAID/DJANGO. This function was
  !> intented to couple Pigasus poisson solver to the hex-mesh.
  !> Output file : boxsplines_connectivity.txt
  !> @param[in]  mesh pointer to the hexagonal mesh
  !> @param[in]  deg  integer designing boxsplines degree
  subroutine sll_s_write_connectivity(mesh, deg)
    type(sll_t_hex_mesh_2d), pointer :: mesh
    sll_int32, intent(in)          :: deg
    sll_int32                      :: out_unit
    character(len=28), parameter   :: name = "boxsplines_connectivity.txt"
    sll_int32                      :: nZ_indices(3*deg*deg)
    sll_int32  :: num_ele
    !sll_int32  :: ele_contained
    sll_int32  :: non_zero
    !sll_int32  :: ierr
    sll_int32  :: i
    !sll_int32  :: s2
    !sll_int32  :: s3
    !sll_int32  :: dist
    sll_int32  :: val

    ! Number of non Zero splines depends on the degree
    non_zero = 3*deg*deg

    ! We open file
    open( file=name, status="replace", form="formatted", newunit=out_unit )

    ! We write total number of cells
    write(out_unit, "(i6)") mesh%num_triangles

    do num_ele = 1,mesh%num_triangles
       ! We write cell ID number
       write(out_unit, "(i6)") sll_f_change_elements_notation(mesh, num_ele)
       ! We write number of non zero
       write(out_unit, "(i6)") non_zero
       ! We write the indices of the non zero splines
       nZ_indices = non_zeros_splines(mesh, num_ele, deg)
       do i=1,non_zero
          val = nZ_indices(i)
          write(out_unit, "(i6)", advance="no") val
          write(out_unit, "(a)", advance="no") ","
       end do
       write(out_unit,"(a)")""
    end do

    close(out_unit)

  end subroutine sll_s_write_connectivity


  !> @brief Generic sub-routine defined for 2D box spline types.
  !> Deallocates the memory associated with the given box spline object.
  !> @param[inout] spline_object.
  subroutine delete_box_spline_2d(spline)
    type(sll_t_box_spline_2d),  intent(inout), pointer :: spline
    sll_int32 :: ierr

    if( .not. associated(spline) ) then
       print *, 'delete_box_spline_2D(): passed spline is not associated'
       STOP
    end if
    call sll_o_delete(spline%mesh)
    SLL_DEALLOCATE(spline%coeffs, ierr)
    SLL_DEALLOCATE(spline, ierr)
  end subroutine delete_box_spline_2d

end module sll_m_box_splines
