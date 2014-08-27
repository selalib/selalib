!**************************************************************
!  This module defines box splines (of arbitrary degree)
!  on a hexagonal mesh subdivided in equilateral triangles
!  Reference :
!     @Condat2006 "Three-directional box splines"
!  Author : 
!     Laura Mendoza (mela@ipp.mpg.de)
!************************************************************** 

module box_splines
#include "sll_working_precision.h"
#include "sll_memory.h"
#include "sll_assert.h"
#include "sll_splines.h"
#include "sll_utilities.h"
use hex_pre_filters
use hex_mesh

implicit none

type box_spline_2d
   type(hex_mesh_2d), pointer  :: mesh
   ! Boundary conditions definition
   sll_int32 SLL_PRIV :: bc_type
   ! Spline coefficients
   sll_real64, dimension(:), pointer :: coeffs
end type box_spline_2d

  
    
contains  ! ****************************************************************


#define MAKE_GET_SLOT_FUNCTION( fname, datatype, slot, ret_type )    \
  function fname( spline_obj ) result(val);                \
    type(datatype), pointer :: spline_obj;                 \
    ret_type :: val;                                       \
    val = spline_obj%slot;                                 \
  end function fname


  function new_box_spline_2d( &
       mesh,         &
       bc_type)
    
    type(box_spline_2d), pointer  :: new_box_spline_2d
    type(hex_mesh_2d),   pointer  :: mesh
    sll_int32,  intent(in)        :: bc_type
    sll_int32                     :: ierr


    SLL_ALLOCATE( new_box_spline_2d, ierr )
    
    new_box_spline_2d%mesh => mesh

    new_box_spline_2d%bc_type = bc_type

    SLL_ALLOCATE( new_box_spline_2d%coeffs(1:mesh%num_pts_tot), ierr )

  end function new_box_spline_2d


  subroutine compute_box_spline_2d( data, deg, spline )
    sll_real64, dimension(:), intent(in), target :: data  ! data to be fit
    sll_int32, intent(in)                        :: deg
    type(box_spline_2d), pointer, intent(in)     :: spline
    sll_int32  :: bc
    sll_int32  :: bc_selector


    if( .not. associated(spline) ) then
       ! FIXME: THROW ERROR
       print *, 'ERROR: compute_box_spline_2d(): ', &
            'uninitialized spline object passed as argument. '
       print *, "Exiting..."
       STOP
    end if
 
    bc = spline%bc_type
    
    ! Treat the bc_selector variable essentially like a bit field, to 
    ! accumulate the information on the different boundary conditions
    ! given. This scheme allows to add more types of boundary conditions
    ! if necessary.
    bc_selector = 0

    ! We make every case explicit to facilitate adding more BC types in
    ! the future.
    if( spline%bc_type .eq. SLL_DIRICHLET ) then 
       bc_selector = bc_selector + 1
    end if
    if( spline%bc_type .eq. SLL_PERIODIC ) then
       bc_selector = bc_selector + 2
    end if
    if( spline%bc_type .eq. SLL_NEUMANN ) then 
       bc_selector = bc_selector + 4
    end if

    select case (bc_selector)
       case ( 1 ) 
          ! boundary condition type is dirichlet
          call compute_box_spline_2d_diri( data, deg, spline )
       case ( 2 )
          ! boundary condition type is periodic
          call compute_box_spline_2d_prdc( data, spline )
       case ( 4 )
          ! boundary condition type is neumann
          call compute_box_spline_2d_neum( data, spline )
       case default
          print *, 'ERROR: compute_box_spline_2d(): ', &
            'did not recognize given boundary condition combination.'
       STOP
    end select

  
  end subroutine compute_box_spline_2d


  subroutine compute_box_spline_2d_diri( data, deg, spline )
    sll_real64, dimension(:), intent(in), target  :: data  ! data to be fit
    type(box_spline_2d), pointer                  :: spline
    sll_int32, intent(in)                         :: deg
    sll_int32  :: num_pts_tot
    sll_int32  :: k1_ref, k2_ref
    sll_int32  :: i,k
    sll_int32  :: nei
    sll_int32  :: num_pts_radius

    num_pts_tot = spline%mesh%num_pts_tot
    ! we will work on a radius of 'deg' cells
    ! we compute the number of total points on that radius
    num_pts_radius = 3*(2*deg)*(2*deg+1) + 1 

    do i = 1, num_pts_tot

       spline%coeffs(i) = real(0,f64)
       k1_ref = spline%mesh%global_to_hex1(i)
       k2_ref = spline%mesh%global_to_hex2(i)

       ! We don't need to fo through all points, just till a certain radius
       ! which depends on the degree of the spline we are evaluating
       do k = 1, num_pts_radius
          nei = spline%mesh%local_hex_to_global(k1_ref, k2_ref, k)
          if ((nei .lt. num_pts_tot).and.(nei .gt. 0)) then
             spline%coeffs(i) = spline%coeffs(i) + data(nei) * & 
                                pre_filter_pfir(spline%mesh, k, deg)
          else
             ! Boundary conditions (BC) to be treated here :
             ! With dirichlet boundary conditions data(out_of_domain) = 0
             spline%coeffs(i) = spline%coeffs(i)
          end if
       end do
    end do
    
  end subroutine compute_box_spline_2d_diri


  subroutine compute_box_spline_2d_prdc( data, spline )
    sll_real64, dimension(:), intent(in), target :: data  ! data to be fit
    type(box_spline_2d), pointer      :: spline
    sll_int32  :: num_pts_tot
    sll_int32  :: i

    print *, ' WARNING : BOUNDARY CONDITIONS PERIODIC NOT &
      & YET IMPLEMENTED'
    num_pts_tot = spline%mesh%num_pts_tot
    do i = 1, num_pts_tot
       spline%coeffs(i) = real(0,f64)*data(i)
    end do

  end subroutine compute_box_spline_2d_prdc


  subroutine compute_box_spline_2d_neum( data, spline )
    sll_real64, dimension(:), intent(in), target :: data  ! data to be fit
    type(box_spline_2d), pointer      :: spline
    sll_int32  :: num_pts_tot
    sll_int32  :: i

    print *, ' WARNING : BOUNDARY CONDITIONS PERIODIC NOT &
      & YET IMPLEMENTED'
    num_pts_tot = spline%mesh%num_pts_tot
    do i = 1, num_pts_tot
       spline%coeffs(i) = real(0,f64)*data(i)
    end do

  end subroutine compute_box_spline_2d_neum

 
  function choose (n, k) result (res)
    sll_int32, intent (in) :: n
    sll_int32, intent (in) :: k
    sll_real64 :: res
    if (k.lt.0) then
       res = 0._f64
    else if (n .lt. k) then
       res = 0._f64
    else 
       res = real(sll_factorial(n),f64) / real((sll_factorial(k) * sll_factorial(n - k)), f64)
    end if
  end function choose


  function chi_gen_val(x1_in,x2_in,deg) result(val)
    ! Computes the value of the box spline (chi) of degree deg
    ! on the point of cartesian coordiantes (x1, x2)
    ! Reference : @Condat and Van De Ville (2006)
    !             "Three directional Box Splines:
    !             Characterization and Efficient Evaluation."
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
       u = x1 - x2/sll_sqrt3
       v = x1 + x2/sll_sqrt3
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
             print *, " K = ", K
          end if
          do L = -deg, CEILING(v)-1
             if ((x1_in.eq.0.).and.(x2_in.eq.0.8)) then
                print *, "    L = ", L
             end if
             do i = 0,min(deg+K, deg+L)
                if ((x1_in.eq.0.).and.(x2_in.eq.0.8)) then
                   print *, "      i = ", i
                end if
                coeff = (-1.0_f64)**(K+L+i)* &
                     choose(deg,i-K)*     &
                     choose(deg,i-L)*     &
                     choose(deg,i)
                if ((x1_in.eq.0.).and.(x2_in.eq.0.8)) then
                   print *, "      coeff = ", coeff
                end if
                do d = 0,deg-1
                   if ((x1_in.eq.0.).and.(x2_in.eq.0.8)) then
                      print *, "          d = ", d
                   end if
                   aux=abs(v-L-u+K)  
                   aux2=(u-K+v-L-aux)/2._f64
                   if(aux2.lt.0.) then
                      aux2 = 0._f64
                   end if
                   val = val + coeff*choose(deg-1+d,d)   &
                        /real(sll_factorial(2*deg-1+d), f64) &
                        /real(sll_factorial(deg -1 -d), f64) &
                        * aux**(deg-1-d) &
                        * aux2**(2*deg-1+d)
                   if ((x1_in.eq.0.).and.(x2_in.eq.0.8)) then
                      print *, "            aux, aux2, val = ", aux, aux2, val
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
       u = x1 - x2/sll_sqrt3
       v = x1 + x2/sll_sqrt3
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


  function change_basis_x1(spline, x1, x2) result(x1_basis)
    ! This function allows to change a point of coordinates (x1, x2)
    ! on the spline basis to the mesh basis
    type(box_spline_2d), pointer    :: spline
    sll_real64, intent(in) :: x1
    sll_real64, intent(in) :: x2 
    sll_real64             :: delta_q
    sll_real64             :: k1_basis
    sll_real64             :: k2_basis
    sll_real64             :: x1_basis
    sll_real64             :: q11, q12
    sll_real64             :: q21, q22
    sll_real64             :: r11, r12
    sll_real64             :: r21, r22

    ! Algorithms basis
    r11 = 0.5_f64
    r12 = -sll_sqrt3 * 0.5_f64
    r21 = 0.5_f64
    r22 =  sll_sqrt3 * 0.5_f64

    ! Getting mesh generator vectors coordinates
    q11 = spline%mesh%r1_x1
    q12 = spline%mesh%r1_x2
    q21 = spline%mesh%r2_x1
    q22 = spline%mesh%r2_x2

    !change of basis :
    delta_q  = q11*q22 - q12*q21
    k1_basis = 1./delta_q*(q22*x1 - q21*x2)
    k2_basis = 1./delta_q*(q11*x2 - q12*x1)
    x1_basis = r11*k1_basis+r21*k2_basis
    
  end function change_basis_x1


  function change_basis_x2(spline, x1, x2) result(x2_basis)
    ! This function allows to change a point of coordinates (x1, x2)
    ! on the spline basis to the mesh basis
    type(box_spline_2d), pointer           :: spline
    sll_real64, intent(in) :: x1
    sll_real64, intent(in) :: x2 
    sll_real64             :: delta_q
    sll_real64             :: k1_basis
    sll_real64             :: k2_basis
    sll_real64             :: x2_basis
    sll_real64             :: q11, q12
    sll_real64             :: q21, q22
    sll_real64             :: r11, r12
    sll_real64             :: r21, r22

    ! Getting spline generator vectors coordinates
    ! Algorithms basis
    r11 = 0.5_f64
    r12 = -sll_sqrt3 * 0.5_f64
    r21 = 0.5_f64
    r22 =  sll_sqrt3 * 0.5_f64

    ! Getting mesh generator vectors coordinates
    q11 = spline%mesh%r1_x1
    q12 = spline%mesh%r1_x2
    q21 = spline%mesh%r2_x1
    q22 = spline%mesh%r2_x2

    !change of basis :
    delta_q  = q11*q22 - q12*q21
    k1_basis = 1./delta_q*(q22*x1 - q21*x2)
    k2_basis = 1./delta_q*(q11*x2 - q12*x1)
    x2_basis = r12*k1_basis+r22*k2_basis
    
  end function change_basis_x2


  function hex_interpolate_value(mesh_geom, x1, x2, spline, deg) result(val)
    ! Interpolates point of cartesian coordinates x1, x2)
    type(hex_mesh_2d), pointer, intent(in) :: mesh_geom
    type(box_spline_2d), pointer           :: spline
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
    k1_asso = cart_to_hex1(mesh_geom, x1, x2)
    k2_asso = cart_to_hex2(mesh_geom, x1, x2)

    ! Then we will do a loop for all the points 
    ! on the envelopping rhomboid of radius=deg
    do ki= 1-deg, deg
       do kj= 1-deg, deg
        
          k1  = k1_asso + ki
          k2  = k2_asso + kj
          distance = cells_to_origin(k1, k2)
          
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
             if (spline%bc_type .eq. SLL_DIRICHLET) then
                val = val
                ind = spline%mesh%hex_to_global(k1_asso, k2_asso)
             else
                print *, "Error : Boundary condition type not yet implemented"
                STOP
             end if
          end if
       end do
    end do
  end function hex_interpolate_value


end module box_splines

