!===========================================================================
!> Quasi-neutrality solver for
!>  4D Vlasov-Poisson hybrid simulation
!>
!> Rk1: Use DJANGO's poisson solver
!>
!> \date 2015-11-09
!> \author V. Grandgirard, L. Mendoza
!---------------------------------------------------------------------------
module QNsolver_DJANGO_VP4D_module
#include "sll_working_precision.h"
#include "sll_memory.h"
#include "sll_assert.h"

  use fdistribu_VP4D_module
  use field2d_VP4D_module
  use mesh_VP4D_module
  ! SeLaLib modules:
  use sll_m_cartesian_meshes
  ! Jorek modules:
  use JRK_MODEL_DEF
  use JRK_MODEL

  implicit none

  !-----------------------------------------------------------------------------
  type, public :: QNsolver_DJANGO_VP4D_t

    !> Flag to know if the mass matrix is also being used
    logical :: use_mass_matrix

    !> JOREK/DJANGO poisson, containing all useful entities for solver
    type(JRK_POISSON_2D) :: poisson

    !> Mesh where the poisson is bieng solved
    type(mesh_VP4D_t), pointer :: mesh4d

    !> RHS term
    type(field2D_VP4D_t) :: QN_RHS

    !> Pointer containing the contribution to the RHS
    sll_real64, dimension(:), allocatable :: ptr_contrib

    !> Pointer cointaining the connectivity between Selalib and DJANGO elements
    sll_int32, dimension(:,:), pointer   :: ptr_connec_table

  end type QNsolver_DJANGO_VP4D_t
  !-----------------------------------------------------------------------------


!===============================================================================
contains
!===============================================================================


  !-----------------------------------------------------------------------------
  !> @details Function to declare a new Django Poisson solver
  !-----------------------------------------------------------------------------
  subroutine new_QNsolver_DJANGO_VP4D( &
      QNsolver, &
      mesh4d,   &
      bound_cond, &
      spline_degree )

    type(QNsolver_DJANGO_VP4D_t)  , intent(inout) :: QNsolver !> solver to init
    type(mesh_VP4D_t),      target, intent(in)    :: mesh4d   !> mesh of model
    type(boundary_conditions_2d_t), intent(in)    :: bound_cond
    type(spline_degree_2d_t)      , intent(in)    :: spline_degree

    ! Local
    sll_int32  :: num_cells1
    sll_int32  :: num_cells2
    sll_int32  :: li_error
    sll_int32, dimension(:,:) ,   pointer :: ptr_connec_table
    class(sll_cartesian_mesh_2d), pointer :: cartesian_mesh2d

    ! Initialization of mesh and extraction of the mesh in the poloidal plane
    QNsolver%mesh4d  => mesh4d
    cartesian_mesh2d => mesh4d%eta1_eta2_mesh2d

    ! Create the connectivity table between selalib and django meshes
    num_cells1 = cartesian_mesh2d%num_cells1
    num_cells2 = cartesian_mesh2d%num_cells2
    SLL_ALLOCATE(ptr_connec_table(num_cells1, num_cells2), li_error)

    ! Initialization of the table, and update in solver.
    call connectivity_loc_to_glob( mesh4d % eta1_eta2_mesh2d, &
         ptr_connec_table )
    QNsolver % ptr_connec_table => ptr_connec_table

    !*** Initialization of the 2D field for RHS term ***
    call new_field2D_VP4D( &
        QNsolver%QN_RHS, &
        "QN_RHS_field", &
        mesh4d, &
        bound_cond, &
        spline_degree )

  end subroutine new_QNsolver_DJANGO_VP4D
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  !> @details Function to declare a new Django Poisson solver
  !-----------------------------------------------------------------------------
  subroutine delete_QNsolver_DJANGO_VP4D( QNsolver )

    type(QNsolver_DJANGO_VP4D_t), intent(inout) :: QNsolver

    call delete_field2D_VP4D( QNsolver%QN_RHS )

  end subroutine delete_QNsolver_DJANGO_VP4D


  !-----------------------------------------------------------------------------
  !> @details Function to declare a new Django Poisson solver
  !-----------------------------------------------------------------------------
  subroutine factorize_QNsolver_DJANGO_VP4D( &
      QNsolver, &
      testcase_jorek_id, &
      use_mass_matrix )

    type(QNsolver_DJANGO_VP4D_t), intent(inout) :: QNsolver !> solver to init
    sll_int32                   , intent(in)    :: testcase_jorek_id
    logical                     , intent(in)    :: use_mass_matrix

    ! Local
    sll_int32  :: n_rows
    sll_int32  :: li_error
    sll_int32,  dimension(2) :: arr_shape

    ! Initialization of the poisson model using the flag for the mass matrix
    QNsolver%use_mass_matrix = use_mass_matrix


    call initialize_model( QNsolver%poisson, &
        use_mass_matrix = QNsolver%use_mass_matrix, &
        testcase_id = testcase_jorek_id )

    ! Assembling model
    call ASSEMBLE_MODEL( QNsolver%poisson )

    ! Getting size of the RHS and allocating it by getting matrix shape
    call GET_SHAPE_MODEL( QNsolver%poisson, arr_shape )
    n_rows = arr_shape(1) ! RHS size is the 1st size
    SLL_ALLOCATE(QNsolver%ptr_contrib(1:n_rows), li_error)

  end subroutine factorize_QNsolver_DJANGO_VP4D
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  subroutine compRHS_QNsolver_DJANGO_VP4D(QNsolver, &
      fdistribu )

    type(QNsolver_DJANGO_VP4D_t), intent(inout) :: QNsolver
    type(fdistribu_VP4D_t)      , intent(in)    :: fdistribu

    !LOCAL 
    sll_int32   :: iloc1, iloc2
    sll_int32   :: Neta1_loc, Neta2_loc
    sll_int32   :: Neta1, Neta2
    sll_int32   :: ieta1, ieta2
    sll_int32   :: num_cells1, num_cells2
    sll_int32   :: num_quad
    sll_int32   :: iquad
    sll_int32   :: i1, i2
    sll_int32   :: i3, i4
    sll_int32   :: Nvx, Nvy
    sll_int32   :: i_elmt
    sll_int32   :: ierr
    sll_real64  :: eta1_coo, eta2_coo
    sll_real64  :: delta_eta1, delta_eta2
    sll_real64  :: lr_weight
    sll_real64  :: val_at_quad
    sll_real64  :: determinant
    sll_real64  :: dv_space, delta_f, intf_dvxdvy
    sll_real64, dimension(1,1:2) :: jrk_pos
    sll_real64, dimension(:,:), allocatable :: quad_points
    sll_real64, dimension(:),   allocatable :: quad_weights
    class(sll_cartesian_mesh_2d), pointer :: mesh2d
    type(sll_arbitrary_degree_spline_interpolator_2d) :: interp2d_QN_RHS !> interpolator for phi
    !---> Used to compute \int rho |J| deta1 deta2/ \int |J| deta1 deta2
    sll_real64 :: jacob_tmp, val_jac
    sll_real64 :: coef_int_deta1, coef_int_deta2
    sll_real64 :: int_J_deta1deta2, int_rho_J_deta1deta2


    mesh2d => QNsolver%mesh4d%eta1_eta2_mesh2d
    num_cells1 = mesh2d%num_cells1
    num_cells2 = mesh2d%num_cells2

    Neta1_loc = size(fdistribu%val4d_seqx3x4,1)
    Neta2_loc = size(fdistribu%val4d_seqx3x4,2)
    Nvx       = size(fdistribu%val4d_seqx3x4,3)
    Nvy       = size(fdistribu%val4d_seqx3x4,4)

    !*** Computation of the RHS locally in (x1,x2) directions ***
    do iloc2 = 1,Neta2_loc
      do iloc1 = 1,Neta1_loc
        intf_dvxdvy = 0._f64
        do i3 = 1,Nvx
          do i4 = 1,Nvy
            dv_space = QNsolver%mesh4d%coef_int_dvx(i3) * &
                QNsolver%mesh4d%coef_int_dvy(i4)
!VG!            delta_f = fdistribu%val4d_seqx3x4(iloc1,iloc2,i3,i4) - &
!VG!                equilibrium%feq_vxvy(i3,i4)
            delta_f = fdistribu%val4d_seqx3x4(iloc1,iloc2,i3,i4)
            intf_dvxdvy = intf_dvxdvy + delta_f*dv_space
          end do
        end do
        QNsolver%QN_RHS%val2d_parx1x2(iloc1,iloc2) = intf_dvxdvy
      end do
    end do

    !--> Compute rho(x,y) = \int f(x,y,vx,vy) dvx dvy
    !-->   (sequential in (x1,x2)=(x,y)
    call remap_parx1x2_to_seqx1x2_field2d_VP4D(QNsolver%QN_RHS)

    !--> Compute \int rho(eta1,eta2) |J| deta1 deta2 / \int |J| deta1 deta2
    Neta1 = size(QNsolver%QN_RHS%val2d_seqx1x2,1)
    Neta2 = size(QNsolver%QN_RHS%val2d_seqx1x2,2)
    int_J_deta1deta2     = 0._f64
    int_rho_J_deta1deta2 = 0._f64
    do i2 = 1,Neta2
      coef_int_deta2 = QNsolver%mesh4d%coef_int_deta2(i2)
      do i1 = 1,Neta1
        jacob_tmp      = QNsolver%mesh4d%transf_eta1eta2_xy%jacobian_at_node(i1,i2)
        val_jac        = abs(jacob_tmp)
        coef_int_deta1 = QNsolver%mesh4d%coef_int_deta1(i1)

        int_J_deta1deta2     = int_J_deta1deta2 + &
            coef_int_deta1*coef_int_deta2*val_jac
        int_rho_J_deta1deta2 = int_rho_J_deta1deta2 + &
            QNsolver%QN_RHS%val2d_seqx1x2(i1,i2) * &
            coef_int_deta1*coef_int_deta2*val_jac
      end do
    end do
    ! TODO : both int_rho_... and int_J... are being multiplied by val_jac
    int_rho_J_deta1deta2 = int_rho_J_deta1deta2/int_J_deta1deta2


    do i2 = 1,Neta2
      do i1 = 1,Neta1
        QNsolver%QN_RHS%val2d_seqx1x2(i1,i2) = &
            QNsolver%QN_RHS%val2d_seqx1x2(i1,i2)-int_rho_J_deta1deta2
      end do
    end do

    QNsolver%ptr_contrib(:) = 0._f64

    ! Initialization and getting the quadrature points:
    num_quad = QNsolver % poisson % space_trial % quadrature % oi_n_points
    SLL_ALLOCATE(quad_points(2, num_quad), ierr)
    SLL_ALLOCATE(quad_weights(num_quad), ierr)
    quad_points(1, :) = QNsolver % poisson % space_trial % quadrature % opr_points(1, :)
    quad_points(2, :) = QNsolver % poisson % space_trial % quadrature % opr_points(2, :)
    quad_weights( : ) = QNsolver % poisson % space_trial % quadrature % opr_weights(:)

    ! Initialization of Phi interpolator (to interpolate elec. pot. in quad points)
    ! Initialization of the interpolators 
    call interp2d_QN_RHS%initialize( &
        mesh2d%num_cells1 + 1, &
        mesh2d%num_cells2 + 1, &
        mesh2d%eta1_min, &
        mesh2d%eta1_max, &
        mesh2d%eta2_min, &
        mesh2d%eta2_max, &
        QNsolver%QN_RHS%bound_cond%left_eta1, &
        QNsolver%QN_RHS%bound_cond%right_eta1, &
        QNsolver%QN_RHS%bound_cond%left_eta2, &
        QNsolver%QN_RHS%bound_cond%right_eta2, &
        QNsolver%QN_RHS%spline_degree%eta1, &
        QNsolver%QN_RHS%spline_degree%eta2 )

    ! Computing the coefficients, given the values at mesh points 
    call interp2d_QN_RHS%compute_interpolants( &
         QNsolver%QN_RHS%val2d_seqx1x2)

    delta_eta1 = mesh2d % delta_eta1
    delta_eta2 = mesh2d % delta_eta2

    do ieta2 = 1, num_cells2
      do ieta1 = 1, num_cells1
         do iquad = 1, num_quad

            eta1_coo = mesh2d % eta1_node(ieta1, ieta2) + quad_points(1, iquad) * delta_eta1
            eta2_coo = mesh2d % eta2_node(ieta1, ieta2) + quad_points(2, iquad) * delta_eta2

            ! Getting the element index in jorek
            call find_jorek_position(QNsolver, eta1_coo, eta2_coo, i_elmt, jrk_pos)

            ! TODO : the jacobian should be computed on the quadrature points, not at the knots!
            determinant = QNsolver%mesh4d%transf_eta1eta2_xy%jacobian(eta1_coo, eta2_coo)
            val_at_quad = interp2d_QN_RHS % interpolate_value(eta1_coo, eta2_coo)
            lr_weight   = val_at_quad * quad_weights(iquad) * determinant &
                 / num_cells1 / num_cells2

            ! adding contribution in the table "lpr_contrib"
            ! routine in jorek/src/fem
            call add_contribution_field_model( QNsolver%poisson, &
                 i_elmt, jrk_pos, lr_weight, QNsolver%ptr_contrib )

         end do
      end do
   end do

   call set_rhs_model( QNsolver%poisson, QNsolver%ptr_contrib, &
        use_mass_matrix = QNsolver%use_mass_matrix )

   call interp2d_QN_RHS%delete( )

  end subroutine compRHS_QNsolver_DJANGO_VP4D
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  subroutine solve_QNsolver_DJANGO_VP4D( QNsolver, &
      fdistribu, &
      Phi )

    type(QNsolver_DJANGO_VP4D_t), intent(inout) :: QNsolver
    type(fdistribu_VP4D_t)      , intent(in)    :: fdistribu
    type(field2d_VP4D_t)        , intent(inout) :: Phi

    sll_int32  :: ieta1, ieta2
    sll_int32  :: num_cells1, num_cells2
    sll_int32  :: i_elmt
    sll_real64 :: eta1_coo
    sll_real64 :: eta2_coo
    sll_real64, dimension(3,1)   :: lpr_fields ! field eval
    sll_real64, dimension(1,1:2) :: jrk_pos
    class(sll_cartesian_mesh_2d), pointer :: mesh2d
    sll_int32, dimension(2) :: shape_rhs
    sll_int32 :: ierr
    sll_real64, dimension(:), allocatable :: temp_rhs

    call compRHS_QNsolver_DJANGO_VP4D( QNsolver, fdistribu )

    call SOLVE_MODEL(QNsolver%poisson)
!    call PROJECT_MODEL(QNsolver%poisson)

    call GET_SHAPE_MODEL(QNsolver%poisson, shape_rhs)
    SLL_ALLOCATE(temp_rhs(shape_rhs(1)), ierr)
    call get_rhs_model( QNsolver%poisson, temp_rhs, QNsolver%use_mass_matrix)

    mesh2d => QNsolver%mesh4d%eta1_eta2_mesh2d
    num_cells1 = mesh2d%num_cells1
    num_cells2 = mesh2d%num_cells2

    do ieta2 = 1,num_cells2
      do ieta1 = 1,num_cells1

         eta1_coo = mesh2d % eta1_node(ieta1, ieta2)
         eta2_coo = mesh2d % eta2_node(ieta1, ieta2)

        call find_jorek_position( QNsolver, eta1_coo, eta2_coo, i_elmt, jrk_pos)

        lpr_fields (:,:) = 0._f64
        call evaluate_field_model (QNsolver%poisson, i_elmt, &
            jrk_pos, lpr_fields)

        Phi%val2d_seqx1x2(ieta1,ieta2) = lpr_fields(1,1)
      end do
    end do

    call DIAGNOSTICS_MODEL(QNsolver %  poisson, 0 )

  end subroutine solve_QNsolver_DJANGO_VP4D
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  !> @details This function returns the connectivity table between the local
  !> and the global position
  !-----------------------------------------------------------------------------
  subroutine connectivity_loc_to_glob(mesh2d, connec_table)
    type(sll_cartesian_mesh_2d), pointer , intent(in)    :: mesh2d
    sll_int32, dimension(:,:),             intent(inout) :: connec_table
    
    ! Local
    sll_int32  :: li_x
    sll_int32  :: li_y
    sll_int32  :: li_i

    SLL_ASSERT(UBOUND(connec_table,1).eq.(mesh2d%num_cells1))
    SLL_ASSERT(UBOUND(connec_table,2).eq.(mesh2d%num_cells2))

    ! TODO : to change, adapt for SL
    li_i = 1
    do li_y = 1, mesh2d % num_cells2
       do li_x = 1, mesh2d % num_cells1
          connec_table(li_x,li_y) = li_i
          li_i = li_i + 1
       end do
    end do

  end subroutine connectivity_loc_to_glob
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  subroutine find_jorek_position( QNsolver, eta1_coo, eta2_coo, i_elmt, jrk_pos )

    type(QNsolver_DJANGO_VP4D_t), intent(in)  :: QNsolver
    sll_real64, intent(in)                     :: eta1_coo
    sll_real64, intent(in)                     :: eta2_coo
    sll_int32,  intent(out)                   :: i_elmt
    sll_real64, dimension(1,1:2), intent(out) :: jrk_pos
    !LOCAL
    sll_int32  :: li_selmt1 
    sll_int32  :: li_selmt2
    sll_int32  :: num_cells
    sll_real64 :: delta_eta1
    sll_real64 :: delta_eta2

    ! Computing the deltas:
    delta_eta1 = QNsolver % mesh4d % eta1_eta2_mesh2d % delta_eta1
    delta_eta2 = QNsolver % mesh4d % eta1_eta2_mesh2d % delta_eta2

    ! Computing the element indices:
    li_selmt1 = INT(eta1_coo/delta_eta1 + 1) 
    li_selmt2 = INT(eta2_coo/delta_eta2 + 1)

    ! Using connectivity table to find elt number in jorek
    i_elmt = QNsolver%ptr_connec_table(li_selmt1, li_selmt2)

    ! finding position in jorek setting
    jrk_pos(1,1) = eta1_coo/delta_eta1 - int(eta1_coo/delta_eta1)
    jrk_pos(1,2) = eta2_coo/delta_eta2 - int(eta2_coo/delta_eta2)

    ! Testing if i_elmt is in the right boudaries
    num_cells = QNsolver % mesh4d % eta1_eta2_mesh2d % num_cells1 * &
         QNsolver % mesh4d % eta1_eta2_mesh2d % num_cells2
    if ((i_elmt .gt. num_cells) .or. (i_elmt .lt. 0)) then
       print *, "ERROR in find_jorek_position():"
       print *, " Arguments are eta1, eta2 =", eta1_coo, eta2_coo
       print *, " Number of cells =", QNsolver % mesh4d % eta1_eta2_mesh2d % num_cells1, QNsolver % mesh4d % eta1_eta2_mesh2d % num_cells2
       print *, " Deltas =", delta_eta1, delta_eta2
       print *, " eta interv =", QNsolver % mesh4d % eta1_eta2_mesh2d % eta1_min, QNsolver % mesh4d % eta1_eta2_mesh2d % eta1_max
       print *, " Computed index of element in each direction =", li_selmt1, li_selmt2
       print *, " Computed RESULT = ", i_elmt
       STOP
    end if
  end subroutine find_jorek_position
  !-----------------------------------------------------------------------------


end module QNsolver_DJANGO_VP4D_module
