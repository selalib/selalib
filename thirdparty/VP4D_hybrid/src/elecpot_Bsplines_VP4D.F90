!===========================================================================
!> For definition of the electrostatic potential and its derivatives for
!>  4D Vlasov-Poisson simulation in generalized coordinates using
!>  B-Spline element
!> 
!> \date 2014-08-20
!> \author V. Grandgirard, A. Back
!---------------------------------------------------------------------------
module elecpot_Bsplines_VP4D_module
#include "sll_working_precision.h"
#include "sll_assert.h"
#include "sll_memory.h"

  use field2d_VP4D_module
  use sll_m_arbitrary_degree_spline_interpolator_2d
  use sll_m_scalar_field_2d_base
  use sll_m_scalar_field_2d

  implicit none

  !===========================================================================
  !> 
  !---------------------------------------------------------------------------
  type, public :: elecpot_Bsplines_VP4D_t

    type(field2d_VP4D_t) :: Phi
    type(field2d_VP4D_t) :: dPhi_deta1
    type(field2d_VP4D_t) :: dPhi_deta2

    !> For computation of the derivative of the electrostatic potential
    type(sll_t_scalar_field_2d_discrete), pointer :: Phi_eta1eta2
    type(sll_t_arbitrary_degree_spline_interpolator_2d) :: interp2d_Phi_eta1eta2

  end type elecpot_Bsplines_VP4D_t
  !---------------------------------------------------------------------------

contains

  !===========================================================================
  !> Electrostatic potential and its derivatives: Allocation 
  !---------------------------------------------------------------------------
  subroutine new_elecpot_Bsplines_VP4D( &
      elecpot, &
      mesh4d, &
      bound_cond, &
      spline_degree )

    type(elecpot_Bsplines_VP4D_t) , intent(inout) :: elecpot
    type(mesh_VP4D_t)             , intent(in)    :: mesh4d
    type(boundary_conditions_2d_t), intent(in)    :: bound_cond
    type(spline_degree_2d_t)      , intent(in)    :: spline_degree

    !-> Local variables
    !--> For cartesian meshes
    class(sll_t_cartesian_mesh_2d), pointer :: cartesian_mesh2d_eta1eta2

    cartesian_mesh2d_eta1eta2 => mesh4d%eta1_eta2_mesh2d

    !*** Initialization of the 2D fields ***
    call new_field2d_VP4D( &
        elecpot%Phi, &
        'Phi', &
        mesh4d, &
        bound_cond, &
        spline_degree )

    call new_field2d_VP4D( &
        elecpot%dPhi_deta1, &
        'dPhi_deta1', &
        mesh4d, &
        bound_cond, &
        spline_degree )

    call new_field2d_VP4D( &
        elecpot%dPhi_deta2, &
        'dPhi_deta2', &
        mesh4d, &
        bound_cond, &
        spline_degree )

    !--> Initialization of the interpolators 
    call elecpot%interp2d_Phi_eta1eta2%init( &
        cartesian_mesh2d_eta1eta2%num_cells1+1, &
        cartesian_mesh2d_eta1eta2%num_cells2+1, &
        cartesian_mesh2d_eta1eta2%eta1_min, &
        cartesian_mesh2d_eta1eta2%eta1_max, &
        cartesian_mesh2d_eta1eta2%eta2_min, &
        cartesian_mesh2d_eta1eta2%eta2_max, &
        bound_cond%left_eta1, &
        bound_cond%right_eta1, &
        bound_cond%left_eta2, &
        bound_cond%right_eta2, &
        spline_degree%eta1, &
        spline_degree%eta2 )

    !--> Initialization of the fields
    elecpot%Phi_eta1eta2 => sll_f_new_scalar_field_2d_discrete( &
        "Phi_eta1eta2", &
        elecpot%interp2d_Phi_eta1eta2, &     
        mesh4d%transf_eta1eta2_xy, &
        bound_cond%left_eta1, &
        bound_cond%right_eta1, &
        bound_cond%left_eta2, &
        bound_cond%right_eta2 )

  end subroutine new_elecpot_Bsplines_VP4D
  !---------------------------------------------------------------------------


  !===========================================================================
  !> Electrostatic potential and its derivatives: Deallocation 
  !---------------------------------------------------------------------------
  subroutine delete_elecpot_Bsplines_VP4D( elecpot )

    type(elecpot_Bsplines_VP4D_t), intent(inout) :: elecpot

    call delete_field2d_VP4D( elecpot%Phi )
    call delete_field2d_VP4D( elecpot%dPhi_deta1 )
    call delete_field2d_VP4D( elecpot%dPhi_deta2 )

    call elecpot%interp2d_Phi_eta1eta2%delete( )

    call sll_o_delete( elecpot%Phi_eta1eta2  )

  end subroutine delete_elecpot_Bsplines_VP4D
  !---------------------------------------------------------------------------


  !===========================================================================
  !>  Compute the derivatives of the electrostatic potential
  !---------------------------------------------------------------------------
  subroutine derivatives_elecpot_Bsplines_VP4D( &
      elecpot, &
      mesh4d )

    type(elecpot_Bsplines_VP4D_t), intent(inout) :: elecpot
    type(mesh_VP4D_t)            , intent(in)    :: mesh4d

    !sll_int32  :: iloc1, iloc2
    sll_int32  :: ieta1, ieta2
    !sll_int32  :: loc2d_sz_x1, loc2d_sz_x2
    sll_int32  :: Neta1, Neta2
    sll_real64 :: eta1_point, eta2_point

    Neta1 = size(mesh4d%eta1_grid)
    Neta2 = size(mesh4d%eta2_grid)

    !*** Compute dPhi_deta1 and dPhi_deta2 (by using Phi%val2d_seqx1x2) ***
    call elecpot%Phi_eta1eta2%set_field_data( &
        elecpot%Phi%val2d_seqx1x2(:,:) )
    call elecpot%Phi_eta1eta2%update_interpolation_coefficients( )
    do ieta2 = 1,Neta2
      eta2_point = mesh4d%eta2_grid(ieta2)
      do ieta1 = 1,Neta1
        eta1_point = mesh4d%eta1_grid(ieta1)
        elecpot%dPhi_deta1%val2d_seqx1x2(ieta1,ieta2) = &
            elecpot%Phi_eta1eta2%first_deriv_eta1_value_at_point( &
            eta1_point, eta2_point )

        elecpot%dPhi_deta2%val2d_seqx1x2(ieta1,ieta2) = &
            elecpot%Phi_eta1eta2%first_deriv_eta2_value_at_point( &
            eta1_point, eta2_point )
      end do
    end do
     
  end subroutine derivatives_elecpot_Bsplines_VP4D
  !---------------------------------------------------------------------------

end module elecpot_Bsplines_VP4D_module
!---------------------------------------------------------------------------
