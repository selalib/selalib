!===========================================================================
!> For definition of the electrostatic potential and its derivatives for
!>  4D drift-kinetic simulation in polar coordinates
!> 
!> \date 2014-09-18
!> \author V. Grandgirard
!---------------------------------------------------------------------------
module elecpot_DK4D_polar_module
#include "sll_working_precision.h"
#include "sll_assert.h"
#include "sll_memory.h"

  use field3d_DK4D_module
  use sll_m_boundary_condition_descriptors
  use sll_m_cartesian_meshes
  use sll_m_cubic_spline_interpolator_1d
  use sll_m_cubic_spline_interpolator_2d
  use sll_m_interpolators_1d_base
  use sll_m_interpolators_2d_base

  implicit none

  !===========================================================================
  !> 
  !---------------------------------------------------------------------------
  type, public :: elecpot_DK4D_polar_t

    type(field3d_DK4D_t) :: Phi
    type(field3d_DK4D_t) :: dPhi_deta1
    type(field3d_DK4D_t) :: dPhi_deta2
    type(field3d_DK4D_t) :: dPhi_deta3

    !> Interpolators for Phi
    class(sll_c_interpolator_2d), pointer :: interp2d_Phi_eta1eta2
    class(sll_c_interpolator_1d), pointer :: interp1d_Phi_eta3

  end type elecpot_DK4D_polar_t
  !---------------------------------------------------------------------------

contains

  !===========================================================================
  !> Electrostatic potential and its derivatives: Allocation 
  !---------------------------------------------------------------------------
  subroutine new_elecpot_DK4D_polar( &
      elecpot, &
      mesh4d, &
      bound_cond, &
      spline_degree )

    type(elecpot_DK4D_polar_t)    , intent(inout) :: elecpot
    type(mesh_DK4D_t)             , intent(in)    :: mesh4d
    type(boundary_conditions_3d_t), intent(in)    :: bound_cond
    type(spline_degree_3d_t)      , intent(in)    :: spline_degree

    !-> Local variables
    !--> For cartesian meshes
    type(sll_t_cartesian_mesh_2d), pointer :: cartesian_mesh2d_eta1eta2
    type(sll_t_cartesian_mesh_1d), pointer :: cartesian_mesh1d_eta3
    type(sll_t_cartesian_mesh_1d), pointer :: cartesian_mesh1d_vpar

    cartesian_mesh2d_eta1eta2 => mesh4d%eta1_eta2_mesh2d
    cartesian_mesh1d_eta3     => mesh4d%eta3_mesh1d
    cartesian_mesh1d_vpar     => mesh4d%vpar_mesh1d

    !*** Initialization of the 3D fields ***
    call new_field3d_DK4D( &
        elecpot%Phi, &
        'Phi', &
        mesh4d, &
        bound_cond, &
        spline_degree )

    call new_field3d_DK4D( &
        elecpot%dPhi_deta1, &
        'dPhi_deta1', &
        mesh4d, &
        bound_cond, &
        spline_degree )

    call new_field3d_DK4D( &
        elecpot%dPhi_deta2, &
        'dPhi_deta2', &
        mesh4d, &
        bound_cond, &
        spline_degree )

    call new_field3d_DK4D( &
        elecpot%dPhi_deta3, &
        'dPhi_deta3', &
        mesh4d, &
        bound_cond, &
        spline_degree )

    !*** Initialization of the interpolators for Phi ***
    elecpot%interp2d_Phi_eta1eta2 => sll_f_new_cubic_spline_interpolator_2d( &
        cartesian_mesh2d_eta1eta2%num_cells1+1, &
        cartesian_mesh2d_eta1eta2%num_cells2+1, &
        cartesian_mesh2d_eta1eta2%eta1_min, &
        cartesian_mesh2d_eta1eta2%eta1_max, &
        cartesian_mesh2d_eta1eta2%eta2_min, &
        cartesian_mesh2d_eta1eta2%eta2_max, &
        sll_p_hermite, &        
        sll_p_periodic, &
        const_eta1_min_slope = 0._f64, & !to prevent problem on the boundary
        const_eta1_max_slope = 0._f64 )

    elecpot%interp1d_Phi_eta3 => sll_f_new_cubic_spline_interpolator_1d( &
        cartesian_mesh1d_eta3%num_cells+1, &
        cartesian_mesh1d_eta3%eta_min, &
        cartesian_mesh1d_eta3%eta_max, &
        sll_p_periodic )

  end subroutine new_elecpot_DK4D_polar
  !---------------------------------------------------------------------------


  !===========================================================================
  !> Electrostatic potential and its derivatives: Deallocation 
  !---------------------------------------------------------------------------
  subroutine delete_elecpot_DK4D_polar( elecpot )

    type(elecpot_DK4D_polar_t), intent(inout) :: elecpot

    call delete_field3d_DK4D( elecpot%Phi )
    call delete_field3d_DK4D( elecpot%dPhi_deta1 )
    call delete_field3d_DK4D( elecpot%dPhi_deta2 )
    call delete_field3d_DK4D( elecpot%dPhi_deta3 )

    call elecpot%interp2d_Phi_eta1eta2%delete( )
!\begin todo: developper la methode delete
!VG!    call elecpot%interp1d_Phi_eta3%delete( ) 
!\end todo

  end subroutine delete_elecpot_DK4D_polar
  !---------------------------------------------------------------------------


  !===========================================================================
  !>  Compute the derivatives of the electrostatic potential
  !---------------------------------------------------------------------------
  subroutine derivatives_elecpot_DK4D_polar( &
      elecpot, &
      mesh4d )

    type(elecpot_DK4D_polar_t), intent(inout) :: elecpot
    type(mesh_DK4D_t)         , intent(in)    :: mesh4d

    sll_int32  :: iloc1, iloc2, iloc3
    sll_int32  :: ieta1, ieta2, ieta3
    sll_int32  :: loc3d_sz_x1, loc3d_sz_x2, loc3d_sz_x3
    sll_int32  :: Neta1, Neta2, Neta3
    sll_int32  :: ierr
    sll_real64 :: eta1_point, eta2_point, eta3_point
    sll_real64, dimension(:), pointer :: elec_pot1d_seqx3_tmp

    Neta1 = size(mesh4d%eta1_grid)
    Neta2 = size(mesh4d%eta2_grid)
    Neta3 = size(mesh4d%eta3_grid) 

    !*** Temporary allocation ***
    SLL_ALLOCATE( elec_pot1d_seqx3_tmp(Neta3), ierr )

    call sll_o_compute_local_sizes( elecpot%Phi%layout3d_seqx1x2, &
      loc3d_sz_x1, &
      loc3d_sz_x2, &
      loc3d_sz_x3)
    
    SLL_ASSERT( loc3d_sz_x1.eq.Neta1 )
    SLL_ASSERT( loc3d_sz_x2.eq.Neta2 )

    !*** Compute dPhi_deta1 and dPhi_deta2 (by using Phi%val3d_seqx1x2) ***
    do iloc3 = 1,loc3d_sz_x3
      call elecpot%interp2d_Phi_eta1eta2%compute_interpolants( &
          elecpot%Phi%val3d_seqx1x2(1:Neta1,1:Neta2,iloc3) )
      do ieta2 = 1,Neta2
        eta2_point = mesh4d%eta2_grid(ieta2)
        do ieta1 = 1,Neta1
          eta1_point = mesh4d%eta1_grid(ieta1)
          elecpot%dPhi_deta1%val3d_seqx1x2(ieta1,ieta2,iloc3) = &
              elecpot%interp2d_Phi_eta1eta2%interpolate_from_interpolant_derivative_eta1( &
              eta1_point, eta2_point )
          elecpot%dPhi_deta2%val3d_seqx1x2(ieta1,ieta2,iloc3) = &
              elecpot%interp2d_Phi_eta1eta2%interpolate_from_interpolant_derivative_eta2( &
              eta1_point, eta2_point )
        end do
      end do
    end do
     
    !*** Compute dPhi_deta3 (by using Phi%val3d_seqx3) ***
    call sll_o_compute_local_sizes( &
        elecpot%Phi%layout3d_seqx3, &
        loc3d_sz_x1, &
        loc3d_sz_x2, &
        loc3d_sz_x3)
    
    SLL_ASSERT( loc3d_sz_x3.eq.Neta3 )

    do iloc2 = 1,loc3d_sz_x2
      do iloc1 = 1,loc3d_sz_x1
        do ieta3 = 1,Neta3
          elec_pot1d_seqx3_tmp(ieta3) = &
              elecpot%Phi%val3d_seqx3(iloc1,iloc2,ieta3)
        end do       
        call elecpot%interp1d_Phi_eta3%compute_interpolants( &
            elec_pot1d_seqx3_tmp )
        do ieta3 = 1,Neta3
          eta3_point = mesh4d%eta3_grid(ieta3)
          elecpot%dPhi_deta3%val3d_seqx3(iloc1,iloc2,ieta3) = &
              elecpot%interp1d_Phi_eta3%interpolate_from_interpolant_derivative_eta1( &
              eta3_point )
        end do
      end do
    end do

    !*** Fill Phi%val3d_seqx3 ***
    call sll_o_apply_remap_3d( &
        elecpot%dPhi_deta1%seqx1x2_to_seqx3, &
        elecpot%dPhi_deta1%val3d_seqx1x2, &
        elecpot%dPhi_deta1%val3d_seqx3)    
    call sll_o_apply_remap_3d( &
        elecpot%dPhi_deta2%seqx1x2_to_seqx3, &
        elecpot%dPhi_deta2%val3d_seqx1x2, &
        elecpot%dPhi_deta2%val3d_seqx3 )    
    call sll_o_apply_remap_3d( &
        elecpot%dPhi_deta3%seqx3_to_seqx1x2, &
        elecpot%dPhi_deta3%val3d_seqx3, &
        elecpot%dPhi_deta3%val3d_seqx1x2 )

    !*** Deallocation ***    
    SLL_DEALLOCATE( elec_pot1d_seqx3_tmp,ierr )

  end subroutine derivatives_elecpot_DK4D_polar
  !---------------------------------------------------------------------------

end module elecpot_DK4D_polar_module
!---------------------------------------------------------------------------
