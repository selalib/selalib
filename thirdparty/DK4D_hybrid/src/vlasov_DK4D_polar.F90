!===========================================================================
!> Solving of the vlasov equation for
!>  4D drift-kinetic simulation in polar coordinates
!>
!> \date 2014-08-21
!> \author V. Grandgirard
!---------------------------------------------------------------------------
module vlasov_DK4D_polar_module
#include "sll_working_precision.h"
#include "sll_memory.h"
#include "sll_assert.h"

  use elecpot_DK4D_polar_module
  use fdistribu_DK4D_module
  use field3d_DK4D_module
  use magnetconf_DK4D_module
  use mesh_DK4D_module
  use sll_m_cubic_spline_interpolator_1d
  use sll_m_cubic_spline_interpolator_2d
  use sll_m_cartesian_meshes
  use sll_m_advection_1d_base
  use sll_m_advection_1d_periodic
  use sll_m_advection_1d_BSL
  use sll_m_advection_2d_base
  use sll_m_advection_2d_BSL
  use sll_m_characteristics_1d_base
  use sll_m_characteristics_1d_explicit_euler
  use sll_m_characteristics_2d_base
  use sll_m_characteristics_2d_explicit_euler
  use sll_m_periodic_interp, only:&
       sll_p_spline
  use sll_m_boundary_condition_descriptors, only: &
       sll_p_set_to_limit
  implicit none

  !===========================================================================
  !> 
  !---------------------------------------------------------------------------
  type, public :: vlasov_DK4D_polar_t

    !> For computation of the interpolation of the distribution function
    class(sll_c_interpolator_2d), pointer :: interp2d_f_eta1eta2
    class(sll_c_interpolator_1d), pointer :: interp1d_f_vpar

    !> For characteristics
    class(sll_c_characteristics_2d_base), pointer :: charac2d_eta1eta2
    class(sll_c_characteristics_1d_base), pointer :: charac1d_vpar

    !> For advection
    class(sll_c_advector_2d), pointer :: advec_x1x2
    class(sll_c_advector_1d), pointer :: advec_x3
    class(sll_c_advector_1d), pointer :: advec_x4

  end type vlasov_DK4D_polar_t
  !---------------------------------------------------------------------------

contains

  !===========================================================================
  !> Vlasov: Allocation 
  !---------------------------------------------------------------------------
  subroutine new_vlasov_DK4D_polar( &
      vlasov, &
      mesh4d, &
      bound_cond, &
      spline_degree )

    type(vlasov_DK4D_polar_t)     , intent(inout) :: vlasov
    type(mesh_DK4D_t)             , intent(in)    :: mesh4d
    type(boundary_conditions_4d_t), intent(in)    :: bound_cond
    type(spline_degree_4d_t)      , intent(in)    :: spline_degree

    !-> Local variables
    sll_int32 :: Neta3, Nvpar
    !--> For cartesian meshes
    type(sll_t_cartesian_mesh_2d), pointer :: logical_mesh2d_eta1eta2
    type(sll_t_cartesian_mesh_1d), pointer :: logical_mesh1d_eta3
    type(sll_t_cartesian_mesh_1d), pointer :: logical_mesh1d_vpar

    Neta3 = size(mesh4d%eta3_grid)
    Nvpar = size(mesh4d%vpar_grid)

    logical_mesh2d_eta1eta2 => mesh4d%eta1_eta2_mesh2d
    logical_mesh1d_eta3     => mesh4d%eta3_mesh1d
    logical_mesh1d_vpar     => mesh4d%vpar_mesh1d

    !*** Allocation of the 2D advection operator in (eta1,eta2) ***
    !--> Allocation of the 2D interpolator for the distribution function
    vlasov%interp2d_f_eta1eta2 => sll_f_new_cubic_spline_interpolator_2d( &
        logical_mesh2d_eta1eta2%num_cells1+1, &
        logical_mesh2d_eta1eta2%num_cells2+1, &
        logical_mesh2d_eta1eta2%eta1_min, &
        logical_mesh2d_eta1eta2%eta1_max, &
        logical_mesh2d_eta1eta2%eta2_min, &
        logical_mesh2d_eta1eta2%eta2_max, &
        sll_p_hermite, &
        sll_p_periodic, &
        const_eta1_min_slope = 0._f64, & !to prevent problem on the boundary
        const_eta1_max_slope = 0._f64 )
    !--> Allocation for the 2D characteristics
    vlasov%charac2d_eta1eta2 => sll_f_new_explicit_euler_2d_charac( &
        logical_mesh2d_eta1eta2%num_cells1+1, &
        logical_mesh2d_eta1eta2%num_cells2+1, &
        bc_type_1 = sll_p_set_to_limit, &
        bc_type_2 = sll_p_periodic, &
        eta1_min  = logical_mesh2d_eta1eta2%eta1_min, &
        eta1_max  = logical_mesh2d_eta1eta2%eta1_max, &
        eta2_min  = logical_mesh2d_eta1eta2%eta2_min, &
        eta2_max  = logical_mesh2d_eta1eta2%eta2_max )
    !--> Allocation for the 2D advection
    vlasov%advec_x1x2 => sll_f_new_advector_2d_bsl( &
        vlasov%interp2d_f_eta1eta2, &
        vlasov%charac2d_eta1eta2, &
        logical_mesh2d_eta1eta2%num_cells1+1, &
        logical_mesh2d_eta1eta2%num_cells2+1, &
        eta1_min = logical_mesh2d_eta1eta2%eta1_min, &
        eta1_max = logical_mesh2d_eta1eta2%eta1_max, &
        eta2_min = logical_mesh2d_eta1eta2%eta2_min, &
        eta2_max = logical_mesh2d_eta1eta2%eta2_max )

    !*** Allocation of the 1D advection operator in eta3 ***
    vlasov%advec_x3 => sll_f_new_periodic_1d_advector( &
        logical_mesh1d_eta3%num_cells, &
        logical_mesh1d_eta3%eta_min, &
        logical_mesh1d_eta3%eta_max, &
        sll_p_spline, & 
        spline_degree%eta3+1 ) 

    !*** Allocation of the 1D advection operator in vpar ***
    !\todo : Why is it a periodic advector in vpar direction ?
    !\todo : Why BSL_advect_1d_constant is not used ?
!VG!    vlasov%advec_x4 => sll_f_new_periodic_1d_advector( &
!VG!        logical_mesh1d_vpar%num_cells, &
!VG!        logical_mesh1d_vpar%eta_min, &
!VG!        logical_mesh1d_vpar%eta_max, &
!VG!        SPLINE, & 
!VG!        spline_degree%vpar+1 ) 
    !--> Allocation of the 1D interpolator
    vlasov%interp1d_f_vpar => sll_f_new_cubic_spline_interpolator_1d( &
        logical_mesh1d_vpar%num_cells+1, &
        logical_mesh1d_vpar%eta_min, &
        logical_mesh1d_vpar%eta_max, &
        sll_p_hermite, &
        slope_left  = 0._f64, &
        slope_right = 0._f64 )

    !--> Allocation for the 1D charac
    vlasov%charac1d_vpar => sll_f_new_charac_1d_explicit_euler( &
        logical_mesh1d_vpar%num_cells+1, &
        bc_type = sll_p_set_to_limit, &
        eta_min = logical_mesh1d_vpar%eta_min, &
        eta_max = logical_mesh1d_vpar%eta_max )
    !--> Allocation for the 1D advection
    vlasov%advec_x4 => sll_f_new_advector_1d_bsl( &
        vlasov%interp1d_f_vpar, &
        vlasov%charac1d_vpar, &
        logical_mesh1d_vpar%num_cells+1, &
        eta_coords = mesh4d%vpar_grid )

  end subroutine new_vlasov_DK4D_polar
  !---------------------------------------------------------------------------


  !===========================================================================
  !> Vlasov: Deallocation 
  !---------------------------------------------------------------------------
  subroutine delete_vlasov_DK4D_polar( vlasov )

    type(vlasov_DK4D_polar_t), intent(inout) :: vlasov

    ! TODO
    !call vlasov%interp2d_f_eta1eta2%free( )
!VG!    call vlasov%interp1d_f_vpar%delete()
!\begin todo: developper la methode delete
!VG!    call vlasov%advec_eta1_eta2%delete( )
!\end todo
    !call vlasov%advec_x3%delete( )
    !call vlasov%advec_x4%delete( ) 

  end subroutine delete_vlasov_DK4D_polar
  !---------------------------------------------------------------------------


  !===========================================================================
  !> Vlasov: Global iteration of the drift-kinetic 4D simulation
  !---------------------------------------------------------------------------
  subroutine iteration_vlasov_DK4D_polar( &
      vlasov, &
      fdistribu, &
      mesh4d, &
      magnetconf, &
      elecpot, &
      dt, &
      advec_vpar_step1_coef, &
      advec_eta3_step1_coef, &
      advec_eta1eta2_step1_coef, &
      advec_vpar_step2_coef, &
      advec_eta3_step2_coef )

    use sll_m_collective

    type(vlasov_DK4D_polar_t) , intent(inout) :: vlasov
    type(fdistribu_DK4D_t)    , intent(inout) :: fdistribu
    type(mesh_DK4D_t)         , intent(in)    :: mesh4d
    type(magnetconf_DK4D_t)   , intent(in)    :: magnetconf
    type(elecpot_DK4D_polar_t), intent(in)    :: elecpot
    sll_real64                , intent(in)    :: dt
    sll_real64                , intent(in)    :: advec_vpar_step1_coef
    sll_real64                , intent(in)    :: advec_eta3_step1_coef
    sll_real64                , intent(in)    :: advec_eta1eta2_step1_coef
    sll_real64      , optional, intent(in)    :: advec_vpar_step2_coef
    sll_real64      , optional, intent(in)    :: advec_eta3_step2_coef

    !--> local variables
    sll_int32 :: my_rank

    my_rank = sll_f_get_collective_rank( sll_v_world_collective )
        
    !*** Advection in vpar direction ***
    if ( advec_vpar_step1_coef.ne.0._F64 ) then
      if ( my_rank == 0) print*,' ==> first vpar advection with dt = ', &
          advec_vpar_step1_coef
      call advec1D_vpar_polar( &
          vlasov, &
          fdistribu, &
          mesh4d, &
          elecpot, &
          advec_vpar_step1_coef*dt )
    end if

    !*** Advection in eta3 direction ***
    if ( advec_eta3_step1_coef.ne.0._F64 ) then
      if ( my_rank == 0) print*,' ==> first eta3 advection with dt = ', &
          advec_eta3_step1_coef
      call advec1D_eta3_polar( &
          vlasov, &
          fdistribu, &
          mesh4d, &
          advec_eta3_step1_coef*dt )
    end if

    !*** Advection in (eta1,eta2) directions ***
    call sll_o_apply_remap_4d( &
        fdistribu%seqx3x4_to_seqx1x2, &
        fdistribu%val4d_seqx3x4,&
        fdistribu%val4d_seqx1x2 )
    if ( advec_eta1eta2_step1_coef.ne.0._F64 ) then
      if ( my_rank == 0) print*,' ==> (eta1,eta2) advection with dt = ', &
          advec_eta1eta2_step1_coef
      call advec2D_eta1eta2_polar( vlasov , &
          fdistribu, &
          mesh4d, &
          magnetconf, &
          elecpot, &
          advec_eta1eta2_step1_coef*dt )
      call sll_o_apply_remap_4d( &
          fdistribu%seqx1x2_to_seqx3x4, &
          fdistribu%val4d_seqx1x2, &
          fdistribu%val4d_seqx3x4 )
    end if

    !*** Advection in eta3 direction ***
    if ( present(advec_eta3_step2_coef) ) then
      if ( advec_eta3_step2_coef.ne.0._F64 ) then
        if ( my_rank == 0) print*,' ==> second eta3 advection with dt = ', &
            advec_eta3_step2_coef
        call advec1D_eta3_polar( &
            vlasov, &
            fdistribu, &
            mesh4d, &
            advec_eta3_step2_coef*dt )
      end if
    end if

    !*** Advection in vpar direction ***
    if ( present(advec_vpar_step2_coef) ) then
      if ( advec_vpar_step2_coef.ne.0._F64 ) then
        if ( my_rank == 0) print*,' ==> second vpar advection with dt = ', &
            advec_vpar_step2_coef
        call advec1D_vpar_polar( &
            vlasov, &
            fdistribu, &
            mesh4d, &
            elecpot, &
            advec_vpar_step2_coef*dt)
      end if
    end if

    !*** Fill fdistribu%val4d_seqx1x2 (to be able to solve QN)
    if ( present(advec_eta3_step2_coef) .or. &
        present(advec_vpar_step2_coef) ) then
      if ( (advec_vpar_step2_coef+advec_eta3_step2_coef).ne.0._F64 ) then
        call sll_o_apply_remap_4d( &
            fdistribu%seqx3x4_to_seqx1x2, &
            fdistribu%val4d_seqx3x4,&
            fdistribu%val4d_seqx1x2 )
      end if
    end if

  end subroutine iteration_vlasov_DK4D_polar
  !---------------------------------------------------------------------------


  !===========================================================================
  ! Vlasov: 1D advection in vpar direction
  !---------------------------------------------------------------------------
  subroutine advec1D_vpar_polar( &
      vlasov , &
      fdistribu, &
      mesh4d, &
      elecpot, &
      deltat_advec )

    type(vlasov_DK4D_polar_t) , intent(inout) :: vlasov
    type(fdistribu_DK4D_t)    , intent(inout) :: fdistribu
    type(mesh_DK4D_t)         , intent(in)    :: mesh4d
    type(elecpot_DK4D_polar_t), intent(in)    :: elecpot
    sll_real64                , intent(in)    :: deltat_advec

    sll_int32  :: ierr
    sll_int32  :: iloc1, iloc2
    sll_int32  :: ieta3, ivpar
    sll_int32  :: loc4d_sz_x1, loc4d_sz_x2
    sll_int32  :: loc4d_sz_x3, loc4d_sz_x4
    sll_int32  :: Neta3, Nvpar
    sll_real64 :: E_eta3
    sll_real64 :: vpar_point, vpar_min, vpar_max, alpha4
    sll_real64, dimension(:), pointer :: f1d_vpar_tmp
    sll_real64, dimension(:), pointer :: f1d_vpar_new_tmp

    SLL_ASSERT( size(fdistribu%val4d_seqx3x4,1).eq.size(elecpot%Phi%val3d_seqx3,1) )
    SLL_ASSERT( size(fdistribu%val4d_seqx3x4,2).eq.size(elecpot%Phi%val3d_seqx3,2) )
    SLL_ASSERT( size(fdistribu%val4d_seqx3x4,3).eq.size(elecpot%Phi%val3d_seqx3,3) )

    Neta3    = size( mesh4d%eta3_grid )
    Nvpar    = size( mesh4d%vpar_grid )
    vpar_min = minval( mesh4d%vpar_grid )
    vpar_max = maxval( mesh4d%vpar_grid )

    SLL_ALLOCATE( f1d_vpar_tmp(Nvpar), ierr )
    SLL_ALLOCATE( f1d_vpar_new_tmp(Nvpar), ierr )

    call sll_o_compute_local_sizes( fdistribu%layout4d_seqx3x4, &
      loc4d_sz_x1, &
      loc4d_sz_x2, &
      loc4d_sz_x3, &
      loc4d_sz_x4 )

    !---> dvpar/dt = E_eta3 with E_eta3 = -dPhi/deta3
    do iloc2 = 1,loc4d_sz_x2
      do iloc1 = 1,loc4d_sz_x1
        do ieta3 = 1,Neta3
          E_eta3 = - elecpot%dPhi_deta3%val3d_seqx3(iloc1,iloc2,ieta3)
          do ivpar = 1,Nvpar
            f1d_vpar_tmp(ivpar) = &
                fdistribu%val4d_seqx3x4(iloc1,iloc2,ieta3,ivpar)
          end do          

!baremettre
!\todo : Comprendre pourquoi cette advection ne marche pas
!VG!          call vlasov%advec_x4%advect_1d_constant( &
!VG!            E_eta3, &
!VG!            deltat_advec, &
!VG!            f1d_vpar_tmp(1:Nvpar), &
!VG!            f1d_vpar_new_tmp(1:Nvpar) )
          call vlasov%interp1d_f_vpar%compute_interpolants(f1d_vpar_tmp)
          do ivpar = 1,Nvpar
            alpha4     = deltat_advec*E_eta3 
            vpar_point = mesh4d%vpar_grid(ivpar) - alpha4
            vpar_point = max(min(vpar_point,vpar_max),vpar_min)
            f1d_vpar_new_tmp(ivpar) = &
                vlasov%interp1d_f_vpar%interpolate_from_interpolant_value(vpar_point)
          end do
!earemettre

          do ivpar = 1,Nvpar
            fdistribu%val4d_seqx3x4(iloc1,iloc2,ieta3,ivpar) = &
                f1d_vpar_new_tmp(ivpar)
          end do
        end do
      end do
    end do

    SLL_DEALLOCATE(f1d_vpar_tmp,ierr)
    SLL_DEALLOCATE(f1d_vpar_new_tmp,ierr)

  end subroutine advec1D_vpar_polar
  !---------------------------------------------------------------------------


  !===========================================================================
  ! Vlasov: 1D advection in eta3 direction
  !---------------------------------------------------------------------------
  subroutine advec1D_eta3_polar( &
      vlasov , &
      fdistribu, &
      mesh4d, &
      deltat_advec )

    type(vlasov_DK4D_polar_t) , intent(inout) :: vlasov
    type(fdistribu_DK4D_t)    , intent(inout) :: fdistribu
    type(mesh_DK4D_t)         , intent(in)    :: mesh4d
    sll_real64                , intent(in)    :: deltat_advec

    sll_int32  :: ierr
    sll_int32  :: iloc1, iloc2
    sll_int32  :: ieta3, ivpar
    sll_int32  :: loc4d_sz_x1, loc4d_sz_x2
    sll_int32  :: loc4d_sz_x3, loc4d_sz_x4
    sll_int32  :: Neta3, Nvpar
    sll_real64 :: vpar_point
    sll_real64, dimension(:), pointer :: f1d_eta3_tmp
    sll_real64, dimension(:), pointer :: f1d_eta3_new_tmp

    Neta3 = size(mesh4d%eta3_grid)
    Nvpar = size(mesh4d%vpar_grid)

    call sll_o_compute_local_sizes( fdistribu%layout4d_seqx3x4, &
        loc4d_sz_x1, &
        loc4d_sz_x2, &
        loc4d_sz_x3, &
        loc4d_sz_x4 )    

    !---> dphi/dt = vpar
    SLL_ALLOCATE(f1d_eta3_tmp(Neta3),ierr)
    SLL_ALLOCATE(f1d_eta3_new_tmp(Neta3),ierr)

    do iloc2 = 1,loc4d_sz_x2
      do iloc1 = 1,loc4d_sz_x1
        do ivpar = 1,Nvpar
          vpar_point = mesh4d%vpar_grid(ivpar)
          do ieta3 = 1,Neta3
            f1d_eta3_tmp(ieta3) = fdistribu%val4d_seqx3x4(iloc1,iloc2,ieta3,ivpar)
          end do
          call vlasov%advec_x3%advect_1d_constant(&
            vpar_point, &
            deltat_advec, &
            f1d_eta3_tmp(1:Neta3), &
            f1d_eta3_new_tmp(1:Neta3) )
          do ieta3 = 1,Neta3
            fdistribu%val4d_seqx3x4(iloc1,iloc2,ieta3,ivpar) = &
                f1d_eta3_new_tmp(ieta3)
          end do
        end do
      end do
    end do

    SLL_DEALLOCATE(f1d_eta3_tmp,ierr)
    SLL_DEALLOCATE(f1d_eta3_new_tmp,ierr)

  end subroutine advec1D_eta3_polar
  !---------------------------------------------------------------------------


  !===========================================================================
  ! Vlasov: 2D advection in (eta1,eta2) directions
  !---------------------------------------------------------------------------
  subroutine advec2D_eta1eta2_polar( &
      vlasov, &
      fdistribu, &
      mesh4d, &
      magnetconf, &
      elecpot, &
      deltat_advec )

    type(vlasov_DK4D_polar_t) , intent(inout) :: vlasov
    type(fdistribu_DK4D_t)    , intent(inout) :: fdistribu
    type(mesh_DK4D_t)         , intent(in)    :: mesh4d
    type(magnetconf_DK4D_t)   , intent(in)    :: magnetconf
    type(elecpot_DK4D_polar_t), intent(in)    :: elecpot
    sll_real64                , intent(in)    :: deltat_advec

    !--> Local variables
    sll_int32  :: ierr
    sll_int32  :: iloc3, iloc4
    sll_int32  :: ieta1, ieta2
    sll_int32  :: Neta1, Neta2
    sll_int32  :: loc4d_sz_x1, loc4d_sz_x2
    sll_int32  :: loc4d_sz_x3, loc4d_sz_x4
    sll_real64, dimension(:,:), pointer :: f2d_eta1eta2_tmp
    sll_real64, dimension(:,:), pointer :: f2d_eta1eta2_new_tmp
    !---> For electric field
    sll_real64 :: eta1_point 
    sll_real64, dimension(:,:), pointer :: efield_eta1_tmp
    sll_real64, dimension(:,:), pointer :: efield_eta2_tmp

    Neta1 = size( mesh4d%eta1_grid )
    Neta2 = size( mesh4d%eta2_grid )

    call sll_o_compute_local_sizes( fdistribu%layout4d_seqx1x2, &
        loc4d_sz_x1, &
        loc4d_sz_x2, &
        loc4d_sz_x3, &
        loc4d_sz_x4 )    

    SLL_ALLOCATE( f2d_eta1eta2_tmp(Neta1,Neta2), ierr )
    SLL_ALLOCATE( f2d_eta1eta2_new_tmp(Neta1,Neta2), ierr )
    SLL_ALLOCATE( efield_eta1_tmp(Neta1,Neta2), ierr )
    SLL_ALLOCATE( efield_eta2_tmp(Neta1,Neta2), ierr )

    do iloc4 = 1,loc4d_sz_x4
      do iloc3 = 1,loc4d_sz_x3
        do ieta2 = 1,Neta2
          do ieta1 = 1,Neta1
            eta1_point = mesh4d%eta1_grid(ieta1)
            f2d_eta1eta2_tmp(ieta1,ieta2) = &
                fdistribu%val4d_seqx1x2(ieta1,ieta2,iloc3,iloc4)
            !\todo : Understand why the signs are opposite in SELALIB simulation
            !\todo : Check what is done in the hybrid version
            !--> vExB.gradr = -1/(B r) dPhi/dtheta
            efield_eta1_tmp(ieta1,ieta2) = &
                - elecpot%dPhi_deta2%val3d_seqx1x2(ieta1,ieta2,iloc3) / &
                eta1_point
            !--> vExB.gradtheta = 1/(B r) dPhi/dr
            efield_eta2_tmp(ieta1,ieta2) = &
                elecpot%dPhi_deta1%val3d_seqx1x2(ieta1,ieta2,iloc3) / &
                eta1_point
          end do
        end do

        call vlasov%advec_x1x2%advect_2d( &
            efield_eta1_tmp, &
            efield_eta2_tmp, &
            deltat_advec, &
            f2d_eta1eta2_tmp(1:Neta1,1:Neta2), &
            f2d_eta1eta2_new_tmp(1:Neta1,1:Neta2) )
        do ieta2 = 1,Neta2
          do ieta1 = 1,Neta1        
            fdistribu%val4d_seqx1x2(ieta1,ieta2,iloc3,iloc4) = &
                f2d_eta1eta2_new_tmp(ieta1,ieta2)
          end do
        end do
      end do
    end do

    SLL_DEALLOCATE(f2d_eta1eta2_tmp,ierr)
    SLL_DEALLOCATE(f2d_eta1eta2_new_tmp,ierr)
    SLL_DEALLOCATE(efield_eta1_tmp,ierr)
    SLL_DEALLOCATE(efield_eta2_tmp,ierr)

  end subroutine advec2D_eta1eta2_polar
  !---------------------------------------------------------------------------

end module vlasov_DK4D_polar_module
!---------------------------------------------------------------------------
