!===========================================================================
!> Solving of the vlasov equation for
!>  4D drift-kinetic hybrid simulation
!>
!> \date 2014-08-21
!> \author V. Grandgirard
!---------------------------------------------------------------------------
module vlasov_DK4D_hybrid_module
#include "sll_working_precision.h"
#include "sll_memory.h"
#include "sll_assert.h"

  use elecpot_DK4D_hybrid_module
  use fdistribu_DK4D_module
  use field3d_DK4D_module
  use magnetconf_DK4D_module
  use mesh_DK4D_module
!eaoter
!VG!  use sll_m_arbitrary_degree_spline_interpolator_1d
!VG!  use sll_m_arbitrary_degree_spline_interpolator_2d
  use sll_m_cubic_spline_interpolator_1d
  use sll_m_cubic_spline_interpolator_2d
!baoter
  use sll_m_cartesian_meshes
  use sll_m_scalar_field_1d
  use sll_m_scalar_field_2d


  implicit none

  !===========================================================================
  !> 
  !---------------------------------------------------------------------------
  type, public :: vlasov_DK4D_hybrid_t

    !> For computation of the interpolation of the distribution function
!VG!    type(sll_t_arbitrary_degree_spline_interpolator_2d) :: interp2d_f_eta1eta2
!VG!    type(sll_arbitrary_degree_spline_interpolator_1d) :: interp1d_f_eta3
!VG!    type(sll_arbitrary_degree_spline_interpolator_1d) :: interp1d_f_vpar
    type(sll_t_cubic_spline_interpolator_2d) :: interp2d_f_eta1eta2
    type(sll_t_cubic_spline_interpolator_1d) :: interp1d_f_eta3
    type(sll_t_cubic_spline_interpolator_1d) :: interp1d_f_vpar

  end type vlasov_DK4D_hybrid_t
  !---------------------------------------------------------------------------

contains

  !===========================================================================
  !> Vlasov: Allocation 
  !---------------------------------------------------------------------------
  subroutine new_vlasov_DK4D_hybrid( vlasov, &
      mesh4d, &
      bound_cond, &
      spline_degree )

    type(vlasov_DK4D_hybrid_t)    , intent(inout) :: vlasov
    type(mesh_DK4D_t)             , intent(in)    :: mesh4d
    type(boundary_conditions_4d_t), intent(in)    :: bound_cond
    type(spline_degree_4d_t)      , intent(in)    :: spline_degree

    !-> Local variables
    sll_int32 :: Neta3, Nvpar
    !--> For cartesian meshes
    type(sll_t_cartesian_mesh_2d), pointer :: logical_mesh2d_eta1eta2
    type(sll_t_cartesian_mesh_1d), pointer :: logical_mesh1d_eta3
    type(sll_t_cartesian_mesh_1d), pointer :: logical_mesh1d_vpar
    !--> For electic field boundary conditions and spline degree
    type(boundary_conditions_3d_t) :: elec_field_bc
    type(spline_degree_3d_t)       :: elec_field_spline_degree

    Neta3 = size(mesh4d%eta3_grid)
    Nvpar = size(mesh4d%vpar_grid)

    logical_mesh2d_eta1eta2 => mesh4d%eta1_eta2_mesh2d
    logical_mesh1d_eta3     => mesh4d%eta3_mesh1d
    logical_mesh1d_vpar     => mesh4d%vpar_mesh1d

    !*** Allocation for the interpolations of the distribution function ***
    !--> Initialization of the interpolators 
    !\todo : Find a solution of be able to use the arbitrary spline
    call vlasov%interp2d_f_eta1eta2%init( &
      logical_mesh2d_eta1eta2%num_cells1+1, &
      logical_mesh2d_eta1eta2%num_cells2+1, &
      logical_mesh2d_eta1eta2%eta1_min, &
      logical_mesh2d_eta1eta2%eta1_max, &
      logical_mesh2d_eta1eta2%eta2_min, &
      logical_mesh2d_eta1eta2%eta2_max, &
       sll_p_hermite, &
       sll_p_periodic )

    call vlasov%interp1d_f_eta3%init( &
       Neta3, &
       logical_mesh1d_eta3%eta_min, &
       logical_mesh1d_eta3%eta_max, &
       sll_p_periodic )

    call vlasov%interp1d_f_vpar%init( &
       Nvpar, &
       logical_mesh1d_vpar%eta_min, &
       logical_mesh1d_vpar%eta_max, &
       sll_p_hermite )

!VG!    call vlasov%interp2d_f_eta1eta2%initialize( &
!VG!      logical_mesh2d_eta1eta2%num_cells1+1, &
!VG!      logical_mesh2d_eta1eta2%num_cells2+1, &
!VG!      logical_mesh2d_eta1eta2%eta1_min, &
!VG!      logical_mesh2d_eta1eta2%eta1_max, &
!VG!      logical_mesh2d_eta1eta2%eta2_min, &
!VG!      logical_mesh2d_eta1eta2%eta2_max, &
!VG!      bound_cond%left_eta1, &
!VG!      bound_cond%right_eta1, &
!VG!      bound_cond%left_eta2, &
!VG!      bound_cond%right_eta2, &
!VG!      spline_degree%eta1, &
!VG!      spline_degree%eta2 )
!VG!
!VG!    call vlasov%interp1d_f_eta3%initialize( &
!VG!       Neta3, &
!VG!       logical_mesh1d_eta3%eta_min, &
!VG!       logical_mesh1d_eta3%eta_max, &
!VG!       bound_cond%left_eta3, &
!VG!       bound_cond%right_eta3, &
!VG!       spline_degree%eta3 )
!VG!
!VG!    call vlasov%interp1d_f_vpar%initialize( &
!VG!       Nvpar, &
!VG!       logical_mesh1d_vpar%eta_min, &
!VG!       logical_mesh1d_vpar%eta_max, &
!VG!       bound_cond%left_vpar, &
!VG!       bound_cond%right_vpar, &
!VG!       spline_degree%vpar )

  end subroutine new_vlasov_DK4D_hybrid
  !---------------------------------------------------------------------------


  !===========================================================================
  !> Vlasov: Deallocation 
  !---------------------------------------------------------------------------
  subroutine delete_vlasov_DK4D_hybrid( vlasov )

    type(vlasov_DK4D_hybrid_t), intent(inout) :: vlasov

    call vlasov%interp2d_f_eta1eta2%delete( )
    call sll_o_delete( vlasov%interp1d_f_eta3 )
    call sll_o_delete( vlasov%interp1d_f_vpar )

  end subroutine delete_vlasov_DK4D_hybrid
  !---------------------------------------------------------------------------


  !===========================================================================
  !> Vlasov: Global iteration of the drift-kinetic 4D simulation
  !---------------------------------------------------------------------------
  subroutine iteration_vlasov_DK4D_hybrid( &
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

    type(vlasov_DK4D_hybrid_t) , intent(inout) :: vlasov
    type(fdistribu_DK4D_t)     , intent(inout) :: fdistribu
    type(mesh_DK4D_t)          , intent(in)    :: mesh4d
    type(magnetconf_DK4D_t)    , intent(in)    :: magnetconf
    type(elecpot_DK4D_hybrid_t), intent(in)    :: elecpot
    sll_real64                 , intent(in)    :: dt
    sll_real64                 , intent(in)    :: advec_vpar_step1_coef
    sll_real64                 , intent(in)    :: advec_eta3_step1_coef
    sll_real64                 , intent(in)    :: advec_eta1eta2_step1_coef
    sll_real64       , optional, intent(in)    :: advec_vpar_step2_coef
    sll_real64       , optional, intent(in)    :: advec_eta3_step2_coef

    !--> local variables
    sll_int32 :: my_rank

    my_rank = sll_f_get_collective_rank( sll_v_world_collective )

    !*** Advection in vpar direction ***
    if ( advec_vpar_step1_coef.ne.0._F64 ) then
      if ( my_rank == 0) print*,' ==> first vpar advection with dt = ', &
          advec_vpar_step1_coef
      call advec1D_vpar_hybrid( &
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
      call advec1D_eta3_hybrid( &
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
      call advec2D_eta1eta2_hybrid( vlasov , &
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
        call advec1D_eta3_hybrid( &
            vlasov, &
            fdistribu, &
            mesh4d, &
            advec_eta3_step1_coef*dt )
      end if
    end if

    !*** Advection in vpar direction ***
    if ( present(advec_vpar_step2_coef) ) then
      if ( advec_vpar_step2_coef.ne.0._F64 ) then
        if ( my_rank == 0) print*,' ==> second vpar advection with dt = ', &
            advec_vpar_step2_coef
        call advec1D_vpar_hybrid( &
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

  end subroutine iteration_vlasov_DK4D_hybrid
  !---------------------------------------------------------------------------


  !===========================================================================
  ! Vlasov: 1D advection in vpar direction
  !---------------------------------------------------------------------------
  subroutine advec1D_vpar_hybrid( vlasov , &
      fdistribu, &
      mesh4d, &
      elecpot, &
      deltat_advec )

    type(vlasov_DK4D_hybrid_t) , intent(inout) :: vlasov
    type(fdistribu_DK4D_t)     , intent(inout) :: fdistribu
    type(mesh_DK4D_t)          , intent(in)    :: mesh4d
    type(elecpot_DK4D_hybrid_t), intent(in)    :: elecpot
    sll_real64                 , intent(in)    :: deltat_advec

    sll_int32  :: ierr
    sll_int32  :: iloc1, iloc2
    sll_int32  :: ieta3, ivpar
    sll_int32  :: loc4d_sz_x1, loc4d_sz_x2
    sll_int32  :: loc4d_sz_x3, loc4d_sz_x4
    sll_int32  :: Neta3, Nvpar
    sll_real64 :: vpar_point, vpar_min, vpar_max, E_eta3, alpha4
    sll_real64, dimension(:), pointer :: f1d_vpar_tmp

    !PN elec_field_eta1_3d is not a member of vlasov
    !SLL_ASSERT( size(fdistribu%val4d_seqx3x4,1).eq.size(vlasov%elec_field_eta1_3d%val3d_seqx3,1) )
    !SLL_ASSERT( size(fdistribu%val4d_seqx3x4,2).eq.size(vlasov%elec_field_eta1_3d%val3d_seqx3,2) )
    !SLL_ASSERT( size(fdistribu%val4d_seqx3x4,3).eq.size(vlasov%elec_field_eta1_3d%val3d_seqx3,3) )

    Neta3    = size( mesh4d%eta3_grid )
    Nvpar    = size( mesh4d%vpar_grid )
    vpar_min = minval( mesh4d%vpar_grid )
    vpar_max = maxval( mesh4d%vpar_grid )

    SLL_ALLOCATE( f1d_vpar_tmp(Nvpar), ierr )

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
          call vlasov%interp1d_f_vpar%compute_interpolants(f1d_vpar_tmp)   
          do ivpar = 1,Nvpar
            ! si change of coordinates in 3D
            ! val_jac  = sim%transf_xy%jacobian(eta1,eta2,eta3) 
            alpha4     = deltat_advec*E_eta3 
            vpar_point = mesh4d%vpar_grid(ivpar) - alpha4
            vpar_point = max(min(vpar_point,vpar_max),vpar_min)

            fdistribu%val4d_seqx3x4(iloc1,iloc2,ieta3,ivpar) = &
              vlasov%interp1d_f_vpar%interpolate_from_interpolant_value(vpar_point)
          end do
        end do
      end do
    end do

    SLL_DEALLOCATE(f1d_vpar_tmp,ierr)

  end subroutine advec1D_vpar_hybrid
  !---------------------------------------------------------------------------


  !===========================================================================
  ! Vlasov: 1D advection in eta3 direction
  !---------------------------------------------------------------------------
  subroutine advec1D_eta3_hybrid( vlasov , &
      fdistribu, &
      mesh4d, &
      deltat_advec )

    type(vlasov_DK4D_hybrid_t), intent(inout) :: vlasov
    type(fdistribu_DK4D_t)    , intent(inout) :: fdistribu
    type(mesh_DK4D_t)         , intent(in)    :: mesh4d
    sll_real64                , intent(in)    :: deltat_advec

    sll_int32  :: ierr
    sll_int32  :: iloc1, iloc2
    sll_int32  :: ieta3, ivpar
    sll_int32  :: loc4d_sz_x1, loc4d_sz_x2
    sll_int32  :: loc4d_sz_x3, loc4d_sz_x4
    sll_int32  :: Neta3, Nvpar
    sll_real64 :: eta3_point, alpha3, vpar_point
    sll_real64, dimension(:), pointer :: f1d_eta3_tmp

    Neta3 = size(mesh4d%eta3_grid)
    Nvpar = size(mesh4d%vpar_grid)

    call sll_o_compute_local_sizes( fdistribu%layout4d_seqx3x4, &
        loc4d_sz_x1, &
        loc4d_sz_x2, &
        loc4d_sz_x3, &
        loc4d_sz_x4 )    

    !---> dphi/dt = vpar
    SLL_ALLOCATE(f1d_eta3_tmp(Neta3),ierr)
    do iloc2 = 1,loc4d_sz_x2
      do iloc1 = 1,loc4d_sz_x1
        do ivpar = 1,Nvpar
          vpar_point = mesh4d%vpar_grid(ivpar)
          do ieta3 = 1,Neta3
            f1d_eta3_tmp(ieta3) = fdistribu%val4d_seqx3x4(iloc1,iloc2,ieta3,ivpar)
          end do
          call vlasov%interp1d_f_eta3%compute_interpolants(f1d_eta3_tmp)  

          do ieta3 = 1,Neta3
            alpha3     = deltat_advec*vpar_point
            eta3_point = mesh4d%eta3_grid(ieta3) - alpha3

            fdistribu%val4d_seqx3x4(iloc1,iloc2,ieta3,ivpar) = &
                vlasov%interp1d_f_eta3%interpolate_from_interpolant_value(eta3_point)
          end do
        end do
      end do
    end do

    SLL_DEALLOCATE(f1d_eta3_tmp,ierr)

  end subroutine advec1D_eta3_hybrid
  !---------------------------------------------------------------------------


  !===========================================================================
  ! Vlasov: 2D advection in (eta1,eta2) directions
  !---------------------------------------------------------------------------
  subroutine advec2D_eta1eta2_hybrid( &
      vlasov, &
      fdistribu, &
      mesh4d, &
      magnetconf, &
      elecpot, &
      deltat_advec )

    type(vlasov_DK4D_hybrid_t) , intent(inout) :: vlasov
    type(fdistribu_DK4D_t)     , intent(inout) :: fdistribu
    type(mesh_DK4D_t)          , intent(in)    :: mesh4d
    type(magnetconf_DK4D_t)    , intent(in)    :: magnetconf
    type(elecpot_DK4D_hybrid_t), intent(in)    :: elecpot
    sll_real64                 , intent(in)    :: deltat_advec

    sll_int32  :: ierr
    sll_int32  :: iloc3, iloc4
    sll_int32  :: ieta1, ieta2
    sll_int32  :: Neta1, Neta2
    sll_int32  :: loc4d_sz_x1, loc4d_sz_x2
    sll_int32  :: loc4d_sz_x3, loc4d_sz_x4
    sll_real64 :: eta1_point, eta2_point 
    sll_real64 :: alpha1, alpha2, E_eta1, E_eta2
    sll_real64 :: val_jac,val_B
    sll_real64, dimension(:,:), pointer :: f2d_eta1eta2_tmp

    Neta1 = size( mesh4d%eta1_grid )
    Neta2 = size( mesh4d%eta2_grid )

    call sll_o_compute_local_sizes( fdistribu%layout4d_seqx1x2, &
        loc4d_sz_x1, &
        loc4d_sz_x2, &
        loc4d_sz_x3, &
        loc4d_sz_x4 )    

    !---> ( deta1/dt )                              ( 0 -1 ) (  d phi/deta1 ) 
    !     (          )=  1 / (B* jac(eta_1,eta2) )  (      ) (              )
    !---> ( deta2/dt )                              ( 1  0 ) (  d phi/deta2 )
    SLL_ALLOCATE( f2d_eta1eta2_tmp(Neta1,Neta2), ierr )

    do iloc4 = 1,loc4d_sz_x4
      do iloc3 = 1,loc4d_sz_x3
        do ieta2 = 1,Neta2
          do ieta1 = 1,Neta1
            f2d_eta1eta2_tmp(ieta1,ieta2) = &
                fdistribu%val4d_seqx1x2(ieta1,ieta2,iloc3,iloc4)            
          end do
        end do

        call vlasov%interp2d_f_eta1eta2%compute_interpolants( f2d_eta1eta2_tmp )

        do ieta2 = 1,Neta2
          do ieta1 = 1,Neta1
            E_eta1  = - elecpot%dPhi_deta1%val3d_seqx1x2(ieta1,ieta2,iloc3)
            E_eta2  = - elecpot%dPhi_deta2%val3d_seqx1x2(ieta1,ieta2,iloc3)
            val_jac = mesh4d%transf_eta1eta2_xy%jacobian(eta1_point,eta2_point)
            val_B   = magnetconf%B_xy(ieta1,ieta2)
            alpha1  =  deltat_advec*E_eta2 / (val_jac*val_B)
            alpha2  =  -deltat_advec*E_eta1 / (val_jac*val_B)
            eta1_point = mesh4d%eta1_grid(ieta1) - alpha1
            eta2_point = mesh4d%eta2_grid(ieta2) - alpha2

            if (eta1_point <= mesh4d%eta1_eta2_mesh2d%eta1_min) then 
              eta1_point = mesh4d%eta1_eta2_mesh2d%eta1_min
            else if ( eta1_point>= mesh4d%eta1_eta2_mesh2d%eta1_max) then 
              eta1_point = mesh4d%eta1_eta2_mesh2d%eta1_max
            end if

            fdistribu%val4d_seqx1x2(ieta1,ieta2,iloc3,iloc4) = &
                vlasov%interp2d_f_eta1eta2%interpolate_from_interpolant_value( eta1_point,eta2_point )

          end do
        end do
      end do
    end do

    SLL_DEALLOCATE(f2d_eta1eta2_tmp,ierr)

  end subroutine advec2D_eta1eta2_hybrid
  !---------------------------------------------------------------------------

end module vlasov_DK4D_hybrid_module
!---------------------------------------------------------------------------
