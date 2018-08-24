!===========================================================================
!> Solving of the vlasov equation for
!>  4D Vlasov-Poisson hybrid simulation
!>
!> \date 2015-03-11
!> \author V. Grandgirard, A. Back
!---------------------------------------------------------------------------
module vlasov_VP4D_module
#include "sll_working_precision.h"
#include "sll_memory.h"
#include "sll_assert.h"

  use elecpot_Bsplines_VP4D_module
  use fdistribu_VP4D_module
  use field2d_VP4D_module
  use mesh_VP4D_module
  use sll_m_arbitrary_degree_spline_interpolator_1d
  use sll_m_cartesian_meshes
  use sll_m_scalar_field_1d
  use sll_m_scalar_field_2d


  implicit none

contains

  !===========================================================================
  !> Vlasov: Global iteration of the drift-kinetic 4D simulation
  !---------------------------------------------------------------------------
  subroutine iteration_vlasov_Bsplines_VP4D( &
      fdistribu, &
      mesh4d, &
      elecpot, &
      dt, &
      advec2D_scheme, &
      advec_vy_step1_coef, &
      advec_vx_step1_coef, &
      advec_eta1eta2_step1_coef, &
      advec_vy_step2_coef, &
      advec_vx_step2_coef )

    type(fdistribu_VP4D_t)       , intent(inout) :: fdistribu
    type(mesh_VP4D_t)            , intent(in)    :: mesh4d
    type(elecpot_Bsplines_VP4D_t), intent(in)    :: elecpot
    sll_real64                   , intent(in)    :: dt
    sll_int32                    , intent(in)    :: advec2D_scheme
    sll_real64                   , intent(in)    :: advec_vy_step1_coef
    sll_real64                   , intent(in)    :: advec_vx_step1_coef
    sll_real64                   , intent(in)    :: advec_eta1eta2_step1_coef
    sll_real64         , optional, intent(in)    :: advec_vy_step2_coef
    sll_real64         , optional, intent(in)    :: advec_vx_step2_coef

    !*** Advection in vy direction ***
    if ( advec_vy_step1_coef.ne.0._F64 ) then
      call advec1D_vy( &
          fdistribu, &
          mesh4d, &
          elecpot, &
          advec_vy_step1_coef*dt )
    end if

    !*** Advection in vx direction ***
    if ( advec_vx_step1_coef.ne.0._F64 ) then
      call advec1D_vx( &
          fdistribu, &
          mesh4d, &
          elecpot, &
          advec_vx_step1_coef*dt )
    end if

    !*** Advection in (eta1,eta2) directions ***
    call sll_o_apply_remap_4d( &
        fdistribu%seqx3x4_to_seqx1x2, &
        fdistribu%val4d_seqx3x4,&
        fdistribu%val4d_seqx1x2 )
    if ( advec_eta1eta2_step1_coef.ne.0._F64 ) then
      call advec2D_eta1eta2( &
          fdistribu, &
          mesh4d, &
          advec_eta1eta2_step1_coef*dt, &
          advec2D_scheme )
      call sll_o_apply_remap_4d( &
          fdistribu%seqx1x2_to_seqx3x4, &
          fdistribu%val4d_seqx1x2, &
          fdistribu%val4d_seqx3x4 )
    end if

    !*** Advection in vx direction ***
    if ( present(advec_vx_step2_coef) ) then
      call advec1D_vx( &
          fdistribu, &
          mesh4d, &
          elecpot, &
          advec_vx_step1_coef*dt )
    end if

    !*** Advection in vy direction ***
    if ( present(advec_vy_step2_coef) ) then
      call advec1D_vy( &
          fdistribu, &
          mesh4d, &
          elecpot, &
          advec_vy_step2_coef*dt)
    end if

    !*** Fill fdistribu%val4d_seqx1x2 (to be able to solve QN)
    if ( present(advec_vx_step2_coef) .or. &
        present(advec_vy_step2_coef) ) then
      if ( (advec_vy_step2_coef+advec_vx_step2_coef).ne.0._F64 ) then
        call sll_o_apply_remap_4d( &
            fdistribu%seqx3x4_to_seqx1x2, &
            fdistribu%val4d_seqx3x4,&
            fdistribu%val4d_seqx1x2 )
      end if
    end if

  end subroutine iteration_vlasov_Bsplines_VP4D
  !---------------------------------------------------------------------------


  !===========================================================================
  ! Vlasov: 1D advection in vy direction
  !---------------------------------------------------------------------------
  subroutine advec1D_vy( &
      fdistribu, &
      mesh4d, &
      elecpot, &
      deltat_advec )

    type(fdistribu_VP4D_t)       , intent(inout) :: fdistribu
    type(mesh_VP4D_t)            , intent(in)    :: mesh4d
    type(elecpot_Bsplines_VP4D_t), intent(in)    :: elecpot
    sll_real64                   , intent(in)    :: deltat_advec

    sll_int32  :: ierr
    sll_int32  :: iloc1, iloc2, ieta1, ieta2
    sll_int32  :: ivx, ivy
    sll_int32  :: Nvx, Nvy
    sll_int32  :: loc4d_sz_x1, loc4d_sz_x2
    sll_int32  :: loc4d_sz_x3, loc4d_sz_x4
    logical    :: outofdomain_vy
    sll_real64 :: Ex_point, Ey_point
    sll_real64 :: vy_point, alpha4
    sll_real64, dimension(:), pointer :: f1d_vy_tmp
    sll_real64, dimension(2,2)        :: inv_JacobMat
    sll_int32 , dimension(4)          :: glob_ind4d

    Nvx    = size( mesh4d%vx_grid )
    Nvy    = size( mesh4d%vy_grid )

    SLL_ALLOCATE( f1d_vy_tmp(Nvy), ierr )

    call sll_o_compute_local_sizes( fdistribu%layout4d_seqx3x4, &
      loc4d_sz_x1, &
      loc4d_sz_x2, &
      loc4d_sz_x3, &
      loc4d_sz_x4 )    

    do iloc2 = 1,loc4d_sz_x2
      do iloc1 = 1,loc4d_sz_x1
        do ivx = 1,Nvx
          do ivy = 1,Nvy
            f1d_vy_tmp(ivy) = &
                fdistribu%val4d_seqx3x4(iloc1,iloc2,ivx,ivy)
          end do
          call fdistribu%interp1d_vy%compute_interpolants(f1d_vy_tmp)   

          do ivy = 1,Nvy
            glob_ind4d(:) = sll_o_local_to_global( &
                fdistribu%layout4d_seqx3x4, &
                (/iloc1,iloc2,ivx,ivy/))
            ieta1        = glob_ind4d(1)
            ieta2        = glob_ind4d(2)
            inv_JacobMat = mesh4d%inv_Jacobian_matrix(ieta1,ieta2,:,:)

            Ex_point = -elecpot%dPhi_deta1%val2d_seqx1x2(ieta1,ieta2)
            Ey_point = -elecpot%dPhi_deta2%val2d_seqx1x2(ieta1,ieta2)

            !\todo: Understand what is done by Aurore
            !\todo: Check the displacement sign
            !\todo: BECAREFUL not the same sign than Aurore
            alpha4     = deltat_advec * &
                ( inv_JacobMat(1,2)*Ex_point + &
                inv_JacobMat(2,2)*Ey_point ) 
            vy_point = mesh4d%vy_grid(ivy) - alpha4

            !*** Treatment of boundary conditions ***
            outofdomain_vy = check_outofdomain_(vy_point, &
                mesh4d%vy_min, mesh4d%vy_max)
            if ( outofdomain_vy ) then
              call modify_position_( vy_point, &
                  mesh4d%vy_min, mesh4d%vy_max, mesh4d%vy_length, &
                  fdistribu%bound_cond%left_vy, &
                  fdistribu%bound_cond%right_vy )
            end if            

            fdistribu%val4d_seqx3x4(iloc1,iloc2,ivx,ivy) = &
              fdistribu%interp1d_vy%interpolate_from_interpolant_value(vy_point)
          end do
        end do
      end do
    end do

    SLL_DEALLOCATE(f1d_vy_tmp,ierr)

  end subroutine advec1D_vy
  !---------------------------------------------------------------------------


  !===========================================================================
  ! Vlasov: 1D advection in vx direction
  !---------------------------------------------------------------------------
  subroutine advec1D_vx( &
      fdistribu, &
      mesh4d, &
      elecpot, &
      deltat_advec )

    type(fdistribu_VP4D_t)       , intent(inout) :: fdistribu
    type(mesh_VP4D_t)            , intent(in)    :: mesh4d
    type(elecpot_Bsplines_VP4D_t), intent(in)    :: elecpot
    sll_real64                   , intent(in)    :: deltat_advec

    sll_int32  :: ierr
    sll_int32  :: iloc1, iloc2, ieta1, ieta2
    sll_int32  :: ivx, ivy
    sll_int32  :: loc4d_sz_x1, loc4d_sz_x2
    sll_int32  :: loc4d_sz_x3, loc4d_sz_x4
    sll_int32  :: Nvx, Nvy
    logical    :: outofdomain_vx
    sll_real64 :: Ex_point, Ey_point
    sll_real64 :: vx_point, alpha3
    sll_real64, dimension(:), pointer :: f1d_vx_tmp
    sll_real64, dimension(2,2)        :: inv_JacobMat
    sll_int32 , dimension(4)          :: glob_ind4d

    Nvx = size(mesh4d%vx_grid)
    Nvy = size(mesh4d%vy_grid)

    SLL_ALLOCATE(f1d_vx_tmp(Nvx),ierr)

    call sll_o_compute_local_sizes( fdistribu%layout4d_seqx3x4, &
        loc4d_sz_x1, &
        loc4d_sz_x2, &
        loc4d_sz_x3, &
        loc4d_sz_x4 )    

    do iloc2 = 1,loc4d_sz_x2
      do iloc1 = 1,loc4d_sz_x1
        do ivy = 1,Nvy
          do ivx = 1,Nvx
            f1d_vx_tmp(ivx) = fdistribu%val4d_seqx3x4(iloc1,iloc2,ivx,ivy)
          end do
          call fdistribu%interp1d_vx%compute_interpolants(f1d_vx_tmp)  

          do ivx = 1,Nvx
            glob_ind4d(:) = sll_o_local_to_global( &
                fdistribu%layout4d_seqx3x4, &
                (/iloc1,iloc2,ivx,ivy/))
            ieta1        = glob_ind4d(1)
            ieta2        = glob_ind4d(2)
            inv_JacobMat = mesh4d%inv_Jacobian_matrix(ieta1,ieta2,:,:)

            Ex_point = -elecpot%dPhi_deta1%val2d_seqx1x2(ieta1,ieta2)
            Ey_point = -elecpot%dPhi_deta2%val2d_seqx1x2(ieta1,ieta2)

            !\todo: Understand what is done by Aurore
            !\todo: Check the displacement sign
            !\todo: BECAREFUL not the same sign than Aurore
            alpha3   = deltat_advec * &
                ( inv_JacobMat(1,1)*Ex_point + &
                inv_JacobMat(2,1)*Ey_point )
            vx_point = mesh4d%vx_grid(ivx) - alpha3

            !*** Treatment of boundary conditions ***
            outofdomain_vx = check_outofdomain_(vx_point, &
                mesh4d%vx_min, mesh4d%vx_max)
            if ( outofdomain_vx ) then
              call modify_position_( vx_point, &
                  mesh4d%vx_min, mesh4d%vx_max, mesh4d%vx_length, &
                  fdistribu%bound_cond%left_vx, &
                  fdistribu%bound_cond%right_vx )
            end if            

            fdistribu%val4d_seqx3x4(iloc1,iloc2,ivx,ivy) = &
                fdistribu%interp1d_vx%interpolate_from_interpolant_value(vx_point)
          end do
        end do
      end do
    end do

    SLL_DEALLOCATE(f1d_vx_tmp,ierr)

  end subroutine advec1D_vx
  !---------------------------------------------------------------------------


  !===========================================================================
  ! Vlasov: 2D advection in (eta1,eta2) directions
  !---------------------------------------------------------------------------
  subroutine advec2D_eta1eta2( &
      fdistribu, &
      mesh4d, &
      deltat_advec, &
      advec2D_scheme )

    type(fdistribu_VP4D_t), intent(inout) :: fdistribu
    type(mesh_VP4D_t)     , intent(in)    :: mesh4d
    sll_real64            , intent(in)    :: deltat_advec
    sll_int32             , intent(in)    :: advec2D_scheme

    sll_int32  :: ierr
    sll_int32  :: iloc3, iloc4, ivx, ivy
    sll_int32  :: ieta1, ieta2
    sll_int32  :: Neta1, Neta2
    sll_int32  :: loc4d_sz_x1, loc4d_sz_x2
    sll_int32  :: loc4d_sz_x3, loc4d_sz_x4
    sll_real64 :: eta1_point, eta2_point
    sll_real64 :: vx_point, vy_point
    sll_real64 :: alpha1, alpha2
    !logical    :: outofdomain_eta1, outofdomain_eta2
    sll_real64, dimension(:,:), pointer :: f2d_eta1eta2_tmp
    sll_int32 , dimension(4)            :: glob_ind4d
    sll_real64, dimension(2,2)          :: inv_JacobMat
    !---> For Runge-Kutta scheme
    sll_real64 :: eta1_RK2, eta2_RK2
    sll_real64 :: alpha1_RK2, alpha2_RK2 
    sll_real64, dimension(2,2) :: inv_JacobMat_RK2

    Neta1 = size( mesh4d%eta1_grid )
    Neta2 = size( mesh4d%eta2_grid )

    call sll_o_compute_local_sizes( fdistribu%layout4d_seqx1x2, &
        loc4d_sz_x1, &
        loc4d_sz_x2, &
        loc4d_sz_x3, &
        loc4d_sz_x4 )    

    !---> ( deta1/dt )=  J^{-1} ( vx )
    !---> ( deta2/dt )          ( vy )
    SLL_ALLOCATE( f2d_eta1eta2_tmp(Neta1,Neta2), ierr )

    do iloc4 = 1,loc4d_sz_x4
      do iloc3 = 1,loc4d_sz_x3
        do ieta2 = 1,Neta2
          do ieta1 = 1,Neta1
            f2d_eta1eta2_tmp(ieta1,ieta2)=&
                fdistribu%val4d_seqx1x2(ieta1,ieta2,iloc3,iloc4)            
          end do
        end do

        call fdistribu%interp2d_eta1eta2%compute_interpolants( f2d_eta1eta2_tmp )

        do ieta2 = 1,Neta2
          do ieta1 = 1,Neta1
            glob_ind4d(:) = sll_o_local_to_global( &
                fdistribu%layout4d_seqx1x2, &
                (/ieta1,ieta2,iloc3,iloc4/))
            ivx      = glob_ind4d(3) 
            ivy      = glob_ind4d(4)

            eta1_point = mesh4d%eta1_grid(ieta1)
            eta2_point = mesh4d%eta2_grid(ieta2)
            vx_point   = mesh4d%vx_grid(ivx)
            vy_point   = mesh4d%vy_grid(ivy)

            inv_JacobMat = mesh4d%inv_Jacobian_matrix(ieta1,ieta2,:,:)

            if (advec2D_scheme .eq. SLL_EULER) then 
              !*** Euler explicit scheme ***
              alpha1 = deltat_advec * &
                  ( inv_JacobMat(1,1)*vx_point + &
                  inv_JacobMat(1,2)*vy_point )
              alpha2 = deltat_advec * &
                  (inv_JacobMat(2,1)*vx_point + &
                  inv_JacobMat(2,2)*vy_point )

            elseif ( advec2D_scheme .eq. SLL_RUNGEKUTTA) then
              !*** Runge-Kutta of second order ***
              alpha1_RK2 = deltat_advec * &
                  ( inv_JacobMat(1,1)*vx_point + &
                  inv_JacobMat(1,2)*vy_point)
              alpha2_RK2 = deltat_advec * &
                  ( inv_JacobMat(2,1)*vx_point + &
                  inv_JacobMat(2,2)*vy_point)
              eta1_RK2 = eta1_point - alpha1_RK2
              eta2_RK2 = eta2_point - alpha2_RK2

              !*** Check eta1 and eta2 position according to ***
              !***  boundary conditions                      ***
              call check_eta1eta2_( mesh4d, fdistribu, &
                  eta1_RK2, eta2_RK2 )

              inv_JacobMat_RK2(:,:) = &
                  mesh4d%transf_eta1eta2_xy%inverse_jacobian_matrix( &
                  eta1_RK2, eta2_RK2 )
              alpha1 = 0.5_f64 * deltat_advec * &
                  ( ( inv_JacobMat(1,1) + inv_JacobMat_RK2(1,1) ) *vx_point &
                  + ( inv_JacobMat(1,2) + inv_JacobMat_RK2(1,2) ) *vy_point )
              alpha2 = 0.5_f64 * deltat_advec * &
                  ( ( inv_JacobMat(2,1) + inv_JacobMat_RK2(2,1) ) *vx_point &
                  + ( inv_JacobMat(2,2) + inv_JacobMat_RK2(2,2) ) *vy_point )
            end if

            eta1_point = eta1_point - alpha1
            eta2_point = eta2_point - alpha2

            !*** Check eta1 and eta2 position according to ***
            !***  boundary conditions                      ***
            call check_eta1eta2_( mesh4d, fdistribu, &
                eta1_point, eta2_point )

            !*** Computation of the interpolation value ***
            fdistribu%val4d_seqx1x2(ieta1,ieta2,iloc3,iloc4) = &
                fdistribu%interp2d_eta1eta2%interpolate_from_interpolant_value( &
                eta1_point,eta2_point )
          end do
        end do
      end do
    end do

    SLL_DEALLOCATE(f2d_eta1eta2_tmp,ierr)

  end subroutine advec2D_eta1eta2
  !---------------------------------------------------------------------------


  !===========================================================================
  !> Check if the particle is out of the domain
  !---------------------------------------------------------------------------
  subroutine check_eta1eta2_( mesh4d, fdistribu, &
      eta1_point, eta2_point )

    type(fdistribu_VP4D_t), intent(in)    :: fdistribu
    type(mesh_VP4D_t)     , intent(in)    :: mesh4d
    sll_real64            , intent(inout) :: eta1_point
    sll_real64            , intent(inout) :: eta2_point

    logical :: outofdomain_eta1, outofdomain_eta2

    outofdomain_eta1 = check_outofdomain_(eta1_point, &
        mesh4d%eta1_min, mesh4d%eta1_max)
    outofdomain_eta2 = check_outofdomain_(eta2_point, &
        mesh4d%eta2_min, mesh4d%eta2_max)

    if ( outofdomain_eta1 ) then
      call modify_position_( eta1_point, &
          mesh4d%eta1_min, mesh4d%eta1_max, mesh4d%eta1_length, &
          fdistribu%bound_cond%left_eta1, &
          fdistribu%bound_cond%right_eta1 )
    end if
    if ( outofdomain_eta2 ) then
      call modify_position_( eta2_point, &
          mesh4d%eta2_min, mesh4d%eta2_max, mesh4d%eta2_length, &
          fdistribu%bound_cond%left_eta2, &
          fdistribu%bound_cond%right_eta2 )
    end if
    
  end subroutine check_eta1eta2_
  !---------------------------------------------------------------------------


  !===========================================================================
  !> Check if the particle is out of the domain
  !---------------------------------------------------------------------------
  function check_outofdomain_( x_point, &
      x_min, x_max ) result(out_of_domain)

    sll_real64, intent(in) :: x_point
    sll_real64, intent(in) :: x_min
    sll_real64, intent(in) :: x_max

    logical :: out_of_domain

    out_of_domain = .false.
    if ( (x_point<x_min) .or. (x_point>x_max) ) then
      out_of_domain = .true.
    end if

  end function check_outofdomain_
  !---------------------------------------------------------------------------


  !===========================================================================
  !> Modify the position of a particle out of domain depending on the
  !>  boundary conditions.
  !>  
  !> Rk: If Dirichlet boundary conditions the particle are fixed at the
  !>     boundary of the domain (this assumes that the distribution function
  !>     is initialized to 0. at the boundary)
  !---------------------------------------------------------------------------
  subroutine modify_position_( x_point, &
      x_min, x_max, x_length, &
      x_bound_cond_left, x_bound_cond_right )

    sll_real64, intent(inout) :: x_point
    sll_real64, intent(in)    :: x_min
    sll_real64, intent(in)    :: x_max
    sll_real64, intent(in)    :: x_length
    sll_int32,  intent(in)    :: x_bound_cond_left
    sll_int32,  intent(in)    :: x_bound_cond_right

    if ( x_point < x_min ) then
      select case ( x_bound_cond_left )
      case ( sll_p_periodic )
        x_point = x_max - mod(x_min-x_point,x_length)
      case ( sll_p_dirichlet )
        x_point = x_min
      case default
        print*, "Left boundary condition not treated"
        stop
      end select
    end if

    if ( x_point > x_max ) then
      select case ( x_bound_cond_right )
      case ( sll_p_periodic )
        x_point = x_min + mod(x_point-x_max,x_length)
      case ( sll_p_dirichlet )
        x_point = x_max
      case default
        print*, "Right boundary condition not treated"
        stop
      end select
    end if

  end subroutine modify_position_
  !---------------------------------------------------------------------------

end module vlasov_VP4D_module
!---------------------------------------------------------------------------
