!> \file sll_lin_trans_pic_2d.F90
!> \namespace sll_lin_trans_pic_2d
!> \brief  
!> The ltpic module provides capabilities for the linearly-transformed pic method. More details to come (this version is not functional yet).
!>
module sll_lin_trans_pic_2d
#include "sll_working_precision.h"
#include "sll_memory.h"
#include "sll_assert.h"
  use sll_flow_base
  
  implicit none
  
  type  ::  lin_trans_pic_2d
     sll_int32                         :: bspline_degree
     sll_int32                         :: num_particles_x
     sll_int32                         :: num_particles_v
     sll_real64                        :: xmin
     sll_real64                        :: xmax
     sll_real64                        :: vmin
     sll_real64                        :: vmax
     sll_int32                         :: qi_stencil_radius ! the quasi-interpolation stencil depends on the degree of the bsplines
     sll_int32                         :: qi_grid_npx
     sll_int32                         :: qi_grid_npv
     sll_real64                        :: qi_grid_xmin
     sll_real64                        :: qi_grid_xmax
     sll_real64                        :: qi_grid_vmin
     sll_real64                        :: qi_grid_vmax
     
     sll_real64, dimension(:,:), pointer :: weights
     sll_real64, dimension(:,:), pointer :: coord_x
     sll_real64, dimension(:,:), pointer :: coord_v
     sll_real64, dimension(:,:), pointer :: deform_matrix_xx
     sll_real64, dimension(:,:), pointer :: deform_matrix_xv
     sll_real64, dimension(:,:), pointer :: deform_matrix_vx
     sll_real64, dimension(:,:), pointer :: deform_matrix_vv

     sll_real64, dimension(:,:), pointer :: qi_coefs
     sll_real64, dimension(:,:), pointer :: data ! e.g., the data to be approximated by the particles
  end type sll_spline_2D

contains


  function new_ltpic_2d( &
            num_particles_x,
            num_particles_v,
            xmin, 
            xmax,
            vmin,
            vmax,
            bspline_degree )  

    type(lin_trans_pic_2d), pointer         :: new_ltpic_2d
    sll_int32,  intent(in)              :: num_particles_x
    sll_int32,  intent(in)              :: num_particles_v        
    sll_real64, intent(in)              :: xmin
    sll_real64, intent(in)              :: xmax
    sll_real64, intent(in)              :: vmin
    sll_real64, intent(in)              :: vmax
    sll_int32,  intent(in)              :: bspline_degree
    sll_real64, dimension(:), pointer   :: qi_coefs_1d
    sll_real64                          :: qi_st_radius
    sll_real64                          :: grid_dx
    sll_real64                          :: grid_dv
    sll_int32                           :: i
    sll_int32                           :: j
                
    SLL_ALLOCATE( new_ltpic_2d, ierr )
    new_ltpic_2d%num_particles_x = num_particles_x
    new_ltpic_2d%num_particles_v = num_particles_v
    new_ltpic_2d%xmin = xmin
    new_ltpic_2d%xmax = xmax
    new_ltpic_2d%vmin = vmin
    new_ltpic_2d%vmax = vmax
    
    grid_dx = (xmax-xmin)/num_particles_x
    grid_dv = (vmax-vmin)/num_particles_v
    
    new_ltpic_2d%bspline_degree = bspline_degree   
    select case (bspline_degree)
    case (1)
       qi_st_radius = 0
       new_ltpic_2d%qi_stencil_radius = qi_st_radius
       SLL_ALLOCATE( qi_coefs_1d           (-qi_st_radius:qi_st_radius)                            ,    ierr )
       SLL_ALLOCATE( new_ltpic_2d%qi_coefs (-qi_st_radius:qi_st_radius, -qi_st_radius:qi_st_radius),    ierr )
       qi_coefs_1d(0) = 1
       do i = -qi_st_radius, qi_st_radius
          do j = -qi_st_radius, qi_st_radius
             new_ltpic_2d%qi_coefs(i,j) = qi_coefs_1d(i) * qi_coefs_1d(j)
          end do
       end do
    case (3)
       qi_st_radius = 1
       new_ltpic_2d%qi_stencil_radius = qi_st_radius
       SLL_ALLOCATE( qi_coefs_1d           (-qi_st_radius:qi_st_radius)                            ,    ierr )
       SLL_ALLOCATE( new_ltpic_2d%qi_coefs (-qi_st_radius:qi_st_radius, -qi_st_radius:qi_st_radius),    ierr )
       qi_coefs_1d(-1) = -1./6
       qi_coefs_1d( 0) =  8./6
       qi_coefs_1d( 1) = -1./6
       do i = -qi_st_radius, qi_st_radius
          do j = -qi_st_radius, qi_st_radius
             new_ltpic_2d%qi_coefs(i,j) = qi_coefs_1d(i) * qi_coefs_1d(j)
          end do
       end do
    case (5)
       qi_st_radius = 4
       new_ltpic_2d%qi_stencil_radius = qi_st_radius
       SLL_ALLOCATE( qi_coefs_1d           (-qi_st_radius:qi_st_radius)                            ,    ierr )
       SLL_ALLOCATE( new_ltpic_2d%qi_coefs (-qi_st_radius:qi_st_radius, -qi_st_radius:qi_st_radius),    ierr )
       qi_coefs_1d(-4) =    1./14400
       qi_coefs_1d(-3) =    13./3600
       qi_coefs_1d(-2) =      7./225
       qi_coefs_1d(-1) = -1469./3600
       qi_coefs_1d( 0) =    503./288
       qi_coefs_1d( 1) = -1469./3600
       qi_coefs_1d( 2) =      7./225
       qi_coefs_1d( 3) =    13./3600
       qi_coefs_1d( 4) =    1./14400
       do i = -qi_st_radius, qi_st_radius
          do j = -qi_st_radius, qi_st_radius
             new_ltpic_2d%qi_coefs(i,j) = qi_coefs_1d(i) * qi_coefs_1d(j)
          end do
       end do
    case default
       print *, 'ERROR: new_ltpic_2d(): not recognized bspline_degree -- value =', bspline_degree
       STOP
    end select
    new_ltpic_2d%qi_grid_npx  = num_particles_x + 2*new_ltpic_2d%qi_stencil_radius
    new_ltpic_2d%qi_grid_npv  = num_particles_v + 2*new_ltpic_2d%qi_stencil_radius
    new_ltpic_2d%qi_grid_xmin = xmin - new_ltpic_2d%qi_stencil_radius*grid_dx
    new_ltpic_2d%qi_grid_xmax = xmax + new_ltpic_2d%qi_stencil_radius*grid_dx
    new_ltpic_2d%qi_grid_vmin = vmin - new_ltpic_2d%qi_stencil_radius*grid_dv
    new_ltpic_2d%qi_grid_vmax = vmax + new_ltpic_2d%qi_stencil_radius*grid_dv
    
    SLL_ALLOCATE( new_ltpic_2d%weights         (num_particles_x,num_particles_v),   ierr )
    SLL_ALLOCATE( new_ltpic_2d%coord_x         (num_particles_x,num_particles_v),   ierr )
    SLL_ALLOCATE( new_ltpic_2d%coord_v         (num_particles_x,num_particles_v),   ierr )
    SLL_ALLOCATE( new_ltpic_2d%deform_matrix_xx(num_particles_x,num_particles_v),   ierr )
    SLL_ALLOCATE( new_ltpic_2d%deform_matrix_xv(num_particles_x,num_particles_v),   ierr )
    SLL_ALLOCATE( new_ltpic_2d%deform_matrix_vx(num_particles_x,num_particles_v),   ierr )
    SLL_ALLOCATE( new_ltpic_2d%deform_matrix_vv(num_particles_x,num_particles_v),   ierr )

    SLL_ALLOCATE( new_ltpic_2d%data            (new_ltpic_2d%qi_grid_npx,new_ltpic_2d%qi_grid_npv),   ierr )

  end function new_ltpic_2d



  subroutine load_ltpic_2d( target_density_initializer, ltpic_object )
    class(scalar_field_2d_initializer_base), intent(in), pointer :: target_density_initializer   ! intial field to be approximated
    type(lin_trans_pic_2d), pointer          :: ltpic_object
    sll_int32                            :: npx
    sll_int32                            :: npv
    sll_int32                            :: k_x
    sll_int32                            :: k_v
    sll_int32                            :: l_x
    sll_int32                            :: l_v
    sll_int32                            :: qi_st_radius
    sll_real64                           :: qi_weight
    sll_real64, dimension(:,:), pointer  :: qi_coefs
         
    call target_density_initializer%f_of_x1x2(ltpic_object%data)   ! writes the target density values on the 'data' array    
    
    npx = ltpic_object%num_particles_x
    npv = ltpic_object%num_particles_v    
    qi_st_radius = ltpic_object%qi_stencil_radius    
    qi_coefs => ltpic_object%qi_coefs

    ! initial particles have no deformation, hence set D_k = (1. 0. // 0. 1.) for every particle index k = (k_x,k_v)
    ltpic_object%deform_matrix_xx = 1.
    ltpic_object%deform_matrix_xv = 0.
    ltpic_object%deform_matrix_vx = 0.
    ltpic_object%deform_matrix_vv = 1.

    ! compute particle coordinates and weights
    do k_x = 1,npx 
       do k_v = 1,npv           
          ! due to the quasi-interpolation stencils, there are more nodes (in the target mesh) than particles -- hence the offset
          ltpic_object%coord_x(k_x,k_v) = target_density_initializer%mesh%x1_at_node( k_x + qi_st_radius, k_v + qi_st_radius )
          ltpic_object%coord_v(k_x,k_v) = target_density_initializer%mesh%x2_at_node( k_x + qi_st_radius, k_v + qi_st_radius )                    
          ! compute the weight with local quasi-interpolation scheme
          qi_weight = 0
          do l_x = -qi_st_radius, qi_st_radius
             do l_v = -qi_st_radius, qi_st_radius
                ! here the offset is for the same reason as above ('data' member has same number of nodes than the target mesh)
                qi_weight = qi_weight + qi_coefs(l_x,l_v) * ltpic_object%data( k_x+l_x + qi_st_radius, k_v+l_v + qi_st_radius )
             end do
          end do
          ltpic_object%weight(k_x,k_v) = qi_weight          
       end do
    end do
  end subroutine load_ltpic_2d




  subroutine deposit_charges( rho )  ! NOT WRITTEN YET
  end subroutine


  subroutine transport_ltpic_2d( flow, ltpic_object )
    class(flow_base), intent(in),  pointer    :: flow
    type(lin_trans_pic_2d), pointer               :: ltpic_object
    sll_int32                                 :: npx
    sll_int32                                 :: npv
    sll_real64                                :: grid_dx
    sll_real64                                :: grid_dv
    sll_real64                                :: inv_2dx
    sll_real64                                :: inv_2dv
    sll_real64                                :: fx
    sll_real64                                :: fv
    sll_real64                                :: fx_plus_h
    sll_real64                                :: fv_plus_h
    sll_real64                                :: fx_minus_h
    sll_real64                                :: fv_minus_h
    sll_real64                                :: fwd_jcbn_mtrx_xx
    sll_real64                                :: fwd_jcbn_mtrx_vx
    sll_real64                                :: fwd_jcbn_mtrx_xv
    sll_real64                                :: fwd_jcbn_mtrx_vv
    sll_real64                                :: bck_jcbn_mtrx_xx
    sll_real64                                :: bck_jcbn_mtrx_vx
    sll_real64                                :: bck_jcbn_mtrx_xv
    sll_real64                                :: bck_jcbn_mtrx_vv
    sll_real64                                :: dfrm_mtrx_xx
    sll_real64                                :: dfrm_mtrx_xv
    sll_real64                                :: dfrm_mtrx_vx
    sll_real64                                :: dfrm_mtrx_vv
    sll_real64                                :: new_dfrm_mtrx_xx
    sll_real64                                :: new_dfrm_mtrx_xv
    sll_real64                                :: new_dfrm_mtrx_vx
    sll_real64                                :: new_dfrm_mtrx_vv

    
    grid_dx = (ltpic_object%xmax-ltpic_object%xmin)/ltpic_object%num_particles_x
    grid_dv = (ltpic_object%vmax-ltpic_object%vmin)/ltpic_object%num_particles_v
    inv_2dx = 1./(2*grid_dx)
    inv_2dv = 1./(2*grid_dv)
    
    ! update particle coordinates and deformation matrices
    npx = ltpic_object%num_particles_x
    npv = ltpic_object%num_particles_v    
    do k_x = 1,npx 
       do k_v = 1,npv
          x = ltpic_object%coord_x(k_x,k_v)
          v = ltpic_object%coord_v(k_x,k_v)

          ! update the particle center 
          call flow%flow_at_xv(x,v, fx,fv)
          ltpic_object%coord_x(k_x,k_v) = fx
          ltpic_object%coord_v(k_x,k_v) = fv

          ! compute the approximate jacobian matrix of the forward flow (use centered finite differences)
          call flow%flow_at_xv( x+grid_dx, v, fx_plus_h,  fv_plus_h  )
          call flow%flow_at_xv( x-grid_dx, v, fx_minus_h, fv_minus_h )
          fwd_jcbn_mtrx_xx = (fx_plus_h-fx_minus_h)*inv_2dx  ! finite diff for d(F_x)/dx
          fwd_jcbn_mtrx_vx = (fv_plus_h-fv_minus_h)*inv_2dx  ! finite diff for d(F_v)/dx
          call flow%flow_at_xv( x, v+grid_dx, fx_plus_h,  fv_plus_h  )
          call flow%flow_at_xv( x, v-grid_dx, fx_minus_h, fv_minus_h )
          fwd_jcbn_mtrx_xv = (fx_plus_h-fx_minus_h)*inv_2dv  ! finite diff for d(F_x)/dv
          fwd_jcbn_mtrx_vv = (fv_plus_h-fv_minus_h)*inv_2dv  ! finite diff for d(F_v)/dv
  
          ! compute the approximate jacobian matrix of the backward flow (with normalized determinant) 
          sqrt_det_jcbn_mtrx = sqrt(abs(fwd_jcbn_mtrx_xx*fwd_jcbn_mtrx_vv - fwd_jcbn_mtrx_xv*fwd_jcbn_mtrx_vx))
          bck_jcbn_mtrx_xx =  sqrt_det_jcbn_mtrx*fwd_jcbn_mtrx_vv
          bck_jcbn_mtrx_xv = -sqrt_det_jcbn_mtrx*fwd_jcbn_mtrx_xv
          bck_jcbn_mtrx_vx = -sqrt_det_jcbn_mtrx*fwd_jcbn_mtrx_vx
          bck_jcbn_mtrx_vv =  sqrt_det_jcbn_mtrx*fwd_jcbn_mtrx_xx

          ! update the particle deformation matrix
          dfrm_mtrx_xx = ltpic_object%deform_matrix_xx
          dfrm_mtrx_xv = ltpic_object%deform_matrix_xv
          dfrm_mtrx_vx = ltpic_object%deform_matrix_vx
          dfrm_mtrx_vv = ltpic_object%deform_matrix_vv
          ltpic_object%deform_matrix_xx = dfrm_mtrx_xx*bck_jcbn_mtrx_xx + dfrm_mtrx_xv*bck_jcbn_mtrx_vx
          ltpic_object%deform_matrix_xv = dfrm_mtrx_xx*bck_jcbn_mtrx_xv + dfrm_mtrx_xv*bck_jcbn_mtrx_vv
          ltpic_object%deform_matrix_vx = dfrm_mtrx_vx*bck_jcbn_mtrx_xx + dfrm_mtrx_vv*bck_jcbn_mtrx_vx
          ltpic_object%deform_matrix_vv = dfrm_mtrx_vx*bck_jcbn_mtrx_xv + dfrm_mtrx_vv*bck_jcbn_mtrx_vv

       end do
    end do
  
  end subroutine
  
end module sll_lin_trans_pic_2d

