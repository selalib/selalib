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
     sll_int32                         :: bspline_cp        ! = (p+1)/2 : the B-spline radius on the cardinal grid (h=1)
     sll_int32                         :: qi_stencil_radius ! the quasi-interpolation stencil depends on the degree of the bsplines
     sll_int32                         :: qi_grid_npx
     sll_int32                         :: qi_grid_npv
     sll_real64                        :: qi_grid_xmin
     sll_real64                        :: qi_grid_xmax
     sll_real64                        :: qi_grid_vmin
     sll_real64                        :: qi_grid_vmax
     sll_real64                        :: nc_poisson_mesh
     sll_real64                        :: xmin_poisson_mesh
     sll_real64                        :: dx_poisson_mesh
     sll_real64                        :: elementary_charge
     sll_int32                         :: bc_type           ! periodic, open domain
     
     sll_real64, dimension(:,:), pointer :: weight
     sll_real64, dimension(:,:), pointer :: coord_x
     sll_real64, dimension(:,:), pointer :: coord_v
     sll_real64, dimension(:,:), pointer :: deform_matrix_xx
     sll_real64, dimension(:,:), pointer :: deform_matrix_xv
     sll_real64, dimension(:,:), pointer :: deform_matrix_vx
     sll_real64, dimension(:,:), pointer :: deform_matrix_vv

     sll_real64, dimension(:,:), pointer :: qi_coefs
     sll_real64, dimension(:,:), pointer :: data ! e.g., the data to be approximated by the particles
  end type sll_spline_2D


#ifdef STDF95
   integer, parameter :: PERIODIC_LTPIC = 0, OPEN_DOMAIN_LTPIC = 1
#else
  enum, bind(C)
     enumerator :: PERIODIC_LTPIC = 0, OPEN_DOMAIN_LTPIC = 1
  end enum
#endif



contains


  function new_ltpic_2d(                                          &
                          num_particles_x,                        &
                          num_particles_v,                        &
                          xmin,                                   &
                          xmax,                                   &
                          vmin,                                   &
                          vmax,                                   &
                          bspline_degree,                         &
                          nc_poisson_mesh,                        &
                          xmin_poisson_mesh,                      &
                          xmax_poisson_mesh,                      &
                          elementary_charge,                      &
                          bc_type                                 &
                          )  

    type(lin_trans_pic_2d), pointer     :: new_ltpic_2d
    sll_int32,  intent(in)              :: num_particles_x
    sll_int32,  intent(in)              :: num_particles_v        
    sll_real64, intent(in)              :: xmin
    sll_real64, intent(in)              :: xmax
    sll_real64, intent(in)              :: vmin
    sll_real64, intent(in)              :: vmax
    sll_int32,  intent(in)              :: bspline_degree
    sll_int32,  intent(in)              :: nc_poisson_mesh
    sll_real64, intent(in)              :: xmin_poisson_mesh
    sll_real64, intent(in)              :: xmax_poisson_mesh
    sll_real64, intent(in)              :: elementary_charge
    sll_real64, dimension(:), pointer   :: qi_coefs_1d
    sll_real64                          :: qi_st_radius
    sll_real64                          :: hx_parts
    sll_real64                          :: hv_parts
    sll_int32                           :: i
    sll_int32                           :: j
                
    SLL_ASSERT( xmin < xmax )        
    SLL_ASSERT( vmin < vmax )        
    SLL_ASSERT( xmin_poisson_mesh < xmax_poisson_mesh )        
    SLL_ASSERT( num_particles_x > 0 )
    SLL_ASSERT( num_particles_v > 0 )
    SLL_ASSERT( nc_poisson_mesh > 0 )
    SLL_ASSERT( (xmax_poisson_mesh-xmin_poisson_mesh == xmax-xmin) .or. (bc_type .ne. PERIODIC_LTPIC) )

    SLL_ALLOCATE( new_ltpic_2d, ierr )
    new_ltpic_2d%num_particles_x = num_particles_x
    new_ltpic_2d%num_particles_v = num_particles_v
    new_ltpic_2d%xmin = xmin
    new_ltpic_2d%xmax = xmax
    new_ltpic_2d%vmin = vmin
    new_ltpic_2d%vmax = vmax
    new_ltpic_2d%nc_poisson_mesh = nc_poisson_mesh
    new_ltpic_2d%xmin_poisson_mesh = xmin_poisson_mesh
    new_ltpic_2d%dx_poisson_mesh = (xmax_poisson_mesh-xmin_poisson_mesh)/num_particles_x
    new_ltpic_2d%elementary_charge = elementary_charge    
    new_ltpic_2d%bc_type = bc_type
    
    hx_parts = (xmax-xmin)/num_particles_x
    hv_parts = (vmax-vmin)/num_particles_v
    
    new_ltpic_2d%bspline_degree = bspline_degree   
    new_ltpic_2d%bspline_cp     = (bspline_degree+1)/2
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
       print *, 'ERROR: new_ltpic_2d(): this bspline degree not implemented in the ltpic module -- value =', bspline_degree
       STOP
    end select
    new_ltpic_2d%qi_grid_npx  = num_particles_x + 2*new_ltpic_2d%qi_stencil_radius
    new_ltpic_2d%qi_grid_npv  = num_particles_v + 2*new_ltpic_2d%qi_stencil_radius
    new_ltpic_2d%qi_grid_xmin = xmin - new_ltpic_2d%qi_stencil_radius*hx_parts
    new_ltpic_2d%qi_grid_xmax = xmax + new_ltpic_2d%qi_stencil_radius*hx_parts
    new_ltpic_2d%qi_grid_vmin = vmin - new_ltpic_2d%qi_stencil_radius*hv_parts
    new_ltpic_2d%qi_grid_vmax = vmax + new_ltpic_2d%qi_stencil_radius*hv_parts
    
    SLL_ALLOCATE( new_ltpic_2d%weight          (num_particles_x,num_particles_v),   ierr )
    SLL_ALLOCATE( new_ltpic_2d%coord_x         (num_particles_x,num_particles_v),   ierr )
    SLL_ALLOCATE( new_ltpic_2d%coord_v         (num_particles_x,num_particles_v),   ierr )
    SLL_ALLOCATE( new_ltpic_2d%deform_matrix_xx(num_particles_x,num_particles_v),   ierr )
    SLL_ALLOCATE( new_ltpic_2d%deform_matrix_xv(num_particles_x,num_particles_v),   ierr )
    SLL_ALLOCATE( new_ltpic_2d%deform_matrix_vx(num_particles_x,num_particles_v),   ierr )
    SLL_ALLOCATE( new_ltpic_2d%deform_matrix_vv(num_particles_x,num_particles_v),   ierr )

    SLL_ALLOCATE( new_ltpic_2d%data            (new_ltpic_2d%qi_grid_npx,new_ltpic_2d%qi_grid_npv),   ierr )

  end function new_ltpic_2d



  ! here we initialize the particles, so that the density approximates the target_density_initializer
  ! note: we don't distinguish here between the periodic and the open_domain cases. In principle the target should be periodic
  !       in the periodic case, but for simplicity we do not use that property here (and don't check it, either).  
  subroutine load_ltpic_2d( target_density_initializer, ltpic_object )
    class(scalar_field_2d_initializer_base), intent(in), pointer :: target_density_initializer   ! intial field to be approximated
    type(lin_trans_pic_2d), pointer      :: ltpic_object

    call target_density_initializer%f_of_x1x2(ltpic_object%data)   ! writes the target density values on the 'data' array        
    call approximate_data( ltpic_object )
  end subroutine load_ltpic_2d


  ! here we (re-)initialize the particles, so that the new density approximates the present one
  ! (as above we don't distinguish between the periodic and the open_domain cases)
  subroutine remap_ltpic_2d( ltpic_object )
     call write_f_on_grid (                                                                                &
                            ltpic_object%qi_grid_xmin,ltpic_object%qi_grid_xmax,ltpic_object%qi_grid_npx,  &
                            ltpic_object%qi_grid_vmin,ltpic_object%qi_grid_vmax,ltpic_object%qi_grid_npv,  &
                            ltpic_object%data, ltpic_object )
     call approximate_data( ltpic_object )
  end subroutine remap_ltpic_2d
  
  
  ! reset the object as a (new) collection of particles with cartesian nodes and weights computed to approximate the stored data
  subroutine approximate_data( ltpic_object )
    type(lin_trans_pic_2d), pointer      :: ltpic_object
    sll_int32                            :: npx
    sll_int32                            :: npv
    sll_int32                            :: k_x
    sll_int32                            :: k_v
    sll_int32                            :: l_x
    sll_int32                            :: l_v
    sll_int32                            :: qi_st_radius
    sll_real64                           :: qi_weight
    sll_real64, dimension(:,:), pointer  :: qi_coefs
    sll_real64                           :: xmin
    sll_real64                           :: vmin
    sll_real64                           :: hx_parts
    sll_real64                           :: hv_parts

    npx           = ltpic_object%num_particles_x
    npv           = ltpic_object%num_particles_v    
    xmin          = ltpic_object%xmin
    hx_parts      = (ltpic_object%xmax-ltpic_object%xmin)/npx
    vmin          = ltpic_object%vmin
    hv_parts      = (ltpic_object%vmax-ltpic_object%vmin)/npv
    qi_st_radius  = ltpic_object%qi_stencil_radius    
    qi_coefs      => ltpic_object%qi_coefs

    ! initial particles have no deformation, hence set D_k = (1. 0. // 0. 1.) for all the particles
    ltpic_object%deform_matrix_xx = 1.
    ltpic_object%deform_matrix_xv = 0.
    ltpic_object%deform_matrix_vx = 0.
    ltpic_object%deform_matrix_vv = 1.

    ! compute particle coordinates and weights
    do k_x = 1,npx 
       do k_v = 1,npv           
          ltpic_object%coord_x(k_x,k_v) = xmin + k_x*hx_parts
          ltpic_object%coord_v(k_x,k_v) = vmin + k_v*hv_parts
          ! compute the weight with local quasi-interpolation scheme
          qi_weight = 0
          do l_x = -qi_st_radius, qi_st_radius
             do l_v = -qi_st_radius, qi_st_radius
                ! due to the quasi-interpolation stencils there are more nodes (in the data grid) than particles -- hence the offset
                qi_weight = qi_weight + qi_coefs(l_x,l_v) * ltpic_object%data( k_x+l_x + qi_st_radius, k_v+l_v + qi_st_radius )
             end do
          end do
          ltpic_object%weight(k_x,k_v) = qi_weight          
       end do
    end do

  end subroutine approximate_data


  ! deposit the charges with a PIC-like (point particle) algorithm, on the poisson grid
  ! { x_i = xmin_poisson_mesh + (i-1) * dx_poisson_mesh    for    i = 1 .. nc_poisson_mesh }
  subroutine deposit_charges( rho, ltpic_object )
    sll_real64, dimension(:), intent(inout)   :: rho
    type(lin_trans_pic_2d), pointer           :: ltpic_object

    sll_float64                               :: x
    sll_float64                               :: v
    sll_float64                               :: inv_hx_pm
    sll_float64                               :: hx_pm
    sll_float64                               :: xmin_pm
    sll_float64                               :: charge
    sll_float64                               :: x_i
    sll_int32                                 :: npx
    sll_int32                                 :: npv
    sll_int32                                 :: k_x
    sll_int32                                 :: k_v
    sll_int32                                 :: i_min
    sll_int32                                 :: i_max
    sll_int32                                 :: i
    sll_int32                                 :: i_per
    sll_int32                                 :: nc_pm
    sll_int32                                 :: cp
    sll_int32                                 :: degree
    logical                                   :: periodic
    
    hx_pm                 = ltpic_object%dx_poisson_mesh
    inv_hx_pm             = 1./hx_pm
    xmin_pm               = ltpic_object%xmin_poisson_mesh
    degree                = ltpic_object%bspline_degree
    cp                    = ltpic_object%bspline_cp
    nc_pm                 = ltpic_object%nc_poisson_mesh
    periodic              = (ltpic_object%bc_type == PERIODIC_LTPIC)

    ! Check that rho vector is associated to the poisson mesh
    SLL_ASSERT( size(rh0)==nc_pm+1 )

    ! loop over the particles
    npx = ltpic_object%num_particles_x
    npv = ltpic_object%num_particles_v    
    do k_x = 1,npx 
       do k_v = 1,npv
          x = ltpic_object%coord_x(k_x,k_v)
          v = ltpic_object%coord_v(k_x,k_v)
          charge = ltpic_object%elementary_charge*ltpic_object%weight(k_x,k_v)
          
          if( charge .neq. 0 ) then
             ! loop over i such that Bspline_h(x_i - x) does not vanish: ie, such that x_i (see above) is in ]x-hx*cp,x+hx*cp[
             i_min = 1 + int(floor  ( inv_hx_pm*(x-xmin_pm)-cp )) + 1
             i_max = 1 + int(ceiling( inv_hx_pm*(x-xmin_pm)+cp )) - 1
             do i = i_min, i_max
                if (periodic) then
                   i_rho = 1 + mod( i-1 + 2*nc_pm, nc_pm )    ! periodic boundary conditions                
                else
                   i_rho = i
                end if
                if( i_rho >= 1 .and. i_rho <= nc_pm ) then
                   x_i = xmin_pm + (i-1) * hx_pm
                   rho(i_rho) = rho(i_rho) + charge*b_spline_1d(degree,inv_hx_pm*(x-x_i))
                end if
             end do
          end if
       end do
    end do
    if (periodic) then
       rho(nc_pm+1) = rho(1)
    end if
  end subroutine deposit_charges


  subroutine write_f_on_grid ( xmin_grid,xmax_grid,npx_grid,vmin_grid,vmax_grid,npv_grid,grid_values,ltpic_object )
    type(lin_trans_pic_2d), pointer             :: ltpic_object
    sll_int32, intent(in)                       :: xmin_grid
    sll_int32, intent(in)                       :: xmax_grid
    sll_int32, intent(in)                       :: npx_grid
    sll_int32, intent(in)                       :: vmin_grid
    sll_int32, intent(in)                       :: vmax_grid
    sll_int32, intent(in)                       :: npv_grid
    sll_real64, dimension(:,:), intent(inout)   :: grid_values
    sll_int32                                   :: npx
    sll_int32                                   :: npv
    sll_int32                                   :: k_x
    sll_int32                                   :: k_v
    sll_int32                                   :: cp
    sll_int32                                   :: degree
    sll_int32                                   :: m_min
    sll_int32                                   :: m_max
    sll_int32                                   :: m
    sll_int32                                   :: ix_min
    sll_int32                                   :: ix_max    
    sll_int32                                   :: ix
    sll_int32                                   :: iv_min
    sll_int32                                   :: iv_max    
    sll_int32                                   :: iv
    sll_real64                                  :: hx_parts
    sll_real64                                  :: hv_parts
    sll_real64                                  :: inv_hx_parts
    sll_real64                                  :: inv_hv_parts
    sll_real64                                  :: x
    sll_real64                                  :: v
    sll_real64                                  :: D_xx
    sll_real64                                  :: D_xv
    sll_real64                                  :: D_vx
    sll_real64                                  :: D_vv
    sll_real64                                  :: part_radius_x
    sll_real64                                  :: part_radius_v
    sll_real64                                  :: x_im
    sll_real64                                  :: v_i
    sll_real64                                  :: hx_grid
    sll_real64                                  :: hv_grid
    sll_real64                                  :: inv_hx_grid
    sll_real64                                  :: inv_hv_grid
    sll_real64                                  :: period_x
    sll_real64                                  :: inv_period_x
    logical                                     :: periodic
    
    cp           = ltpic_object%bspline_cp
    degree       = ltpic_object%bspline_degree
    npx          = ltpic_object%num_particles_x
    npv          = ltpic_object%num_particles_v    
    hx_parts     = (ltpic_object%xmax-ltpic_object%xmin)/npx
    hv_parts     = (ltpic_object%vmax-ltpic_object%vmin)/npv
    hx_grid      = (xmax_grid-xmin_grid)/npx_grid
    hv_grid      = (vmax_grid-vmin_grid)/npv_grid
    inv_hx_parts = 1./hx_parts
    inv_hv_parts = 1./hv_parts
    inv_hx_grid  = 1./hx_grid
    inv_hv_grid  = 1./hv_grid

    periodic     = (ltpic_object%bc_type == PERIODIC_LTPIC)
    if( periodic ) then
       period_x     = ltpic_object%xmax-ltpic_object%xmin
       inv_period_x = 1./period_x
    else
       ! values don't matter in this case
       period_x     = 1
       inv_period_x = 1
    end if
    
    ! loop over the particles first, then over the relevant grid points
    do k_x = 1,npx 
       do k_v = 1,npv
          x       = ltpic_object%coord_x(k_x,k_v)
          v       = ltpic_object%coord_v(k_x,k_v)
          D_xx    = ltpic_object%deform_matrix_xx(k_x,k_v)
          D_xv    = ltpic_object%deform_matrix_xv(k_x,k_v)
          D_vx    = ltpic_object%deform_matrix_vx(k_x,k_v)
          D_vv    = ltpic_object%deform_matrix_vv(k_x,k_v)
          weight  = ltpic_object%weight(k_x,k_v)
          
          part_radius_x = cp * ( abs( hx_parts*D_vv ) + abs( hv_parts*D_xv ) )
          part_radius_v = cp * ( abs( hx_parts*D_vx ) + abs( hv_parts*D_xx ) )
          
          if( periodic ) then
             ! loop over every m such that there may be an (x_{i,m},v_i) in the support of the particle (more details in the doc)
             m_min = int(floor  ( inv_period_x*(xmin_grid-x-part_radius_x) )) + 1
             m_max = int(ceiling( inv_period_x*(xmax_grid-x+part_radius_x) )) - 1
          else
             ! no m-loop in the non-periodic case
             m_min = 0
             m_max = 0
          end if                    
          do m = m_min, m_max
             ! loop over every i=(ix,iv) such that (x_{i,m},v_i) may be in the support of the particle (more details in the doc)
             ix_min = 1 + int(floor  ( inv_hx_grid*(x-xmin_grid+m*period_x-part_radius_x) )) + 1
             ix_max = 1 + int(ceiling( inv_hx_grid*(x-xmin_grid+m*period_x+part_radius_x) )) - 1
             iv_min = 1 + int(floor  ( inv_hv_grid*(v-vmin_grid-part_radius_v) )) + 1
             iv_max = 1 + int(ceiling( inv_hv_grid*(v-vmin_grid+part_radius_v) )) - 1          
             do ix = ix_min, ix_max
                if( ix >= 1 .and. ix <= npx_grid ) then
                   x_im = xmin_grid + (ix-1) * hx_grid - m*period_x
                   do iv = iv_min, iv_max
                      if( iv >= 1 .and. iv <= npv_grid ) then
                          v_i = vmin_grid + (iv-1) * hv_grid
                          grid_values(ix,iv) = grid_values(ix,iv)                                                                  &
                                 + weight * inv_hx_parts * b_spline_1d(degree,inv_hx_parts*(D_xx*(x_im-x) + Dxv*(v_i-v)))          &
                                          * inv_hv_parts * b_spline_1d(degree,inv_hv_parts*(D_vx*(x_im-x) + Dvv*(v_i-v)))
                      end if
                  end do
              end if
          end do
       end do
    end do
  end subroutine write_f_on_grid


  ! computes the phase-space density carried by one (weighted) particle, 
  ! i.e., computes w_k * \varphi_h(D_k(x-x_k,v-v_k))   where  \varphi_h(x,v) = 1/h**2 * bspline_1d(x/h) * bspline_1d(v/h)
  function particle_value ( x,v,w_k,x_k,v_k,Dxx_k,Dxv_k,Dvx_k,Dvv_k,inv_h_x,inv_h_v )
    sll_real64                                :: particle_value
    sll_real64, intent(in)                    :: x
    sll_real64, intent(in)                    :: v
    sll_real64, intent(in)                    :: w_k
    sll_real64, intent(in)                    :: x_k
    sll_real64, intent(in)                    :: v_k
    sll_real64, intent(in)                    :: Dxx_k
    sll_real64, intent(in)                    :: Dxv_k
    sll_real64, intent(in)                    :: Dvx_k
    sll_real64, intent(in)                    :: Dvv_k
    sll_real64, intent(in)                    :: inv_h_x
    sll_real64, intent(in)                    :: inv_h_v

    particle_value = w_k * inv_h_x * b_spline_1d(degree,inv_h_x*(Dxx*(x-x_k) + Dxv*(v-v_k)))            &
                         * inv_h_v * b_spline_1d(degree,inv_h_v*(Dvx*(x-x_k) + Dvv*(v-v_k)))

  end function particle_value


  function b_spline_1d (degree, x)
    sll_real64                 :: b_spline_1d
    sll_int32,  intent(in)     :: degree
    sll_real64, intent(in)     :: x
    sll_real64                 :: x_aux

    b_spline_1d = -10000
    select case (degree)
    case (1)
       ! affine spline
       if ( x <= -1 .or. x >= 1 ) then
          b_spline_1d = 0
          return 
       end if
       if ( x < 0 ) then
          b_spline_1d = x+1
          return
       end if
       b_spline_1d = 1-x
       return
    case (3)
       ! cubic spline
       x_aux = x+2
       if ( x_aux <= 0 .or. x_aux >= 4 ) then
          b_spline_1d = 0
          return 
       end if
       if ( x_aux < 1 ) then
          b_spline_1d = 0.16666666666666666*(x_aux**3)
          return 
       end if
       if ( x_aux < 2 ) then
          b_spline_1d = 0.6666666666666666 - 2.*x_aux + 2.*(x_aux**2) - 0.5*(x_aux**3)
          return 
       end if
       if ( x_aux < 3 ) then
          b_spline_1d = -7.333333333333333 + 10.*x_aux - 4.*(x_aux**2) + 0.5*(x_aux**3)
          return 
       end if
       b_spline_1d = 10.66666666666666 - 8.*x_aux + 2.*(x_aux**2) - 0.16666666666666666*(x_aux**3)
       return        
    case(5)
       ! quintic spline
       x_aux = x+3
       if ( x_aux <= 0 .or. x_aux >= 6 ) then
          b_spline_1d = 0
          return 
       end if
       if ( x_aux < 1 ) then
          b_spline_1d = 0.00833333333333333*(x_aux**5)
          return 
       end if
       if ( x_aux < 2 ) then
          b_spline_1d =-0.0416666666667 *(x_aux**5) + 0.25 *(x_aux**4) - 0.5 *(x_aux**3) + 0.5 *(x_aux**2) - 0.25 *(x_aux) + 0.05
          return 
       end if
       if ( x_aux < 3 ) then
          b_spline_1d = 0.0833333333333 *(x_aux**5) - 1.0 *(x_aux**4) + 4.5 *(x_aux**3) - 9.5 *(x_aux**2) + 9.75 *(x_aux) - 3.95
          return 
       end if
       if ( x_aux < 4 ) then
          b_spline_1d =-0.0833333333333 *(x_aux**5) + 1.5 *(x_aux**4) - 10.5 *(x_aux**3) + 35.5 *(x_aux**2) - 57.75 *(x_aux) + 36.55
          return 
       end if        
       if ( x_aux < 5 ) then
          b_spline_1d = 0.0416666666667 *(x_aux**5) - 1.0 *(x_aux**4) + 9.5 *(x_aux**3) - 44.5 *(x_aux**2) + 102.25 *(x_aux) - 91.45
          return 
       end if
       b_spline = -0.008333333333333333 * ((x_aux-6.)**5)
       return
    
    case default
       print *, 'ERROR: b_spline_1d(): this degree not implemented in the ltpic module -- value =', degree
       STOP
    end select
  end function b_spline_1d

  subroutine transport_ltpic_2d ( flow, ltpic_object )
    class(flow_base), intent(in),  pointer    :: flow
    type(lin_trans_pic_2d), pointer           :: ltpic_object
    sll_int32                                 :: npx
    sll_int32                                 :: npv
    sll_real64                                :: hx_parts
    sll_real64                                :: hv_parts
    sll_real64                                :: inv_2hx
    sll_real64                                :: inv_2hv
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

    
    hx_parts = (ltpic_object%xmax-ltpic_object%xmin)/ltpic_object%num_particles_x
    hv_parts = (ltpic_object%vmax-ltpic_object%vmin)/ltpic_object%num_particles_v
    inv_2hx = 1./(2*hx_parts)
    inv_2hv = 1./(2*hv_parts)
    
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
          call flow%flow_at_xv( x+hx_parts, v, fx_plus_h,  fv_plus_h  )
          call flow%flow_at_xv( x-hx_parts, v, fx_minus_h, fv_minus_h )
          fwd_jcbn_mtrx_xx = (fx_plus_h-fx_minus_h)*inv_2hx  ! finite diff for d(F_x)/dx
          fwd_jcbn_mtrx_vx = (fv_plus_h-fv_minus_h)*inv_2hx  ! finite diff for d(F_v)/dx
          call flow%flow_at_xv( x, v+hx_parts, fx_plus_h,  fv_plus_h  )
          call flow%flow_at_xv( x, v-hx_parts, fx_minus_h, fv_minus_h )
          fwd_jcbn_mtrx_xv = (fx_plus_h-fx_minus_h)*inv_2hv  ! finite diff for d(F_x)/dv
          fwd_jcbn_mtrx_vv = (fv_plus_h-fv_minus_h)*inv_2hv  ! finite diff for d(F_v)/dv
  
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
  
  end subroutine transport_ltpic_2d
  
end module sll_lin_trans_pic_2d

