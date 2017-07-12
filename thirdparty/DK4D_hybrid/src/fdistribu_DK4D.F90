!===========================================================================
!> Distribution function for 4D drift-kinetic hybrid simulation
!>  f(x,y,eta3,vpar)
!>
!> \date 2014-08-20
!> \author V. Grandgirard
!---------------------------------------------------------------------------
module fdistribu_DK4D_module
#include "sll_working_precision.h"
#include "sll_assert.h"
#include "sll_memory.h"

  use equilibrium_DK4D_module
  use input_DK4D_module
  use mesh_DK4D_module
  use sll_m_arbitrary_degree_spline_interpolator_1d
  use sll_m_arbitrary_degree_spline_interpolator_2d
  use sll_m_collective, only : sll_v_world_collective, &
      sll_f_get_collective_size, sll_f_get_collective_rank, &
      sll_s_collective_barrier
  use sll_m_remapper
  use sll_m_utilities, only : sll_f_is_even
  use utils_DK4D_module

  implicit none

  !===========================================================================
  !> 
  !---------------------------------------------------------------------------
  type, public :: fdistribu_DK4D_t

    !> Boundary conditions
    type(boundary_conditions_4d_t) :: bound_cond

    !> Sequential in (x1,x2) and parallel in (x3,x4)
    type(sll_t_layout_4d)               , pointer :: layout4d_seqx1x2
    sll_real64, dimension(:,:,:,:), pointer :: val4d_seqx1x2 

    !> Parallel in (x1,x2) and sequential in (x3,x4) 
    type(sll_t_layout_4d)               , pointer :: layout4d_seqx3x4
    sll_real64, dimension(:,:,:,:), pointer :: val4d_seqx3x4

    !> For remapping
    type(sll_t_remap_plan_4D_real64), pointer :: seqx1x2_to_seqx3x4
    type(sll_t_remap_plan_4D_real64), pointer :: seqx3x4_to_seqx1x2

    !> For interpolations
    type(spline_degree_4d_t)          :: spline_degree
    type(sll_t_arbitrary_degree_spline_interpolator_2d)     :: interp2d_f_eta1eta2
    type(sll_t_arbitrary_degree_spline_interpolator_1d) :: interp1d_f_eta3
    type(sll_t_arbitrary_degree_spline_interpolator_1d) :: interp1d_f_vpar

  end type fdistribu_DK4D_t
  !---------------------------------------------------------------------------

contains

  !===========================================================================
  !> Distribution function: Allocation 
  !---------------------------------------------------------------------------
  subroutine new_fdistribu_DK4D( fdistribu, mesh4d, &
      bound_cond, spline_degree )

    type(fdistribu_DK4D_t)        , intent(inout) :: fdistribu
    type(mesh_DK4D_t)             , intent(in)    :: mesh4d
    type(boundary_conditions_4d_t), intent(in)    :: bound_cond
    type(spline_degree_4d_t)      , intent(in)    :: spline_degree


    !-> Local variables 
    sll_int32 :: ierr
    !--> For MPI parallelization
    sll_int32 :: world_size
    sll_int32 :: my_rank
    sll_int32 :: power2
    sll_int32 :: nproc_x1, nproc_x2, nproc_x3, nproc_x4
    !--> For parallel distribution
    sll_int32 :: npoints_x1, npoints_x2, npoints_x3, npoints_x4
    sll_int32 :: loc4d_sz_x1, loc4d_sz_x2, loc4d_sz_x3, loc4d_sz_x4

    !**********************************************************
    !*** Initialization of the parallel distribution        ***
    !***  layout for sequential operations in x3 and x4.    ***
    !***  Make an even split for x1 and x2, or as close as  ***
    !***  even if the power of 2 is odd. This should        ***
    !***  be packaged in some sort of routine and set up    ***
    !***  at initialization time.                           ***
    !**********************************************************
    world_size = sll_f_get_collective_size(sll_v_world_collective)
    my_rank    = sll_f_get_collective_rank(sll_v_world_collective)
    power2     = int(log(real(world_size))/log(2.0))

    !**********************************************************************
    !*** Initialization of parallel layout of f4d in (x1,x2) directions ***
    !***  (x1,x2) : parallelized layout                                 ***
    !***  (x3,x4) : sequential                                          ***
    !**********************************************************************
    !--> special case N = 1, so power2 = 0
    if ( power2 == 0 ) then
      nproc_x1 = 1
      nproc_x2 = 1
      nproc_x3 = 1
      nproc_x4 = 1
    end if
    
    if ( sll_f_is_even(power2) ) then
      nproc_x1 = 2**(power2/2)
      nproc_x2 = 2**(power2/2)
      nproc_x3 = 1
      nproc_x4 = 1
    else 
      nproc_x1 = 2**((power2-1)/2)
      nproc_x2 = 2**((power2+1)/2)
      nproc_x3 = 1
      nproc_x4 = 1
    end if

    fdistribu%layout4d_seqx3x4  => sll_f_new_layout_4d( sll_v_world_collective )
    npoints_x1 = mesh4d%eta1_eta2_mesh2d%num_cells1 + 1
    npoints_x2 = mesh4d%eta1_eta2_mesh2d%num_cells2 + 1
    npoints_x3 = mesh4d%eta3_mesh1d%num_cells + 1
    npoints_x4 = mesh4d%vpar_mesh1d%num_cells + 1

    call sll_o_initialize_layout_with_distributed_array( &
      npoints_x1, &
      npoints_x2, &
      npoints_x3, &
      npoints_x4, &
      nproc_x1, &
      nproc_x2, &
      nproc_x3, &
      nproc_x4, &
      fdistribu%layout4d_seqx3x4 )
        
    call sll_o_compute_local_sizes( fdistribu%layout4d_seqx3x4, &
      loc4d_sz_x1, &
      loc4d_sz_x2, &
      loc4d_sz_x3, &
      loc4d_sz_x4 )    
    SLL_ALLOCATE( fdistribu%val4d_seqx3x4(loc4d_sz_x1,loc4d_sz_x2,loc4d_sz_x3,loc4d_sz_x4),ierr )

    !**********************************************************************
    !*** Initialization of parallel layout of f4d in (x3,x4) directions ***
    !***  (x1,x2) : sequential                                          ***
    !***  (x3,x4) : parallelized layout                                 ***
    !**********************************************************************
    SLL_ASSERT( npoints_x3.ge.power2 )
    nproc_x3 = nproc_x1*nproc_x2
    nproc_x1 = 1
    nproc_x2 = 1
    nproc_x4 = 1

    fdistribu%layout4d_seqx1x2  => sll_f_new_layout_4d( sll_v_world_collective )
    call sll_o_initialize_layout_with_distributed_array( &
      npoints_x1, & 
      npoints_x2, & 
      npoints_x3, &
      npoints_x4, &
      nproc_x1, &
      nproc_x2, &
      nproc_x3, &
      nproc_x4, &
      fdistribu%layout4d_seqx1x2 )
    
    ! Allocate the array needed to store the local chunk 
    ! of the distribution function data. First compute the 
    ! local sizes. Since the remap operations
    ! are out-of-place, we will allocate two different arrays, 
    ! one for each layout.
    call sll_o_compute_local_sizes( fdistribu%layout4d_seqx1x2, &
      loc4d_sz_x1, &
      loc4d_sz_x2, &
      loc4d_sz_x3, &
      loc4d_sz_x4 )
    SLL_ALLOCATE( fdistribu%val4d_seqx1x2(loc4d_sz_x1,loc4d_sz_x2,loc4d_sz_x3,loc4d_sz_x4), ierr )

    !****************************************
    !*** Initialization of the remaping   ***
    !****************************************
    fdistribu%seqx3x4_to_seqx1x2 => &
      sll_o_new_remap_plan( fdistribu%layout4d_seqx3x4, &
      fdistribu%layout4d_seqx1x2, fdistribu%val4d_seqx3x4 )
    
    fdistribu%seqx1x2_to_seqx3x4 => &
      sll_o_new_remap_plan( fdistribu%layout4d_seqx1x2, &
      fdistribu%layout4d_seqx3x4, fdistribu%val4d_seqx1x2 )

    !********************************************
    !*** Initialization of the interpolations ***
    !********************************************
    !--> Initialize  boundary conditions
    fdistribu%bound_cond%left_eta1  = bound_cond%left_eta1
    fdistribu%bound_cond%right_eta1 = bound_cond%right_eta1
    fdistribu%bound_cond%left_eta2  = bound_cond%left_eta2
    fdistribu%bound_cond%right_eta2 = bound_cond%right_eta2
    fdistribu%bound_cond%left_eta3  = bound_cond%left_eta3
    fdistribu%bound_cond%right_eta3 = bound_cond%right_eta3
    fdistribu%bound_cond%left_vpar  = bound_cond%left_vpar
    fdistribu%bound_cond%right_vpar = bound_cond%right_vpar

    !--> Initialize the spline degree
    fdistribu%spline_degree%eta1 = spline_degree%eta1
    fdistribu%spline_degree%eta2 = spline_degree%eta2
    fdistribu%spline_degree%eta3 = spline_degree%eta3
    fdistribu%spline_degree%vpar = spline_degree%vpar

    !--> Initialize the different interpolators
    call fdistribu%interp2d_f_eta1eta2%init( &
      mesh4d%eta1_eta2_mesh2d%num_cells1+1, &
      mesh4d%eta1_eta2_mesh2d%num_cells2+1, &
      mesh4d%eta1_eta2_mesh2d%eta1_min, &
      mesh4d%eta1_eta2_mesh2d%eta1_max, &
      mesh4d%eta1_eta2_mesh2d%eta2_min, &
      mesh4d%eta1_eta2_mesh2d%eta2_max, &
      fdistribu%bound_cond%left_eta1, &
      fdistribu%bound_cond%right_eta1, &
      fdistribu%bound_cond%left_eta2, &
      fdistribu%bound_cond%right_eta2, &
      fdistribu%spline_degree%eta1, &
      fdistribu%spline_degree%eta2 )

    call fdistribu%interp1d_f_eta3%init( &
       mesh4d%eta3_mesh1d%num_cells + 1, &
       mesh4d%eta3_mesh1d%eta_min, &
       mesh4d%eta3_mesh1d%eta_max, &
       fdistribu%bound_cond%left_eta3, &
       fdistribu%bound_cond%right_eta3, &
       fdistribu%spline_degree%eta3 )

    call fdistribu%interp1d_f_vpar%init( &
       mesh4d%vpar_mesh1d%num_cells + 1, &
       mesh4d%vpar_mesh1d%eta_min, &
       mesh4d%vpar_mesh1d%eta_max, &
       fdistribu%bound_cond%left_vpar, &
       fdistribu%bound_cond%right_vpar, &
       fdistribu%spline_degree%vpar )

  end subroutine new_fdistribu_DK4D
  !---------------------------------------------------------------------------


  !===========================================================================
  !> Distribution function: Deallocation 
  !---------------------------------------------------------------------------
  subroutine delete_fdistribu_DK4D( fdistribu )

    type(fdistribu_DK4D_t), intent(inout) :: fdistribu

    sll_int32 :: ierr

    SLL_DEALLOCATE( fdistribu%val4d_seqx1x2, ierr )
    SLL_DEALLOCATE( fdistribu%val4d_seqx3x4, ierr )
    call sll_o_delete( fdistribu%layout4d_seqx1x2 )
    call sll_o_delete( fdistribu%layout4d_seqx3x4 )

  end subroutine delete_fdistribu_DK4D
  !---------------------------------------------------------------------------


  !===========================================================================
  !> Distribution function: Initialization
  !---------------------------------------------------------------------------
  subroutine init_fdistribu_DK4D( fdistribu, &
      input_data, mesh4d, equilibrium )

    use sll_m_constants, only: sll_p_pi
    
    type(fdistribu_DK4D_t)  , intent(inout) :: fdistribu
    type(input_DK4D_t)      , intent(in)    :: input_data
    type(mesh_DK4D_t)       , intent(in)    :: mesh4d
    type(equilibrium_DK4D_t), intent(in)    :: equilibrium

    !-> Local variables
    !--> For g(x,y) definition 
    sll_int32  :: ierr
    sll_int32  :: ix, iy, Nx, Ny
    sll_real64 :: x_point, y_point
    sll_real64 :: r_peak, r_min, r_max, r_point
    sll_real64, dimension(:), pointer :: r_grid
    sll_real64, dimension(:,:), pointer :: g_xy
    !--> For parallel remaping
    sll_int32  :: i1, i2, i3, i4
    sll_int32  :: iloc1, iloc2
    sll_int32  :: loc4d_sz_x1, loc4d_sz_x2, loc4d_sz_x3, loc4d_sz_x4
    sll_int32, dimension(1:4) :: glob_ind4d
    !--> For f(x,y,eta3,vpar) definition
    sll_int32  :: Nphi, Nvpar
    sll_real64 :: Ltheta, Lphi, ktheta, kphi
    sll_real64 :: theta_j, phi_k
    
    !*** Definition of g(x,y) which is used to force    ***
    !***  the pertubation equal to 0 for rmin and rmax  ***
    !---> BECAREFUL: This definition is not valid for the generalized case
    !> \todo Check if this term is useful even for circular geometry 
    Nx = size(mesh4d%xgrid_2d,1)
    Ny = size(mesh4d%xgrid_2d,2)
    SLL_ALLOCATE( g_xy(Nx,Ny), ierr )
    if ( mesh4d%geometry_type == CIRCULAR_GEOMETRY ) then
      r_grid => mesh4d%rgrid_2d(:,1)
      r_min  = r_grid(1)
      r_max  = r_grid(size(r_grid))
      r_peak = equilibrium%r_peak
      do iy = 1, Ny
        do ix = 1, Nx
          x_point = mesh4d%xgrid_2d(ix,iy)
          y_point = mesh4d%ygrid_2d(ix,iy)
          r_point = sqrt(x_point**2+y_point**2)
          g_xy(ix,iy) = exp(-(r_point-r_peak)**2/ &
              (4._f64*input_data%deltarn/input_data%deltarTi)) 
        end do
      end do
    else
      g_xy(:,:) = 1._f64
    end if

    !*** Initialization of the distribution function f4d_seqx3x4 ***
    call sll_o_compute_local_sizes( fdistribu%layout4d_seqx3x4, &
      loc4d_sz_x1, &
      loc4d_sz_x2, &
      loc4d_sz_x3, &
      loc4d_sz_x4 )

    Ltheta = 2._f64*sll_p_pi
    ktheta = 2._f64*sll_p_pi*real(input_data%mmode)/Ltheta

    Nphi   = size(mesh4d%eta3_grid)
    Lphi   = abs(mesh4d%eta3_grid(Nphi)-mesh4d%eta3_grid(1))
    kphi   = 2._f64*sll_p_pi*real(input_data%nmode)/Lphi

    Nvpar  = size(mesh4d%vpar_grid)
    do i4 = 1,Nvpar
      do i3 = 1,Nphi
        do iloc2 = 1,loc4d_sz_x2
          do iloc1 = 1,loc4d_sz_x1
            glob_ind4d(:) = sll_o_local_to_global( &
                fdistribu%layout4d_seqx3x4, &
                (/iloc1,iloc2,i3,i4/) )
            i1 = glob_ind4d(1)
            i2 = glob_ind4d(2)
            theta_j = modulo(atan2( &
                mesh4d%ygrid_2d(i1,i2), &
                mesh4d%xgrid_2d(i1,i2)), 2._f64*sll_p_pi)
            phi_k = mesh4d%eta3_grid(i3) 
            fdistribu%val4d_seqx3x4(iloc1,iloc2,i3,i4) = &
                equilibrium%feq_xyvpar(i1,i2,i4) * &
                ( 1._f64 + input_data%eps_perturb * &
                cos(ktheta*theta_j + kphi*phi_k) * g_xy(i1,i2) )
          end do
        end do
      end do
    end do

!VG!!=>MB!    k_x2  = 2._f64*sll_p_pi/(sim%m_x2%eta_max - sim%m_x2%eta_min)
!VG!!=>MB!    k_x3  = 2._f64*sll_p_pi/(sim%m_x3%eta_max - sim%m_x3%eta_min)
!VG!!=>MB!
!VG!!=>MB!    rpeak = x1_min+sim%rho_peak*(x1_max-x1_min)
!VG!!=>MB!
!VG!!=>MB!    do iloc4 = 1,loc4d_sz_x4
!VG!!=>MB!      do iloc3 = 1,loc4d_sz_x3
!VG!!=>MB!        do iloc2 = 1,loc4d_sz_x2
!VG!!=>MB!          do iloc1 = 1,loc4d_sz_x1
!VG!!=>MB!            glob_ind(:) = sll_o_local_to_global(layout, &
!VG!!=>MB!              (/iloc1,iloc2,iloc3,iloc4/))
!VG!!=>MB!            i1 = glob_ind(1)
!VG!!=>MB!            i2 = glob_ind(2)
!VG!!=>MB!            i3 = glob_ind(3)
!VG!!=>MB!            i4 = glob_ind(4)
!VG!!=>MB!            tmp_mode = cos(real(sim%nmode,f64)*k_x3*x3_node(i3)&
!VG!!=>MB!               +real(sim%mmode,f64)*k_x2*x2_node(i2))
!VG!!=>MB!            tmp = exp(-(x1_node(i1)-rpeak)**2/(4._f64*sim%deltarn/sim%deltarTi)) 
!VG!!=>MB!            f4d(iloc1,iloc2,iloc3,iloc4) = &
!VG!!=>MB!              (1._f64+tmp_mode*sim%eps_perturb*tmp)*sim%feq_x1x4(i1,i4)          
!VG!!=>MB!          end do
!VG!!=>MB!        end do
!VG!!=>MB!      end do
!VG!!=>MB!    end do

    SLL_DEALLOCATE( g_xy, ierr )

    call sll_o_apply_remap_4d( fdistribu%seqx3x4_to_seqx1x2, &
        fdistribu%val4d_seqx3x4, fdistribu%val4d_seqx1x2 )

  end subroutine init_fdistribu_DK4D
  !---------------------------------------------------------------------------


  !===========================================================================
  !> Distribution function: Initialization
  !---------------------------------------------------------------------------
  subroutine print_fdistribu_DK4D( &
      fdistribu, &
      ix1_diag, &
      ix2_diag, &
      ix3_diag, &
      ix4_diag, &
      idiag_num )

    use sll_m_collective
    use sll_m_hdf5_io_serial, only: &
        sll_s_hdf5_ser_file_create, sll_s_hdf5_ser_file_close, &
        sll_o_hdf5_ser_write_array, sll_t_hdf5_ser_handle

    type(fdistribu_DK4D_t), intent(in) :: fdistribu
    sll_int32             , intent(in) :: ix1_diag
    sll_int32             , intent(in) :: ix2_diag
    sll_int32             , intent(in) :: ix3_diag
    sll_int32             , intent(in) :: ix4_diag
    sll_int32   , optional, intent(in) :: idiag_num 

    !-> Local variables
    !--> For parallelisation
    sll_int32 :: my_rank
    !--> For  HDF5 saving
    integer   :: file_err
    type(sll_t_hdf5_ser_handle) :: handle    !< file handle
    character(len=80)    :: filename_HDF5
    character(20) , save :: numfmt = "'_d',i5.5"
    !--> For cross-section saving 
    sll_int32 :: ix1_diag_loc, ix2_diag_loc
    sll_int32 :: ix1_min_loc, ix1_max_loc
    sll_int32 :: ix2_min_loc, ix2_max_loc
    sll_int32 :: ix3_diag_loc, ix4_diag_loc
    sll_int32 :: ieta3_min_loc, ieta3_max_loc
    sll_int32 :: ivpar_min_loc, ivpar_max_loc
    sll_int32, dimension(1:4) :: indx_loc
    sll_int32, dimension(1:4) :: indx_diag

    !*** Construction of the name for the output HDF5 file ***
    if ( present(idiag_num) ) then
      write(filename_HDF5,'(A,'//numfmt//',A)') &
          'fdistribu', idiag_num, ".h5"
    else
      write(filename_HDF5,'(A)') "fdistribu.h5"      
    end if

    
    !*** HDF5 saving ***
    my_rank = sll_f_get_collective_rank(sll_v_world_collective)

    !--> HDF5 saving of fdistribu_x1x2
    ieta3_min_loc = sll_o_get_layout_k_min( fdistribu%layout4d_seqx1x2, my_rank )
    ieta3_max_loc = sll_o_get_layout_k_max( fdistribu%layout4d_seqx1x2, my_rank )
    ivpar_min_loc = sll_o_get_layout_l_min( fdistribu%layout4d_seqx1x2, my_rank )
    ivpar_max_loc = sll_o_get_layout_l_max( fdistribu%layout4d_seqx1x2, my_rank )

    if ( (ieta3_min_loc.le.ix3_diag) .and. (ix3_diag.le.ieta3_max_loc) .and. &
        (ivpar_min_loc.le.ix4_diag) .and. (ix4_diag.le.ivpar_max_loc) ) then

      indx_loc = sll_o_global_to_local( fdistribu%layout4d_seqx1x2, &
          (/1,1,ix3_diag,ix4_diag/) )
      ix3_diag_loc = indx_loc(3)
      ix4_diag_loc = indx_loc(4)

      call sll_s_hdf5_ser_file_create( trim(filename_HDF5), handle, file_err )
      !--> Saving of fdistribu_indx_diag
      indx_diag(1) = ix1_diag
      indx_diag(2) = ix2_diag
      indx_diag(3) = ix3_diag
      indx_diag(4) = ix4_diag
      call sll_o_hdf5_ser_write_array( &
          handle,  &
          indx_diag, &
          'fdistribu_indx_diag', file_err )
      !--> Saving of fdistribu_x1x2
      call sll_o_hdf5_ser_write_array( &
          handle, &
          fdistribu%val4d_seqx1x2(:,:,ix3_diag_loc,ix4_diag_loc), &
          'fdistribu_x1x2', file_err )
      call sll_s_hdf5_ser_file_close( handle, file_err )

    end if

    call sll_s_collective_barrier(sll_v_world_collective)

    !--> HDF5 saving of fdistribu_x3x4
    ix1_min_loc = sll_o_get_layout_i_min( fdistribu%layout4d_seqx3x4, my_rank )
    ix1_max_loc = sll_o_get_layout_i_max( fdistribu%layout4d_seqx3x4, my_rank )
    ix2_min_loc = sll_o_get_layout_j_min( fdistribu%layout4d_seqx3x4, my_rank )
    ix2_max_loc = sll_o_get_layout_j_max( fdistribu%layout4d_seqx3x4, my_rank )

    if ( (ix1_min_loc.le.ix1_diag) .and. (ix1_diag.le.ix1_max_loc) .and. &
        (ix2_min_loc.le.ix2_diag) .and. (ix2_diag.le.ix2_max_loc) ) then

      indx_loc = sll_o_global_to_local( fdistribu%layout4d_seqx3x4, &
          (/ix1_diag,ix2_diag,1,1/) )
      ix1_diag_loc = indx_loc(1)
      ix2_diag_loc = indx_loc(2)

      call sll_s_hdf5_ser_file_create( trim(filename_HDF5), handle, file_err )
      !--> Saving of fdistribu_x1x2
      call sll_o_hdf5_ser_write_array( &
          handle, &
          fdistribu%val4d_seqx3x4(ix1_diag_loc,ix2_diag_loc,:,:), &
          'fdistribu_x3x4', file_err )
      call sll_s_hdf5_ser_file_close( handle, file_err )

    end if

  end subroutine print_fdistribu_DK4D
  !---------------------------------------------------------------------------


end module fdistribu_DK4D_module
!---------------------------------------------------------------------------
