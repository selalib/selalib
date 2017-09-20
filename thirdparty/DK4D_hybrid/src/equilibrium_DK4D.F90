!===========================================================================
!> Initialization of the equilibrium for 4D drift-kinetic hybrid
!>  simulation
!>
!> \date 2014-08-19
!> \author V. Grandgirard
!---------------------------------------------------------------------------
module equilibrium_DK4D_module
#include "sll_working_precision.h"
#include "sll_memory.h"

  use input_DK4D_module
  use mesh_DK4D_module
  use sll_m_boundary_condition_descriptors

  implicit none

  !===========================================================================
  !> 
  !---------------------------------------------------------------------------
  type, public :: equilibrium_DK4D_t

    !> Profile parameters
    sll_real64 :: r_peak
    sll_real64 :: inv_Ln
    sll_real64 :: inv_LTi
    sll_real64 :: inv_LTe
    sll_real64 :: deltarn
    sll_real64 :: deltarTi
    sll_real64 :: deltarTe

    !> Density and temperature profiles
    sll_real64, dimension(:,:), pointer :: n0_xy
    sll_real64, dimension(:,:), pointer :: Ti_xy
    sll_real64, dimension(:,:), pointer :: Te_xy

    sll_real64, dimension(:)  , pointer :: n0_r
    sll_real64, dimension(:)  , pointer :: Ti_r
    sll_real64, dimension(:)  , pointer :: Te_r

    !> Equilibrium distribution function
    sll_real64, dimension(:,:,:), pointer :: feq_xyvpar

  end type equilibrium_DK4D_t
  !---------------------------------------------------------------------------

  contains


  !===========================================================================
  !> Equilibrium: Allocation 
  !---------------------------------------------------------------------------
  subroutine new_equilibrium_DK4D( equilibrium, mesh4d ) 

    type(equilibrium_DK4D_t), intent(inout) :: equilibrium
    type(mesh_DK4D_t)       , intent(in)    :: mesh4d

    sll_int32 :: Nx, Ny, Nvpar
    sll_int32 :: Nr
    sll_int32 :: ierr

    !*** Array allocations ***
    !-> For density and temperature profiles
    Nx = size(mesh4d%xgrid_2d,1)
    Ny = size(mesh4d%xgrid_2d,2)
    SLL_ALLOCATE( equilibrium%n0_xy(Nx,Ny), ierr )
    SLL_ALLOCATE( equilibrium%Ti_xy(Nx,Ny), ierr )
    SLL_ALLOCATE( equilibrium%Te_xy(Nx,Ny), ierr )
    
    Nr = Nx
    SLL_ALLOCATE( equilibrium%n0_r(Nr), ierr )
    SLL_ALLOCATE( equilibrium%Ti_r(Nr), ierr )
    SLL_ALLOCATE( equilibrium%Te_r(Nr), ierr )
    
    !-> For equilibrium distribution function
    Nvpar = size(mesh4d%vpar_grid)
    SLL_ALLOCATE( equilibrium%feq_xyvpar(Nx,Ny,Nvpar), ierr )

  end subroutine new_equilibrium_DK4D
  !---------------------------------------------------------------------------


  !===========================================================================
  !> Equilibrium: Deallocation 
  !---------------------------------------------------------------------------
  subroutine delete_equilibrium_DK4D( equilibrium )

    type(equilibrium_DK4D_t), intent(inout) :: equilibrium
    
    sll_int32 :: ierr

    SLL_DEALLOCATE( equilibrium%n0_xy, ierr )
    SLL_DEALLOCATE( equilibrium%Ti_xy, ierr )
    SLL_DEALLOCATE( equilibrium%Te_xy, ierr )

    SLL_DEALLOCATE( equilibrium%n0_r, ierr )
    SLL_DEALLOCATE( equilibrium%Ti_r, ierr )
    SLL_DEALLOCATE( equilibrium%Te_r, ierr )

    SLL_DEALLOCATE( equilibrium%feq_xyvpar, ierr )

  end subroutine delete_equilibrium_DK4D
  !---------------------------------------------------------------------------


  !===========================================================================
  !> Equilibrium: Initialization
  !---------------------------------------------------------------------------
  subroutine init_equilibrium_DK4D( equilibrium, &
      input_data, mesh4d )

    type(equilibrium_DK4D_t), intent(inout) :: equilibrium
    type(input_DK4D_t)      , intent(in)    :: input_data
    type(mesh_DK4D_t)       , intent(in)    :: mesh4d

    sll_real64 :: r_min, r_max, Lr, R0

    !*** Initialization of the density and temperature profiles in (x,y) ***
    R0      = mesh4d%major_radius
    r_min   = minval(mesh4d%rgrid_2d(:,1))
    r_max   = maxval(mesh4d%rgrid_2d(:,1))
    Lr      = abs(r_max-r_min)

    equilibrium%r_peak  = r_min + input_data%rho_peak * Lr

    equilibrium%inv_Ln  = input_data%kappan/R0
    equilibrium%inv_LTi = input_data%kappaTi/R0
    equilibrium%inv_LTe = input_data%kappaTe/R0

    equilibrium%deltarn  = input_data%deltarn*Lr
    equilibrium%deltarTi = input_data%deltarTi*Lr
    equilibrium%deltarTe = input_data%deltarTe*Lr

    call initialize_nT_profiles( equilibrium, mesh4d, &
        equilibrium%r_peak, &
        equilibrium%inv_Ln, &
        equilibrium%deltarn, &
        equilibrium%inv_LTi, &
        equilibrium%deltarTi, &
        equilibrium%inv_LTe, &
        equilibrium%deltarTe )

    !*** Initialization of the equilibrium distribution function ***
    call compute_fequilibrium_xy( &
        mesh4d%xgrid_2d, &
        mesh4d%ygrid_2d, &
        mesh4d%vpar_grid, &
        equilibrium%n0_xy, &
        equilibrium%Ti_xy, &
        equilibrium%feq_xyvpar )

  end subroutine init_equilibrium_DK4D


  !===========================================================================
  !> Equilibrium: Print in HDF5
  !---------------------------------------------------------------------------
  subroutine print_equilibrium_DK4D( equilibrium, filename )

    use sll_m_collective
    use sll_m_hdf5_io_serial, only: sll_s_hdf5_ser_file_create, &
         sll_o_hdf5_ser_write_array, &
         sll_s_hdf5_ser_file_close, &
         sll_t_hdf5_ser_handle

    type(equilibrium_DK4D_t), intent(in) :: equilibrium
    character(len=*)        , intent(in) :: filename

    sll_int32 :: my_rank

    !--> For  HDF5 saving
    integer   :: file_err
    type(sll_t_hdf5_ser_handle) :: handle    !< file handle

    my_rank = sll_f_get_collective_rank(sll_v_world_collective)
    if (my_rank.eq.0) then
      call sll_s_hdf5_ser_file_create( trim(filename), handle, file_err )
      call sll_o_hdf5_ser_write_array( handle, equilibrium%n0_r, &
          'n0_r', file_err )
      call sll_o_hdf5_ser_write_array( handle, equilibrium%Ti_r, &
          'Ti_r', file_err )
      call sll_o_hdf5_ser_write_array( handle, equilibrium%Te_r, &
          'Te_r', file_err )
      call sll_o_hdf5_ser_write_array( handle, equilibrium%n0_xy, &
          'n0_xy', file_err )
      call sll_o_hdf5_ser_write_array( handle, equilibrium%Ti_xy, &
          'Ti_xy', file_err )
      call sll_o_hdf5_ser_write_array( handle, equilibrium%Te_xy, &
          'Te_xy', file_err )
      call sll_o_hdf5_ser_write_array( handle, equilibrium%feq_xyvpar, &
          'feq_xyvpar', file_err )
      call sll_s_hdf5_ser_file_close( handle, file_err )
    end if

  end subroutine print_equilibrium_DK4D


  !===========================================================================
  !> Initialization of the density and temperature profiles 
  !---------------------------------------------------------------------------
  subroutine initialize_nT_profiles( equilibrium, mesh4d, &
      r_peak, &
      inv_Ln, deltarn, &
      inv_LTi, deltarTi, &
      inv_LTe, deltarTe )

    use sll_m_arbitrary_degree_spline_interpolator_1d

    type(equilibrium_DK4D_t), intent(inout) :: equilibrium
    type(mesh_DK4D_t)       , intent(in)    :: mesh4d
    sll_real64, intent(in) :: r_peak
    sll_real64, intent(in) :: inv_Ln
    sll_real64, intent(in) :: deltarn
    sll_real64, intent(in) :: inv_LTi
    sll_real64, intent(in) :: deltarTi
    sll_real64, intent(in) :: inv_LTe
    sll_real64, intent(in) :: deltarTe

    !-> local variables 
    !---> for profile definition
    sll_real64, dimension(3) :: params_n0
    sll_real64, dimension(3) :: params_Te
    sll_real64, dimension(3) :: params_Ti
    !--> for (x,y) profiles
    sll_int32  :: ix, iy, Nx, Ny
    sll_real64 :: x_point, y_point
    !--> for radial profiles
    sll_int32  :: ir, Nr
    sll_real64 :: r_point
    sll_real64, dimension(:), pointer :: r_grid
    !--> for normalization of n0(r)
    sll_real64 :: r_min, r_max
    sll_real64 :: n0_norm
    !--> for interpolation
    type(sll_t_arbitrary_degree_spline_interpolator_1d) :: interp1d_n0_r
    type(sll_t_arbitrary_degree_spline_interpolator_1d) :: interp1d_Ti_r
    type(sll_t_arbitrary_degree_spline_interpolator_1d) :: interp1d_Te_r
    
    !*** parameters for radial density profiles ***
    params_n0(1) = inv_Ln
    params_n0(2) = deltarn
    params_n0(3) = r_peak

    !--> parameter for radial temperature ions profiles
    params_Ti(1) = inv_LTi
    params_Ti(2) = deltarTi
    params_Ti(3) = r_peak

    !--> parameters for radial temperature electrons profiles
    params_Te(1) = inv_LTe
    params_Te(2) = deltarTe
    params_Te(3) = r_peak

    !*** Initialization of the radial profiles ****
    equilibrium%n0_r(:) = 0.0_f64
    equilibrium%Ti_r(:) = 0.0_f64
    equilibrium%Te_r(:) = 0.0_f64
    r_grid => mesh4d%rgrid_2d(:,1)
    Nr = size(r_grid)
    do ir = 1,Nr
      r_point  = r_grid(ir)
      equilibrium%n0_r(ir) = compute_profile_r( r_point, params_n0 )
      equilibrium%Ti_r(ir) = compute_profile_r( r_point, params_Ti ) 
      equilibrium%Te_r(ir) = compute_profile_r( r_point, params_Te )
    end do
    !--> Normalization of n0_r
    r_min   = r_grid(1)
    r_max   = r_grid(Nr)
    n0_norm = 0.5_f64 * (equilibrium%n0_r(1)*r_min + &
        equilibrium%n0_r(Nr)*r_max)
    do ir = 2,Nr-1
      r_point = r_grid(ir)
      n0_norm = n0_norm + equilibrium%n0_r(ir)*r_point
    end do
    n0_norm          = n0_norm/real(Nr-1,f64)
    equilibrium%n0_r = equilibrium%n0_r/n0_norm

    !*** Initialization of the profiles in (x,y) ****
    !--> Initialization of the interpolators
    call interp1d_n0_r%init( &
        Nr, &
        r_min, &
        r_max, &
        bc_left=sll_p_dirichlet, &
        bc_right=sll_p_dirichlet, &
        spline_degree=3 )
    call interp1d_Ti_r%init( &
        Nr, &
        r_min, &
        r_max, &
        bc_left=sll_p_hermite, &
        bc_right=sll_p_hermite, &
        spline_degree=3 )
    call interp1d_Te_r%init( &
        Nr, &
        r_min, &
        r_max, &
        bc_left=sll_p_hermite, &
        bc_right=sll_p_hermite, &
        spline_degree=3 )

    call interp1d_n0_r%compute_interpolants( &
        equilibrium%n0_r )
    call interp1d_Ti_r%compute_interpolants( &
        equilibrium%Ti_r )
    call interp1d_Te_r%compute_interpolants( &
        equilibrium%Te_r )

    !--> Computation of the profiles by interpolation
    Nx = size(mesh4d%xgrid_2d,1)
    Ny = size(mesh4d%xgrid_2d,2)
    do ix = 1, Nx
      do iy = 1, Ny 
        x_point = mesh4d%xgrid_2d(ix,iy)
        y_point = mesh4d%ygrid_2d(ix,iy)
        r_point = sqrt(x_point**2+y_point**2)
        r_point = min(max(r_point,r_min),r_max)
        equilibrium%n0_xy(ix,iy) = &
            interp1d_n0_r%interpolate_from_interpolant_value( r_point )
        equilibrium%Ti_xy(ix,iy) = &
            interp1d_Ti_r%interpolate_from_interpolant_value( r_point )
        equilibrium%Te_xy(ix,iy) = &
            interp1d_Te_r%interpolate_from_interpolant_value( r_point )

      end do
    end do

    !--> Destruction of the interpolators
    call sll_o_delete( interp1d_n0_r )
    call sll_o_delete( interp1d_Ti_r )
    call sll_o_delete( interp1d_Te_r )

  end subroutine initialize_nT_profiles
  !---------------------------------------------------------------------------


  !===========================================================================
  !> Initialization of the profiles solution of the following equation
  !>
  !>       1   d T(r)    = - kappa cosh^(-2)( (r- rx)/delta ) 
  !>      T(r) dr
  !>
  !>  kappa = params(1)
  !>  delta = params(2)
  !>  rx    = params(3)
  !---------------------------------------------------------------------------
  function compute_profile_r( r, params_profil )

    sll_real64, intent(in) :: r
    sll_real64, dimension(:), intent(in) :: params_profil

    sll_real64 :: compute_profile_r
    sll_real64 :: kappa
    sll_real64 :: delta
    sll_real64 :: rx

    kappa = params_profil(1)
    delta = params_profil(2)
    rx    = params_profil(3)
    compute_profile_r = exp(-kappa*delta*tanh((r-rx)/delta))

  end function compute_profile_r
  !---------------------------------------------------------------------------


  !===========================================================================
  ! Initialization of the 3D array for the equilibrium
  !  distribution function feq(x,y,vpar)
  !  feq(x,y,vpar) = n0(x,y)/(2*pi*Ti(x,y))**(1/2) * 
  !                    exp(-0.5*vpar**2/Ti(x,y))
  !---------------------------------------------------------------------------
  subroutine compute_fequilibrium_xy( xgrid_2d, ygrid_2d, &
      vpar_grid, n0_xy, Ti_xy, feq_xyvpar )

    use sll_m_constants, only : sll_p_pi
    
    sll_real64, dimension(:,:)  , intent(in)  :: xgrid_2d
    sll_real64, dimension(:,:)  , intent(in)  :: ygrid_2d
    sll_real64, dimension(:)    , intent(in)  :: vpar_grid
    sll_real64, dimension(:,:)  , intent(in)  :: n0_xy
    sll_real64, dimension(:,:)  , intent(in)  :: Ti_xy
    sll_real64, dimension(:,:,:), intent(out) :: feq_xyvpar

    sll_int32 :: Npt1, Npt2, Nvpar
    sll_int32 :: ix, iy, ivpar

    Npt1  = size(xgrid_2d,1)
    Npt2  = size(xgrid_2d,2)
    Nvpar = size(vpar_grid,1)

    do ivpar = 1,Nvpar
      do iy = 1,Npt2
        do ix = 1,Npt1
          feq_xyvpar(ix,iy,ivpar) = n0_xy(ix,iy) * &
              exp(-0.5_f64*vpar_grid(ivpar)**2/Ti_xy(ix,iy)) / &
              sqrt(2._f64*sll_p_pi*Ti_xy(ix,iy))
        end do
      end do
    end do

  end subroutine compute_fequilibrium_xy
  !---------------------------------------------------------------------------

end module equilibrium_DK4D_module
!---------------------------------------------------------------------------
