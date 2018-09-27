module sll_m_diagnostics
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#include "sll_assert.h"

  use sll_m_working_precision, only: f64

  use sll_m_constants, only: sll_p_twopi

  use sll_m_spline_2d, only: sll_t_spline_2d

  use sll_m_electric_field, only: sll_t_electric_field

  use sll_m_polar_mapping_analytical, only: sll_c_polar_mapping_analytical

  use sll_m_polar_mapping_iga, only: sll_t_polar_mapping_iga

  use sll_m_simulation_state, only: sll_t_simulation_state

  use sll_m_gauss_legendre_integration, only: &
    sll_f_gauss_legendre_points, &
    sll_f_gauss_legendre_weights

  use sll_m_hdf5_io_serial, only: &
    sll_t_hdf5_ser_handle, &
    sll_o_hdf5_ser_write_array

  implicit none

  public :: sll_t_diagnostics

  private
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  ! Working precision
  integer, parameter :: wp = f64

  type :: sll_t_diagnostics

    integer :: Nk1, Nk2, Nq1, Nq2, nx1, nx2

    type(sll_t_simulation_state), pointer :: sim_state => null()

    real(wp), allocatable :: tau_eta1(:)
    real(wp), allocatable :: tau_eta2(:)

    real(wp), allocatable :: quad_points_eta1(:,:)
    real(wp), allocatable :: quad_points_eta2(:,:)
    real(wp), allocatable :: quad_weights_eta1(:,:)
    real(wp), allocatable :: quad_weights_eta2(:,:)

    real(wp), allocatable :: phi_quad_eq(:,:,:,:)
    real(wp), allocatable :: volume(:,:,:,:)

    real(wp), allocatable :: x1_grid(:)
    real(wp), allocatable :: x2_grid(:)
    real(wp), allocatable :: x1x2_inverse_grid(:,:,:)

  contains

    procedure :: init                        => s_diagnostics__init
    procedure :: write_scalar_data           => s_diagnostics__write_scalar_data
    procedure :: write_on_interpolation_grid => s_diagnostics__write_on_interpolation_grid
    procedure :: write_on_cartesian_grid     => s_diagnostics__write_on_cartesian_grid
    procedure :: free                        => s_diagnostics__free

  end type sll_t_diagnostics

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
contains
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine s_diagnostics__init( &
    self            , &
    ncells1         , &
    ncells2         , &
    p1              , &
    p2              , &
    nx1             , &
    nx2             , &
    tau_eta1        , &
    tau_eta2        , &
    breaks_eta1     , &
    breaks_eta2     , &
    mapping_discrete, &
    mapping_analytic, &
    sim_state )
    class(sll_t_diagnostics)             , intent(inout) :: self
    integer                              , intent(in   ) :: ncells1
    integer                              , intent(in   ) :: ncells2
    integer                              , intent(in   ) :: p1
    integer                              , intent(in   ) :: p2
    integer                              , intent(in   ) :: nx1
    integer                              , intent(in   ) :: nx2
    real(wp)                             , intent(in   ) :: tau_eta1(:)
    real(wp)                             , intent(in   ) :: tau_eta2(:)
    real(wp)                             , intent(in   ) :: breaks_eta1(:)
    real(wp)                             , intent(in   ) :: breaks_eta2(:)
    type(sll_t_polar_mapping_iga)        , intent(in   ) :: mapping_discrete
    class(sll_c_polar_mapping_analytical), intent(in   ) :: mapping_analytic
    type(sll_t_simulation_state), target , intent(in   ) :: sim_state

    integer  :: i1, i2, k1, k2, q1, q2
    real(wp) :: eta(2), x(2)

    ! Cartesian grid
    self % nx1 = nx1
    self % nx2 = nx2

    ! For quadrature points
    self % Nk1 = ncells1
    self % Nk2 = ncells2
    self % Nq1 = 1 + p1
    self % Nq2 = 1 + p2

    allocate( self % tau_eta1( size( tau_eta1 ) ), source = tau_eta1 )
    allocate( self % tau_eta2( size( tau_eta2 ) ), source = tau_eta2 )

    allocate( self % quad_points_eta1 ( self % Nq1, self % Nk1 ) )
    allocate( self % quad_points_eta2 ( self % Nq2, self % Nk2 ) )
    allocate( self % quad_weights_eta1( self % Nq1, self % Nk1 ) )
    allocate( self % quad_weights_eta2( self % Nq2, self % Nk2 ) )

    allocate( self % volume     ( self % Nq1, self % Nq2, self % Nk1, self % Nk2 ) )
    allocate( self % phi_quad_eq( self % Nq1, self % Nq2, self % Nk1, self % Nk2 ) )

    allocate( self % x1_grid( nx1 ) )
    allocate( self % x2_grid( nx2 ) )
    allocate( self % x1x2_inverse_grid( 2, nx1, nx2 ) )

    ! Quadrature points and weights along s
    do k1 = 1, self % Nk1
      self % quad_points_eta1 (:,k1) = sll_f_gauss_legendre_points ( self % Nq1, breaks_eta1(k1), breaks_eta1(k1+1) )
      self % quad_weights_eta1(:,k1) = sll_f_gauss_legendre_weights( self % Nq1, breaks_eta1(k1), breaks_eta1(k1+1) )
    end do

    ! Quadrature points and weights along theta
    do k2 = 1, self % Nk2
      self % quad_points_eta2 (:,k2) = sll_f_gauss_legendre_points ( self % Nq2, breaks_eta2(k2), breaks_eta2(k2+1) )
      self % quad_weights_eta2(:,k2) = sll_f_gauss_legendre_weights( self % Nq2, breaks_eta2(k2), breaks_eta2(k2+1) )
    end do

    self % sim_state => sim_state

    do k2 = 1, self % Nk2
      do k1 = 1, self % Nk1
        do q2 = 1, self % Nq2
          do q1 = 1, self % Nq1
            eta(1) = self % quad_points_eta1(q1,k1)
            eta(2) = self % quad_points_eta2(q2,k2)
            self % volume(q1,q2,k1,k2) = abs( mapping_discrete % jdet( eta ) ) * self % quad_weights_eta1(q1,k1) * self % quad_weights_eta2(q2,k2)
            self % phi_quad_eq(q1,q2,k1,k2) = self % sim_state % spline_2d_phi % eval( eta(1), eta(2) )
          end do
        end do
      end do
    end do

    ! Write Cartesian grid
    do i2 = 1, nx2
      do i1 = 1, nx1
        self % x1_grid(i1) = 2.0_wp * real( i1-1, wp ) / real( nx1-1, wp ) - 1.0_wp
        self % x2_grid(i2) = 2.0_wp * real( i2-1, wp ) / real( nx2-1, wp ) - 1.0_wp
        x = (/ self % x1_grid(i1), self % x2_grid(i2) /)
        eta(1) = sqrt( x(1)**2 + x(2)**2 )
        eta(2) = modulo( atan2( x(2), x(1) ), sll_p_twopi )
        self % x1x2_inverse_grid(:,i1,i2) = mapping_analytic % eval_inverse( x, eta, tol=1.0e-14_wp, maxiter=100 )
      end do
    end do

  end subroutine s_diagnostics__init

  !-----------------------------------------------------------------------------
  subroutine s_diagnostics__write_scalar_data( self, file_unit, time )
    class(sll_t_diagnostics), intent(in) :: self
    integer                 , intent(in) :: file_unit
    real(wp)                , intent(in) :: time

    integer :: k1, k2, q1, q2

    real(wp) :: mass, energy, l2_norm_phi
    real(wp) :: eta(2), El(2)

    mass        = 0.0_wp
    energy      = 0.0_wp
    l2_norm_phi = 0.0_wp

    associate( Nk1 => self % Nk1, &
               Nk2 => self % Nk2, &
               Nq1 => self % Nq1, &
               Nq2 => self % Nq2, &
               quad_points_eta1 => self % quad_points_eta1, &
               quad_points_eta2 => self % quad_points_eta2, &
               volume           => self % volume          , &
               phi_quad_eq      => self % phi_quad_eq     , &
               spline_2d_rho    => self % sim_state % spline_2d_rho, &
               spline_2d_phi    => self % sim_state % spline_2d_phi, &
               electric_field   => self % sim_state % electric_field )

      do k2 = 1, Nk2
        do k1 = 1, Nk1
          do q2 = 1, Nq2
            do q1 = 1, Nq1

              eta(1) = quad_points_eta1(q1,k1)
              eta(2) = quad_points_eta2(q2,k2)

              mass = mass + volume(q1,q2,k1,k2) * spline_2d_rho % eval( eta(1), eta(2) )

              El     = electric_field % eval( eta )
              energy = energy + volume(q1,q2,k1,k2) * ( El(1)**2 + El(2)**2 )

              l2_norm_phi = l2_norm_phi + volume(q1,q2,k1,k2) * ( spline_2d_phi % eval( eta(1), eta(2) ) - phi_quad_eq(q1,q2,k1,k2) )**2

            end do
          end do
        end do
      end do

    end associate

    l2_norm_phi = sqrt( l2_norm_phi )

    write( file_unit, '(4g24.15)' ) time, mass, energy, l2_norm_phi
    flush( file_unit )

  end subroutine s_diagnostics__write_scalar_data

  !-----------------------------------------------------------------------------
  subroutine s_diagnostics__write_on_interpolation_grid( self, file_id, iteration )
    class(sll_t_diagnostics)   , intent(in) :: self
    type(sll_t_hdf5_ser_handle), intent(in) :: file_id
    integer                    , intent(in) :: iteration

    integer  :: i1, i2, h5_error
    real(wp) :: eta(2), E(2)
    real(wp), allocatable :: phi(:,:), Ex(:,:), Ey(:,:)

    character(len=32) :: attr_name

    associate( ntau1 => size( self % tau_eta1 ), ntau2 => size( self % tau_eta2 ) )

      allocate( phi( ntau1, ntau2+1 ) )
      allocate( Ex ( ntau1, ntau2+1 ) )
      allocate( Ey ( ntau1, ntau2+1 ) )

      ! Write phi, Ex and Ey on interpolation grid
      do i2 = 1, ntau2
        do i1 = 1, ntau1
          eta(1) = self % tau_eta1(i1)
          eta(2) = self % tau_eta2(i2)
          phi(i1,i2) = self % sim_state % spline_2d_phi % eval( eta(1), eta(2) )
          E = self % sim_state % electric_field % eval( eta )
          Ex(i1,i2) = E(1)
          Ey(i1,i2) = E(2)
        end do
      end do

      ! Apply periodicity along theta
      phi(:,ntau2+1) = phi(:,1)
      Ex (:,ntau2+1) = Ex (:,1)
      Ey (:,ntau2+1) = Ey (:,1)

    end associate

    ! Write data to HDF5 file
    write( attr_name, '(a,i0)' ) "/rho_", iteration
    call sll_o_hdf5_ser_write_array( file_id, self % sim_state % rho, trim(attr_name), h5_error )
    write( attr_name, '(a,i0)' ) "/phi_", iteration
    call sll_o_hdf5_ser_write_array( file_id, phi, trim(attr_name), h5_error )
    write( attr_name, '(a,i0)' ) "/Ex_" , iteration
    call sll_o_hdf5_ser_write_array( file_id, Ex, trim(attr_name), h5_error )
    write( attr_name, '(a,i0)' ) "/Ey_" , iteration
    call sll_o_hdf5_ser_write_array( file_id, Ey, trim(attr_name), h5_error )

    deallocate( phi )
    deallocate( Ex  )
    deallocate( Ey  )

  end subroutine s_diagnostics__write_on_interpolation_grid

  !-----------------------------------------------------------------------------
  subroutine s_diagnostics__write_on_cartesian_grid( self, file_id, iteration )
    class(sll_t_diagnostics)   , intent(in) :: self
    type(sll_t_hdf5_ser_handle), intent(in) :: file_id
    integer                    , intent(in) :: iteration

    integer  :: i1, i2, h5_error
    real(wp) :: eta(2), E(2)
    real(wp), allocatable :: Ex_cart(:,:), Ey_cart(:,:)

    character(len=32) :: attr_name

    if ( iteration == 0 ) then
      write( attr_name, '(a)' ) "/x1_cart"
      call sll_o_hdf5_ser_write_array( file_id, self % x1_grid, trim(attr_name), h5_error )
      write( attr_name, '(a)' ) "/x2_cart"
      call sll_o_hdf5_ser_write_array( file_id, self % x2_grid, trim(attr_name), h5_error )
    end if

    associate( nx1 => self % nx1, nx2 => self % nx2 )

      allocate( Ex_cart( nx1, nx2 ) )
      allocate( Ey_cart( nx1, nx2 ) )

      ! Write electric field on Cartesian grid
      do i2 = 1, nx2
        do i1 = 1, nx1
          E = self % sim_state % electric_field % eval( self % x1x2_inverse_grid(:,i1,i2) )
          Ex_cart(i1,i2) = E(1)
          Ey_cart(i1,i2) = E(2)
        end do
      end do

    end associate

    ! Write data to HDF5 file
    write( attr_name, '(a,i0)' ) "/Ex_cart_", iteration
    call sll_o_hdf5_ser_write_array( file_id, Ex_cart, trim(attr_name), h5_error )
    write( attr_name, '(a,i0)' ) "/Ey_cart_", iteration
    call sll_o_hdf5_ser_write_array( file_id, Ey_cart, trim(attr_name), h5_error )

    deallocate( Ex_cart )
    deallocate( Ey_cart )

  end subroutine s_diagnostics__write_on_cartesian_grid

  !-----------------------------------------------------------------------------
  subroutine s_diagnostics__free( self )
    class(sll_t_diagnostics), intent(inout) :: self

    deallocate( self % tau_eta1 )
    deallocate( self % tau_eta2 )

    deallocate( self % quad_points_eta1 )
    deallocate( self % quad_points_eta2 )
    deallocate( self % quad_weights_eta1 )
    deallocate( self % quad_weights_eta2 )

    deallocate( self % phi_quad_eq )
    deallocate( self % volume      )

    deallocate( self % x1_grid )
    deallocate( self % x2_grid )
    deallocate( self % x1x2_inverse_grid )

    nullify( self % sim_state )

  end subroutine s_diagnostics__free

end module sll_m_diagnostics
