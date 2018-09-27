module sll_m_diagnostics
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#include "sll_assert.h"

  use sll_m_working_precision, only: f64

  use sll_m_spline_2d, only: sll_t_spline_2d

  use sll_m_electric_field, only: sll_t_electric_field

  use sll_m_polar_mapping_iga, only: sll_t_polar_mapping_iga

  use sll_m_simulation_state, only: sll_t_simulation_state

  use sll_m_gauss_legendre_integration, only: &
    sll_f_gauss_legendre_points, &
    sll_f_gauss_legendre_weights

  implicit none

  public :: sll_t_diagnostics

  private
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  ! Working precision
  integer, parameter :: wp = f64

  type :: sll_t_diagnostics

    integer :: Nk1, Nk2, Nq1, Nq2, file_unit

    type(sll_t_simulation_state), pointer :: sim_state => null()

    real(wp), allocatable :: quad_points_eta1(:,:)
    real(wp), allocatable :: quad_points_eta2(:,:)
    real(wp), allocatable :: quad_weights_eta1(:,:)
    real(wp), allocatable :: quad_weights_eta2(:,:)

    real(wp), allocatable :: phi_quad_eq(:,:,:,:)
    real(wp), allocatable :: volume(:,:,:,:)

  contains

    procedure :: init       => s_diagnostics__init
    procedure :: write_data => s_diagnostics__write_data
    procedure :: free       => s_diagnostics__free

  end type sll_t_diagnostics

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
contains
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine s_diagnostics__init( &
    self            , &
    file_unit       , &
    ncells1         , &
    ncells2         , &
    p1              , &
    p2              , &
    breaks_eta1     , &
    breaks_eta2     , &
    mapping_discrete, &
    sim_state )
    class(sll_t_diagnostics)            , intent(inout) :: self
    integer                             , intent(in   ) :: file_unit
    integer                             , intent(in   ) :: ncells1
    integer                             , intent(in   ) :: ncells2
    integer                             , intent(in   ) :: p1
    integer                             , intent(in   ) :: p2
    real(wp)                            , intent(in   ) :: breaks_eta1(:)
    real(wp)                            , intent(in   ) :: breaks_eta2(:)
    type(sll_t_polar_mapping_iga)       , intent(in   ) :: mapping_discrete
    type(sll_t_simulation_state), target, intent(in   ) :: sim_state

    integer  :: k1, k2, q1, q2
    real(wp) :: eta(2)

    self % file_unit = file_unit

    ! For quadrature points
    self % Nk1 = ncells1
    self % Nk2 = ncells2
    self % Nq1 = 1 + p1
    self % Nq2 = 1 + p2

    self % sim_state => sim_state

    allocate( self % quad_points_eta1 ( self % Nq1, self % Nk1 ) )
    allocate( self % quad_points_eta2 ( self % Nq2, self % Nk2 ) )
    allocate( self % quad_weights_eta1( self % Nq1, self % Nk1 ) )
    allocate( self % quad_weights_eta2( self % Nq2, self % Nk2 ) )

    allocate( self % volume     ( self % Nq1, self % Nq2, self % Nk1, self % Nk2 ) )
    allocate( self % phi_quad_eq( self % Nq1, self % Nq2, self % Nk1, self % Nk2 ) )

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

  end subroutine s_diagnostics__init

  !-----------------------------------------------------------------------------
  subroutine s_diagnostics__write_data( self, time )
    class(sll_t_diagnostics), intent(in) :: self
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

              l2_norm_phi = l2_norm_phi + volume(q1,q2,k1,k2) * &
                            ( spline_2d_phi % eval( eta(1), eta(2) ) - phi_quad_eq(q1,q2,k1,k2) )**2

            end do
          end do
        end do
      end do

    end associate

    l2_norm_phi = sqrt( l2_norm_phi )

    write( self % file_unit, '(4g24.15)' ) time, mass, energy, l2_norm_phi
    flush( self % file_unit )

  end subroutine s_diagnostics__write_data

  !-----------------------------------------------------------------------------
  subroutine s_diagnostics__free( self )
    class(sll_t_diagnostics), intent(inout) :: self

    deallocate( self % quad_points_eta1 )
    deallocate( self % quad_points_eta2 )
    deallocate( self % quad_weights_eta1 )
    deallocate( self % quad_weights_eta2 )

    deallocate( self % phi_quad_eq )
    deallocate( self % volume      )

    nullify( self % sim_state )

  end subroutine s_diagnostics__free

end module sll_m_diagnostics
