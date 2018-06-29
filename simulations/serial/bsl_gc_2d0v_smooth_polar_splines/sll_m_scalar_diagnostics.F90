module sll_m_scalar_diagnostics
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#include "sll_assert.h"

  use sll_m_working_precision, only: f64

  use sll_m_spline_2d, only: sll_t_spline_2d

  use sll_m_electric_field, only: sll_t_electric_field

  implicit none

  public :: sll_t_scalar_diagnostics

  private
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  ! Working precision
  integer, parameter :: wp = f64

  type :: sll_t_scalar_diagnostics

    integer :: Nk1, Nk2, Nq1, Nq2, file_unit

    type(sll_t_spline_2d)     , pointer :: spline_2d_rho  => null()
    type(sll_t_spline_2d)     , pointer :: spline_2d_phi  => null()
    type(sll_t_electric_field), pointer :: electric_field => null()

    real(wp), pointer :: quad_points_eta1(:,:) => null()
    real(wp), pointer :: quad_points_eta2(:,:) => null()
    real(wp), pointer :: phi_quad_eq(:,:,:,:)  => null()
    real(wp), pointer :: volume(:,:,:,:)       => null()

  contains

    procedure :: init       => s_scalar_diagnostics__init
    procedure :: write_data => s_scalar_diagnostics__write_data
    procedure :: free       => s_scalar_diagnostics__free

  end type sll_t_scalar_diagnostics

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
contains
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine s_scalar_diagnostics__init( &
    self            , &
    file_unit       , &
    spline_2d_rho   , &
    spline_2d_phi   , &
    electric_field  , &
    quad_points_eta1, &
    quad_points_eta2, &
    phi_quad_eq     , &
    volume )
    class(sll_t_scalar_diagnostics)   , intent(inout) :: self
    integer                           , intent(in   ) :: file_unit
    type(sll_t_spline_2d)     , target, intent(in   ) :: spline_2d_rho
    type(sll_t_spline_2d)     , target, intent(in   ) :: spline_2d_phi
    type(sll_t_electric_field), target, intent(in   ) :: electric_field 
    real(wp)                  , target, intent(in   ) :: quad_points_eta1(:,:)
    real(wp)                  , target, intent(in   ) :: quad_points_eta2(:,:)
    real(wp)                  , target, intent(in   ) :: phi_quad_eq(:,:,:,:)
    real(wp)                  , target, intent(in   ) :: volume(:,:,:,:)

    self % file_unit = file_unit

    self % spline_2d_rho => spline_2d_rho
    self % spline_2d_phi => spline_2d_phi

    self % electric_field => electric_field

    self % quad_points_eta1 => quad_points_eta1
    self % quad_points_eta2 => quad_points_eta2

    self % phi_quad_eq => phi_quad_eq
    self % volume      => volume

    self % Nk1 = size( self % quad_points_eta1, 2 )
    self % Nk2 = size( self % quad_points_eta2, 2 )
    self % Nq1 = size( self % quad_points_eta1, 1 )
    self % Nq2 = size( self % quad_points_eta2, 1 )

  end subroutine s_scalar_diagnostics__init

  !-----------------------------------------------------------------------------
  subroutine s_scalar_diagnostics__write_data( self, time )
    class(sll_t_scalar_diagnostics), intent(in) :: self
    real(wp)                       , intent(in) :: time

    integer :: k1, k2, q1, q2

    real(wp) :: mass, energy, l2_norm_phi
    real(wp) :: eta(2), El(2)

    mass        = 0.0_wp
    energy      = 0.0_wp
    l2_norm_phi = 0.0_wp

    associate( Nk1 => self%Nk1, Nk2 => self%Nk2, Nq1 => self%Nq1, Nq2 => self%Nq2 )

      do k2 = 1, Nk2
        do k1 = 1, Nk1
          do q2 = 1, Nq2
            do q1 = 1, Nq1

              eta(1) = self % quad_points_eta1(q1,k1)
              eta(2) = self % quad_points_eta2(q2,k2)

              mass = mass + self % volume(q1,q2,k1,k2) * self % spline_2d_rho % eval( eta(1), eta(2) )

              El     = self % electric_field % eval( eta )
              energy = energy + self % volume(q1,q2,k1,k2) * ( El(1)**2 + El(2)**2 )

              l2_norm_phi = l2_norm_phi + self % volume(q1,q2,k1,k2) * &
                            ( self % spline_2d_phi % eval( eta(1), eta(2) ) - self % phi_quad_eq(q1,q2,k1,k2) )**2

            end do
          end do
        end do
      end do

    end associate

    l2_norm_phi = sqrt( l2_norm_phi )

    write( self % file_unit, '(4g24.15)' ) time, mass, energy, l2_norm_phi
    flush( self % file_unit )

  end subroutine s_scalar_diagnostics__write_data

  !-----------------------------------------------------------------------------
  subroutine s_scalar_diagnostics__free( self )
    class(sll_t_scalar_diagnostics), intent(inout) :: self

    nullify( self % spline_2d_rho  )
    nullify( self % spline_2d_phi  )
    nullify( self % electric_field )

    nullify( self % quad_points_eta1 )
    nullify( self % quad_points_eta2 )

    nullify( self % phi_quad_eq )
    nullify( self % volume      )

  end subroutine s_scalar_diagnostics__free

end module sll_m_scalar_diagnostics
