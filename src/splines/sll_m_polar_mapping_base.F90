module sll_m_polar_mapping_base
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#include "sll_assert.h"
#include "sll_errors.h"

  use sll_m_working_precision, only: f64

  use sll_m_constants, only: sll_p_twopi

  use sll_m_hdf5_io_serial, only: sll_t_hdf5_ser_handle

  implicit none

  public :: sll_c_polar_mapping

  private
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  !> Working precision
  integer, parameter :: wp = f64

  !> Abstract  type, polar mapping
  type, abstract :: sll_c_polar_mapping

  contains

    ! Deferred procedures
    procedure(i_fun_eval)      , deferred :: eval
    procedure(i_fun_jmat)      , deferred :: jmat ! Jacobian matrix
    procedure(i_sub_store_data), deferred :: store_data

    ! Non-deferred procedures
    procedure :: jdet         => f_polar_mapping__jdet ! Jacobian determinant
    procedure :: eval_inverse => f_polar_mapping__eval_inverse ! inverse mapping (Newton method)

  end type sll_c_polar_mapping

  ! Interfaces for deferred procedures
  abstract interface

    SLL_PURE function i_fun_eval( self, eta ) result( x )
      import sll_c_polar_mapping, wp
      class(sll_c_polar_mapping), intent(in) :: self
      real(wp)                  , intent(in) :: eta(2)
      real(wp) :: x(2)
    end function i_fun_eval

    SLL_PURE function i_fun_jmat( self, eta ) result( jmat )
      import sll_c_polar_mapping, wp
      class(sll_c_polar_mapping), intent(in) :: self
      real(wp)                  , intent(in) :: eta(2)
      real(wp) :: jmat(2,2)
    end function i_fun_jmat

    subroutine i_sub_store_data( self, n1, n2, file_id )
      import sll_c_polar_mapping, sll_t_hdf5_ser_handle
      class(sll_c_polar_mapping) , intent(in) :: self
      integer                    , intent(in) :: n1
      integer                    , intent(in) :: n2
      type(sll_t_hdf5_ser_handle), intent(in) :: file_id
    end subroutine i_sub_store_data

  end interface

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
contains
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  SLL_PURE function f_polar_mapping__jdet( self, eta ) result( jdet )
    class(sll_c_polar_mapping), intent(in) :: self
    real(wp)                  , intent(in) :: eta(2)
    real(wp) :: jdet

    real(wp) :: jmat(2,2)

    jmat = self % jmat( eta )

    jdet = jmat(1,1) * jmat(2,2) - jmat(1,2) * jmat(2,1)

  end function f_polar_mapping__jdet

  !-----------------------------------------------------------------------------
  function f_polar_mapping__eval_inverse( self, x, eta0, tol, maxiter ) result( eta )
    class(sll_c_polar_mapping), intent(in) :: self
    real(wp)                  , intent(in) :: x(2)
    real(wp)                  , intent(in) :: eta0(2)
    real(wp)                  , intent(in) :: tol
    integer                   , intent(in) :: maxiter
    real(wp) :: eta(2)

    real(wp) :: jdet
    real(wp) :: eta_new(2), x_new(2), temp(2), delx(2)
    real(wp) :: jmat(2,2)
    integer  :: i

    ! To handle error
    character(len=*), parameter :: this_fun_name = "sll_c_polar_mapping % eval_inverse"
    character(len=64) :: err_msg

    ! First guess:
    ! - use circle inverse mapping for pole of singularity
    ! - use given coordinates for all other points in domain
    if ( abs( eta0(1) ) < 1.0e-14_wp ) then
      eta_new(1) = sqrt( x(1)**2 + x(2)**2 )
      eta_new(2) = modulo( atan2( x(2), x(1) ), sll_p_twopi )
    else
      eta_new = eta0
    end if

    ! Newton method
    do i = 1, maxiter

      ! Function to be inverted
      temp = self % eval( eta_new ) - x

      ! Jacobian matrix
      jmat = self % jmat( eta_new )

      ! Jacobian determinant (needed for inverse of Jacobian matrix)
      jdet = self % jdet( eta_new )

      ! Apply inverse of Jacobian matrix
      delx(1) = (   jmat(2,2) * temp(1) - jmat(1,2) * temp(2) ) / jdet
      delx(2) = ( - jmat(2,1) * temp(1) + jmat(1,1) * temp(2) ) / jdet

      ! New guess
      eta = eta_new - delx

      ! Apply periodicity for eta2
      eta(2) = modulo( eta(2), sll_p_twopi )

      ! Deal with points out of domain
      if ( eta(1) < 0.0_wp ) then; eta(1) = 1.0e-14_wp;         end if
      if ( eta(1) > 1.0_wp ) then; eta(1) = 1.0_wp    ; return; end if

      ! Check convergence using Euclidean distance
      x_new = self % eval( eta )
      if ( norm2( x_new - x ) < tol ) return

      eta_new = eta

    end do

    ! Newton method did not converge
    write(err_msg,'(a,i0,a)') "Newton method did not converge after ", maxiter, " iterations"
    SLL_ERROR( this_fun_name, err_msg )

  end function f_polar_mapping__eval_inverse

end module sll_m_polar_mapping_base
