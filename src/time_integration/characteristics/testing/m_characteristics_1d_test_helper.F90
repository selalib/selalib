!**************************************************************
!  Copyright INRIA
!  Authors : 
!     CALVI project team
!  
!  This code SeLaLib (for Semi-Lagrangian-Library) 
!  is a parallel library for simulating the plasma turbulence 
!  in a tokamak.
!  
!  This software is governed by the CeCILL-B license 
!  under French law and abiding by the rules of distribution 
!  of free software.  You can  use, modify and redistribute 
!  the software under the terms of the CeCILL-B license as 
!  circulated by CEA, CNRS and INRIA at the following URL
!  "http://www.cecill.info". 
!**************************************************************

! MODULE: m_characteristics_1d_test_helper
!
! DESCRIPTION:
!> @ingroup characteristics
!> @brief   helper module for testing 1D characteristics integrators
!> @author  Yaman Güçlü
!> @detail  this module contains test-case definitions (flow field and 
!>          analytical solution) as well as a subroutine that tests the numerics

module m_characteristics_1d_test_helper
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  use sll_m_working_precision, only: &
    f64

  use sll_m_constants, only: &
    sll_p_pi

  use sll_m_characteristics_1d_base, only: &
    sll_c_characteristics_1d_base

  implicit none

  public :: &
    sll_c_ch1d_test_case, &
    sll_t_ch1d_constant_flow,  &
    sll_t_ch1d_stationary_compression_wave, &
    sll_s_test_characteristics_1d
  
  private
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  integer, parameter :: wp = f64

  !-----------------------------------------------------------------------------
  type, abstract :: sll_c_ch1d_test_case
  contains
    procedure(i_get_domain), deferred :: eta_min
    procedure(i_get_domain), deferred :: eta_max
    procedure(i_flow_field), deferred :: flow_field
    procedure(i_exact_soln), deferred :: exact_soln
  end type

  abstract interface

    function i_get_domain( self ) result( eta )
      import :: sll_c_ch1d_test_case, wp
      class(sll_c_ch1d_test_case), intent(in) :: self
      real(wp) :: eta
    end function i_get_domain

    function i_flow_field( self, eta ) result( a )
      import :: sll_c_ch1d_test_case, wp
      class(sll_c_ch1d_test_case), intent(in) :: self
      real(wp),                    intent(in) :: eta
      real(wp) :: a
    end function i_flow_field

    function i_exact_soln( self, eta0, t ) result( eta )
      import :: sll_c_ch1d_test_case, wp
      class(sll_c_ch1d_test_case), intent(in) :: self
      real(wp),                    intent(in) :: eta0
      real(wp),                    intent(in) :: t
      real(wp) :: eta
    end function i_exact_soln

  end interface

  !-----------------------------------------------------------------------------
  type, extends( sll_c_ch1d_test_case ) :: sll_t_ch1d_constant_flow
    real(wp), private :: eta_domain(2) = [-1.7_wp, 3.3_wp]
    real(wp), private :: a = 1.0_wp
  contains
    procedure :: eta_min    => f_ch1d_cf__eta_min
    procedure :: eta_max    => f_ch1d_cf__eta_max
    procedure :: flow_field => f_ch1d_cf__flow_field
    procedure :: exact_soln => f_ch1d_cf__exact_soln
  end type sll_t_ch1d_constant_flow

  !-----------------------------------------------------------------------------
  type, extends( sll_c_ch1d_test_case ) :: sll_t_ch1d_stationary_compression_wave
    real(wp), private :: eta_domain(2) = [0.0_wp, 2.0_wp*sll_p_pi]
    real(wp), private :: c = 2.0_wp          !> free parameter c>1
    real(wp), private :: b = sqrt( 3.0_wp )  !> b = sqrt( -1+c^2 )
  contains
    procedure :: eta_min    => f_ch1d_scw__eta_min
    procedure :: eta_max    => f_ch1d_scw__eta_max
    procedure :: flow_field => f_ch1d_scw__flow_field
    procedure :: exact_soln => f_ch1d_scw__exact_soln
  end type sll_t_ch1d_stationary_compression_wave

!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
contains
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

  !-----------------------------------------------------------------------------
  ! COMPUTE CHARACTERISTIC WITH GIVEN ALGORITHM AND TEST-CASE
  !-----------------------------------------------------------------------------
  function sll_s_test_characteristics_1d( engine, test_case, Ncells, dt ) &
    result( error )

    class(sll_c_ch1d_test_case),          intent(in) :: test_case
    class(sll_c_characteristics_1d_base), intent(in) :: engine
    integer,                              intent(in) :: Ncells
    real(wp),                             intent(in) :: dt
    real(wp) :: error

    integer               :: i
    integer               :: Npts
    real(wp)              :: eta_min
    real(wp)              :: eta_max
    real(wp), allocatable :: A(:)
    real(wp), allocatable :: eta(:)
    real(wp), allocatable :: output(:)
    real(wp), allocatable :: output_ex(:)

    ! Allocate arrays
    Npts = Ncells+1
    allocate(      A   (Npts) )
    allocate(    eta   (Npts) ) 
    allocate( output   (Npts) )
    allocate( output_ex(Npts) )

    ! Computational grid = initial conditions for characteristics
    eta_min = test_case % eta_min()
    eta_max = test_case % eta_max()
    do i=1,Npts
      eta(i) = eta_min + (eta_max-eta_min) * real(i-1,wp)/real(Ncells,wp)
      A(i) = test_case % flow_field( eta(i) )
      output_ex(i) = test_case % exact_soln( eta(i), dt )
    end do

    ! Numerical solution 
    call engine%compute_characteristics( &
      A      =   A, &
      dt     = -dt, &
      input  = eta, &
      output = output )

    ! Compute L-infinity norm of error - mod is used because domain is periodic
    error = maxval( mod( output-output_ex, eta_max-eta_min ) )

    ! Deallocate arrays
    deallocate( A )
    deallocate( eta )
    deallocate( output )
    deallocate( output_ex )

  end function sll_s_test_characteristics_1d

  !-----------------------------------------------------------------------------
  ! CONSTANT ADVECTION
  !-----------------------------------------------------------------------------
  function f_ch1d_cf__eta_min( self ) result( eta_min )
    class(sll_t_ch1d_constant_flow), intent(in) :: self
    real(wp) :: eta_min
    eta_min = self%eta_domain(1)
  end function f_ch1d_cf__eta_min

  function f_ch1d_cf__eta_max( self ) result( eta_max )
    class(sll_t_ch1d_constant_flow), intent(in) :: self
    real(wp) :: eta_max
    eta_max = self%eta_domain(2)
  end function f_ch1d_cf__eta_max

  function f_ch1d_cf__flow_field( self, eta ) result( a )
    class(sll_t_ch1d_constant_flow), intent(in) :: self
    real(wp),                        intent(in) :: eta
    real(wp) :: a
    a = self%a
  end function f_ch1d_cf__flow_field

  function f_ch1d_cf__exact_soln( self, eta0, t ) result( eta )
    class(sll_t_ch1d_constant_flow), intent(in) :: self
    real(wp),                        intent(in) :: eta0
    real(wp),                        intent(in) :: t
    real(wp) :: eta

    ! Exact solution
    eta = eta0 + self%a * t

    ! Apply periodic boundary conditions
    eta = eta - self%eta_min()
    eta = mod( eta, self%eta_max()-self%eta_min() )
    eta = eta + self%eta_min()
  end function f_ch1d_cf__exact_soln

  !-----------------------------------------------------------------------------
  ! STATIONARY COMPRESSION WAVE
  !-----------------------------------------------------------------------------
  function f_ch1d_scw__eta_min( self ) result( eta_min )
    class(sll_t_ch1d_stationary_compression_wave), intent(in) :: self
    real(wp) :: eta_min
    eta_min = self%eta_domain(1)
  end function f_ch1d_scw__eta_min

  function f_ch1d_scw__eta_max( self ) result( eta_max )
    class(sll_t_ch1d_stationary_compression_wave), intent(in) :: self
    real(wp) :: eta_max
    eta_max = self%eta_domain(2)
  end function f_ch1d_scw__eta_max

  function f_ch1d_scw__flow_field( self, eta ) result( a )
    class(sll_t_ch1d_stationary_compression_wave), intent(in) :: self
    real(wp),                                      intent(in) :: eta
    real(wp) :: a
    a = self%c + sin( eta )
  end function f_ch1d_scw__flow_field

  function f_ch1d_scw__exact_soln( self, eta0, t ) result( eta )
    class(sll_t_ch1d_stationary_compression_wave), intent(in) :: self
    real(wp),                                      intent(in) :: eta0
    real(wp),                                      intent(in) :: t
    real(wp) :: eta

    ! Exact solution (periodic on [0,2pi])
    associate( c => self%c, b => self%b )
      eta = 2.0_wp*atan( b/c*tan( atan( (1.0_wp+c*tan( 0.5_wp*eta0 ))/b ) &
        +0.5_wp*b*t )-1.0_wp/c )
    end associate
  end function f_ch1d_scw__exact_soln

end module m_characteristics_1d_test_helper
