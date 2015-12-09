!------------------------------------------------------------------------------
! SELALIB
!------------------------------------------------------------------------------
! MODULE:   sll_m_ode_integrator_base
!
! DESCRIPTION:
!> @ingroup ode_integrators
!> @authors Yaman Güçlü    - <yaman.guclu@gmail.com>
!> @authors Marco Restelli - <marco.restelli@gmail.com>
!> @brief   Abstract types for: 1) generic ODE system, and 2) ODE integrator.
!> @todo    Add detailed description
!------------------------------------------------------------------------------
module sll_m_ode_integrator_base

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#include "sll_working_precision.h"

  use sll_m_vector_space_base, only: &
    sll_vector_space_base

  implicit none

  public :: &
    sll_ode_base, &
    sll_ode_integrator_base

  private
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  !----------------------------------------------------------------------------
  ! Abstract type: sll_ode_base
  !----------------------------------------------------------------------------
  !> @brief   ODE system
  !> @details Abstract type for ODE systems, to be fed to an ODE integrator.
  !>          ODE systems are in the form \f$ \frac{dy}{dt} = f(t,y) \f$, where
  !>          *y* is the state vector, of <class(sll_vector_space)>, and *t* is
  !>          time, a double precision real number.
  !>          The user of the library should provide a derived ODE class, and 
  !>          possibly a derived vector space for its data.
  !>
  type, abstract :: sll_ode_base
  contains
    procedure( i_rhs ), deferred :: rhs   !< evaluate *f(t,y)*
!    procedure                    :: solve !< solve an implicit problem
  end type sll_ode_base

  !----------------------------------------------------------------------------
  ! Abstract type: sll_ode_integrator_base
  !----------------------------------------------------------------------------
  !> @brief   Base class for standard ODE integrators
  !> @details Abstract type for ODE integrators that do not use specific
  !>          information about the structure of the ODE system (e.g., no 
  !>          splitting and no partitioning is contemplated at this time).
  !>          The user of the library should know about the exposed interface 
  !>          of this class in order to use any integrator.  Additionally,
  !>          specific integrators will have other options and methods.
  !>
  type, abstract :: sll_ode_integrator_base

    !> Pointer to ODE object
    class( sll_ode_base ), pointer :: ode => null()

    !> Storage for intermediate stage computation
    class( sll_vector_space_base ), allocatable :: work(:)

  contains
    procedure( i_init  ), deferred :: init
    procedure( i_step  ), deferred :: step
    procedure( i_clean ), deferred :: clean

  end type sll_ode_integrator_base

  !----------------------------------------------------------------------------
  ! Abstract interface: i_rhs
  !----------------------------------------------------------------------------
  !> @brief        Compute the time derivative of the state vector
  !> @param[inout] self  ODE system, caller
  !> @param[in]    t     Time
  !> @param[in]    y     State vector
  !> @param[inout] ydot  Time derivative of *y*, also a vector
  !>
  abstract interface
   subroutine i_rhs( self, t, y, ydot )
    use sll_m_working_precision
    import sll_ode_base, sll_vector_space_base
    class( sll_ode_base )         , intent( inout ) :: self
    sll_real64                    , intent( in    ) :: t
    class( sll_vector_space_base ), intent( in    ) :: y
    class( sll_vector_space_base ), intent( inout ) :: ydot
   end subroutine i_rhs
  end interface

  !----------------------------------------------------------------------------
  ! Abstract interface: i_init
  !----------------------------------------------------------------------------
  !> @brief        Initialize the time integrator
  !> @note         The actual argument of *ode* must be a pointer or a
  !>               non-pointer target.  The latter option is supported by
  !>               GFortran since version 4.6, but not by IFort as of version
  !>               15.0.  See Sec. 20.14.2 (Automatic pointer targetting) in
  !>               "Modern Fortran Explained" [Metcalf-Reid-Cohen].
  !> @param[out]   self  ODE integrator, caller
  !> @param[in]    ode   ODE system, to which a pointer is stored
  !> @param[in]    t0    Initial time
  !> @param[inout] y0    Initial conditions (state vector)
  !>
  abstract interface
   subroutine i_init( self, ode, t0, y0 )
    use sll_m_working_precision
    import sll_ode_integrator_base, sll_ode_base, sll_vector_space_base
    class( sll_ode_integrator_base ), intent(   out ) :: self
    class( sll_ode_base ), pointer  , intent( in    ) :: ode
    sll_real64                      , intent( in    ) :: t0
    class( sll_vector_space_base )  , intent( inout ) :: y0
   end subroutine i_init
  end interface

  !----------------------------------------------------------------------------
  ! Abstract interface: i_step
  !----------------------------------------------------------------------------
  !> @brief        Advance the solution by one time step
  !> @param[inout] self  ODE integrator, caller
  !> @param[in]    t     Current time
  !> @param[in]    y     Current solution (state vector)
  !> @param[in]    h     Time step size
  !> @param[inout] ynew  New solution, at time *t+h*
  !>
  abstract interface
   subroutine i_step( self, t, y, h, ynew )
    use sll_m_working_precision
    import sll_ode_integrator_base, sll_ode_base, sll_vector_space_base
    class( sll_ode_integrator_base ), intent( inout ) :: self
    sll_real64                      , intent( in    ) :: t
    class( sll_vector_space_base )  , intent( in    ) :: y
    sll_real64                      , intent( in    ) :: h
    class( sll_vector_space_base )  , intent( inout ) :: ynew
   end subroutine i_step
  end interface

  !----------------------------------------------------------------------------
  ! Abstract interface: i_clean
  !----------------------------------------------------------------------------
  !> @brief        Clean up the time integrator
  !> @param[inout] self  ODE integrator, caller
  !>
  abstract interface
   subroutine i_clean( self )
    import sll_ode_integrator_base
    class( sll_ode_integrator_base ), intent( inout ) :: self
   end subroutine i_clean
  end interface

!==============================================================================
!contains
!==============================================================================

!  subroutine solve( self, t, y, ynew, b, sigma )
!    class( sll_ode_base ), intent( inout ) :: self
!    sll_real64           , intent( in    ) :: t
!  end subroutine solve

end module sll_m_ode_integrator_base
