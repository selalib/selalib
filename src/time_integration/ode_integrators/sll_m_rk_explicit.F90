!----------------------
! A note about notation
!----------------------
!  rkXY:
!       e = explicit
!       i = implicit (fully)
!       d = diagonally implicit
!       m = embedded
!       l = low-storage
!----------------------

module sll_m_rk_explicit

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#include "sll_assert.h"
#include "sll_working_precision.h"

  use sll_m_ode_integrator_base, only: &
    sll_ode_base, &
    sll_ode_integrator_base

  use sll_m_vector_space_base, only: &
    sll_vector_space_base

  implicit none

  public :: &
    sll_rk1e_fwd_euler, &
    sll_rk2e_heun, &
    sll_rk2e_midpoint, &
    sll_rk2e_ralston, &
    sll_rk3e_heun3, &
    sll_rk4e_classic

  private
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  
  !----------------------------------------------------------------------------
  ! Explicit Euler method
  !----------------------------------------------------------------------------
  type, extends( sll_ode_integrator_base ) :: sll_rk1e_fwd_euler

  contains
    procedure :: init  =>  init__rk1e_fwd_euler
    procedure :: step  =>  step__rk1e_fwd_euler
    procedure :: clean => clean__rk1e_fwd_euler

  end type sll_rk1e_fwd_euler

  !----------------------------------------------------------------------------
  ! Explicit midpoint method
  !----------------------------------------------------------------------------
  type, extends( sll_ode_integrator_base ) :: sll_rk2e_midpoint

  contains
    procedure :: init  =>  init__rk2e_midpoint
    procedure :: step  =>  step__rk2e_midpoint
    procedure :: clean => clean__rk2e_midpoint

  end type sll_rk2e_midpoint

  !----------------------------------------------------------------------------
  ! Explicit trapezoidal rule (Heun's method)
  !----------------------------------------------------------------------------
  type, extends( sll_ode_integrator_base ) :: sll_rk2e_heun

  contains
    procedure :: init  =>  init__rk2e_heun
    procedure :: step  =>  step__rk2e_heun
    procedure :: clean => clean__rk2e_heun

  end type sll_rk2e_heun

  !----------------------------------------------------------------------------
  ! Ralston's method
  !----------------------------------------------------------------------------
  type, extends( sll_ode_integrator_base ) :: sll_rk2e_ralston

  contains
    procedure :: init  =>  init__rk2e_ralston
    procedure :: step  =>  step__rk2e_ralston
    procedure :: clean => clean__rk2e_ralston

  end type sll_rk2e_ralston

  !----------------------------------------------------------------------------
  ! 3rd-order method by Heun
  !----------------------------------------------------------------------------
  type, extends( sll_ode_integrator_base ) :: sll_rk3e_heun3

  contains
    procedure :: init  =>  init__rk3e_heun3
    procedure :: step  =>  step__rk3e_heun3
    procedure :: clean => clean__rk3e_heun3

  end type sll_rk3e_heun3

  !----------------------------------------------------------------------------
  ! The 'classic' 4th-order RK
  !----------------------------------------------------------------------------
  type, extends( sll_ode_integrator_base ) :: sll_rk4e_classic

  contains
    procedure :: init  =>  init__rk4e_classic
    procedure :: step  =>  step__rk4e_classic
    procedure :: clean => clean__rk4e_classic

  end type sll_rk4e_classic

!==============================================================================
contains
!==============================================================================

  !----------------------------------------------------------------------------
  ! Explicit Euler method
  !----------------------------------------------------------------------------
  subroutine init__rk1e_fwd_euler( self, ode, t0, y0 )
    class( sll_rk1e_fwd_euler  )  , intent(   out ) :: self
    class( sll_ode_base ), pointer, intent( in    ) :: ode
    sll_real64                    , intent( in    ) :: t0
    class( sll_vector_space_base ), intent( inout ) :: y0

    self%ode => ode
    SLL_ASSERT(t0>=0.0)
    SLL_ASSERT(storage_size(y0)>0)

  end subroutine init__rk1e_fwd_euler

  !----------------------------------------------------------------------------
  subroutine step__rk1e_fwd_euler( self, t, y, h, ynew )
    class( sll_rk1e_fwd_euler  )  , intent( inout ) :: self
    sll_real64                    , intent( in    ) :: t
    class( sll_vector_space_base ), intent( in    ) :: y
    sll_real64                    , intent( in    ) :: h
    class( sll_vector_space_base ), intent( inout ) :: ynew

    ! ynew = y + h * f(t,y)
    call self%ode%rhs( t, y, ynew ) ! ynew  = f(t,y)
    call ynew%scal( h )             ! ynew *= h
    call ynew%incr( y )             ! ynew += y
    
  end subroutine step__rk1e_fwd_euler

  !----------------------------------------------------------------------------
  subroutine clean__rk1e_fwd_euler( self )
    class( sll_rk1e_fwd_euler ), intent( inout ) :: self
    ! Do nothing, because method does not require local storage
    SLL_ASSERT(storage_size(self)>0)
  end subroutine clean__rk1e_fwd_euler

  !----------------------------------------------------------------------------
  ! Explicit midpoint method
  !----------------------------------------------------------------------------
  subroutine init__rk2e_midpoint( self, ode, t0, y0 )
    class( sll_rk2e_midpoint )    , intent(   out ) :: self
    class( sll_ode_base ), pointer, intent( in    ) :: ode
    sll_real64                    , intent( in    ) :: t0
    class( sll_vector_space_base ), intent( inout ) :: y0

    self%ode => ode                 ! Store pointer to ODE system
    call y0%source( self%work, 1 )  ! Allocate temporary storage
    SLL_ASSERT(t0>=0.0)

  end subroutine init__rk2e_midpoint

  !----------------------------------------------------------------------------
  subroutine step__rk2e_midpoint( self, t, y, h, ynew )
    class( sll_rk2e_midpoint )    , intent( inout ) :: self
    sll_real64                    , intent( in    ) :: t
    class( sll_vector_space_base ), intent( in    ) :: y
    sll_real64                    , intent( in    ) :: h
    class( sll_vector_space_base ), intent( inout ) :: ynew

    ! Half time step
    sll_real64 :: h2
    h2 = h*0.5_f64
    
    ! yi   = y + h/2 * f(t    ,y )
    ! ynew = y + h   * f(t+h/2,yi)

    ! Compute 1st stage derivative:  k1 = f(t,y)
    call self%ode%rhs( t, y, self%work(1) )

    ! Compute 2nd stage solution:    y2 = y + h/2 * k1
    call ynew%mult_add( h2, self%work(1), y )       ! (use ynew as temporary)

    ! Compute 2nd stage derivative:  k2 = f(t+h/2,y2) 
    call self%ode%rhs( t+h2, ynew, self%work(1) )   ! (overwrite k1)

    ! Compute new solution:  ynew = y + h * k2
    call ynew%mult_add( h, self%work(1), y )

  end subroutine step__rk2e_midpoint

  !----------------------------------------------------------------------------
  subroutine clean__rk2e_midpoint( self )
    class( sll_rk2e_midpoint), intent( inout ) :: self
    deallocate( self%work )
  end subroutine clean__rk2e_midpoint

  !----------------------------------------------------------------------------
  ! Explicit trapezoidal rule (Heun's method)
  !----------------------------------------------------------------------------
  subroutine init__rk2e_heun( self, ode, t0, y0 )
    class( sll_rk2e_heun )        , intent(   out ) :: self
    class( sll_ode_base ), pointer, intent( in    ) :: ode
    sll_real64                    , intent( in    ) :: t0
    class( sll_vector_space_base ), intent( inout ) :: y0

    self%ode => ode                 ! Store pointer to ODE system
    call y0%source( self%work, 2 )  ! Allocate temporary storage
    SLL_ASSERT(t0>=0.0)

  end subroutine init__rk2e_heun

  !----------------------------------------------------------------------------
  subroutine step__rk2e_heun( self, t, y, h, ynew )
    class( sll_rk2e_heun )        , intent( inout ) :: self
    sll_real64                    , intent( in    ) :: t
    class( sll_vector_space_base ), intent( in    ) :: y
    sll_real64                    , intent( in    ) :: h
    class( sll_vector_space_base ), intent( inout ) :: ynew

    ! Compute 1st stage derivative:  k1 = f(t,y)
    call self%ode%rhs( t, y, self%work(1) )

    ! Compute 2nd stage solution:    y2 = y + h * k1
    call ynew%mult_add( h, self%work(1), y )       ! (use ynew as temporary)

    ! Compute 2nd stage derivative:  k2 = f(t+h,y2)
    call self%ode%rhs( t+h, ynew, self%work(2) )

    ! Compute average time derivative, multiplied by 2:
    ! 2*f_avg = k1+k2
    call self%work(1)%incr( self%work(2) )         ! (overwrite k1)

    ! Compute new solution:  ynew = y + h/2 * (k1+k2)
    call ynew%mult_add( h*0.5_f64, self%work(1), y )

  end subroutine step__rk2e_heun

  !----------------------------------------------------------------------------
  subroutine clean__rk2e_heun( self )
    class( sll_rk2e_heun ), intent( inout ) :: self
    deallocate( self%work )
  end subroutine clean__rk2e_heun

  !----------------------------------------------------------------------------
  ! Ralston's method
  !----------------------------------------------------------------------------
  subroutine init__rk2e_ralston( self, ode, t0, y0 )
    class( sll_rk2e_ralston )     , intent(   out ) :: self
    class( sll_ode_base ), pointer, intent( in    ) :: ode
    sll_real64                    , intent( in    ) :: t0
    class( sll_vector_space_base ), intent( inout ) :: y0

    self%ode => ode                 ! Store pointer to ODE system
    call y0%source( self%work, 2 )  ! Allocate temporary storage

    SLL_ASSERT(t0>=0.0)
  end subroutine init__rk2e_ralston

  !----------------------------------------------------------------------------
  subroutine step__rk2e_ralston( self, t, y, h, ynew )
    class( sll_rk2e_ralston )     , intent( inout ) :: self
    sll_real64                    , intent( in    ) :: t
    class( sll_vector_space_base ), intent( in    ) :: y
    sll_real64                    , intent( in    ) :: h
    class( sll_vector_space_base ), intent( inout ) :: ynew

    sll_real64, parameter :: two_thirds = 2.0_f64/3.0_f64

    ! Compute 1st stage derivative:  k1 = f(t,y)
    call self%ode%rhs( t, y, self%work(1) )

    ! Compute 2nd stage solution:    y2 = y + 2/3*h * k1
    call ynew%mult_add( h*two_thirds, self%work(1), y )   ! (use ynew as temp.)

    ! Compute 2nd stage derivative:  k2 = f(t+2/3*h,y2)
    call self%ode%rhs( t+h*two_thirds, ynew, self%work(2) )

    ! Compute average time derivative, multiplied by 4:
    ! 4*f_avg = k1+3*k2
    call self%work(1)%incr_mult( 3.0_f64, self%work(2) )  ! (overwrite k1)

    ! Compute new solution:  ynew = y + h/4 * (k1+3*k2)
    call ynew%mult_add( h*0.25_f64, self%work(1), y )

  end subroutine step__rk2e_ralston

  !----------------------------------------------------------------------------
  subroutine clean__rk2e_ralston( self )
    class( sll_rk2e_ralston ), intent( inout ) :: self
    deallocate( self%work )
  end subroutine clean__rk2e_ralston

  !----------------------------------------------------------------------------
  ! 3rd-order method by Heun
  !----------------------------------------------------------------------------
  subroutine init__rk3e_heun3( self, ode, t0, y0 )
    class( sll_rk3e_heun3 )       , intent(   out ) :: self
    class( sll_ode_base ), pointer, intent( in    ) :: ode
    sll_real64                    , intent( in    ) :: t0
    class( sll_vector_space_base ), intent( inout ) :: y0

    self%ode => ode                 ! Store pointer to ODE system
    call y0%source( self%work, 2 )  ! Allocate temporary storage

    SLL_ASSERT(t0>=0.0)

  end subroutine init__rk3e_heun3

  !----------------------------------------------------------------------------
  subroutine step__rk3e_heun3( self, t, y, h, ynew )
    class( sll_rk3e_heun3 )       , intent( inout ) :: self
    sll_real64                    , intent( in    ) :: t
    class( sll_vector_space_base ), intent( in    ) :: y
    sll_real64                    , intent( in    ) :: h
    class( sll_vector_space_base ), intent( inout ) :: ynew

    ! Compute 1st stage derivative:  k1 = f(t,y)
    call self%ode%rhs( t, y, self%work(1) )

    ! Compute 2nd stage solution:    y2 = y + h/3 * k1
    call ynew%mult_add( h/3.0_f64, self%work(1), y )  ! (use ynew as temporary)

    ! Compute 2nd stage derivative:  k2 = f(t+h/3,y2)
    call self%ode%rhs( t+h/3.0_f64, ynew, self%work(2) )

    ! Compute 3rd stage solution:    y3 = y + 2/3*h *k2
    call ynew%mult_add( h/1.5_f64, self%work(2), y )  ! (use ynew as temporary)

    ! Compute 3rd stage derivative:  k3 = f(t+2/3*h,y3)
    call self%ode%rhs( t+h/1.5_f64, ynew, self%work(2) )  ! (overwrite k2)

    ! Compute average time derivative, multiplied by 4:
    ! 4*f_avg = k1+3*k3
    call self%work(1)%incr_mult( 3.0_f64, self%work(2) )  ! (overwrite k1)

    ! Compute new solution:  ynew = y + h/4 * (k1+3*k2)
    call ynew%mult_add( h*0.25_f64, self%work(1), y )

  end subroutine step__rk3e_heun3

  !----------------------------------------------------------------------------
  subroutine clean__rk3e_heun3( self )
    class( sll_rk3e_heun3 ), intent( inout ) :: self
    deallocate( self%work )
  end subroutine clean__rk3e_heun3

  !----------------------------------------------------------------------------
  ! The 'classic' 4th-order RK
  !----------------------------------------------------------------------------
  subroutine init__rk4e_classic( self, ode, t0, y0 )
    class( sll_rk4e_classic )     , intent(   out ) :: self
    class( sll_ode_base ), pointer, intent( in    ) :: ode
    sll_real64                    , intent( in    ) :: t0
    class( sll_vector_space_base ), intent( inout ) :: y0

    self%ode => ode                 ! Store pointer to ODE system
    call y0%source( self%work, 2 )  ! Allocate temporary storage
    SLL_ASSERT(t0>=0.0)

  end subroutine init__rk4e_classic

  !----------------------------------------------------------------------------
  subroutine step__rk4e_classic( self, t, y, h, ynew )
    class( sll_rk4e_classic )     , intent( inout ) :: self
    sll_real64                    , intent( in    ) :: t
    class( sll_vector_space_base ), intent( in    ) :: y
    sll_real64                    , intent( in    ) :: h
    class( sll_vector_space_base ), intent( inout ) :: ynew

    ! Butcher's table
    sll_real64, parameter :: f3=1.0_f64/3.0_f64, f6=1.0_f64/6.0_f64
    sll_real64, parameter :: &
      c2=0.5_f64, a21=0.5_f64,                                   &
      c3=0.5_f64,              a32=0.5_f64,                      &
      c4=1.0_f64,                           a43=1.0_f64,         &
                   b1=f6     ,  b2=f3     ,  b3=f3     , b4=f6

    ! Stage 1
    call self%ode%rhs ( t, y, self%work(2) )
    call ynew%mult_add( h*b1, self%work(2), y )

    ! Stage 2
    call self%work(1)%mult_add( h*a21, self%work(2), y )
    call self%ode%rhs( t+h*c2, self%work(1), self%work(2) )
    call ynew%incr_mult( h*b2, self%work(2) )

    ! Stage 3
    call self%work(1)%mult_add( h*a32, self%work(2), y )
    call self%ode%rhs( t+h*c3, self%work(1), self%work(2) )
    call ynew%incr_mult( h*b3, self%work(2) )

    ! Stage 4
    call self%work(1)%mult_add( h*a43, self%work(2), y )
    call self%ode%rhs( t+h*c4, self%work(1), self%work(2) )
    call ynew%incr_mult( h*b4, self%work(2) )

  end subroutine step__rk4e_classic

  !----------------------------------------------------------------------------
  subroutine clean__rk4e_classic( self )
    class( sll_rk4e_classic ), intent( inout ) :: self
    deallocate( self%work )
  end subroutine clean__rk4e_classic

end module sll_m_rk_explicit
