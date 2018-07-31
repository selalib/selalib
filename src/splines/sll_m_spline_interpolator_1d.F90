!> @ingroup splines
!> @brief   Interpolator for 1D splines of arbitrary degree,
!>          on uniform and non-uniform grids
!> @author  Yaman Güçlü  - IPP Garching
!> @author  Edoardo Zoni - IPP Garching

module sll_m_spline_interpolator_1d
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#include "sll_assert.h"
#include "sll_errors.h"

  use sll_m_working_precision, only: f64

  use sll_m_boundary_condition_descriptors, only: &
    sll_p_periodic, &
    sll_p_hermite , &
    sll_p_greville

  use sll_m_bsplines_base, only: sll_c_bsplines
  use sll_m_spline_1d, only: sll_t_spline_1d

  use sll_m_spline_matrix, only: &
    sll_c_spline_matrix, &
    sll_s_spline_matrix_new

  implicit none

  public :: &
    sll_t_spline_interpolator_1d, &
    sll_s_spline_1d_compute_num_cells

  private
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  !> Working precision
  integer, parameter :: wp = f64

  !> Allowed boundary conditions
  integer, parameter :: &
    allowed_bcs(1:3) = [sll_p_periodic, sll_p_hermite, sll_p_greville]

  !-----------------------------------------------------------------------------
  !> 1D spline interpolator
  !-----------------------------------------------------------------------------
  type :: sll_t_spline_interpolator_1d

    ! Private attributes
    class(sll_c_bsplines),          pointer, private :: bspl => null()
    integer                                , private ::  bc_xmin
    integer                                , private ::  bc_xmax
    integer                                , private :: nbc_xmin
    integer                                , private :: nbc_xmax
    integer                                , private :: odd
    integer                                , private :: offset
    real(wp)                               , private :: dx
    real(wp),                   allocatable, private :: tau(:)
    class(sll_c_spline_matrix), allocatable, private :: matrix

  contains

    procedure :: init                => s_spline_interpolator_1d__init
    procedure :: free                => s_spline_interpolator_1d__free
    procedure :: get_interp_points   => s_spline_interpolator_1d__get_interp_points
    procedure :: compute_interpolant => s_spline_interpolator_1d__compute_interpolant

  end type sll_t_spline_interpolator_1d

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
contains
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  !-----------------------------------------------------------------------------
  !> @brief      Calculate number of cells from number of interpolation points
  !> @details    Important for parallelization: for given spline degree and BCs,
  !>             calculate the number of grid cells that yields the desired
  !>             number of interpolation points
  !>
  !> @param[in]  degree   spline degree
  !> @param[in]  bc_xmin  boundary condition type at left  boundary (x=xmin)
  !> @param[in]  bc_xmax  boundary condition type at right boundary (x=xmax)
  !> @param[in]  nipts    desired number of interpolation points
  !> @param[out] ncells   calculated number of cells in domain
  !-----------------------------------------------------------------------------
  subroutine sll_s_spline_1d_compute_num_cells( &
      degree , &
      bc_xmin, &
      bc_xmax, &
      nipts  , &
      ncells )

    integer, intent(in   ) :: degree
    integer, intent(in   ) :: bc_xmin
    integer, intent(in   ) :: bc_xmax
    integer, intent(in   ) :: nipts
    integer, intent(  out) :: ncells

    ! Sanity checks
    SLL_ASSERT( degree > 0  )
    SLL_ASSERT( any( bc_xmin == allowed_bcs ) )
    SLL_ASSERT( any( bc_xmax == allowed_bcs ) )

    if (any( [bc_xmin,bc_xmax]==sll_p_periodic ) .and. bc_xmin /= bc_xmax) then
      SLL_ERROR( "sll_s_spline_1d_compute_num_cells", "Incompatible BCs" )
    end if

    if (bc_xmin == sll_p_periodic) then
      ncells = nipts
    else
      associate( nbc_xmin => merge( degree/2, 0, bc_xmin == sll_p_hermite ), &
                 nbc_xmax => merge( degree/2, 0, bc_xmax == sll_p_hermite ) )
        ncells = nipts + nbc_xmin + nbc_xmax - degree
      end associate
    end if

  end subroutine sll_s_spline_1d_compute_num_cells

  !-----------------------------------------------------------------------------
  !> @brief      Initialize a 1D spline interpolator object
  !> @param[out] self     1D spline interpolator
  !> @param[in]  bspl     B-splines (basis)
  !> @param[in]  bc_xmin  boundary condition at xmin
  !> @param[in]  bc_xmax  boundary condition at xmax
  !-----------------------------------------------------------------------------
  subroutine s_spline_interpolator_1d__init( self, bspl, bc_xmin, bc_xmax )

    class(sll_t_spline_interpolator_1d), intent(  out) :: self
    class(sll_c_bsplines),       target, intent(in   ) :: bspl
    integer                            , intent(in   ) :: bc_xmin
    integer                            , intent(in   ) :: bc_xmax

    integer :: kl, ku
    character(len=32) :: matrix_type

    ! Sanity checks
    SLL_ASSERT( any( bc_xmin == allowed_bcs ) )
    SLL_ASSERT( any( bc_xmax == allowed_bcs ) )
    if (bspl % periodic) then
      SLL_ASSERT( bc_xmin == sll_p_periodic )
      SLL_ASSERT( bc_xmax == sll_p_periodic )
    end if

    ! Check that these are not radial B-splines on polar grid
    SLL_ASSERT( .not. bspl % radial )

    ! Save pointer to B-splines
    ! (later needed to verify 1D spline input to 'compute_interpolant')
    self % bspl => bspl

    ! Save other data
    self % bc_xmin  = bc_xmin
    self % bc_xmax  = bc_xmax
    self % nbc_xmin = merge ( bspl%degree/2, 0, bc_xmin == sll_p_hermite )
    self % nbc_xmax = merge ( bspl%degree/2, 0, bc_xmax == sll_p_hermite )
    self % odd      = modulo( bspl%degree  , 2 )
    self % offset   = merge ( bspl%degree/2, 0, bspl%periodic )

    ! Save average cell size for normalization of derivatives
    self%dx = (bspl%xmax - bspl%xmin) / real(bspl%ncells,f64)

    ! Compute interpolation points and number of diagonals in linear system
    if (bspl % uniform) then
      call s_compute_interpolation_points_uniform( self, self%tau )
      call s_compute_num_diags_uniform( self, kl, ku )
    else
      call s_compute_interpolation_points_non_uniform( self, self%tau )
      call s_compute_num_diags_non_uniform( self, kl, ku )
    end if

    ! Special case: linear spline
    ! No need for matrix assembly
    if (self % bspl % degree == 1) return

    ! Determine matrix storage type and allocate matrix
    associate( nbasis => self % bspl % nbasis )

      if (bc_xmin == sll_p_periodic) then
        if (kl+1+ku >= nbasis) then
          matrix_type = "dense"
        else
          matrix_type = "periodic_banded"
        end if
      else
        matrix_type = "banded"
      end if

      call sll_s_spline_matrix_new( self%matrix, matrix_type, nbasis, kl, ku )

    end associate ! nbasis

    ! Fill in entries A_ij of linear system A*x = b for interpolation
    call s_build_system( self, self%matrix )

    ! Factorize matrix A to speed-up solution for any given b
    call self % matrix % factorize()

  end subroutine s_spline_interpolator_1d__init

  !-----------------------------------------------------------------------------
  !> @brief        Destroy local objects and free allocated memory
  !> @param[inout] self  1D spline interpolator
  !-----------------------------------------------------------------------------
  subroutine s_spline_interpolator_1d__free( self )

    class(sll_t_spline_interpolator_1d), intent(inout) :: self

    nullify   ( self % bspl )
    deallocate( self % tau  )

    if (allocated( self % matrix )) then
      call self % matrix % free()
      deallocate( self % matrix )
    end if

  end subroutine s_spline_interpolator_1d__free

  !-----------------------------------------------------------------------------
  !> @brief      Get coordinates of interpolation points (1D grid)
  !> @param[in]  self  1D spline interpolator
  !> @param[out] tau   x coordinates of interpolation points
  !-----------------------------------------------------------------------------
  subroutine s_spline_interpolator_1d__get_interp_points( self, tau )

    class(sll_t_spline_interpolator_1d), intent(in   ) :: self
    real(wp),               allocatable, intent(  out) :: tau(:)

    SLL_ASSERT( allocated( self%tau ) )
    allocate( tau(size(self%tau)), source=self%tau )

  end subroutine s_spline_interpolator_1d__get_interp_points

  !-----------------------------------------------------------------------------
  !> @brief        Compute interpolating 1D spline
  !> @details      Compute coefficients of 1D spline that interpolates function
  !>               values on grid. If Hermite BCs are used, function derivatives
  !>               at appropriate boundary are also needed.
  !>
  !> @param[in]    self        1D spline interpolator
  !> @param[inout] spline      1D spline
  !> @param[in]    gtau        function values of interpolation points
  !> @param[in]    derivs_xmin (optional) array with boundary conditions at xmin
  !> @param[in]    derivs_xmax (optional) array with boundary conditions at xmax
  !-----------------------------------------------------------------------------
  subroutine s_spline_interpolator_1d__compute_interpolant( self, &
      spline, gtau, derivs_xmin, derivs_xmax )

    class(sll_t_spline_interpolator_1d), intent(in   ) :: self
    type (sll_t_spline_1d)             , intent(inout) :: spline
    real(wp)                           , intent(in   ) :: gtau(:)
    real(wp),                  optional, intent(in   ) :: derivs_xmin(:)
    real(wp),                  optional, intent(in   ) :: derivs_xmax(:)

    character(len=*), parameter :: &
      this_sub_name = "sll_t_spline_interpolator_1d % compute_interpolant"

    integer :: i

    SLL_ASSERT( size(gtau) == self%bspl%nbasis - self%nbc_xmin - self%nbc_xmax )
    SLL_ASSERT( spline % belongs_to_space( self % bspl ) )

    associate ( nbasis   =>   self % bspl % nbasis, &
                degree   =>   self % bspl % degree, &
                nbc_xmin =>   self % nbc_xmin     , &
                nbc_xmax =>   self % nbc_xmax     , &
                g        =>   self % offset       , &
                bcoef    => spline % bcoef )

      ! Special case: linear spline
      if (degree == 1) then
        bcoef(1:nbasis) = gtau(1:nbasis)
        ! Periodic only: "wrap around" coefficients onto extended array
        if (self%bspl%periodic) then
          bcoef(nbasis+1) = bcoef(1)
        end if
        return
      end if

      ! Hermite boundary conditions at xmin, if any
      ! NOTE: For consistency with the linear system, the i-th derivative
      !       provided by the user must be multiplied by dx^i
      if ( self%bc_xmin == sll_p_hermite ) then
        if ( present( derivs_xmin ) ) then
          bcoef(1:nbc_xmin) = &
                     [(derivs_xmin(i)*self%dx**(i+self%odd-1), i=nbc_xmin,1,-1)]
        else
          SLL_ERROR(this_sub_name,"Hermite BC at xmin requires derivatives")
        end if
      end if

      ! Interpolation points
      bcoef(nbc_xmin+1+g:nbasis-nbc_xmax+g) = gtau(:)

      ! Hermite boundary conditions at xmax, if any
      ! NOTE: For consistency with the linear system, the i-th derivative
      !       provided by the user must be multiplied by dx^i
      if ( self%bc_xmax == sll_p_hermite ) then
        if ( present( derivs_xmax ) ) then
          bcoef(nbasis-nbc_xmax+1:nbasis) = &
                        [(derivs_xmax(i)*self%dx**(i+self%odd-1), i=1,nbc_xmax)]
        else
          SLL_ERROR(this_sub_name,"Hermite BC at xmax requires derivatives")
        end if
      end if

      ! Solve linear system and compute coefficients
      call self % matrix % solve_inplace( bcoef(1+g:nbasis+g) )

      ! Periodic only: "wrap around" coefficients onto extended array
      if (self%bc_xmin == sll_p_periodic) then
        bcoef(1:g)                      = bcoef(nbasis+1:nbasis+g)
        bcoef(nbasis+1+g:nbasis+degree) = bcoef(1+g:degree)
      end if

    end associate

  end subroutine s_spline_interpolator_1d__compute_interpolant

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!
  !!                            PRIVATE SUBROUTINES
  !!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  ! TODO: move subroutine to init()
  !-----------------------------------------------------------------------------
  !> @brief        Private subroutine for assembling and factorizing linear
  !>               system for any kind of boundary conditions at xmin and xmax
  !> @param[in]    self   1D spline interpolator object
  !> @param[inout] matrix generic 'spline' matrix (dense/banded/periodic-banded)
  !-----------------------------------------------------------------------------
  subroutine s_build_system( self, matrix )

    class(sll_t_spline_interpolator_1d), intent(in   ) :: self
    class(sll_c_spline_matrix)         , intent(inout) :: matrix

    integer  :: i,j,d,s
    integer  :: j0,d0
    integer  :: order
    real(wp) :: x
    real(wp) :: values(self%bspl%degree+1)
    real(wp), allocatable :: derivs(:,:)
#ifdef __PGI
    real(wp), allocatable :: h(:)
#endif

    ! NEW
    integer :: jmin

    associate( nbasis   => self % bspl % nbasis, &
               degree   => self % bspl % degree, &
               nbc_xmin => self % nbc_xmin     , &
               nbc_xmax => self % nbc_xmax )

      if ( any( [self%bc_xmin,self%bc_xmax] == sll_p_hermite ) ) then
        allocate ( derivs (0:degree/2, 1:degree+1) )
      end if

      ! Hermite boundary conditions at xmin, if any
      if ( self%bc_xmin == sll_p_hermite ) then
        x = self%bspl%xmin
        call self % bspl % eval_basis_and_n_derivs( x, nbc_xmin, derivs, jmin )

        ! In order to improve the condition number of the matrix, we normalize all
        ! derivatives by multiplying the i-th derivative by dx^i

#ifdef __PGI
        allocate(h(ubound(derivs,1)))
        h = [(self%dx**i, i=1, ubound(derivs,1))]
        do j = lbound(derivs,2), ubound(derivs,2)
          derivs(1:,j) = derivs(1:,j) * h(1:)
        end do
        deallocate(h)
#else
        associate( h => [(self%dx**i, i=1, ubound(derivs,1))] )
          do j = lbound(derivs,2), ubound(derivs,2)
            derivs(1:,j) = derivs(1:,j) * h(1:)
          end do
        end associate
#endif

        do i = 1, nbc_xmin
          ! iterate only to deg as last bspline is 0
          order = nbc_xmin-i+self%odd
          do j = 1, degree
            call matrix % set_element( i, j, derivs(order,j) )
          end do
        end do

      end if

      ! Interpolation points
      do i = nbc_xmin+1, nbasis-nbc_xmax
        x = self%tau(i-nbc_xmin)
        call self % bspl % eval_basis( x, values, jmin )
        do s = 1, degree+1
          j = modulo( jmin-self%offset+s-2, nbasis ) + 1
          call matrix % set_element( i, j, values(s) )
        end do
      end do

      ! Hermite boundary conditions at xmax, if any
      if ( self%bc_xmax == sll_p_hermite ) then
        x = self%bspl%xmax
        call self % bspl % eval_basis_and_n_derivs( x, nbc_xmax, derivs, jmin )

        ! In order to improve the condition number of the matrix, we normalize all
        ! derivatives by multiplying the i-th derivative by dx^i
#ifdef __PGI
        allocate(h(ubound(derivs,1)))
        h = [(self%dx**i, i=1, ubound(derivs,1))]
        do j = lbound(derivs,2), ubound(derivs,2)
          derivs(1:,j) = derivs(1:,j) * h(1:)
        end do
        deallocate(h)
#else
        associate( h => [(self%dx**i, i=1, ubound(derivs,1))] )
          do j = lbound(derivs,2), ubound(derivs,2)
            derivs(1:,j) = derivs(1:,j) * h(1:)
          end do
        end associate
#endif

        do i = nbasis-nbc_xmax+1, nbasis
          order = i-(nbasis-nbc_xmax+1)+self%odd
          j0 = nbasis-degree
          d0 = 1
          do s = 1, degree
            j = j0 + s
            d = d0 + s
            call matrix % set_element( i, j, derivs(order,d) )
          end do
        end do

      end if

      if ( allocated( derivs ) ) deallocate( derivs )

    end associate ! nbasis, degree

  end subroutine s_build_system

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!
  !!                       PRIVATE SUBROUTINES, UNIFORM
  !!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !-----------------------------------------------------------------------------
  subroutine s_compute_interpolation_points_uniform( self, tau )
    class(sll_t_spline_interpolator_1d), intent(in   ) :: self
    real(wp),               allocatable, intent(  out) :: tau(:)

    integer :: i, ntau, isum
    integer, allocatable :: iknots(:)

    associate( nbasis   => self % bspl % nbasis, &
               ncells   => self % bspl % ncells, &
               degree   => self % bspl % degree, &
               xmin     => self % bspl % xmin  , &
               xmax     => self % bspl % xmax  , &
               dx       => self % dx           , &
               bc_xmin  => self %  bc_xmin     , &
               bc_xmax  => self %  bc_xmax     , &
               nbc_xmin => self % nbc_xmin     , &
               nbc_xmax => self % nbc_xmax )

      ! Determine size of tau and allocate tau
      ntau = nbasis - nbc_xmin - nbc_xmax
      allocate( tau(1:ntau) )

      if ( bc_xmin == sll_p_periodic ) then

      ! Periodic case:
      !   . for odd  degree, interpolation points are breakpoints (last excluded)
      !   . for even degree, interpolation points are cell centers
      associate( shift => merge( 0.5_wp, 0.0_wp, self%odd == 0 ) )
        tau = [(xmin + (real(i-1,wp)+shift)*dx, i=1,ntau)]
      end associate

      else

      ! Non-periodic case: create array of temporary knots (integer shifts only)
      ! in order to compute interpolation points using Greville-style averaging:
      ! tau(i) = xmin + average(knots(i+1-degree:i)) * dx
      allocate( iknots (2-degree:ntau) )

      ! Additional knots near x=xmin
      associate( r => 2-degree, s => -nbc_xmin )
        select case (bc_xmin)
          case (sll_p_greville); iknots(r:s) = 0
          case (sll_p_hermite ); iknots(r:s) = [(i,i=r-s-1,-1)]
        end select
      end associate

      ! Knots inside the domain
      associate( r => -nbc_xmin+1, s => -nbc_xmin+1+ncells )
        iknots(r:s) = [(i,i=0,ncells)]
      end associate

      ! Additional knots near x=xmax
      associate( r => -nbc_xmin+1+ncells+1, s => ntau )
        select case (bc_xmax)
          case (sll_p_greville); iknots(r:s) = ncells
          case (sll_p_hermite ); iknots(r:s) = [(i,i=ncells+1,ncells+1+s-r)]
        end select
      end associate

      ! Compute interpolation points using Greville-style averaging
      associate( inv_deg => 1.0_wp / real( degree, wp ) )
        do i = 1, ntau
          isum = sum( iknots(i+1-degree:i) )
          if (modulo( isum, degree ) == 0) then
            tau(i) = xmin + real(isum/degree,wp) * dx
          else
            tau(i) = xmin + real(isum,wp) * inv_deg * dx
          end if
        end do
      end associate

      ! Non-periodic case, odd degree: fix round-off issues
      if ( self%odd == 1 ) then
        tau(1)    = xmin
        tau(ntau) = xmax
      end if

      deallocate( iknots )

      end if  ! (bc_xmin == sll_p_periodic)

    end associate

  end subroutine s_compute_interpolation_points_uniform

  !-----------------------------------------------------------------------------
  subroutine s_compute_num_diags_uniform( self, kl, ku )
    class(sll_t_spline_interpolator_1d), intent(in   ) :: self
    integer                            , intent(  out) :: kl
    integer                            , intent(  out) :: ku

    associate( degree => self % bspl % degree )

      ! FIXME: In Hermite case ku and kl computed in general case when derivatives
      !        of B-splines do not vanish at boundary
      ! TODO: reduce kl and ku to take advantage of uniform grid
      select case( self % bc_xmin )
        case ( sll_p_periodic ); ku = (degree+1)/2
        case ( sll_p_hermite  ); ku = max( (degree+1)/2, degree-1 )
        case ( sll_p_greville ); ku = degree
      end select

      select case( self % bc_xmax )
        case ( sll_p_periodic ); kl = (degree+1)/2
        case ( sll_p_hermite  ); kl = max( (degree+1)/2, degree-1 )
        case ( sll_p_greville ); kl = degree
      end select

    end associate

  end subroutine s_compute_num_diags_uniform

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!
  !!                       PRIVATE SUBROUTINES, NON-UNIFORM
  !!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !-----------------------------------------------------------------------------
  subroutine s_compute_interpolation_points_non_uniform( self, tau )
    class(sll_t_spline_interpolator_1d), intent(in   ) :: self
    real(wp),               allocatable, intent(  out) :: tau(:)

   integer :: i, ntau
   real(wp), allocatable :: temp_knots(:)

   associate( nbasis   => self % bspl % nbasis, &
              ncells   => self % bspl % ncells, &
              degree   => self % bspl % degree, &
              xmin     => self % bspl % xmin  , &
              xmax     => self % bspl % xmax  , &
              dx       => self % dx           , &
              bc_xmin  => self %  bc_xmin     , &
              bc_xmax  => self %  bc_xmax     , &
              nbc_xmin => self % nbc_xmin     , &
              nbc_xmax => self % nbc_xmax )

   associate( breaks   => self % bspl % knots(1:ncells+1) )

      ! Determine size of tau and allocate tau
      ntau = nbasis - nbc_xmin - nbc_xmax
      allocate( tau(1:ntau) )

      ! Array of temporary knots needed to compute interpolation points
      ! using Greville-style averaging: tau(i) = average(temp_knots(i+1-degree:i))
      allocate( temp_knots( 2-degree : ntau ) )

      if ( bc_xmin == sll_p_periodic ) then

        associate( k => degree/2 )
          temp_knots(:) = self%bspl%knots(2-degree+k:ntau+k)
        end associate

      else

        associate( r => 2-degree, s => -nbc_xmin )
          select case (bc_xmin)
          case (sll_p_greville); temp_knots(r:s) = breaks(1)
          case (sll_p_hermite ); temp_knots(r:s) = 2.0_wp*breaks(1) - breaks(2+s-r:2:-1)
          end select
        end associate

        associate( r => -nbc_xmin+1, s => -nbc_xmin+1+ncells )
          temp_knots(r:s) = breaks(:)
        end associate

        associate( r => -nbc_xmin+1+ncells+1, s => ntau )
          select case (bc_xmax)
          case (sll_p_greville)
            temp_knots(r:s) = breaks(ncells+1)
          case (sll_p_hermite )
            temp_knots(r:s) = 2.0_wp*breaks(ncells+1) - breaks(ncells:ncells+r-s:-1)
          end select
        end associate

      end if

      ! Compute interpolation points using Greville-style averaging
      associate( inv_deg => 1.0_wp / real( degree, wp ) )
        do i = 1, ntau
          tau(i) = sum( temp_knots(i+1-degree:i) ) * inv_deg
        end do
      end associate

      ! Periodic case: apply periodic BCs to interpolation points
      if ( bc_xmin == sll_p_periodic ) then
        tau(:) = modulo( tau(:)-xmin, xmax-xmin ) + xmin

      ! Non-periodic case, odd degree: fix round-off issues
      else if ( self%odd == 1 ) then
        tau(1)    = xmin
        tau(ntau) = xmax
      end if

    end associate ! breaks
    end associate ! all other variables

  end subroutine s_compute_interpolation_points_non_uniform

  !-----------------------------------------------------------------------------
  subroutine s_compute_num_diags_non_uniform( self, kl, ku )
    class(sll_t_spline_interpolator_1d), intent(in   ) :: self
    integer                            , intent(  out) :: kl
    integer                            , intent(  out) :: ku

    integer  :: i,j,s,d,icell
    real(wp) :: x

    ku = 0
    kl = 0

    associate( nbasis   => self % bspl % nbasis, &
               degree   => self % bspl % degree, &
               offset   => self % offset       , &
               nbc_xmin => self % nbc_xmin     , &
               nbc_xmax => self % nbc_xmax )

      if (self%bc_xmin == sll_p_periodic) then

        do i = 1, nbasis
          x = self % tau(i)
          icell = self % bspl % find_cell( x )
          do s = 1, degree+1
            j = modulo(icell-offset-2+s,nbasis)+1
            d = j-i
            if (d >  nbasis/2) then; d = d-nbasis; else &
            if (d < -nbasis/2) then; d = d+nbasis; end if
            ku = max( ku,  d )
            kl = max( kl, -d )
          end do
        end do

      else

        do i = nbc_xmin+1, nbasis-nbc_xmax
          x = self % tau(i-nbc_xmin)
          icell = self % bspl % find_cell( x )
          do s = 1, degree+1
            j = icell-1+s
            d = j-i
            ku = max( ku,  d )
            kl = max( kl, -d )
          end do
        end do
        ! FIXME: In Hermite case ku and kl computed in general case where
        !        derivatives of B-splines do not vanish at boundary
        ku = ku + nbc_xmin
        kl = kl + nbc_xmax

      end if

    end associate

  end subroutine s_compute_num_diags_non_uniform

end module sll_m_spline_interpolator_1d
