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

!> @ingroup fft
!> @brief
!> Functions for Fast Fourier Transform
!> @details
!> These functions do not depend on external library
module sll_m_fft
#ifndef DOXYGEN_SHOULD_SKIP_THIS
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#include "sll_working_precision.h"
#include "sll_assert.h"
#include "sll_memory.h"
#include "sll_errors.h"
#include "sll_fftw.h"

  use sll_m_constants, only: &
    sll_p_pi

  use sll_m_utilities, only : &
    sll_f_is_power_of_two

  implicit none 

  public :: &
    sll_t_fft, &
    sll_p_fft_forward, &
    sll_p_fft_backward, &
    sll_p_fft_measure, &
    sll_p_fft_patient, &
    sll_p_fft_estimate, &
    sll_p_fft_exhaustive, &
    sll_p_fft_wisdom_only, &
    sll_s_fft_print_defaultlib, &
    sll_f_fft_allocate_aligned_complex, &
    sll_f_fft_allocate_aligned_real, &
    sll_s_fft_init_r2r_1d, &
    sll_s_fft_init_c2r_1d, &
    sll_s_fft_init_r2c_1d, &
    sll_s_fft_init_c2c_1d, &
    sll_s_fft_init_r2c_2d, &
    sll_s_fft_init_c2r_2d, &
    sll_s_fft_init_c2c_2d, &
    sll_s_fft_exec_r2r_1d, &
    sll_s_fft_exec_c2r_1d, &
    sll_s_fft_exec_r2c_1d, &
    sll_s_fft_exec_c2c_1d, &
    sll_s_fft_exec_r2c_2d, &
    sll_s_fft_exec_c2r_2d, &
    sll_s_fft_exec_c2c_2d, &
    sll_s_fft_set_mode_c2r_1d, &
    sll_f_fft_get_mode_r2c_1d, &
    sll_s_fft_free


  private
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++    


  !> Derived type for ftt plan
  type :: sll_t_fft
     ! twiddle factors complex case
     sll_comp64, allocatable, private :: t(:) 
     ! twiddles factors real case
     sll_real64, allocatable, private :: twiddles(:)
     ! twiddles factors real case
     sll_real64, allocatable, private :: twiddles_n(:)

     logical                          :: normalized !< Boolean telling whether or not values of the FFT should be normalized by \a problem_shape
     sll_int32                        :: direction
     sll_int32                        :: problem_rank
     sll_int32, allocatable           :: problem_shape(:)
     sll_int32, allocatable, private  :: scramble_index(:)
     sll_int32, private                :: transform_type !< Type of the transform. Use for assertion to make sure execution is called of the same type as fft object was initialized for.
  end type sll_t_fft



  !> Set a forward fft
  integer, parameter :: sll_p_fft_forward = -1
  !> Set a backward fft
  integer, parameter :: sll_p_fft_backward = 1

  ! Flags for initialization of the plan: These options are only used in FFTW interface and all set to -1 here
  integer, parameter :: sll_p_fft_measure = -1 !< FFTW planning-rigor flag FFTW_MEASURE (optimized plan) NOTE: planner overwrites the input array during planning  [value 0]
  integer, parameter :: sll_p_fft_patient = -1 !< FFTW planning-rigor flag FFTW_PATIENT (more optimizaton than MEASURE) NOTE: planner overwrites the input array during planning  [value 32]
  integer, parameter :: sll_p_fft_estimate = -1 !< FFTW planning-rigor flag FFTW_ESTIMATE (simple heuristic for planer)   [value 64]. This is our default value
  integer, parameter :: sll_p_fft_exhaustive = -1 !< FFTW planning-rigor flag FFTW_EXHAUSTIVE (more optimization than PATIENT) NOTE: planner overwrites the input array during planning  [value 8]
  integer, parameter :: sll_p_fft_wisdom_only = -1 ! 2097152 !< FFTW planning-rigor flag FFTW_WISDOM_ONLY (planner only initialized if wisdom is available)  [value 2097152]


  ! Flags for the various types of transform (to make sure same type of init and execute functions are used)
  integer, parameter :: p_fftw_c2c_1d = 0
  integer, parameter :: p_fftw_r2r_1d = 1
  integer, parameter :: p_fftw_r2c_1d = 2
  integer, parameter :: p_fftw_c2r_1d = 3
  integer, parameter :: p_fftw_c2c_2d = 4
  integer, parameter :: p_fftw_r2r_2d = 5
  integer, parameter :: p_fftw_r2c_2d = 6
  integer, parameter :: p_fftw_c2r_2d = 7

contains


  !> Debug function
  subroutine sll_s_fft_print_defaultlib()
    print *, 'The library used is SLLFFT'
  end subroutine


  !> Function to allocate an aligned complex array
  function sll_f_fft_allocate_aligned_complex(n) result(data)!, ptr_data, data)
    sll_int32,                 intent(in)  :: n                !< Size of the pointer
    complex(f64), pointer     :: data(:) !< Array to be allocated
   
    allocate(data(n))
    SLL_WARNING('sll_f_fft_allocate_aligned_complex', 'Aligned allocation not implemented for SLLFFT. Usual allocation.')

  end function sll_f_fft_allocate_aligned_complex


  !> Function to allocate an aligned real array
  function sll_f_fft_allocate_aligned_real(n) result(data)!, ptr_data, data)
    sll_int32,                 intent(in)  :: n                !< Size of the pointer
    real(f64), pointer                :: data(:) !< Array to be allocated
    
   
    allocate(data(n))
    SLL_WARNING('sll_f_fft_allocate_aligned_real', 'Aligned allocation not implemented for SLLFFT. Usual allocation.')

  end function sll_f_fft_allocate_aligned_real


  !> Function to reconstruct the complex FFT mode from the data of a r2r transform
  function sll_f_fft_get_mode_r2c_1d(plan,data,k) result(mode)
    type(sll_t_fft),      intent(in)   :: plan !< FFT plan
    sll_real64,           intent(in)   :: data(0:) !< real data produced by r2r transform
    sll_int32,            intent(in)   :: k    !< mode to be extracted
    sll_comp64                         :: mode !< Complex value of kth mode

    sll_int32                   :: n_2, n


    n = plan%problem_shape(1)
    n_2 = n/2 !ishft(n,-1)

      if( k .eq. 0 ) then
        mode = cmplx(data(0),0.0_f64,kind=f64)
      else if( k .eq. n_2 ) then
        mode = cmplx(data(n_2),0.0_f64,kind=f64)
      else if( k .gt. n_2 ) then
        !mode = complex( data(k-n_2) , -data(n-k+n_2) )
        mode = cmplx( data(n-k) , -data(k) ,kind=f64)
      else
        mode = cmplx( data(k) , data(n-k) ,kind=f64)
      endif
  end function


  !> Function to set a complex mode to the real representation of r2r.
  subroutine sll_s_fft_set_mode_c2r_1d(plan,data,new_value,k)
    type(sll_t_fft), intent(in)    :: plan !< FFT planner object
    sll_real64,      intent(out)   :: data(0:) !< Real array to be set
    sll_comp64,      intent(in)    :: new_value !< Complex value of the kth mode
    sll_int32,       intent(in)    :: k !< mode to be set

    sll_int32 :: n_2, n!, index_mode

    n = plan%problem_shape(1)
    n_2 = n/2 !ishft(n,-1)

      if( k .eq. 0 ) then
        data(0) = real(new_value,kind=f64)
      else if( k .eq. n_2 ) then
        data(n_2) = real(new_value,kind=f64)
      else if( k .gt. n_2 ) then
        data(n-k) = real(new_value,kind=f64)
        data(k) = -aimag(new_value)
      else
        data(k) = real(new_value,kind=f64)
        data(n-k) = aimag(new_value)
      endif
  end subroutine 

! THIS IS THE ORDERING THAT SLLFFT R2R PRODUCES WITHOUT THE REORDERING TO FFTW
!!$
!!$  function sll_f_fft_get_mode_r2c_1d(plan,data,k) result(mode)
!!$    type(sll_t_fft), pointer :: plan
!!$    sll_real64, dimension(0:)   :: data
!!$    sll_int32                   :: k, n_2, n
!!$    sll_comp64                  :: mode
!!$
!!$    n = plan%problem_shape(1)
!!$    n_2 = n/2 !ishft(n,-1)
!!$
!!$      if( k .eq. 0 ) then
!!$        mode = cmplx(data(0),0.0_f64,kind=f64)
!!$      else if( k .eq. n_2 ) then
!!$        mode = cmplx(data(1),0.0_f64,kind=f64)
!!$      else if( k .gt. n_2 ) then
!!$        mode = cmplx( data(2*(n-k)) , -data(2*(n-k)+1),kind=f64 )
!!$      else
!!$        mode = cmplx( data(2*k) , data(2*k+1) ,kind=f64)
!!$      endif
!!$  end function
!!$
!!$
!!$  subroutine sll_s_fft_set_mode_c2r_1d(plan,data,new_value,k)
!!$    type(sll_t_fft), pointer :: plan
!!$    sll_real64, dimension(0:)   :: data
!!$    sll_int32                   :: k, n_2, n
!!$    sll_comp64                  :: new_value
!!$
!!$    n = plan%problem_shape(1)
!!$    n_2 = n/2 !ishft(n,-1)
!!$
!!$      if( k .eq. 0 ) then
!!$        data(0) = real(new_value,kind=f64)
!!$      else if( k .eq. n_2 ) then
!!$        data(1) = real(new_value,kind=f64)
!!$      else if( k .gt. n_2 ) then
!!$        data(2*(n-k)) = real(new_value,kind=f64)
!!$        data(2*(n-k)+1) = -aimag(new_value)
!!$      else
!!$        data(2*k) = real(new_value,kind=f64)
!!$        data(2*k+1) = aimag(new_value)
!!$      endif
!!$  end subroutine



! COMPLEX
! ------
! - 1D -
! ------
  !> Create new 1d complex to complex plan
  subroutine sll_s_fft_init_c2c_1d(plan, nx,array_in,array_out,direction,normalized, aligned, optimization) 
    type(sll_t_fft)                         :: plan !< FFT planner object 
    sll_int32,           intent(in)         :: nx !< Number of points
    sll_comp64,          intent(inout)      :: array_in(:) !< (Typical) input array (gets overwritten for certain options)
    sll_comp64,          intent(inout)      :: array_out(:) !< (Typical) output array (gets overwritten for certain options)
    sll_int32,           intent(in)         :: direction  !< Direction of the FFT (\a sll_p_fft_forward or \a sll_p_fft_backward)
    logical, optional,   intent(in)         :: normalized !< Flag to decide if FFT should be normalized by 1/N (default: \a FALSE)
    logical, optional,   intent(in)         :: aligned    !< Flag only used by FFTW.
    sll_int32, optional, intent(in)         :: optimization !< Flag only used by FFTW.

    sll_int32 :: ierr, i

    SLL_ASSERT(size(array_in).ge.nx)
    SLL_ASSERT(size(array_out).ge.nx)
    plan%transform_type = p_fftw_c2c_1d
    plan%direction = direction
    if( present(normalized) ) then
       plan%normalized = normalized
    else
       plan%normalized = .false.
    end if

    plan%problem_rank = 1
    SLL_ALLOCATE(plan%problem_shape(1),ierr)
    plan%problem_shape = (/ nx /)

    SLL_ALLOCATE(plan%scramble_index(0:nx-1),ierr)
    ! For the moment the mode are stored in the natural order
    ! 0,1,2,...,n-1
    do i=0,nx-1
      plan%scramble_index(i) = i
    enddo
    SLL_ALLOCATE(plan%t(1:nx/2),ierr)
    call compute_twiddles(nx,plan%t)
    if ( direction == sll_p_fft_forward ) then
       plan%t = conjg(plan%t)
    end if
    call bit_reverse_complex(nx/2,plan%t)
  end subroutine sll_s_fft_init_c2c_1d

  !> Compute fast Fourier transform in complex to complex mode.
  subroutine sll_s_fft_exec_c2c_1d(plan,array_in,array_out)
    type(sll_t_fft), intent(in)            :: plan !< FFT planner object
    sll_comp64,      intent(inout)         :: array_in(:) !< Complex data to be Fourier transformed
    sll_comp64,      intent(inout)         :: array_out(:) !< Fourier coefficients on output

    sll_real64 :: factor

    SLL_ASSERT( plan%transform_type == p_fftw_c2c_1d )

    if( loc(array_in(1)) .ne. loc(array_out(1))) then ! out-place transform
       array_out = array_in
    endif

    call fft_dit_nr(array_out,plan%problem_shape(1),plan%t,plan%direction)
    call bit_reverse_complex(plan%problem_shape(1),array_out)


    if ( plan%normalized .EQV. .true.) then
       factor = 1.0_f64/real(plan%problem_shape(1),kind=f64)
       array_out = factor*array_out
    endif

  end subroutine


  ! --------------------
  ! - 2D               -
  ! --------------------
  !> Create new 2d complex to complex plan
  subroutine sll_s_fft_init_c2c_2d(plan, nx,ny,array_in,array_out,direction,normalized, aligned, optimization)
    type(sll_t_fft),       intent(out)    :: plan !< initialized planner object
    sll_int32,             intent(in)     :: nx !< Number of points along first dimension
    sll_int32,             intent(in)     :: ny !< Number of points along second dimension
    sll_comp64,            intent(inout)  :: array_in(0:,0:) !< (Typical) input array (gets overwritten for certain options)
    sll_comp64,            intent(inout)  :: array_out(0:,0:) !<(Typical) output array (gets overwritten for certain options)
    sll_int32,             intent(in)     :: direction  !< Direction of the FFT (\a sll_p_fft_forward or \a sll_p_fft_backward)
    logical,   optional,   intent(in)     :: normalized !< Flag to decide if FFT should be normalized by 1/N (default: \a FALSE)
    logical,   optional,   intent(in)     :: aligned    !< Flag for FFTW (unused)
    sll_int32, optional,  intent(in)      :: optimization !< Planning-rigor flag for FFTW (unused)
    sll_int32                                        :: ierr

    ! This does not look good.
    ! 1. Error checking like this should be permanent, not with assertions.
    SLL_ASSERT(size(array_in,dim=1).ge.nx)
    SLL_ASSERT(size(array_in,dim=2).ge.ny)
    SLL_ASSERT(size(array_out,dim=1).ge.nx)
    SLL_ASSERT(size(array_out,dim=2).ge.ny)

    plan%transform_type = p_fftw_c2c_2d
    plan%direction = direction
    if( present(normalized) ) then
       plan%normalized = normalized
    else
       plan%normalized = .false.
    end if

    plan%problem_rank = 2
    SLL_ALLOCATE(plan%problem_shape(2),ierr)
    plan%problem_shape = (/ nx , ny /)

    SLL_ALLOCATE(plan%t(1:nx/2 + ny/2),ierr)

    call compute_twiddles(nx,plan%t(1:nx/2))
    if ( direction == sll_p_fft_forward ) then
       plan%t(1:nx/2) = conjg(plan%t(1:nx/2))
    end if
    call bit_reverse_complex(nx/2,plan%t(1:nx/2))


    call compute_twiddles(ny,plan%t(nx/2+1:nx/2 + ny/2))
    if ( direction == sll_p_fft_forward ) then
       plan%t(ny/2+1:ny/2 + ny/2) = conjg(plan%t(nx/2+1:nx/2 + ny/2))
    end if
    call bit_reverse_complex(ny/2,plan%t(nx/2+1:nx/2 + ny/2))


  end subroutine sll_s_fft_init_c2c_2d


!> Compute fast Fourier transform in complex to complex mode.
  subroutine sll_s_fft_exec_c2c_2d(plan,array_in,array_out)
    type(sll_t_fft),  intent(in)     :: plan !< FFT planner object
    sll_comp64,       intent(inout)  :: array_in(0:,0:) !< Complex data to be Fourier transformed 
    sll_comp64,       intent(inout)  :: array_out(0:,0:) !< Fourier coefficients on output

    sll_int32                                       :: i, nx, ny
    sll_int32, dimension(2)                         :: fft_shape
    sll_real64 :: factor

    SLL_ASSERT( plan%transform_type == p_fftw_c2c_2d )

    if( loc(array_in(1,1)) .ne. loc(array_out(1,1))) then ! out-place transform
       array_out = array_in  ! copy source
    endif
    fft_shape(1:2) = plan%problem_shape(1:2)
    nx = fft_shape(1)
    ny = fft_shape(2)

    do i=0,ny-1
       call fft_dit_nr(array_out(0:nx-1,i),nx,plan%t(1:nx/2),plan%direction)
       call bit_reverse_complex(nx,array_out(0:nx-1,i))
    enddo
    
    do i=0,nx-1
       call fft_dit_nr( &
            array_out(i,0:ny-1), &
            ny, &
            plan%t(nx/2+1:nx/2+ny/2), &
            plan%direction)
       call bit_reverse_complex(ny,array_out(i,0:ny-1))
    enddo

    if ( plan%normalized .EQV. .true.) then
      factor = 1.0_f64/real(nx*ny,kind=f64)
      array_out = factor*array_out
    endif

  end subroutine



! REAL
!> Create new 1d real to real plan
  subroutine sll_s_fft_init_r2r_1d(plan, nx,array_in,array_out,direction,normalized, aligned, optimization)
    type(sll_t_fft),     intent(out)   :: plan !< FFT planner object
    sll_int32,           intent(in)    :: nx !< Number of points
    sll_real64,          intent(inout) :: array_in(:) !< (Typical) input array (gets overwritten for certain options)
    sll_real64,          intent(inout) :: array_out(:) !< (Typical) output array (gets overwritten for certain options)
    sll_int32,           intent(in)    :: direction  !< Direction of the FFT (\a sll_p_fft_forward or \a sll_p_fft_backward)
    logical,   optional, intent(in)    :: normalized !< Flag to decide if FFT should be normalized by 1/N (default: \a FALSE)
    logical,   optional, intent(in)    :: aligned    !< Flag for FFTW (unused)
    sll_int32, optional, intent(in)    :: optimization !< Planning-rigor flag for FFTW (unsed)

    sll_int32                                    :: ierr, i

    SLL_ASSERT(size(array_in).eq.nx)
    SLL_ASSERT(size(array_out).eq.nx)
    plan%transform_type = p_fftw_r2r_1d
    plan%direction = direction
    if( present(normalized) ) then
       plan%normalized = normalized
    else
       plan%normalized = .false.
    end if

    plan%problem_rank = 1
    SLL_ALLOCATE(plan%problem_shape(1),ierr)
    plan%problem_shape = (/ nx /)

    SLL_ALLOCATE(plan%scramble_index(0:nx/2),ierr)
    ! The mode are stored in the order 0,n/2,1,2,3,...,n/2-1
    ! with the representation r_0,r_n/2,r_1,i_1,...,
    ! The modes 0 and n/2 are purely real.
    plan%scramble_index(0) = 0
    plan%scramble_index(1) = nx/2
    do i=1,nx/2-1
      plan%scramble_index(i+1) = i
    enddo

    SLL_ALLOCATE(plan%twiddles(0:nx/2-1),ierr)
    SLL_ALLOCATE(plan%twiddles_n(0:nx-1),ierr)
    call compute_twiddles_real_array( nx, plan%twiddles_n )
    call compute_twiddles_real_array( nx/2, plan%twiddles(0:nx/2-1) )
    if( direction .eq. sll_p_fft_forward ) then
      call conjg_in_pairs(nx/2,plan%twiddles)
    endif
    call bit_reverse_in_pairs( nx/4, plan%twiddles(0:nx/2-1))
  end subroutine sll_s_fft_init_r2r_1d

  !> Compute fast Fourier transform in real to real mode.
  subroutine sll_s_fft_exec_r2r_1d(plan,array_in,array_out)
    type(sll_t_fft),          intent(in)    :: plan !< FFT planner object
    sll_real64,               intent(inout) :: array_in(:) !< Real data to be Fourier transformed
    sll_real64,               intent(inout) :: array_out(:) !< Fourier coefficients in real form (sin/cos coefficients)

    sll_int32 :: nx, k
    sll_real64 :: factor
    sll_real64 :: tmp(plan%problem_shape(1))

    SLL_ASSERT( plan%transform_type == p_fftw_r2r_1d )

    nx = plan%problem_shape(1)

    if( loc(array_in(1)) .ne. loc(array_out(1))) then ! out-place transform
       array_out = array_in
    endif

    ! Change from FFTW ordering back to SLLFFT ordering
    if (plan%direction .EQ. sll_p_fft_backward) then
       tmp = array_out
       do k=1,nx/2-1
          array_out(2*k+1) = tmp(k+1)
          array_out(2*k+2) = tmp(nx-k+1)
       end do
       array_out(1) = tmp(1)
       array_out(2) = tmp(nx/2+1)
    end if

    
    call real_data_fft_dit( array_out, nx, plan%twiddles, plan%twiddles_n, plan%direction )
    
    ! Change to FFTW ordering
    if (plan%direction .EQ. sll_p_fft_forward) then
       tmp = array_out      
       do k=1,nx/2-1
          array_out(k+1) = tmp(2*k+1)
          array_out(nx-k+1) = tmp(2*k+1+1)
       end do
       array_out(nx/2+1) = tmp(2)
    end if

    if( plan%normalized .EQV. .TRUE. ) then
      factor = 1.0_f64/real(plan%problem_shape(1),kind=f64)
      array_out = factor*array_out
    endif

  end subroutine
! END REAL


! REAL TO COMPLEX
  !> Create new 1d real to complex plan for forward FFT
  subroutine sll_s_fft_init_r2c_1d(plan, nx,array_in,array_out, normalized, aligned, optimization)
    type(sll_t_fft),     intent(out)   :: plan !< FFT planner object
    sll_int32,           intent(in)    :: nx !< Number of points
    sll_real64,          intent(inout) :: array_in(:) !< (Typical) input array (gets overwritten for certain options)
    sll_comp64,          intent(out)   :: array_out(:) !< (Typical) output array (gets overwritten for certain options)
    logical,   optional, intent(in)    :: normalized !< Flag to decide if FFT should be normalized by 1/N (default: \a FALSE)
    logical,   optional, intent(in)    :: aligned    !< Flag for FFTW (unused)
    sll_int32, optional, intent(in)    :: optimization !< Planning-rigor flag for FFTW. (unused)

    sll_int32                                    :: ierr

    SLL_ASSERT(size(array_in).eq.nx)
    SLL_ASSERT(size(array_out).eq.nx/2+1)
    plan%transform_type = p_fftw_r2c_1d
    plan%direction = sll_p_fft_forward
    if( present(normalized) ) then
       plan%normalized = normalized
    else
       plan%normalized = .false.
    end if
    plan%problem_rank = 1
    SLL_ALLOCATE(plan%problem_shape(1),ierr)
    plan%problem_shape = (/ nx /)

    SLL_ALLOCATE(plan%twiddles(0:nx/2-1),ierr)
    SLL_ALLOCATE(plan%twiddles_n(0:nx-1),ierr)
    call compute_twiddles_real_array( nx, plan%twiddles_n )
    call compute_twiddles_real_array( nx/2, plan%twiddles(0:nx/2-1) )
    call conjg_in_pairs(nx/2,plan%twiddles)
    call bit_reverse_in_pairs( nx/4, plan%twiddles(0:nx/2-1))
  end subroutine sll_s_fft_init_r2c_1d



  !> Compute fast Fourier transform in real to complex mode.
  subroutine sll_s_fft_exec_r2c_1d(plan,array_in,array_out)
    type(sll_t_fft), intent(in)            :: plan !< FFT planner object
    sll_real64,      intent(inout)         :: array_in(0:) !< Real input data to be Fourier transformed
    sll_comp64,      intent(out)           :: array_out(0:) !< Complex Fourier mode (only first half due to symmetry)

    sll_real64 :: factor
    sll_int32 :: i, nx

    SLL_ASSERT( plan%transform_type == p_fftw_r2c_1d )

    nx = plan%problem_shape(1)

    call real_data_fft_dit( array_in, nx , plan%twiddles, plan%twiddles_n, plan%direction )
    !mode k=0
    array_out(0) = cmplx(array_in(0),0.0_f64,kind=f64)
    !mode k=n/2
    array_out(nx/2) = cmplx(array_in(1),0.0_f64,kind=f64)
    !mode k=1 to k= n-2
    do i=1,nx/2-1
        array_out(i) = cmplx(array_in(2*i),array_in(2*i+1),kind=f64)
        !array_out(plan%N-i) = cmplx(array_in(2*i),-array_in(2*i+1),kind=f64)
    enddo

    if( plan%normalized .EQV. .TRUE. ) then
      factor = 1.0_f64/real(plan%problem_shape(1),kind=f64)
      array_out = factor*array_out
    endif
  end subroutine
! END REAL TO COMPLEX

! COMPLEX TO REAL
  !> Create new 1d complex to real plan for backward FFT
  subroutine sll_s_fft_init_c2r_1d(plan, nx,array_in,array_out, normalized, aligned, optimization)
    type(sll_t_fft),     intent(out)      :: plan !< FFT planner object
    sll_int32,           intent(in)       :: nx !< Number of points
    sll_comp64,          intent(in)       :: array_in(:)  !< (Typical) input array (gets overwritten for certain options)
    sll_real64,          intent(out)      :: array_out(:) !< (Typical) output array (gets overwritten for certain options)
    logical,   optional, intent(in)       :: normalized !< Flag to decide if FFT should be normalized by 1/N (default: \a FALSE)
    logical,   optional, intent(in)       :: aligned    !< Flag for FFTW (unused)
    sll_int32, optional, intent(in)       :: optimization !< Planning-rigor flag for FFTW. (unused)
    
    sll_int32 :: ierr

    SLL_ASSERT(size(array_in).eq.nx/2+1)
    SLL_ASSERT(size(array_out).eq.nx)
    plan%transform_type = p_fftw_c2r_1d
    plan%direction = sll_p_fft_backward
    if( present(normalized) ) then
       plan%normalized = normalized
    else
       plan%normalized = .false.
    end if
    plan%problem_rank = 1
    SLL_ALLOCATE(plan%problem_shape(1),ierr)
    plan%problem_shape = (/ nx /)

    SLL_ALLOCATE(plan%twiddles(0:nx/2-1),ierr)
    SLL_ALLOCATE(plan%twiddles_n(0:nx-1),ierr)
    call compute_twiddles_real_array( nx, plan%twiddles_n )
    call compute_twiddles_real_array( nx/2, plan%twiddles(0:nx/2-1) )
    call bit_reverse_in_pairs( nx/4, plan%twiddles(0:nx/2-1))
  end subroutine sll_s_fft_init_c2r_1d


  !> Compute fast Fourier transform in complex to real mode.
  subroutine sll_s_fft_exec_c2r_1d(plan,array_in,array_out)
    type(sll_t_fft), intent(in)    :: plan !< FFT planner objece
    sll_comp64,      intent(inout) :: array_in(0:) !< Complex Fourier coefficient to be transformed back
    sll_real64,      intent(inout) :: array_out(0:) !< Real result of Fourier transform

    sll_int32                                     :: nx, i
    sll_real64 :: factor

    SLL_ASSERT( plan%transform_type == p_fftw_c2r_1d )

    nx = plan%problem_shape(1)

    !mode k=0
    array_out(0) = real(array_in(0),kind=f64)
    !mode k=n/2
    array_out(1) = real(array_in(nx/2),kind=f64)
    !mode k=1 to k= n-2
    do i=1,nx/2-1
      array_out(2*i) = real(array_in(i),kind=f64)
      array_out(2*i+1) = aimag(array_in(i))
    enddo
    call real_data_fft_dit( array_out, nx , plan%twiddles, plan%twiddles_n, plan%direction )

    if( plan%normalized .EQV. .TRUE. ) then
      factor = 1.0_f64/real(plan%problem_shape(1),kind=f64)
      array_out = factor*array_out
    endif

  end subroutine
! END COMPLEX TO REAL


! REAL TO COMPLEX 2D
  !> Create new 2d complex to real plan for forward FFT
  subroutine sll_s_fft_init_r2c_2d(plan, nx,ny,array_in,array_out,normalized, aligned, optimization)
    type(sll_t_fft),     intent(out)   :: plan !< FFT planner object
    sll_int32,           intent(in)    :: nx !< Number of points along first dimension
    sll_int32,           intent(in)    :: ny !< Number of points along second dimension
    sll_real64,          intent(inout) :: array_in(:,:) !< (Typical) input array (gets overwritten for certain options)
    sll_comp64,          intent(out)   :: array_out(:,:) !< (Typical) output array (gets overwritten for certain options)
    logical,   optional, intent(in)    :: normalized !< Flag to decide if FFT should be normalized by 1/N (default: \a FALSE)
    logical,   optional, intent(in)    :: aligned    !< Flag for FFTW (unused)
    sll_int32, optional, intent(in)    :: optimization !< Planning-rigor flag for FFTW (unused)


    sll_int32                                    :: ierr

    SLL_ASSERT(size(array_in,dim=1).eq.nx)
    SLL_ASSERT(size(array_in,dim=2).eq.ny)
    SLL_ASSERT(size(array_out,dim=1).eq.nx/2+1)
    SLL_ASSERT(size(array_out,dim=2).eq.ny)
    plan%transform_type = p_fftw_r2c_2d
    plan%direction = sll_p_fft_forward
    if( present(normalized) ) then
       plan%normalized = normalized
    else
       plan%normalized = .false.
    end if
    plan%problem_rank = 2
    SLL_ALLOCATE(plan%problem_shape(2),ierr)
    plan%problem_shape = (/ nx , ny /)

    SLL_ALLOCATE(plan%t(1:ny/2),ierr)
    call compute_twiddles(ny,plan%t)
    plan%t = conjg(plan%t)
    call bit_reverse_complex(ny/2,plan%t)

    SLL_ALLOCATE(plan%twiddles(0:nx/2-1),ierr)
    SLL_ALLOCATE(plan%twiddles_n(0:nx-1),ierr)
    call compute_twiddles_real_array( nx, plan%twiddles_n )
    call compute_twiddles_real_array( nx/2, plan%twiddles(0:nx/2-1) )
    call conjg_in_pairs(nx/2,plan%twiddles)
    call bit_reverse_in_pairs( nx/4, plan%twiddles(0:nx/2-1))
  end subroutine sll_s_fft_init_r2c_2d

 !> Compute fast Fourier transform in real to complex mode.
  subroutine sll_s_fft_exec_r2c_2d(plan,array_in,array_out)
    type(sll_t_fft), intent(in)           :: plan      !< FFT planner object
    sll_real64,      intent(inout)        :: array_in(:,:)  !< Real input data to be Fourier transformed
    sll_comp64,      intent(out)          :: array_out(:,:) !< Complex Fourier coefficients (only half part along first dimension due to symmetry)

    sll_int32                                       :: nx, i, ny, k
    sll_real64 :: factor

    SLL_ASSERT( plan%transform_type == p_fftw_r2c_2d )

    nx = plan%problem_shape(1)
    ny = plan%problem_shape(2)

    do k=0,ny-1
      call real_data_fft_dit( array_in(0:nx-1,k), nx , plan%twiddles, plan%twiddles_n, plan%direction )
      array_out(0,k) = cmplx(array_in(0,k),0.0_f64,kind=f64)
      array_out(nx/2,k) = cmplx(array_in(1,k),0.0_f64,kind=f64)
      do i=1,nx/2-1
          array_out(i,k) = cmplx(array_in(2*i,k),array_in(2*i+1,k),kind=f64)
      enddo
    enddo
    do k=0,nx/2
      call fft_dit_nr(array_out(k,0:ny-1),ny,plan%t,plan%direction)
      call bit_reverse_complex(ny,array_out(k,0:ny-1))
    enddo

    if( plan%normalized .EQV. .TRUE. ) then
      factor = 1.0_f64/real(nx*ny,kind=f64)
      array_out = factor*array_out
    endif
  end subroutine
! END REAL TO COMPLEX 2D

! COMPLEX TO REAL 2D
 !> Create new 2d real to complex plan for backward FFT
  subroutine sll_s_fft_init_c2r_2d(plan, nx,ny,array_in,array_out,normalized, aligned, optimization)
    type(sll_t_fft),     intent(out)      :: plan !< FFT planner object
    sll_int32,           intent(in)       :: nx !< Number of points
    sll_int32,           intent(in)      :: ny !< Number of points along second dimension
    sll_comp64,          intent(in)       :: array_in(:,:)  !< (Typical) input array (gets overwritten for certain options)
    sll_real64,          intent(out)      :: array_out(:,:) !< (Typical) output array (gets overwritten for certain options)
    logical,   optional, intent(in)       :: normalized !< Flag to decide if FFT should be normalized by 1/N (default: \a FALSE)
    logical,   optional, intent(in)       :: aligned    !< Flag for FFTW (unused)
    sll_int32, optional, intent(in)       :: optimization !< Planning-rigor flag for FFTW (unused). 

    sll_int32                                    :: ierr

    SLL_ASSERT(size(array_in,dim=1).eq.nx/2+1)
    SLL_ASSERT(size(array_in,dim=2).eq.ny)
    SLL_ASSERT(size(array_out,dim=1).eq.nx)
    SLL_ASSERT(size(array_out,dim=2).eq.ny)

    plan%transform_type = p_fftw_c2r_2d
    plan%direction = sll_p_fft_backward
    if( present(normalized) ) then
       plan%normalized = normalized
    else
       plan%normalized = .false.
    end if
    plan%problem_rank = 2
    SLL_ALLOCATE(plan%problem_shape(2),ierr)
    plan%problem_shape = (/ nx , ny /)

    SLL_ALLOCATE(plan%t(1:ny/2),ierr)
    call compute_twiddles(ny,plan%t)
    call bit_reverse_complex(ny/2,plan%t)

    SLL_ALLOCATE(plan%twiddles(0:nx/2-1),ierr)
    SLL_ALLOCATE(plan%twiddles_n(0:nx-1),ierr)
    call compute_twiddles_real_array( nx, plan%twiddles_n )
    call compute_twiddles_real_array( nx/2, plan%twiddles(0:nx/2-1) )
    call bit_reverse_in_pairs( nx/4, plan%twiddles(0:nx/2-1))
  end subroutine sll_s_fft_init_c2r_2d


  !> Compute fast Fourier transform in complex to real mode.
  subroutine sll_s_fft_exec_c2r_2d(plan,array_in,array_out)
    type(sll_t_fft), intent(in)    :: plan !< FFT planner objece
    sll_comp64,      intent(inout) :: array_in(0:, 0:) !< Complex Fourier coefficient to be transformed back
    sll_real64,      intent(inout) :: array_out(0:, 0:) !< Real result of Fourier transform
    
    sll_int32                                       :: nx, i, ny, k, j
    sll_real64 :: factor

    SLL_ASSERT( plan%transform_type == p_fftw_c2r_2d )

    nx = plan%problem_shape(1)
    ny = plan%problem_shape(2)

    do j=0,nx/2
      call fft_dit_nr(array_in(j,0:ny-1),ny,plan%t,plan%direction)
      call bit_reverse_complex(ny,array_in(j,0:ny-1))
    enddo
    do i=0,ny-1
      array_out(0,i) = real(array_in(0,i),kind=f64)
      array_out(1,i) = real(array_in(nx/2,i),kind=f64)
      do k=1,nx/2-1
        array_out(2*k,i) = real(array_in(k,i),kind=f64)
        array_out(2*k+1,i) = aimag(array_in(k,i))
      enddo
      call real_data_fft_dit( array_out(0:nx-1,i), nx , plan%twiddles, plan%twiddles_n, plan%direction )
    enddo


    if( plan%normalized .EQV. .TRUE. ) then
      factor = 1.0_f64/real(nx*ny,kind=f64)
      array_out = factor*array_out
    endif
  end subroutine
! END COMPLEX TO REAL 2D

! SELALIB SELALIB SELALIB SELALIB SELALIB SELALIB SELALIB SELALIB SELALIB SELALIB
! SELALIB SELALIB SELALIB SELALIB SELALIB SELALIB SELALIB SELALIB SELALIB SELALIB
! SELALIB SELALIB SELALIB SELALIB SELALIB SELALIB SELALIB SELALIB SELALIB SELALIB





  !> Deallocate the fft plan
  subroutine sll_s_fft_free(plan)
   type(sll_t_fft), intent( inout ) :: plan
   sll_int32 :: ierr

      if(allocated(plan%t)) then
        deallocate(plan%t)
      endif
      if(allocated(plan%twiddles)) then
        deallocate(plan%twiddles)
      endif
      if(allocated(plan%twiddles_n)) then
        deallocate(plan%twiddles_n)
      endif
      if(allocated(plan%problem_shape)) then
        deallocate(plan%problem_shape)
      endif
  end subroutine




!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! The following is only concerning the sll implementation of the fft !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  ! We choose the convention in which the direction of the FFT is determined
  ! by the conjugation of the twiddle factors. If
  !
  !                  omega_N = -j2*pi/N
  !
  ! then we call this the FORWARD transform.

  ! Compute the twiddle factors by Singleton's method. N is a factor of 2.
  ! N represents the full size of the transform for which the twiddles are
  ! needed, however, because of redundancies in the values of the twiddle
  ! factors, only half of them need to be stored.
  subroutine compute_twiddles( n, t )
    intrinsic                                 :: exp, real
    sll_int32                                 :: n ! size of data for FFT
    sll_comp64, dimension(1:n/2), intent(out) :: t
    sll_int32                                 :: k
    sll_real64                                :: theta
    if ( n .eq. 0 ) then
       print *, "ERROR: Zero array size passed to compute_twiddle_factors()."
       return
    else if (size(t) .lt. n/2) then
       print *, "ERROR: compute_twiddles(), array is of insufficient size."
       return
    else
       theta   = 2.0_f64*sll_p_pi/real(n,kind=f64)      ! basic angular interval
       ! By whatever means we use to compute the twiddles, some sanity
       ! checks are in order:
       ! t(1)     = (1,0)
       ! t(n/8+1) = (sqrt(2)/2, sqrt(2)/2)
       ! t(n/4+1) = (0,1)
       ! t(n/2+1) = (0,-1) ... but this one is not stored
       do k = 0,n/2-1
          t(k+1) = exp((0.0_f64,1.0_f64)*theta*cmplx(k,0.0,kind=f64))
       end do
       ! might as well fix this by hand since the result isn't exact otherwise:
       t(n/4+1) = (0.0_f64,1.0_f64)
    end if
  end subroutine compute_twiddles

  ! Need the following macros to access a real array as if it were made of
  ! complex numbers. Please notice an important detail:
  ! The validity of these macros depends on the type of indexing of the array.
  ! For example, for a zero-indexed array, the following pair should be used:
#define CREAL0(array, i) array(2*(i))
#define CIMAG0(array, i) array(2*(i)+1)
  ! For a one-indexed array, the following pair is appropriate
#define CREAL1(array, i) array(2*(i)-1)
#define CIMAG1(array, i) array(2*(i))
  ! In this module we use both types of indexings, so one should be aware
  ! of the macro pair to use.

  subroutine compute_twiddles_real_array( n, t )
    intrinsic                                 :: cos, sin, real
    sll_int32                                 :: n ! size of data for FFT
    sll_real64, dimension(0:n-1), intent(out) :: t
    sll_int32                                 :: k
    sll_real64                                :: theta
    SLL_ASSERT(sll_f_is_power_of_two(int(n,i64)))
    theta   = 2.0_f64*sll_p_pi/real(n,kind=f64)      ! basic angular interval
    ! By whatever means we use to compute the twiddles, some sanity
    ! checks are in order:
    ! t(0)   = 1, t(1)      = 0
    ! t(n/8) = (sqrt(2)/2, sqrt(2)/2)
    ! t(n/4) = (0,1)
    ! t(n/2) = (0,-1) ... but this one is not stored
    do k = 0,n/2-1
       CREAL0(t,k) = cos(real(k,kind=f64)*theta)
       CIMAG0(t,k) = sin(real(k,kind=f64)*theta)
    end do
    ! might as well fix this by hand since the result isn't exact otherwise:
    CREAL0(t,0) = 1.0_f64
    CIMAG0(t,0) = 0.0_f64

    if( n >= 4 ) then
    CREAL0(t,n/4) = 0.0_f64
    CIMAG0(t,n/4) = 1.0_f64
    endif

    if( n >= 8 ) then
    CREAL0(t,n/8) = sqrt(2.0_f64)/2.0_f64
    CIMAG0(t,n/8) = sqrt(2.0_f64)/2.0_f64
    endif
  end subroutine compute_twiddles_real_array

  !***************************************************************************
  ! The following function uses Meyer's algorithm (also known as the
  ! Gold-Rader algorithm according to "Rodriguez. J. IEEE Transactions on
  ! Acoustics, Speech and Signal Porcessing. Vol. 37. No. 8. August 1989.
  ! p. 1298.") to yield the in-place, bit-reversed ordering of the array 'a'.
  ! Use with n = power of 2.
  !
  !***************************************************************************
#define MAKE_BIT_REVERSE_FUNCTION( function_name, data_type )   \
  subroutine function_name( n, a );                             \
    intrinsic ishft, iand;                                      \
    data_type, dimension(:), intent(inout) :: a;                \
    integer, intent(in)                    :: n;                \
    integer                                :: i;                \
    integer                                :: j;                \
    integer                                :: k;                \
    data_type                              :: tmp;              \
    SLL_ASSERT(sll_f_is_power_of_two(int(n,i64)));                    \
    j = 0;                                                      \
    k = n;                                                      \
    do i=0,n-2;                                                 \
       if( i<j ) then;                                          \
          tmp    = a(j+1);                                      \
          a(j+1) = a(i+1);                                      \
          a(i+1) = tmp;                                         \
       end if;                                                  \
       k = ishft(n,-1);                                         \
       do while ( k<=j );                                       \
          j = j-k;                                              \
          k = ishft(k,-1);                                      \
       end do;                                                  \
       j = j+k;                                                 \
    end do;                                                     \
  end subroutine function_name

  MAKE_BIT_REVERSE_FUNCTION( bit_reverse_complex, sll_comp64 )
!PN DEFINED BUT NOT USED
! MAKE_BIT_REVERSE_FUNCTION( bit_reverse_integer32, sll_int32 )
! MAKE_BIT_REVERSE_FUNCTION( bit_reverse_integer64, sll_int64 )
! MAKE_BIT_REVERSE_FUNCTION( bit_reverse_real, sll_real64 )

  ! ugly special case to bit-reverse a complex array that is represented
  ! by an array of reals. This is truly awful...
  subroutine bit_reverse_in_pairs( num_pairs, a )
    intrinsic ishft, iand
    integer, intent(in)                    :: num_pairs
    sll_real64, dimension(0:2*num_pairs-1), intent(inout) :: a
    integer                                :: i
    integer                                :: j
    integer                                :: k
    sll_real64, dimension(1:2)             :: tmp
    SLL_ASSERT(sll_f_is_power_of_two(int(num_pairs,i64)))
    j = 0
    k = num_pairs
#define ARRAY(i) a(2*(i):2*(i)+1)
    do i=0,num_pairs-2
       if( i<j ) then
!          write (*,'(a,i4,a,i4)') 'will swap ',i,' and ',j
          tmp      = ARRAY(j)
          ARRAY(j) = ARRAY(i)
          ARRAY(i) = tmp
       end if
       k = ishft(num_pairs,-1)
       do while ( k<=j )
          j = j-k
          k = ishft(k,-1)
       end do
       j = j+k
    end do
  end subroutine bit_reverse_in_pairs
#undef ARRAY

  !Conjugate complex array that is represented by an array of reals.
  subroutine conjg_in_pairs(n,array)
    sll_real64, dimension(:) :: array
    sll_int32 :: n, i

    SLL_ASSERT( size(array) .eq. n )

    do i=2,n,2
      array(i) = -array(i)
    enddo
  end subroutine

  ! *************************************************************************
  ! Decimation in time FFT, natural order input, bit-reversed output (=NR).
  ! Size of the data must be a power of 2. Twiddle array must be in
  ! bit-reversed order. This implementation is 'cache-oblivious'. This is
  ! only a placeholder for an FFT really, it is not very efficient since it:
  ! - is just a simple radix-2
  ! - is not parallelized
  ! - no hardware acceleration
  ! But it is good to test our intarfaces.
  !
  ! Implementation notes:
  !
  ! The following function implements the Cooley-Tukey algorithm:
  !
  ! Y_r^(n-1)  ---------------> Y_r^(n)=Y_r^(n-1) + omega_N^r*Z_r^(n-1) = X_r
  !             *           *
  !                *     *
  !                   *
  !                *     *
  !             *           *
  ! Z_r^(n-1)  ---------------> Y_(r+N/2)^(n-1) - omega_N^r*Z_r^(n-1)=X_(r+N/2)
  !
  ! In terms of the twiddles used at each level, the size two problems use only
  ! the omega_2^1 twiddle, the size-4 problems use omega_4^1 and omega_4^2,
  ! the size-8 problems use omega_8^1, omega_8^2, omega_8^3 and omega_8^4. This
  ! is the expected progression of the twiddle indices as we move deeper into
  ! the recursions.
  !
  ! *************************************************************************

  subroutine fft_dit_nr(data, n, twiddles, sign)
    sll_comp64, dimension(:), intent(inout) :: data
    sll_int32, intent(in)                   :: sign
    sll_comp64, dimension(:), intent(in)    :: twiddles
    sll_int32, intent(in)                   :: n
    SLL_ASSERT(sll_f_is_power_of_two(int(n,i64)))
    SLL_ASSERT(size(data) .ge. n)
    call fft_dit_nr_aux(data, n, twiddles, 0, sign)
  end subroutine fft_dit_nr

  ! Decimation-in-time, natural-order input, bit-reversed order output:
  recursive subroutine fft_dit_nr_aux( dat, size, twiddles, &
                                       twiddle_index, sign )
    intrinsic ishft, conjg
    integer, intent(in)                             :: size
    !sll_comp64, dimension(0:size-1), intent(inout)  :: dat
    sll_comp64, dimension(0:),       intent(inout)  :: dat
    ! It is more convenient when the twiddles are 0-indexed
    sll_comp64, dimension(0:), intent(in)           :: twiddles
    integer, intent(in)                             :: twiddle_index
    integer, optional, intent(in)                   :: sign
    integer                                         :: half
    integer                                         :: j
    integer                                         :: t_index_new
    sll_comp64                                      :: omega
    sll_comp64                                      :: tmp

    half = ishft(size,-1) ! any evidence this is faster?
    !if ( sign == sll_p_fft_backward ) then
    !   omega = twiddles(twiddle_index)
    !else if ( sign == sll_p_fft_forward ) then
    !   omega = conjg(twiddles(twiddle_index))
    !else
    !   stop 'ERROR in =fft_dit_nr_aux= argument sign invalid'
    !end if
    omega = twiddles(twiddle_index)
    ! Do the butterflies for this stage
    do j=0,half-1
       tmp         = dat(j+half)*omega
       dat(j+half) = dat(j) - tmp
       dat(j)      = dat(j) + tmp
    end do
    ! Spawn two recursive calls for the two new problems that have been created:
    ! the upper and lower halves of the data
    if( half > 1) then
       t_index_new = ishft(twiddle_index,1)
       call fft_dit_nr_aux( dat(0:half-1),  &
                            half,           &
                            twiddles,       &
                            t_index_new,    &
                            sign )
       call fft_dit_nr_aux( dat(half:2*half-1), &
                            half,               &
                            twiddles,           &
                            t_index_new+1,      &
                            sign )
    else
       return
    end if
  end subroutine fft_dit_nr_aux

!PN DEFINED BUT NOT USED
!PN  subroutine fft_dit_rn( data, sign )
!PN    sll_comp64, dimension(:), intent(inout) :: data
!PN    sll_int32, intent(in)                   :: sign
!PN    sll_comp64, dimension(:), pointer       :: twiddles
!PN    integer                                 :: n
!PN    sll_int32                               :: ierr
!PN    n = size(data) ! bad
!PN    SLL_ASSERT(sll_f_is_power_of_two(int(n,i64)))
!PN    SLL_ALLOCATE(twiddles(n/2),ierr)
!PN    call compute_twiddles(n,twiddles)
!PN    ! This algorithm uses the twiddles in natural order. The '1'
!PN    ! argument is because fft_dit_rn_aux internally 1-indexes its
!PN    ! arrays, so we are just indicating the first twiddle factor.
!PN    call fft_dit_rn_aux(data, n, twiddles, 1, sign)
!PN  end subroutine fft_dit_rn

  recursive subroutine fft_dit_rn_aux( data,           &
                                       data_size,      &
                                       twiddles,       &
                                       twiddle_stride, &
                                       sign )
    intrinsic ishft, conjg
    sll_int32, intent(in)                               :: data_size
    ! looks like the most natural indexing for this algorithm is
    ! the 1-based indexing...
    sll_comp64, dimension(1:data_size), intent(inout)   :: data
    sll_comp64, dimension(1:data_size/2), intent(in)    :: twiddles
    integer, intent(in)                                 :: twiddle_stride
    sll_int32                                           :: new_stride
    sll_int32                                           :: jtwiddle
    integer, intent(in)                                 :: sign
    integer                                             :: half
    integer                                             :: j
    sll_comp64                                          :: omega
    sll_comp64                                          :: tmp
    jtwiddle = 1
    half = ishft(data_size,-1) ! is this be faster than a division in Fortran?
    if( half > 1) then
       new_stride = twiddle_stride*2
       call fft_dit_rn_aux(data(1:half), half, twiddles, new_stride, sign)
       call fft_dit_rn_aux(data(half+1:data_size),half,twiddles,new_stride,sign)
    end if
    ! Do the butterflies for this stage
    do j=1,half
       !if ( sign == sll_p_fft_forward ) then
       !   omega = conjg(twiddles(jtwiddle))
       !else if ( sign == sll_p_fft_backward ) then
       !   omega = twiddles(jtwiddle)
       !else
       !  stop 'ERROR in =fft_dit_rn_aux= argument sign invalid'
       !end if
       omega = twiddles(jtwiddle)
       tmp          = data(j+half)*omega
       data(j+half) = data(j) - tmp
       data(j)      = data(j) + tmp
       jtwiddle     = jtwiddle + twiddle_stride
    end do
  end subroutine fft_dit_rn_aux



  ! notes:
  ! - the real data is size n
  ! - but we apply a complex fft, treating the n-sized array as n/2
  !   complex numbers. This means that the twiddle array should be
  !   n/4-size (number of complex nums), but n/2 if seen as an array
  !   of reals.
  ! - this is badly named, it is not that the this is an FFT for real
  !   values, but that the argument itself is real and the FFT interprets
  !   it as complex. Change this confusing function name.
  ! Proposed interface for this version of the FFT function:
  ! - Treat it as a function for FFT of complex data, even though the
  !   argument is a real array. This implies:
  ! - The size of the problem is the number of complex numbers, not the
  !   number of elements in the real array. Thus, if the real array is
  !   length 'n', this corresponds to a problem size 'n/2', as this is the
  !   number of complex numbers that are represented.
  recursive subroutine fft_dit_nr_real_array_aux( samples,       &
                                                  num_complex,   &
                                                  twiddles,      &
                                                  twiddle_index, &
                                                  sign )
    intrinsic ishft, conjg
    ! num_complex is the number of complex numbers represented in the array.
    integer, intent(in)                                    :: num_complex
    sll_real64, dimension(0:2*num_complex-1), intent(inout):: samples
    ! It is more convenient when the twiddles are 0-indexed
    ! Declaring the size of the twiddles as in the following line, gives
    ! a runtime error as the array is accessed out of bounds.
    !    sll_real64, dimension(0:num_complex-1), intent(in)     :: twiddles
    sll_real64, dimension(0:), intent(in)     :: twiddles
    integer, intent(in)                                    :: twiddle_index
    integer, intent(in)                                    :: sign
    ! half represents half of the complex problem, not half of the array
    ! size.
    integer                                                :: half
    integer                                                :: j
    integer                                                :: t_index_new
    sll_real64                                             :: omega_re
    sll_real64                                             :: omega_im
    sll_real64                                             :: tmp_re
    sll_real64                                             :: tmp_im
    SLL_ASSERT(num_complex .le. size(samples))
    half = ishft(num_complex,-1) ! would this be faster than a division
                                 ! in Fortran?
    ! select the value of the twiddle factor for this stage depending on
    ! the direction of the transform
    !if ( sign == sll_p_fft_forward ) then
    !   omega_re =  CREAL0(twiddles, twiddle_index)
    !   omega_im = -CIMAG0(twiddles, twiddle_index)
    !else if ( sign == sll_p_fft_backward ) then
    !   omega_re =  CREAL0(twiddles, twiddle_index)
    !   omega_im =  CIMAG0(twiddles, twiddle_index)
    !else
    !   stop 'ERROR in =fft_dit_nr_real_array_aux= argument sign invalid'
    !end if
       omega_re =  CREAL0(twiddles, twiddle_index)
       omega_im =  CIMAG0(twiddles, twiddle_index)
    ! Do the butterflies for this stage
    do j=0,half-1
       ! multiply samples(j+half) and omega, store in tmp
       tmp_re = &
            CREAL0(samples,j+half)*omega_re - CIMAG0(samples,j+half)*omega_im
       tmp_im = &
            CREAL0(samples,j+half)*omega_im + CIMAG0(samples,j+half)*omega_re
       ! subtract tmp from samples(j), store in samples(j+half)
       CREAL0(samples,j+half) = CREAL0(samples,j) - tmp_re
       CIMAG0(samples,j+half) = CIMAG0(samples,j) - tmp_im
       ! add samples(j) and tmp, store in samples(j)
       CREAL0(samples,j) = CREAL0(samples,j) + tmp_re
       CIMAG0(samples,j) = CIMAG0(samples,j) + tmp_im
    end do
    ! Spawn two recursive calls for the two new problems that have been created:
    ! the upper and lower halves of the samples.
    if( half > 1) then
       t_index_new = ishft(twiddle_index,1)
       call fft_dit_nr_real_array_aux( samples(0:2*half-1),      &
                                       half,                     &
                                       twiddles,                 &
                                       t_index_new,              &
                                       sign )
       call fft_dit_nr_real_array_aux( samples(2*half:4*half-1), &
                                       half,                     &
                                       twiddles,                 &
                                       t_index_new+1,            &
                                       sign )
    else
       return
    end if
  end subroutine fft_dit_nr_real_array_aux

  recursive subroutine fft_dit_rn_real_array_aux( data,           &
                                                  num_complex,    &
                                                  twiddles,       &
                                                  twiddle_stride, &
                                                  sign )
    intrinsic ishft, conjg
    ! num_complex is the number of complex numbers represented in the array.
    integer, intent(in)                                    :: num_complex
    sll_real64, dimension(1:2*num_complex), intent(inout)  :: data
    ! For this reverse-to-natural order algorithm, it is more convenient
    ! to use the 1-based indexing. Note that this is a num_complex-sized FFT,
    ! which means that we need num_complex/2 twiddles, but since they are
    ! stored as reals, we need that the twiddle array be of size num_complex.
    sll_real64, dimension(1:num_complex), intent(in)       :: twiddles
    sll_int32, intent(in)                                  :: twiddle_stride
    sll_int32                                              :: new_stride
    sll_int32                                              :: jtwiddle
    integer, intent(in)                                    :: sign
    integer                                                :: half
    integer                                                :: j
    sll_real64                                             :: omega_re
    sll_real64                                             :: omega_im
    sll_real64                                             :: tmp_re
    sll_real64                                             :: tmp_im
    SLL_ASSERT(num_complex .le. size(data))
    jtwiddle = 1
    ! note that half represents half the complex problem size, not half
    ! the array size.
    half = ishft(num_complex,-1) ! would this be faster than a division
                                 ! in Fortran?
    if ( half>1 ) then
       new_stride = twiddle_stride*2
       call fft_dit_rn_real_array_aux( data(1:num_complex), &
                                       half,                &
                                       twiddles,            &
                                       new_stride,          &
                                       sign )
       call fft_dit_rn_real_array_aux( data(num_complex+1:2*num_complex), &
                                       half,                              &
                                       twiddles,                          &
                                       new_stride,                        &
                                       sign )
    end if
    ! Do the butterflies for this stage
    do j=1,half
       ! select the value of the twiddle factor for this stage depending on
       ! the direction of the transform
       !if( sign == sll_p_fft_forward ) then
       !   omega_re =  CREAL1(twiddles, jtwiddle)
       !   omega_im = -CIMAG1(twiddles, jtwiddle)
       !else if ( sign == sll_p_fft_backward ) then
       !   omega_re =  CREAL1(twiddles, jtwiddle)
       !   omega_im =  CIMAG1(twiddles, jtwiddle)
       !else
       !   stop 'ERROR in =fft_dit_rn_real_array_aux= argument sign invalid'
       !end if
       omega_re =  CREAL1(twiddles, jtwiddle)
       omega_im =  CIMAG1(twiddles, jtwiddle)
       ! multiply data(j+half) and omega, store in tmp
       tmp_re = CREAL1(data,j+half)*omega_re - CIMAG1(data,j+half)*omega_im
       tmp_im = CREAL1(data,j+half)*omega_im + CIMAG1(data,j+half)*omega_re
       ! subtract tmp from data(j), store in data(j+half)
       CREAL1(data,j+half) = CREAL1(data,j) - tmp_re
       CIMAG1(data,j+half) = CIMAG1(data,j) - tmp_im
       ! add data(j) and tmp, store in data(j)
       CREAL1(data,j) = CREAL1(data,j) + tmp_re
       CIMAG1(data,j) = CIMAG1(data,j) + tmp_im
       jtwiddle = jtwiddle + twiddle_stride
    end do
  end subroutine fft_dit_rn_real_array_aux


  ! Important stuff that should be checked by the function or preferably
  ! at the plan-level:
  ! - twiddles is needed by the call to the complex FFT that operates on
  !   real-formatted data. Therefore, if the real data to be processed is
  !   size 'n', then, in the first step this data will be interpreted as
  !   a complex array of size 'n/2'. This means that the twiddles array is
  !   composed of 'n/4' complex numbers or 'n/2' reals
  ! - This function should be allowed to operate in an 'out-of-place' way.
  !   In case that an 'in-place' transform is required, the data array
  !   should have at least size n+2 real elements. The reason for this is
  !   that the FFT of real data has the symmetry property:
  !
  !                      X_r = X_(N-r)*
  !
  !   which means that it requires storage for only N/2+1 independent modes.
  !   Two of these modes (0 and N/2) are real-valued. The remaining N/2-1 are
  !   complex valued. In terms of required memory:
  !
  !             2 real modes        =>   2   reals
  !             N/2-1 complex modes =>   N-2 reals
  !                                    +___________
  !                                       N  reals
  !
  !   Thus, in principle we could store all the information in the original
  !   real array. This has several inconveniences. Firstly, it would require
  !   packing the data in 'special' ways (the Numerical Recipes approach),
  !   like putting the real modes share a single complex-valued field. This
  !   requires special packing/unpacking and 'if' statements when it comes
  !   down to reading the modes. Another way could be to require the input
  !   data to be padded to allow storage for the extra modes. This is also
  !   very problematic, as it requires, from the moment of the allocation,
  !   some foresight that a given array will be the subject of an FFT
  !   operation. This may not be so bad... In any case, it seems that
  !   special logic will be required to deliver a requested Fourier mode,
  !   since in some cases one could read a value from an array but in other
  !   cases one needs the conjugate... For now, the present implementation
  !   chooses to pack the data in the minimum space possible, to the 'weird'
  !   packing and unpacking and use special logic to access the Fourier modes.
  subroutine real_data_fft_dit( data, n, twiddles, twiddles_n, sign )
    sll_int32, intent(in)                      :: n
    sll_real64, dimension(0:),   intent(inout) :: data
    sll_real64, dimension(0:n/2-1), intent(in) :: twiddles
    sll_real64, dimension(0:n-1),   intent(in) :: twiddles_n
    sll_int32, intent(in)                      :: sign
    sll_int32                                  :: n_2
    sll_int32                                  :: i
    sll_real64                                 :: hi_re
    sll_real64                                 :: hi_im
    sll_real64                                 :: hn_2mi_re
    sll_real64                                 :: hn_2mi_im
    sll_real64                                 :: tmp_re
    sll_real64                                 :: tmp_im
    sll_real64                                 :: tmp2_re
    sll_real64                                 :: tmp2_im
    sll_real64                                 :: tmp3_re
    sll_real64                                 :: tmp3_im
    sll_real64                                 :: omega_re
    sll_real64                                 :: omega_im
    sll_real64                                 :: s
    SLL_ASSERT( sll_f_is_power_of_two(int(n,i64)) )
    SLL_ASSERT( size(data) .ge. n )
    SLL_ASSERT( size(twiddles) .eq. n/2 )
    n_2 = n/2
    ! First step in computing the FFT of the real array is to launch
    ! a complex FFT on the data, interpreting the real-array data as
    ! complex numbers. In other words, if the data is:
    !
    ! [d1, d2, d3, d4, ..., dn]
    !
    ! We reinterpret each array element as the real or imaginary part
    ! of a sequence of complex numbers totalling in number to half
    ! the number of elements in the original array.
    !
    ! [re1, im1, re2, im2, ..., ren/2, imn/2]
    !
    ! The next part of the algorithm is to construct the modes of the
    ! real FFT from the results of this complex transform. The inverse
    ! transform is computed by essentially inverting the order of this
    ! process and changing some signs...
    if( sign .eq. sll_p_fft_forward ) then
       ! we use the following as the 'switch' to flip signs between
       ! FORWARD and INVERSE transforms.
       s = -1.0_f64
       ! The following call has to change to refer to its natural wrapper...
       call fft_dit_nr_real_array_aux( data(0:n-1), &
                                       n_2,         &
                                       twiddles,    &
                                       0,           &
                                       sll_p_fft_forward )
       ! but our _nr_ algorithm bit reverses the result, so, until we have
       ! some way to index the data correctly we have to do this:
       call bit_reverse_in_pairs( n_2, data(0:n-1) )
    else if (sign .eq. sll_p_fft_backward) then
       s =  1.0_f64
    else
      stop 'ERROR IN =REAL_DATA_FFT_DIT= invalid argument sign'
    end if
    do i=1,n_2/2 ! the index 'i' corresponds to indexing H
       ! sll_p_fft_forward case: We intend to mix the odd/even components that we
       ! have computed into complex numbers H_n. These Complex numbers will
       ! give the modes as in:
       !
       ! F_i = 1/2*(H_i + H_(N/2-i)^*) - i/2*(H_i-H_(N/2-i)^*)*exp(-j*2*pi*i/N)
       !
       ! which is the answer we are after.
       !
       ! sll_p_fft_backward case: The process of decomposing the H_n's into the
       ! corresponding even and odd terms:
       !
       ! (even terms:) F_n^e = (F_n + F_(N/2-n)^*)
       !
       ! (odd terms: ) F_n^o = (F_n - F_(N/2-n)^*)exp(j*2*pi*i/N)
       !
       ! followed by a complex transform to calculate
       !
       !   H_n = F_n^(1) + i*F_n^(2)
       hi_re              =  CREAL0(data,i)
       hi_im              =  CIMAG0(data,i)
       hn_2mi_re          =  CREAL0(data,n_2-i)
       hn_2mi_im          = -CIMAG0(data,n_2-i)
       omega_re           =  CREAL0(twiddles_n,i)
       omega_im           =  s*CIMAG0(twiddles_n,i) ! conjugated for FORWARD
       ! Compute tmp =  1/2*(H_i + H_(N/2-i)^*)
       tmp_re             =  0.5_f64*(hi_re + hn_2mi_re)
       tmp_im             =  0.5_f64*(hi_im + hn_2mi_im)
       ! Compute tmp2 = i/2*(H_n - H_(N/2-n)^*), the sign depends on the
       ! direction of the FFT.
       tmp2_re            =  s*0.5_f64*(hi_im - hn_2mi_im)
       tmp2_im            = -s*0.5_f64*(hi_re - hn_2mi_re)
       ! Multiply tmp2 by the twiddle factor and add to tmp.
       tmp3_re            =  tmp2_re*omega_re - tmp2_im*omega_im
       tmp3_im            =  tmp2_re*omega_im + tmp2_im*omega_re
       CREAL0(data,i)     =  tmp_re - tmp3_re
       CIMAG0(data,i)     =  tmp_im - tmp3_im
       CREAL0(data,n_2-i) =  tmp_re + tmp3_re
       CIMAG0(data,n_2-i) = -tmp_im - tmp3_im
    end do
    if ( sign .eq. sll_p_fft_forward ) then
       ! Set the first and N/2 values independently and pack them in the
       ! memory space provided by the original data.
       tmp_re = data(0)
       data(0) = tmp_re + data(1)   ! mode 0   is real
       data(1) = tmp_re - data(1)   ! mode N/2 is real
    else if ( sign .eq. sll_p_fft_backward ) then
       ! Unpack the modes.
       tmp_re  = data(0)
       data(0) = 0.5_f64*(tmp_re + data(1))
       data(1) = 0.5_f64*(tmp_re - data(1))
       ! The following call has to change to refere to its natural wrapper...
       call fft_dit_nr_real_array_aux( data(0:n-1), &
                                       n_2,         &
                                       twiddles,    &
                                       0,           &
                                       sll_p_fft_backward )
       ! but our _nr_ algorithm bit reverses the result, so, until we have
       ! some way to index the data correctly we have to do this:
       call bit_reverse_in_pairs( n_2, data(0:n-1) )
       data = data*2.0_f64
    end if
  end subroutine real_data_fft_dit

#undef CREAL0
#undef CIMAG0
#undef CREAL1
#undef CIMAG1

#endif /* DOXYGEN_SHOULD_SKIP_THIS */
end module sll_m_fft
