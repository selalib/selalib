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
!> @author ???? and Katharina Kormann (IPP)
!> @brief FFT interface for FFTW
module sll_m_fft
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#include "sll_working_precision.h"
#include "sll_assert.h"
#include "sll_memory.h"
#include "sll_errors.h"
#include "sll_fftw.h"

  use sll_m_fftw3
  use, intrinsic :: iso_c_binding

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
  
  !> Type for FFT plan in SeLaLib
  type sll_t_fft
    fftw_plan, private               :: fftw       !< FFTW plan (specific for FFTW version)
    logical                          :: normalized !< Boolean telling whether or not values of the FFT should be normalized by \a problem_shape
    sll_int32                        :: library    !< Specifies the FFT library that is behind the interface
    sll_int32                        :: direction  !< Direction of the FFT, either \a sll_p_fft_forward (negative sign in exponential) or \a sll_p_fft_backward (positive sign in exponential) 
    sll_int32                        :: problem_rank !< Dimension of FFT
    sll_int32, dimension(:), pointer :: problem_shape !< Array of size \a problem_rank specifying the number of points along each dimension

    sll_comp64, allocatable, private :: scratch(:) !< Scratch data to simulate r2r for not FFTW_F2003
    sll_int32, private               :: transform_type !< Type of the transform. Use for assertion to make sure execution is called of the same type as fft object was initialized for.

  end type sll_t_fft

  
  ! Flags for direction (values as in fftw3.f03)
  integer, parameter :: sll_p_fft_forward = FFTW_FORWARD !< Flag to specify forward FFT (negative sign) [value -1]
  integer, parameter :: sll_p_fft_backward = FFTW_BACKWARD !< Flag to specify backward FFT (positive sign) [value 1]

  ! Flags for initialization of the plan (values as in fftw3.f03)
  integer(C_INT), parameter :: sll_p_fft_measure = FFTW_MEASURE !< FFTW planning-rigor flag FFTW_MEASURE (optimized plan) NOTE: planner overwrites the input array during planning  [value 0]
  integer(C_INT), parameter :: sll_p_fft_patient = FFTW_PATIENT !< FFTW planning-rigor flag FFTW_PATIENT (more optimizaton than MEASURE) NOTE: planner overwrites the input array during planning  [value 32]
  integer(C_INT), parameter :: sll_p_fft_estimate = FFTW_ESTIMATE !< FFTW planning-rigor flag FFTW_ESTIMATE (simple heuristic for planer)   [value 64]. This is our default value
  integer(C_INT), parameter :: sll_p_fft_exhaustive = FFTW_EXHAUSTIVE !< FFTW planning-rigor flag FFTW_EXHAUSTIVE (more optimization than PATIENT) NOTE: planner overwrites the input array during planning  [value 8]
  integer(C_INT), parameter :: sll_p_fft_wisdom_only = FFTW_WISDOM_ONLY ! 2097152 !< FFTW planning-rigor flag FFTW_WISDOM_ONLY (planner only initialized if wisdom is available)  [value 2097152]
  
! Flags to pass when we create a new plan
! We can define 31 different flags.
! The value assigned to the flag can only be a power of two.
! See section "How-to manipulate flags ?" for more information.

  integer, parameter :: FFTW_MOD = 1000000000 !< Flag for FFTW

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

  !> Function to print the FFT library behind the interface.
  subroutine sll_s_fft_print_defaultlib()
    print *, 'The library used is FFTW'
  end subroutine


  !> Function to allocate an aligned complex array
  function sll_f_fft_allocate_aligned_complex(n) result(data)!, ptr_data, data)
    sll_int32,                 intent(in)  :: n                !< Size of the pointer
    complex(C_DOUBLE_COMPLEX), pointer     :: data(:) !< Array to be allocated
    
    type(C_PTR) :: ptr_data         ! C pointer needed for allocation
   
#ifdef FFTW_F2003
    ptr_data = fftw_alloc_complex(int(n,kind=C_SIZE_T))
    call c_f_pointer(ptr_data, data, [n])
#else
    allocate(data(n))
    SLL_WARNING('sll_f_fft_allocate_aligned_complex', 'Aligned allocation only supported by F2003. Usual allocation.')
#endif

  end function sll_f_fft_allocate_aligned_complex


  !> Function to allocate an aligned real array
  function sll_f_fft_allocate_aligned_real(n) result(data)!, ptr_data, data)
    sll_int32,                 intent(in)  :: n                !< Size of the pointer
    real(C_DOUBLE), pointer                :: data(:) !< Array to be allocated
    
    type(C_PTR) :: ptr_data         ! C pointer needed for allocation
   
#ifdef FFTW_F2003
    ptr_data = fftw_alloc_real(int(n,kind=C_SIZE_T))
    call c_f_pointer(ptr_data, data, [n])
#else
    allocate(data(n))
    SLL_WARNING('sll_f_fft_allocate_aligned_real','Aligned allocation only supported by F2003. Usual allocation.')
#endif

  end function sll_f_fft_allocate_aligned_real



  !> Function to reconstruct the complex FFT mode from the data of a r2r transform
  function sll_f_fft_get_mode_r2c_1d(plan,data,k) result(mode)
    type(sll_t_fft),      intent(in)   :: plan !< FFT plan
    sll_real64, dimension(0:), intent(in)   :: data !< real data produced by r2r transform
    sll_int32, intent(in)                   :: k    !< mode to be extracted
    sll_comp64                              :: mode !< Complex value of kth mode

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
    type(sll_t_fft),      intent(in)    :: plan !< FFT planner object
    sll_real64, dimension(0:), intent(out)   :: data !< Real array to be set
    sll_comp64, intent(in)                   :: new_value !< Complex value of the kth mode
    sll_int32, intent(in)                    :: k !< mode to be set

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

! COMPLEX
  !> Create new 1d complex to complex plan
  subroutine sll_s_fft_init_c2c_1d(plan, nx,array_in,array_out,direction,normalized, aligned, optimization)
    type(sll_t_fft)                         :: plan !< FFT planner object 
    sll_int32, intent(in)                        :: nx !< Number of points
    sll_comp64, dimension(:), intent(inout)      :: array_in !< (Typical) input array (gets overwritten for certain options)
    sll_comp64, dimension(:), intent(inout)      :: array_out !< (Typical) output array (gets overwritten for certain options)
    sll_int32, intent(in)                        :: direction  !< Direction of the FFT (\a sll_p_fft_forward or \a sll_p_fft_backward)
    logical, optional,   intent(in)              :: normalized !< Flag to decide if FFT should be normalized by 1/N (default: \a FALSE)
    logical, optional,   intent(in)              :: aligned    !< Flag to decide if FFT routine can assume data alignment (default: \a FALSE). Note that you need to call an aligned initialization if you want to set this option to \a TRUE. 
    sll_int32, optional, intent(in)              :: optimization !< Planning-rigor flag for FFTW. Possible values \a sll_p_fft_estimate, \a sll_p_fft_measure, \a sll_p_fft_patient, \a sll_p_fft_exhaustive, \a sll_p_fft_wisdom_only. (default: \a sll_p_fft_estimate).

    ! local variables
    sll_int32                                    :: ierr
    sll_int32                                    :: flag_fftw

    plan%library = FFTW_MOD
    plan%transform_type = p_fftw_c2c_1d
    plan%direction = direction
    if( present(normalized) ) then
       plan%normalized = normalized
    else
       plan%normalized = .false.
    end if
    ! Set the information about the algorithm to compute the plan. The default is FFTW_ESTIMATE
    if ( present(optimization) ) then
       flag_fftw = optimization
    else
       flag_fftw = FFTW_ESTIMATE
    end if
    if ( present(aligned) ) then
       if (aligned .EQV. .false.) then
          flag_fftw = flag_fftw + FFTW_UNALIGNED
       end if
    else
       flag_fftw = flag_fftw + FFTW_UNALIGNED
    end if
    plan%problem_rank = 1
    SLL_ALLOCATE(plan%problem_shape(1),ierr)
    plan%problem_shape = (/ nx  /)

#ifdef FFTW_F2003
    plan%fftw = fftw_plan_dft_1d(nx,array_in,array_out,direction,flag_fftw)
#else
    call dfftw_plan_dft_1d(plan%fftw,nx,array_in,array_out,direction,flag_fftw)
#endif

  end subroutine sll_s_fft_init_c2c_1d

  !> Compute fast Fourier transform in complex to complex mode.
  subroutine sll_s_fft_exec_c2c_1d(plan,array_in,array_out)
    type(sll_t_fft),     intent(in)            :: plan !< FFT planner object
    sll_comp64, dimension(:), intent(inout)         :: array_in !< Complex data to be Fourier transformed
    sll_comp64, dimension(:), intent(inout)         :: array_out !< Fourier coefficients on output

    sll_real64 :: factor

    SLL_ASSERT( plan%transform_type == p_fftw_c2c_1d )

    call fftw_execute_dft(plan%fftw, array_in, array_out)

    if ( plan%normalized .EQV. .true.) then
       factor = 1.0_f64/real(plan%problem_shape(1),kind=f64)
       array_out = factor*array_out
    endif
  end subroutine 
! END COMPLEX

! COMPLEX 2D
  !> Create new 2d complex to complex plan
  subroutine sll_s_fft_init_c2c_2d(plan, nx, ny, array_in, array_out, direction, normalized, aligned, optimization)
    sll_int32, intent(in)                        :: nx !< Number of points along first dimension
    sll_int32, intent(in)                        :: ny !< Number of points along second dimension
    sll_comp64, dimension(0:,0:), intent(inout)  :: array_in !< (Typical) input array (gets overwritten for certain options)
    sll_comp64, dimension(0:,0:), intent(inout)  :: array_out !<(Typical) output array (gets overwritten for certain options)
    sll_int32, intent(in)                        :: direction  !< Direction of the FFT (\a sll_p_fft_forward or \a sll_p_fft_backward)
    logical, optional,   intent(in)              :: normalized !< Flag to decide if FFT should be normalized by 1/N (default: \a FALSE)
    logical, optional,   intent(in)              :: aligned    !< Flag to decide if FFT routine can assume data alignment (default: \a FALSE). Not that you need to call an aligned initialization if you want to set this option to \a TRUE.
    sll_int32, optional, intent(in)              :: optimization !< Planning-rigor flag for FFTW. Possible values \a sll_p_fft_estimate, \a sll_p_fft_measure, \a sll_p_fft_patient, \a sll_p_fft_exhaustive, \a sll_p_fft_wisdom_only. (default: \a sll_p_fft_estimate). Note that you need to 
    type(sll_t_fft), intent(out)             :: plan !< initialized planner object
 
    sll_int32                                     :: ierr
    sll_int32                                     :: flag_fftw

    plan%library = FFTW_MOD
    plan%transform_type = p_fftw_c2c_2d
    plan%direction = direction
    if( present(normalized) ) then
       plan%normalized = normalized
    else
       plan%normalized = .false.
    end if
    ! Set the information about the algorithm to compute the plan. The default is FFTW_ESTIMATE
    if ( present(optimization) ) then
       flag_fftw = optimization
    else
       flag_fftw = FFTW_ESTIMATE
    end if
    if ( present(aligned) ) then
       if (aligned .EQV. .false.) then
          flag_fftw = flag_fftw + FFTW_UNALIGNED
       end if
    else
       flag_fftw = flag_fftw + FFTW_UNALIGNED
    end if
    plan%problem_rank = 2
    SLL_ALLOCATE(plan%problem_shape(2),ierr)
    plan%problem_shape = (/ nx , ny /)
  
    !We must switch the dimension. It's a fftw convention. 

#ifdef FFTW_F2003
    plan%fftw = fftw_plan_dft_2d(ny,nx,array_in,array_out,direction,flag_fftw)!FFTW_ESTIMATE + FFTW_UNALIGNED)
#else
    call dfftw_plan_dft_2d(plan%fftw,ny,nx,array_in,array_out,direction,flag_fftw)!FFTW_ESTIMATE + FFTW_UNALIGNED)
#endif

  end subroutine sll_s_fft_init_c2c_2d


  !> Compute fast Fourier transform in complex to complex mode.
  subroutine sll_s_fft_exec_c2c_2d(plan,array_in,array_out)
    type(sll_t_fft),  intent(in)            :: plan !< FFT planner object
    sll_comp64, dimension(0:,0:), intent(inout)  :: array_in !< Complex data to be Fourier transformed 
    sll_comp64, dimension(0:,0:), intent(inout)  :: array_out !< Fourier coefficients on output

    sll_real64                                   :: factor

    call fftw_execute_dft(plan%fftw, array_in, array_out)

    if ( plan%normalized .EQV. .true.) then
      factor = 1.0_f64/real(plan%problem_shape(1)*plan%problem_shape(2),kind=f64)
      array_out = factor*array_out
    endif

  end subroutine 
! END COMPLEX 2D

! REAL
  !> Create new 1d real to real plan
  subroutine sll_s_fft_init_r2r_1d(plan,nx,array_in,array_out,direction,normalized, aligned, optimization)
    type(sll_t_fft)                         :: plan !< FFT planner object
    sll_int32, intent(in)                        :: nx !< Number of points
    sll_real64, dimension(:), intent(inout)      :: array_in !< (Typical) input array (gets overwritten for certain options)
    sll_real64, dimension(:), intent(inout)      :: array_out !< (Typical) output array (gets overwritten for certain options)
    sll_int32, intent(in)                        :: direction  !< Direction of the FFT (\a sll_p_fft_forward or \a sll_p_fft_backward)
    logical, optional,   intent(in)              :: normalized !< Flag to decide if FFT should be normalized by 1/N (default: \a FALSE)
    logical, optional,   intent(in)              :: aligned    !< Flag to decide if FFT routine can assume data alignment (default: \a FALSE). Not that you need to call an aligned initialization if you want to set this option to \a TRUE.
    sll_int32, optional, intent(in)              :: optimization !< Planning-rigor flag for FFTW. Possible values \a sll_p_fft_estimate, \a sll_p_fft_measure, \a sll_p_fft_patient, \a sll_p_fft_exhaustive, \a sll_p_fft_wisdom_only. (default: \a sll_p_fft_estimate). Note that you need to 
    
    sll_int32 :: ierr
    sll_int32 :: flag_fftw

    plan%library = FFTW_MOD
    plan%transform_type = p_fftw_r2r_1d
    plan%direction = direction
    if( present(normalized) ) then
       plan%normalized = normalized
    else
       plan%normalized = .false.
    end if
    ! Set the information about the algorithm to compute the plan. The default is FFTW_ESTIMATE
    if ( present(optimization) ) then
       flag_fftw = optimization
    else
       flag_fftw = FFTW_ESTIMATE
    end if
    if ( present(aligned) ) then
       if (aligned .EQV. .false.) then
          flag_fftw = flag_fftw + FFTW_UNALIGNED
       end if
    else
       flag_fftw = flag_fftw + FFTW_UNALIGNED
    end if
    plan%problem_rank = 1
    SLL_ALLOCATE(plan%problem_shape(1),ierr)
    plan%problem_shape = (/ nx /)

    if(direction .eq. sll_p_fft_forward) then
#ifdef FFTW_F2003
      plan%fftw = fftw_plan_r2r_1d(nx,array_in,array_out,FFTW_R2HC,flag_fftw)
#else      
      ! This would be the call to actually create a r2r plan
      !call dfftw_plan_r2r_1d(plan%fftw,nx,array_in,array_out,FFTW_R2HC,flag_fftw)

      allocate(plan%scratch(nx/2+1))
      call dfftw_plan_dft_r2c_1d(plan%fftw,nx,array_in,plan%scratch,flag_fftw)
      SLL_WARNING('sll_s_fft_init_r2r_1d','R2HC not supported without FFTW_F2003. Use implementation based on r2c/c2r transform.')
#endif
    else if(direction .eq. sll_p_fft_backward) then
#ifdef FFTW_F2003
      plan%fftw = fftw_plan_r2r_1d(nx,array_in,array_out,FFTW_HC2R,flag_fftw)
#else
      ! This would be the code to actually create an r2r plan
      !call dfftw_plan_r2r_1d(plan%fftw,nx,array_in,array_out,FFTW_HC2R,flag_fftw)

      allocate(plan%scratch(nx/2+1))
      call dfftw_plan_dft_c2r_1d(plan%fftw,nx,plan%scratch,array_out,flag_fftw)
      SLL_WARNING('sll_s_fft_init_r2r_1d','HC2R not supported without FFTW_F2003. Use implementation based on r2c/c2r transform.')
#endif
    endif
  end subroutine sll_s_fft_init_r2r_1d

  !> Compute fast Fourier transform in real to real mode.
  subroutine sll_s_fft_exec_r2r_1d(plan,array_in,array_out)

    type(sll_t_fft),     intent(inout) :: plan !< FFT planner object
    sll_real64, dimension(:), intent(inout) :: array_in !< Real data to be Fourier transformed
    sll_real64, dimension(:), intent(inout) :: array_out !< Fourier coefficients in real form (sin/cos coefficients)

    sll_real64                              :: factor

    
    SLL_ASSERT( plan%transform_type == p_fftw_r2r_1d )

#ifdef FFTW_F2003
    call fftw_execute_r2r(plan%fftw, array_in, array_out)
    
    if( plan%normalized .EQV. .TRUE. ) then
      factor = 1.0_f64/real(plan%problem_shape(1),kind=f64)
      array_out = factor*array_out
    endif

#else

    if (plan%direction == sll_p_fft_forward) then
       call sll_s_fft_exec_r2c_1d(plan, array_in, plan%scratch)
       array_out(1) = real(plan%scratch(1), f64)
       array_out(plan%problem_shape(1)/2+1) =  real(plan%scratch((plan%problem_shape(1)/2+1)),f64)
       do j=1, plan%problem_shape(1)/2-1
          call sll_s_fft_set_mode_c2r_1d(plan, array_out, plan%scratch(j+1), j)
       end do
    elseif (plan%direction == sll_p_fft_backward) then
       do j=0, plan%problem_shape(1)/2
          plan%scratch(j+1) = sll_f_fft_get_mode_r2c_1d(plan, array_out, j)
       end do
       call sll_s_fft_exec_c2r_1d(plan, plan%scratch ,array_out)
    end if

#endif

  end subroutine


! R2C
  !> Create new 1d real to complex plan for forward FFT
  subroutine sll_s_fft_init_r2c_1d(plan, nx,array_in,array_out, normalized, aligned, optimization)
    type(sll_t_fft), intent(out)            :: plan !< FFT planner object
    sll_int32, intent(in)                        :: nx !< Number of points
    sll_real64, dimension(:), intent(inout)      :: array_in !< (Typical) input array (gets overwritten for certain options)
    sll_comp64, dimension(:), intent(out)        :: array_out !< (Typical) output array (gets overwritten for certain options)
    logical, optional,   intent(in)              :: normalized !< Flag to decide if FFT should be normalized by 1/N (default: \a FALSE)
    logical, optional,   intent(in)              :: aligned    !< Flag to decide if FFT routine can assume data alignment (default: \a FALSE). Not that you need to call an aligned initialization if you want to set this option to \a TRUE.
    sll_int32, optional, intent(in)              :: optimization !< Planning-rigor flag for FFTW. Possible values \a sll_p_fft_estimate, \a sll_p_fft_measure, \a sll_p_fft_patient, \a sll_p_fft_exhaustive, \a sll_p_fft_wisdom_only. (default: \a sll_p_fft_estimate). Note that you need to 

    sll_int32 :: ierr
    sll_int32 :: flag_fftw

    plan%library = FFTW_MOD
    plan%transform_type = p_fftw_r2c_1d
    plan%direction = 0
    if( present(normalized) ) then
       plan%normalized = normalized
    else
       plan%normalized = .false.
    end if
    ! Set the information about the algorithm to compute the plan. The default is FFTW_ESTIMATE
    if ( present(optimization) ) then
       flag_fftw = optimization
    else
       flag_fftw = FFTW_ESTIMATE
    end if
    if ( present(aligned) ) then
       if (aligned .EQV. .false.) then
          flag_fftw = flag_fftw + FFTW_UNALIGNED
       end if
    else
       flag_fftw = flag_fftw + FFTW_UNALIGNED
    end if
    plan%problem_rank = 1
    SLL_ALLOCATE(plan%problem_shape(1),ierr)
    plan%problem_shape = (/ nx /)

#ifdef FFTW_F2003
    plan%fftw = fftw_plan_dft_r2c_1d(nx,array_in,array_out, flag_fftw)
#else
    call dfftw_plan_dft_r2c_1d(plan%fftw,nx,array_in,array_out,flag_fftw)
#endif
  end subroutine sll_s_fft_init_r2c_1d
  
  !> Compute fast Fourier transform in real to complex mode.
  subroutine sll_s_fft_exec_r2c_1d(plan,array_in,array_out)
    type(sll_t_fft),     intent(in)            :: plan !< FFT planner object
    sll_real64, dimension(:), intent(inout)         :: array_in !< Real input data to be Fourier transformed
    sll_comp64, dimension(:), intent(out)           :: array_out !< Complex Fourier mode (only first half due to symmetry)

    sll_real64                                      :: factor 

    SLL_ASSERT( plan%transform_type == p_fftw_r2c_1d )

    call fftw_execute_dft_r2c(plan%fftw, array_in, array_out)

    if( plan%normalized .EQV. .TRUE. ) then
      factor = 1.0_f64/real(plan%problem_shape(1),kind=f64)
      array_out = factor*array_out
    endif
  end subroutine

  !> Create new 2d complex to real plan for forward FFT
  subroutine sll_s_fft_init_r2c_2d(plan,nx,ny,array_in,array_out,normalized, aligned, optimization)
    type(sll_t_fft), intent(out)             :: plan !< FFT planner object
    sll_int32, intent(in)                        :: nx !< Number of points along first dimension
    sll_int32, intent(in)                        :: ny !< Number of points along second dimension
    sll_real64, dimension(:,:), intent(inout)    :: array_in !< (Typical) input array (gets overwritten for certain options)
    sll_comp64, dimension(:,:), intent(out)      :: array_out !< (Typical) output array (gets overwritten for certain options)
    logical, optional,   intent(in)              :: normalized !< Flag to decide if FFT should be normalized by 1/N (default: \a FALSE)
    logical, optional,   intent(in)              :: aligned    !< Flag to decide if FFT routine can assume data alignment (default: \a FALSE). Not that you need to call an aligned initialization if you want to set this option to \a TRUE.
    sll_int32, optional, intent(in)               :: optimization !< Planning-rigor flag for FFTW. Possible values \a sll_p_fft_estimate, \a sll_p_fft_measure, \a sll_p_fft_patient, \a sll_p_fft_exhaustive, \a sll_p_fft_wisdom_only. (default: \a sll_p_fft_estimate). Note that you need to 

    sll_int32 :: ierr    
    sll_int32 :: flag_fftw

    plan%library = FFTW_MOD
    plan%transform_type = p_fftw_r2c_2d
    plan%direction = 0
    if( present(normalized) ) then
       plan%normalized = normalized
    else
       plan%normalized = .false.
    end if
    ! Set the information about the algorithm to compute the plan. The default is FFTW_ESTIMATE
    if ( present(optimization) ) then
       flag_fftw = optimization
    else
       flag_fftw = FFTW_ESTIMATE
    end if
    if ( present(aligned) ) then
       if (aligned .EQV. .false.) then
          flag_fftw = flag_fftw + FFTW_UNALIGNED
       end if
    else
       flag_fftw = flag_fftw + FFTW_UNALIGNED
    end if
      
    plan%problem_rank = 2
    SLL_ALLOCATE(plan%problem_shape(2),ierr)
    plan%problem_shape = (/ nx , ny /)

#ifdef FFTW_F2003
    plan%fftw = fftw_plan_dft_r2c_2d(ny,nx,array_in,array_out,FFTW_ESTIMATE + FFTW_UNALIGNED)
#else
    call dfftw_plan_dft_r2c_2d(plan%fftw,nx,ny,array_in,array_out,FFTW_ESTIMATE)
#endif

  end subroutine sll_s_fft_init_r2c_2d

  !> Compute fast Fourier transform in real to complex mode.
  subroutine sll_s_fft_exec_r2c_2d(plan,array_in,array_out)
    type(sll_t_fft),       intent(in)           :: plan      !< FFT planner object
    sll_real64, dimension(:,:), intent(inout)         :: array_in  !< Real input data to be Fourier transformed
    sll_comp64, dimension(:,:), intent(out)           :: array_out !< Complex Fourier coefficients (only half part along first dimension due to symmetry)
    sll_int32     :: nx, ny
    sll_real64    :: factor


    SLL_ASSERT( plan%transform_type == p_fftw_r2c_2d )

    nx = plan%problem_shape(1)
    ny = plan%problem_shape(2)

    if( (size(array_in,dim=1).ne.nx) .and. (size(array_in,dim=2).ne.ny) ) then
      print * ,'Error in file sll_m_fft.F90'
      print * ,'      in subroutine fft_apply_r2c_2d'
      print * ,'      array_in size problem'
      stop ''
    else if( size(array_in,dim=1).ne.nx/2+1 .and. size(array_in,dim=2).ne.ny ) then
      print * ,'Error in file sll_m_fft.F90'
      print * ,'      in subroutine fft_apply_r2c_2d'
      print * ,'      array_out size problem'
      stop ''
    endif


    call fftw_execute_dft_r2c(plan%fftw, array_in, array_out)

    if( plan%normalized .EQV. .TRUE. ) then
      factor = 1.0_f64/real(nx*ny,kind=f64)
      array_out = factor*array_out
    endif
  end subroutine

!END R2C

! C2R
  !> Create new 1d complex to real plan for backward FFT
  subroutine sll_s_fft_init_c2r_1d(plan,nx,array_in,array_out, normalized, aligned, optimization)
    type(sll_t_fft)                         :: plan !< FFT planner object
    sll_int32, intent(in)                        :: nx !< Number of points
    sll_comp64, dimension(:)                     :: array_in  !< (Typical) input array (gets overwritten for certain options)
    sll_real64, dimension(:)                     :: array_out !< (Typical) output array (gets overwritten for certain options)
    logical, optional,   intent(in)              :: normalized !< Flag to decide if FFT should be normalized by 1/N (default: \a FALSE)
    logical, optional,   intent(in)              :: aligned    !< Flag to decide if FFT routine can assume data alignment (default: \a FALSE). Not that you need to call an aligned initialization if you want to set this option to \a TRUE.
    sll_int32, optional, intent(in)              :: optimization !< Planning-rigor flag for FFTW. Possible values \a sll_p_fft_estimate, \a sll_p_fft_measure, \a sll_p_fft_patient, \a sll_p_fft_exhaustive, \a sll_p_fft_wisdom_only. (default: \a sll_p_fft_estimate). Note that you need to 
    
    sll_int32 :: ierr
    sll_int32 :: flag_fftw

    plan%library = FFTW_MOD
    plan%transform_type = p_fftw_c2r_1d
    plan%direction = 0
    if( present(normalized) ) then
       plan%normalized = normalized
    else
       plan%normalized = .false.
    end if
    ! Set the information about the algorithm to compute the plan. The default is FFTW_ESTIMATE
    if ( present(optimization) ) then
       flag_fftw = optimization
    else
       flag_fftw = FFTW_ESTIMATE
    end if
    if ( present(aligned) ) then
       if (aligned .EQV. .false.) then
          flag_fftw = flag_fftw + FFTW_UNALIGNED
       end if
    else
       flag_fftw = flag_fftw + FFTW_UNALIGNED
    end if
    plan%problem_rank = 1
    SLL_ALLOCATE(plan%problem_shape(1),ierr)
    plan%problem_shape = (/ nx /)

#ifdef FFTW_F2003
    plan%fftw = fftw_plan_dft_c2r_1d(nx,array_in,array_out, flag_fftw)
#else
    call dfftw_plan_dft_c2r_1d(plan%fftw,nx,array_in,array_out,flag_fftw)
#endif

  end subroutine sll_s_fft_init_c2r_1d


  !> Compute fast Fourier transform in complex to real mode.
  subroutine sll_s_fft_exec_c2r_1d(plan,array_in,array_out)
    type(sll_t_fft),         intent(in)    :: plan !< FFT planner objece
    sll_comp64, dimension(:),    intent(inout) :: array_in !< Complex Fourier coefficient to be transformed back
    sll_real64, dimension(:),    intent(inout) :: array_out !< Real result of Fourier transform

    sll_real64                  :: factor

    SLL_ASSERT( plan%transform_type == p_fftw_c2r_1d )

    call fftw_execute_dft_c2r(plan%fftw, array_in, array_out)

    if( plan%normalized .EQV. .TRUE. ) then
      factor = 1.0_f64/real(plan%problem_shape(1),kind=f64)
      array_out = factor*array_out
    endif
  end subroutine


  !> Create new 2d real to complex plan for backward FFT
  !function sll_f_fft_new_plan_c2r_2d(nx,ny,array_in,array_out,normalized, aligned, optimization) result(plan)
  subroutine sll_s_fft_init_c2r_2d(plan,nx,ny,array_in,array_out,normalized, aligned, optimization)
    type(sll_t_fft), intent(out)            :: plan !< FFT planner object
    sll_int32, intent(in)                        :: nx !< Number of point along first dimension
    sll_int32, intent(in)                        :: ny !< Number of points along second dimension
    sll_comp64, dimension(:,:), intent(inout)    :: array_in !< (Typical) input array (gets overwritten for certain options)
    sll_real64, dimension(:,:), intent(out)      :: array_out !< (Typical) output array (gets overwritten for certain options)
    logical, optional,   intent(in)              :: normalized !< Flag to decide if FFT should be normalized by 1/N (default: \a FALSE)
    logical, optional,   intent(in)              :: aligned    !< Flag to decide if FFT routine can assume data alignment (default: \a FALSE). Not that you need to call an aligned initialization if you want to set this option to \a TRUE.
    sll_int32, optional, intent(in)              :: optimization !< Planning-rigor flag for FFTW. Possible values \a sll_p_fft_estimate, \a sll_p_fft_measure, \a sll_p_fft_patient, \a sll_p_fft_exhaustive, \a sll_p_fft_wisdom_only. (default: \a sll_p_fft_estimate). Note that you need to 

    sll_int32 :: ierr
    sll_int32 :: flag_fftw

    plan%library = FFTW_MOD
    plan%transform_type = p_fftw_c2r_2d
    plan%direction = 0
    if( present(normalized) ) then
       plan%normalized = normalized
    else
       plan%normalized = .false.
    end if
    ! Set the information about the algorithm to compute the plan. The default is FFTW_ESTIMATE
    if ( present(optimization) ) then
       flag_fftw = optimization
    else
       flag_fftw = FFTW_ESTIMATE
    end if
    if ( present(aligned) ) then
       if (aligned .EQV. .false.) then
          flag_fftw = flag_fftw + FFTW_UNALIGNED
       end if
    else
       flag_fftw = flag_fftw + FFTW_UNALIGNED
    end if
      
    plan%problem_rank = 2
    SLL_ALLOCATE(plan%problem_shape(2),ierr)
    plan%problem_shape = (/ nx , ny /)
#ifdef FFTW_F2003
    plan%fftw = fftw_plan_dft_c2r_2d(ny,nx,array_in,array_out,FFTW_ESTIMATE + FFTW_UNALIGNED)
#else
    call dfftw_plan_dft_c2r_2d(plan%fftw,nx,ny,array_in,array_out,FFTW_ESTIMATE)
#endif
  end subroutine sll_s_fft_init_c2r_2d
  
  !> Compute fast Fourier transform in complex to real mode.
  subroutine sll_s_fft_exec_c2r_2d(plan,array_in,array_out)
    type(sll_t_fft),         intent(in)          :: plan      !< FFT planner object
    sll_comp64, dimension(1:,1:), intent(inout)       :: array_in  !< Complex Fourier coefficient to be transformed back
    sll_real64, dimension(1:,1:), intent(out)         :: array_out !< Real output of Fourier transform

    sll_int32                                         :: nx, ny
    sll_real64                                        :: factor

    SLL_ASSERT( plan%transform_type == p_fftw_c2r_2d )

    nx = plan%problem_shape(1)
    ny = plan%problem_shape(2)
    call fftw_execute_dft_c2r(plan%fftw, array_in(1:nx/2+1,1:ny), array_out(1:nx,1:ny) )

    if( plan%normalized .EQV. .TRUE. ) then
      factor = 1.0_f64/real(nx*ny,kind=f64)
      array_out = factor*array_out
    endif
  end subroutine
!END C2R
! FFTW FFTW FFTW FFTW FFTW FFTW FFTW FFTW FFTW FFTW FFTW FFTW FFTW FFTW
! FFTW FFTW FFTW FFTW FFTW FFTW FFTW FFTW FFTW FFTW FFTW FFTW FFTW FFTW
! FFTW FFTW FFTW FFTW FFTW FFTW FFTW FFTW FFTW FFTW FFTW FFTW FFTW FFTW

  !> Delete a plan
  subroutine sll_s_fft_free(plan)
   type(sll_t_fft) :: plan !< FFT planner object

   sll_int32 :: ierr

    call fftw_destroy_plan(plan%fftw)
    if(associated(plan%problem_shape)) then
      SLL_DEALLOCATE(plan%problem_shape,ierr)
    endif
    
  end subroutine

end module sll_m_fft
