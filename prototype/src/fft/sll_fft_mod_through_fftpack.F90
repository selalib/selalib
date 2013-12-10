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

module sll_fft
#include "sll_working_precision.h"
#include "sll_utilities.h"
#include "sll_assert.h"
#include "sll_memory.h"

  use sll_fft_utils

  implicit none 

  type sll_fft_plan
    sll_real64, dimension(:), pointer :: twiddles => null()  ! twiddles factors
    sll_int32                         :: style
    sll_int32                         :: library
    sll_int32                         :: direction
    sll_int32                         :: problem_rank
    sll_int32, dimension(:), pointer  :: problem_shape
  end type sll_fft_plan

  interface fft_new_plan
    module procedure fftpack_new_plan_c2c_1d, fftpack_new_plan_r2r_1d, &
                    ! à implémenter
                     fftpack_new_plan_r2c_1d, fftpack_new_plan_c2r_1d , &
                     fftpack_new_plan_c2r_2d, fftpack_new_plan_r2c_2d, &
                     fftpack_new_plan_c2c_2d
  end interface

  interface fft_apply_plan
    module procedure fftpack_apply_plan_c2c_1d, fftpack_apply_plan_r2r_1d, &
                     fftpack_apply_plan_r2c_1d, fftpack_apply_plan_c2r_1d , &
                     fftpack_apply_plan_c2r_2d, fftpack_apply_plan_r2c_2d , &
                     fftpack_apply_plan_c2c_2d
  end interface 
  
  integer, parameter :: FFT_FORWARD = -1
  integer, parameter :: FFT_INVERSE = 1


! Flags to pass when we create a new plan
! We can define 31 different flags.
! The value assigned to the flag can only be a power of two.
! See section "How-to manipulate flags ?" for more information.
  integer, parameter :: FFT_NORMALIZE_FORWARD = 2**0
  integer, parameter :: FFT_NORMALIZE_INVERSE = 2**0
  integer, parameter :: FFT_NORMALIZE         = 2**0
  integer, parameter :: FFT_ONLY_FIRST_DIRECTION  = 2**2
  integer, parameter :: FFT_ONLY_SECOND_DIRECTION = 2**3
  integer, parameter :: FFT_ONLY_THIRD_DIRECTION  = 2**4

  integer, parameter :: FFTPACK_MOD = 100


  interface fft_get_mode
     module procedure fft_get_mode_complx_1d, fft_get_mode_real_1d
  end interface

  interface fft_set_mode
     module procedure fft_set_mode_complx_1d, fft_set_mode_real_1d
  end interface

contains



  subroutine print_defaultfftlib()
    print *, 'The library used is FFTPACK'
  end subroutine


  function fft_get_mode_complx_1d(plan,array,k) result(mode)
    type(sll_fft_plan), pointer :: plan
    sll_comp64, dimension(0:)   :: array
    sll_int32                   :: k
    sll_comp64                  :: mode
    mode = array(k)
  end function

  function fft_get_mode_real_1d(plan,data,k) result(mode)
    type(sll_fft_plan), pointer :: plan
    sll_real64, dimension(0:)   :: data
    sll_int32                   :: k, n_2, n
    sll_comp64                  :: mode

    n = plan%problem_shape(1)
    n_2 = n/2 !ishft(n,-1)

    if( k .eq. 0 ) then
      mode = cmplx(data(0),0.0_f64,kind=f64)
    else if( k .eq. n_2 ) then
      mode = cmplx(data(n-1),0.0_f64,kind=f64)
    else if( k .gt. n_2 ) then
      mode = cmplx( data(2*(n-k)-1) , -data(2*(n-k)) ,kind=f64)
    else
      mode = cmplx( data(2*k-1) , data(2*k) ,kind=f64)
    endif
  end function

  subroutine fft_set_mode_complx_1d(plan,array,new_value,k)
    type(sll_fft_plan), pointer :: plan
    sll_comp64, dimension(0:)   :: array
    sll_int32                   :: k
    sll_comp64                  :: new_value
    array(k) = new_value
  end subroutine

  subroutine fft_set_mode_real_1d(plan,data,new_value,k)
    type(sll_fft_plan), pointer :: plan
    sll_real64, dimension(0:)   :: data
    sll_int32                   :: k, n_2, n, index_mode
    sll_comp64                  :: new_value

    n = plan%problem_shape(1)
    n_2 = n/2 !ishft(n,-1)

    if( k .eq. 0 ) then
      data(0) = real(new_value,kind=f64)
    else if( k .eq. n_2 ) then
      data(n-1) = real(new_value,kind=f64)
    else if( k .gt. n_2 ) then
      data(2*(n-k)-1) = real(new_value,kind=f64)
      data(2*(n-k)) = -dimag(new_value)
    else
      data(2*k-1) = real(new_value,kind=f64)
      data(2*k) = dimag(new_value)
    endif
  end subroutine 

  ! return the index mode of ith stored mode
  ! In the complex output case the mode are stored in the natural order
  !     X_0,X_1,...,X_N-1
  ! In the real output case the order is
  !     r_o,r_1,i_1,r_2,i_2,...,r_N/2-1,i_N/2-1,r_N/2
  ! where X_k is the complex number (r_k,i_k).
  ! X_o and X_N/2 are purely real.
  function fft_ith_stored_mode(plan,i)
    type(sll_fft_plan), pointer :: plan
    sll_int32                   :: i, fft_ith_stored_mode
    fft_ith_stored_mode = i
  end function fft_ith_stored_mode

! FFTPACK FFTPACK FFTPACK FFTPACK FFTPACK FFTPACK FFTPACK FFTPACK FFTPACK
! FFTPACK FFTPACK FFTPACK FFTPACK FFTPACK FFTPACK FFTPACK FFTPACK FFTPACK
! FFTPACK FFTPACK FFTPACK FFTPACK FFTPACK FFTPACK FFTPACK FFTPACK FFTPACK
! COMPLEX
  function fftpack_new_plan_c2c_1d(nx,array_in,array_out,direction,flags) result(plan)
    sll_int32, intent(in)                        :: nx
    sll_comp64, dimension(:)                     :: array_in, array_out
    sll_int32, intent(in)                        :: direction
    sll_int32, optional,  intent(in)             :: flags
    type(sll_fft_plan), pointer                  :: plan
    sll_int32                                    :: ierr

    SLL_ALLOCATE(plan,ierr)
    plan%library = FFTPACK_MOD 
    plan%direction = direction
    if( present(flags) )then
      plan%style = flags
    else
      plan%style = 0_f32
    endif
    plan%problem_rank = 1
    SLL_ALLOCATE(plan%problem_shape(1),ierr)
    plan%problem_shape = (/ nx  /)
    allocate(plan%twiddles(4*nx + 15))
    call zffti(nx,plan%twiddles)
  end function

  subroutine fftpack_apply_plan_c2c_1d(plan,array_in,array_out)
#ifdef STDF95
    type(sll_fft_plan), pointer                   :: plan
#else
    type(sll_fft_plan), pointer, intent(in)       :: plan
#endif
    sll_comp64, dimension(:), intent(inout)       :: array_in, array_out
    sll_int32  :: nx
    sll_real64 :: factor

    if( loc(array_in) .ne. loc(array_out)) then ! out-place transform
       array_out = array_in
    endif

    nx = plan%problem_shape(1)
 
    if( plan%direction .eq. FFT_FORWARD ) then
      call zfftf( nx , array_out ,plan%twiddles )
    else if( plan%direction .eq. FFT_INVERSE ) then
      call zfftb( nx, array_out , plan%twiddles )
    endif

    if( fft_is_present_flag(plan%style,FFT_NORMALIZE) ) then
      factor = 1.0_f64/real(nx,kind=f64)
      array_out = factor*array_out
    endif
  end subroutine
! END COMPLEX

! REAL
  function fftpack_new_plan_r2r_1d(nx,array_in,array_out,direction,flags) result(plan)
    sll_int32, intent(in)                        :: nx
    sll_real64, dimension(:)                     :: array_in, array_out
    sll_int32, intent(in)                        :: direction
    sll_int32, optional,  intent(in)             :: flags
    type(sll_fft_plan), pointer                  :: plan
    sll_int32                                    :: ierr

    SLL_ALLOCATE(plan,ierr)
    plan%library = FFTPACK_MOD 
    plan%direction = direction
    if( present(flags) )then
      plan%style = flags
    else
      plan%style = 0_f32
    endif
    plan%problem_rank = 1
    SLL_ALLOCATE(plan%problem_shape(1),ierr)
    plan%problem_shape = (/ nx  /)

    SLL_ALLOCATE(plan%twiddles(2*nx + 15),ierr)
    call dffti(nx,plan%twiddles)
  end function

  subroutine fftpack_apply_plan_r2r_1d(plan,array_in,array_out)
#ifdef STDF95
    type(sll_fft_plan), pointer                     :: plan
#else
    type(sll_fft_plan), pointer, intent(in)         :: plan
#endif
    sll_real64, dimension(:), intent(inout)         :: array_in, array_out
    sll_int32  :: nx
    sll_real64 :: factor

    if( loc(array_in) .ne. loc(array_out)) then ! out-place transform
       array_out = array_in
    endif   

    nx = plan%problem_shape(1)

    if( plan%direction .eq. FFT_FORWARD ) then
      call dfftf( nx , array_out ,plan%twiddles )
    else if( plan%direction .eq. FFT_INVERSE ) then
      call dfftb( nx, array_out , plan%twiddles )
    endif

    if( fft_is_present_flag(plan%style,FFT_NORMALIZE) ) then
      factor = 1.0_f64/real(nx,kind=f64)
      array_out = factor*array_out
    endif
  end subroutine
! END REAL

  function fftpack_new_plan_r2c_1d(nx,array_in,array_out,flags) result(plan)
    sll_int32, intent(in)                        :: nx
    sll_real64, dimension(:)                     :: array_in
    sll_comp64, dimension(:)                     :: array_out
    sll_int32, optional,  intent(in)             :: flags
    type(sll_fft_plan), pointer                  :: plan
    sll_int32                                    :: ierr

    print*,"function fftpack_new_plan_r2c_1d is not implemented"
  end function

  subroutine fftpack_apply_plan_r2c_1d(plan,array_in,array_out)
    type(sll_fft_plan), pointer                     :: plan
    sll_real64, dimension(0:), intent(inout)        :: array_in
    sll_comp64, dimension(0:), intent(out)          :: array_out
    sll_int32                                       :: nx, i
    sll_real64 :: factor

    print*,"subroutine fftpack_apply_plan_r2c_1d is not implemented"
  end subroutine

  function fftpack_new_plan_c2r_1d(nx,array_in,array_out,flags) result(plan)
    sll_int32, intent(in)                        :: nx
    sll_real64, dimension(:)                     :: array_out
    sll_comp64, dimension(:)                     :: array_in
    sll_int32, optional,  intent(in)             :: flags
    type(sll_fft_plan), pointer                  :: plan
    sll_int32                                    :: ierr

    print*,"function fftpack_new_plan_c2r_1d is not implemented"
  end function

  subroutine fftpack_apply_plan_c2r_1d(plan,array_in,array_out)
    type(sll_fft_plan), pointer                     :: plan
    sll_real64, dimension(0:), intent(inout)        :: array_out
    sll_comp64, dimension(0:), intent(out)          :: array_in
    sll_int32                                       :: nx, i
    sll_real64 :: factor

    print*,"subroutine fftpack_apply_plan_c2r_1d is not implemented"
  end subroutine

  function fftpack_new_plan_c2c_2d(NX,NY,array_in,array_out,direction,flags) &
    result(plan)

    sll_int32, intent(in)                            :: NX,NY
    sll_comp64, dimension(0:,0:), target, intent(in) :: array_in, array_out
    sll_int32, intent(in)                            :: direction
    sll_int32, optional,  intent(in)                 :: flags
    type(sll_fft_plan), pointer                      :: plan
    sll_int32                                        :: ierr    
    !true if dft in the two directions, false otherwise.    
    logical                                          :: two_direction
   
    print*,"function fftpack_new_plan_c2c_2d is not implemented"
  end function fftpack_new_plan_c2c_2d

  subroutine fftpack_apply_plan_c2c_2d(plan,array_in,array_out)
    type(sll_fft_plan), pointer                     :: plan
    sll_comp64, dimension(0:,0:), intent(inout)     :: array_in, array_out
    sll_int32                                       :: i, nx, ny
    sll_int32, dimension(2)                         :: fft_shape
    logical                                         :: two_direction

    print*,"subroutine fftpack_apply_plan_c2c_2d is not implemented"
  end subroutine fftpack_apply_plan_c2c_2d

  function fftpack_new_plan_c2r_2d(nx,ny,array_in,array_out,flags) result(plan)
    sll_int32, intent(in)                        :: nx,ny
    sll_comp64, dimension(:,:), intent(in)       :: array_in
    sll_real64, dimension(:,:), intent(in)       :: array_out
    sll_int32, optional,  intent(in)             :: flags
    type(sll_fft_plan), pointer                  :: plan
    sll_int32                                    :: ierr

   print*,"function fftpack_new_plan_c2r_2d is no implemented"
  end function fftpack_new_plan_c2r_2d

  subroutine fftpack_apply_plan_c2r_2d(plan,array_in,array_out)
    type(sll_fft_plan), pointer                     :: plan
    sll_comp64, dimension(0:,0:), intent(inout)     :: array_in
    sll_real64, dimension(0:,0:), intent(out)       :: array_out
    sll_int32                                       :: nx, i, ny, k, j
    sll_real64 :: factor

   print*,"subroutine fftpack_apply_plan_c2r_2d is not implemented"
  end subroutine fftpack_apply_plan_c2r_2d

  function fftpack_new_plan_r2c_2d(nx,ny,array_in,array_out,flags) result(plan)
    sll_int32, intent(in)                        :: nx,ny
    sll_real64, dimension(:,:), intent(in)       :: array_in
    sll_comp64, dimension(:,:), intent(in)       :: array_out
    sll_int32, optional,  intent(in)             :: flags
    type(sll_fft_plan), pointer                  :: plan
    sll_int32                                    :: ierr

   print*,"function fftpack_new_plan_r2c_2d is not implemented"
  end function

  subroutine fftpack_apply_plan_r2c_2d(plan,array_in,array_out)
    type(sll_fft_plan), pointer                     :: plan
    sll_real64, dimension(0:,0:), intent(inout)     :: array_in
    sll_comp64, dimension(0:,0:), intent(out)       :: array_out
    sll_int32                                       :: nx, i, ny, k
    sll_real64 :: factor

   print*,"subroutine fftpack_apply_plan_r2c_2d is not implemented"
  end subroutine fftpack_apply_plan_r2c_2d


! FFTPACK FFTPACK FFTPACK FFTPACK FFTPACK FFTPACK FFTPACK FFTPACK FFTPACK
! FFTPACK FFTPACK FFTPACK FFTPACK FFTPACK FFTPACK FFTPACK FFTPACK FFTPACK
! FFTPACK FFTPACK FFTPACK FFTPACK FFTPACK FFTPACK FFTPACK FFTPACK FFTPACK


  subroutine fft_delete_plan(plan)
   type(sll_fft_plan), pointer :: plan

    if( .not. associated(plan) ) then
      print * , '  Error in fft_delete_plan subroutine'
      print * , '  you try to delete a plan not associated'
      stop 
    endif

    if(associated(plan%twiddles)) then
      deallocate(plan%twiddles)
      plan%twiddles => null()
    endif
    plan => null()
  end subroutine

end module sll_fft
