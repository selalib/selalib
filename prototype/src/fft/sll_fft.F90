!------------------------------------------------------------------------------
! SELALIB
!------------------------------------------------------------------------------
! MODULE: sll_fft
!
! DESCRIPTION:
!> @file sll_fft.F90
!> @namespace sll_fft
!> @author EDWIN C. GOLCHER & SAMUEL DE SANTIS
!> @brief Interface around fftpack, fftw and the interne selalib fft.
!> @details 
!>
!> 
!>
!> \section how How to use sll_fft module?
!>
!> First, initialize the plan with sll_new_fft function
!> \code plan => sll_new_fft(n,data_type,flags)\endcode
!> You can call sll_new_fft without the flags argument
!> \code plan => sll_new_fft(n,data_type)\endcode
!> In this case, sll_apply_fft computes an unnormalized FFT.
!> 
!> data_type can take two values : FFT_COMPLEX or FFT_REAL.
!> sll_fft module provides only 64bit fast fourier transform.
!> 
!> flags can take the values : FFT_NORMALIZE_FORWARD, FFT_NORMALIZE_INVERSE
!> You can combine flags with "+", by example
!> \code sll_new_fft(n,data_type,FFT_NORMALIZE_FORWARD + FFT_NORMALIZE_INVERSE) \endcode
!> 
!> Second, apply the fft on the data with sll_apply_fft
!> \code plan => sll_apply_fft(plan,data,direction)\endcode
!> 
!>
!> \warning the output of sll_apply_fft is only in-place way and it is scrambled.
!!          Thus, if you want the mode k (i.e. X_k) you must call sll_get_mode(k).
!!          If you want know which mode is in position i in array data call
!!          sll_get_index(i)
!>
!> \section example Examples:
!>
!> For complexe data :
!> \code
!> sll_comp64, dimension(0,n-1) :: data    ! n is the size of the problem
!> type(sll_fft_plan), pointer  :: plan
!> sll_int32                    :: flags
!> sll_int32                    :: direction
!>
!> ** INIT DATA **
!>
!> plan => sll_new_fft(n,FFT_COMPLEX,flags)
!> plan => sll_apply_fft(plan,data,direction)
!>
!> plan => sll_delete_fft(plan)
!> \endcode
!>
!> By example if, the input data is \f$(x_0,x_1,x_2,x_3)\f$ the output is \f$(X_0,X_2,X_1,X_3)\f$.
!> Thus, sll_get_index(1) returns 2 (cause data[1]=X_2) and sll_get_mode(1) returns X_1.
!> 
!> 
!> 
!> For real data : 
!> \code
!> sll_real64, dimension(0,n-1) :: data    ! n is the size of the problem
!> type(sll_fft_plan), pointer  :: plan
!> sll_int32                    :: flags
!> sll_int32                    :: direction
!>
!> ** INIT DATA **
!>
!> plan => sll_new_fft(n,FFT_REAL,flags)
!> plan => sll_apply_fft(plan,data,direction)
!>
!> plan => sll_delete_fft(plan)
!> \endcode
!>
! \warning let p = sll_get_index(i), if p is even data(i) is the real part of X_p, else if p is odd data(i) is the imaginary part of X_p
!>
!>
!>
!> \section what What sll_fft really computes
!>
!> The forward (FFT_FORWARD) DFT of a 1d complex array x of size n computes an array X, where:
!>
!> \f[ X_k = \sum_{i=0}^{n-1} x_i e^{-2\pi i j k/n}. \f]
!>
!> The backward (FFT_INVERSE) DFT computes:
!>
!> \f[ x_i = \sum_{k=0}^{n-1} X_k e^{2\pi k j i/n}. \f]
!>
!> For the real transform, we have
!> \f$ (x_0,x_1,\dots,x_{n-1}) \rightarrow
!!     (r_0,r_{n/2},r_1,i_1,\dots,r_{n/2-1},i_{n/2-1})\f$
!> which must be interpreted as the complex array
!> \f[ \begin{pmatrix} r_0 &,& 0
!!                     \\ r_1 &,& i_1
!!                     \\ \vdots  & & \vdots 
!!                     \\ r_{n/2-1} &,& i_{n/2-1}
!!                     \\ r_{n/2} &,& 0
!!                     \\ r_{n/2-1} &,& -i_{n/2-1}
!!                     \\ \vdots    & & \vdots
!!                     \\ r_1 &,& -i_1 
!! \end{pmatrix}\f] 
!> \warning Note that ffw use \f$(r_0,r_1,\dots,r_{n/2-1},r_{n/2},i_{n/2-1},\dots,i_1)\f$
!!          convention whereas fftpack use \f$(r_0,r_1,i_1,\dots,r_{n/2-1},i_{n/2-1},r_{n/2})\f$
!> 
!>
!------------------------------------------------------------------------------

module sll_fft
  use numeric_constants
#include "conf.h"
#ifndef _NOFFTW
  !use FFTW3
  use, intrinsic :: iso_c_binding
  !include 'fftw3.f03'
#endif
#ifndef _NOFFTPACK
  !use fftpack
#endif
#include "sll_working_precision.h"
#include "misc_utils.h"  
#include "sll_assert.h"
#include "sll_memory.h"
  implicit none
#ifndef _NOFFTW
  !use FFTW3
  !use, intrinsic :: iso_c_binding
  include 'fftw3.f03'
#endif

  ! Basic 1D FFT descriptor:
  ! - N: number of samples
  ! - style : transform style
  ! - index : given a wavenumber on 1:N, return the corresponding array 
  !           index (also on 1:N)
  ! - mode  : given an array index, return the corresponding wavenumber
  type sll_fft_plan
    sll_int32 :: N
    sll_int32 :: style
    sll_int32 :: data_type
    sll_int32,  dimension(:), pointer :: index
    sll_int32,  dimension(:), pointer :: mode
    sll_comp64, dimension(:), pointer :: t          ! twiddle factors complex case
    sll_real64, dimension(:), pointer :: twiddles   ! twiddles factors real case 
    sll_real64, dimension(:), pointer :: twiddles_n ! twiddles factors real case 
    sll_real32, dimension(:), pointer :: wsave ! for use fftpack
    sll_real64, dimension(:), pointer :: dwsave ! for use fftpack
    sll_int32                         :: mod !type of library used
    sll_real64                        :: normalization_factor_forward = 1.0_f64
    sll_real64                        :: normalization_factor_backward = 1.0_f64
  end type sll_fft_plan

  ! We choose the convention in which the direction of the FFT is determined
  ! by the conjugation of the twiddle factors. If 
  !
  !                  omega_N = -j2*pi/N
  ! 
  ! then we call this the FORWARD transform. The following enumeration will
  ! allow us to always use a descriptive qualifier instead of using an
  ! integer (or something else). Since we are using the 'right' sign
  ! convention in the passed parameter, in principle we could use this
  ! number directly and avoid some branching (i.e.: if-statements). This
  ! would be more efficient but would tie the functions to this particular
  ! enumeration, reducing flexibility. However, this is a flexibility we
  ! may choose never to exercise.

  enum, bind(C)
     enumerator :: FFT_FORWARD = -1, FFT_INVERSE = 1 
  end enum

  enum, bind(C)
     enumerator :: FFT_REAL = 0, FFT_COMPLEX = 1
  end enum

  enum, bind(C)
     enumerator :: FFT_NORMALIZE_FORWARD = 1, FFT_NORMALIZE_INVERSE = 2,&
                   FFT_NEGATIVE_PHASE = 4
  end enum

  enum, bind(C)
     enumerator :: FFTPACK_MOD = 100, FFTW_MOD = 1000000000, SLLFFT_MOD = 0
  end enum

  interface bit_reverse
    module procedure bit_reverse_complex, bit_reverse_integer32, &
                     bit_reverse_integer64
  end interface

  interface sll_apply_fft
    module procedure sll_apply_fft_complex32,  sll_apply_fft_real64, sll_apply_fft_complex64
  end interface

  interface sll_get_mode
    module procedure sll_get_mode_real, sll_get_mode_complex, sll_get_mode_complex32
  end interface
contains

  ! data_type = 0 for real
  !           = 1 for complexe
  ! example of call sll_new_fft(N, FFT_NORMALIZE_FORWARD , FFT_COMPLEX)
  function sll_new_fft(n,data_type,style)
    type(sll_fft_plan), pointer     :: sll_new_fft
    sll_int32, intent(in)           :: n, data_type
    sll_int32, optional, intent(in) :: style
    sll_int32                       :: ierr, s, mod, n_2, i
   
    ! By default style = 0
    if(present(style)) then
      if( (style>=FFTPACK_MOD) .and. (style<FFTW_MOD) ) then
        s = style - FFTPACK_MOD
        mod = FFTPACK_MOD
      else if(style>=FFTW_MOD) then
        s = style - FFTW_MOD
        mod = FFTW_MOD
      else
        s = style
        mod = SLLFFT_MOD
      endif
    else
      s = 0
      mod = SLLFFT_MOD
    endif

    if( mod .eq. FFTPACK_MOD ) then
#ifndef _NOFFTPACK
      sll_new_fft => fftpack_new_fft(n,data_type,s)
      return
#else
  stop 'FFTPACK NOT INSTALLED ==sll_new_fft=='
#endif
    endif

    if( mod .eq. FFTW_MOD ) then
#ifndef _NOFFTW
      sll_new_fft => fftw_new_fft(n,data_type,s)
      return
#else
  stop 'FFTW NOT INSTALLED'
#endif
    endif

    SLL_ALLOCATE(sll_new_fft,ierr)
    sll_new_fft%N = n
    sll_new_fft%style = s
    sll_new_fft%data_type = data_type
    sll_new_fft%mod = mod
    !SLL_ALLOCATE(sll_new_fft%index(0:n-1),ierr)

    !sll_new_fft%index => compute_index(n)

    n_2 = n/2
    if(data_type .eq. FFT_COMPLEX) then

      if( s .eq. 1 ) then
        sll_new_fft%normalization_factor_forward = 1.0_f64/(real(n,kind=f64))
        sll_new_fft%normalization_factor_backward = 1.0_f64
      else if( s .eq. 2 ) then
        sll_new_fft%normalization_factor_backward = 1.0_f64/(real(n,kind=f64))
        sll_new_fft%normalization_factor_forward = 1.0_f64
      else if( s .eq. 3 ) then
        sll_new_fft%normalization_factor_forward = 1.0_f64/(real(n,kind=f64))
        sll_new_fft%normalization_factor_backward = 1.0_f64/(real(n,kind=f64))
      endif

      SLL_ALLOCATE(sll_new_fft%mode(0:n-1),ierr)
      do i=0,n-1
        sll_new_fft%mode(i) = compute_mode(i,n)
      enddo

      SLL_ALLOCATE(sll_new_fft%t(1:n_2),ierr)
      call compute_twiddles(n,sll_new_fft%t)
      call bit_reverse(n_2,sll_new_fft%t)
    else if(data_type .eq. FFT_REAL) then

      !do i=0,n-1
      !  sll_new_fft%mode(i) = sll_get_mode_real(i,n)
      !enddo

      if( s .eq. 1 ) then
        sll_new_fft%normalization_factor_forward = 2.0_f64/(real(n,kind=f64))
        sll_new_fft%normalization_factor_backward = 1.0_f64
      else if( s .eq. 2 ) then
        sll_new_fft%normalization_factor_backward = 2.0_f64/(real(n,kind=f64))
        sll_new_fft%normalization_factor_forward = 1.0_f64
      else if( s .eq. 3 ) then
        sll_new_fft%normalization_factor_forward = 2.0_f64/(real(n,kind=f64))
        sll_new_fft%normalization_factor_backward = 2.0_f64/(real(n,kind=f64))
      endif

      SLL_ALLOCATE(sll_new_fft%twiddles(0:n_2-1),ierr)
      SLL_ALLOCATE(sll_new_fft%twiddles_n(0:n-1),ierr)
      call compute_twiddles_real_array( n, sll_new_fft%twiddles_n )
      call compute_twiddles_real_array( n_2, sll_new_fft%twiddles(0:n_2-1) )
      call bit_reverse_in_pairs( n/4, sll_new_fft%twiddles(0:n_2-1))
      !call bit_reverse_in_pairs( n_2/2, sll_new_fft%twiddles )
      !call bit_reverse_in_pairs( n_2, sll_new_fft%twiddles_n )
    else
      stop 'Error in ==sll_new_fft== unknown data_type'
    endif
  end function sll_new_fft

  function sll_get_mode_complex32(k,plan,data) result(mode)
    sll_comp32                                :: mode
    sll_int32, intent(in)                     :: k
    type(sll_fft_plan)                        :: plan
    sll_comp32, dimension(0:) , intent(in)    :: data
    sll_int32                                 :: n_2, n
  
    n = plan%N
    n_2 = n/2

    SLL_ASSERT( (k .ge. 0) .and. (k .lt. n) )
    SLL_ASSERT( size(data) .eq. n )
   
    if(plan%mod .eq. FFTW_MOD) then
      mode = data(k)
    else if(plan%mod .eq. FFTPACK_MOD) then
      mode = data(k)
    else if(plan%mod .eq. SLLFFT_MOD) then
      mode = data(plan%mode(k))
    else
      stop 'ERROR IN =SLL_GET_MODE_COMPLEX32= plan%mod invalid'
    endif
  end function sll_get_mode_complex32
  
  function sll_get_mode_complex(k,plan,data) result(mode)
    sll_comp64                                :: mode
    sll_int32, intent(in)                     :: k
    type(sll_fft_plan)                        :: plan
    sll_comp64, dimension(0:) , intent(in)    :: data
    sll_int32                                 :: n_2, n
  
    n = plan%N
    n_2 = n/2

    SLL_ASSERT( (k .ge. 0) .and. (k .lt. n) )
    SLL_ASSERT( size(data) .eq. n )
   
    if(plan%mod .eq. FFTW_MOD) then
      mode = data(k)
    else if(plan%mod .eq. FFTPACK_MOD) then
      mode = data(k)
    else if(plan%mod .eq. SLLFFT_MOD) then
      mode = data(plan%mode(k))
    else
      stop 'ERROR IN =SLL_GET_MODE_COMPLEX= plan%mod invalid'
    endif
  end function sll_get_mode_complex

  function sll_get_mode_real(k,plan,data) result(mode)
    sll_comp64                                :: mode
    sll_int32, intent(in)                     :: k
    type(sll_fft_plan) , intent(in)           :: plan
    sll_real64, dimension(0:) , intent(in) :: data
    sll_int32                                 :: n_2, n
   
    n = plan%N
    n_2 = n/2

    SLL_ASSERT( (k .ge. 0) .and. (k .lt. n) )
    SLL_ASSERT( size(data) .eq. n )


    if(plan%mod .eq. FFTPACK_MOD) then
      if( k .eq. 0 ) then
        mode = complex(data(0),0.0_f64)
      else if( k .eq. n_2 ) then
        mode = complex(data(n-1),0.0_f64)
      else if( k .gt. n_2 ) then
        mode = complex( data(2*(k-n_2)-1) , -data(2*(k-n_2)) )
      else
        mode = complex( data(2*k-1) , data(2*k) )
      endif
    else if(plan%mod .eq. SLLFFT_MOD) then
      if( k .eq. 0 ) then
        mode = complex(data(0),0.0_f64)
      else if( k .eq. n_2 ) then
        mode = complex(data(1),0.0_f64)
      else if( k .gt. n_2 ) then
        mode = complex( data(2*(k-n_2)) , -data(2*(k-n_2)+1) )
      else
        mode = complex( data(2*k) , data(2*k+1) )
      endif
    else if(plan%mod .eq. FFTW_MOD) then
      if( k .eq. 0 ) then
        mode = complex(data(0),0.0_f64)
      else if( k .eq. n_2 ) then
        mode = complex(data(n_2),0.0_f64)
      else if( k .gt. n_2 ) then
        mode = complex( data(k-n_2) , -data(n-k+n_2) )
      else
        mode = complex( data(k) , data(n-k) )
      endif
    else
      stop 'ERROR IN =SLL_GET_MODE_REAL= plan%mod invalid'
    endif
  end function sll_get_mode_real


! FFTW FFTW FFTW FFTW FFTW FFTW FFTW FFTW FFTW FFTW FFTW FFTW FFTW FFTW
#ifndef _NOFFTW
  function fftw_new_fft(n,data_type,style)
    type(sll_fft_plan), pointer     :: fftw_new_fft
    sll_int32, intent(in)           :: n, data_type
    sll_int32, optional, intent(in) :: style
    sll_int32                       :: ierr, s

    !
    s = 0 !FFTW_ESTIMATE
    if(present(style)) s = style

    SLL_ALLOCATE(fftw_new_fft,ierr)

    fftw_new_fft%N = n
    fftw_new_fft%style = s
    fftw_new_fft%data_type = data_type
    fftw_new_fft%mod = FFTW_MOD
    if( s .eq. 1 ) then
      fftw_new_fft%normalization_factor_forward = 1.0_f64/(real(n,kind=f64))
      fftw_new_fft%normalization_factor_backward = 1.0_f64
    else if( s .eq. 2 ) then
      fftw_new_fft%normalization_factor_backward = 1.0_f64/(real(n,kind=f64))
      fftw_new_fft%normalization_factor_forward = 1.0_f64
    else if( s .eq. 3 ) then
      fftw_new_fft%normalization_factor_forward = 1.0_f64/(real(n,kind=f64))
      fftw_new_fft%normalization_factor_backward = 1.0_f64/(real(n,kind=f64))
    endif
  end function fftw_new_fft

  function fftw_apply_fft_complex64(fft_plan,in,direction)
    type(sll_fft_plan), pointer                     :: fftw_apply_fft_complex64
    type(sll_fft_plan), pointer, intent(in)         :: fft_plan
    integer, intent(in)                             :: direction
    sll_comp64, dimension(0:fft_plan%N-1), intent(inout) :: in
    type(C_PTR) :: plan

    plan = fftw_plan_dft_1d(fft_plan%N,in,in,direction,FFTW_ESTIMATE)
    call fftw_execute_dft(plan, in, in)
    call fftw_destroy_plan(plan)
    fftw_apply_fft_complex64 => fft_plan
  end function fftw_apply_fft_complex64

  function fftw_apply_fft_real64(fft_plan,in,direction)
    type(sll_fft_plan), pointer                     :: fftw_apply_fft_real64
    type(sll_fft_plan), pointer, intent(in)         :: fft_plan
    integer, intent(in)                             :: direction
    sll_real64, dimension(0:fft_plan%N-1), intent(inout) :: in
    type(C_PTR) :: plan

    if(direction .eq. FFT_INVERSE) then
      plan = fftw_plan_r2r_1d(fft_plan%n,in,in,FFTW_HC2R,FFTW_ESTIMATE)
      call fftw_execute_r2r(plan,in, in)
      call fftw_destroy_plan(plan)
      in = fft_plan%normalization_factor_backward*in
      fftw_apply_fft_real64 => fft_plan
    else if( direction .eq. FFT_FORWARD ) then
      plan = fftw_plan_r2r_1d(fft_plan%n,in,in,FFTW_R2HC,FFTW_ESTIMATE)
      call fftw_execute_r2r(plan,in, in)
      call fftw_destroy_plan(plan)
      in = fft_plan%normalization_factor_forward*in
      fftw_apply_fft_real64 => fft_plan
    else
      stop 'ERROR IN =FFTW_APPLY_FFT_REAL64= argument direction invalid'
    endif
  end function fftw_apply_fft_real64
#endif
! FFTW FFTW FFTW FFTW FFTW FFTW FFTW FFTW FFTW FFTW FFTW FFTW FFTW FFTW

! FFTPACK FFTPACK FFTPACK FFTPACK FFTPACK FFTPACK FFTPACK FFTPACK FFTPACK
#ifndef _NOFFTPACK
  function fftpack_new_fft(n,data_type,style)
    type(sll_fft_plan), pointer :: fftpack_new_fft
    sll_int32, intent(in)       :: n, data_type
    sll_int32, optional, intent(in) :: style
    sll_int32                   :: ierr, s

    ! By default style = 0
    s = 0
    if(present(style)) s = style
       
    SLL_ALLOCATE(fftpack_new_fft,ierr)
  
    fftpack_new_fft%N = n
    fftpack_new_fft%style = s
    fftpack_new_fft%mod = FFTPACK_MOD
    fftpack_new_fft%data_type = data_type

    if( s .eq. 1 ) then
      fftpack_new_fft%normalization_factor_forward = 1.0_f64/(real(n,kind=f64))
      fftpack_new_fft%normalization_factor_backward = 1.0_f64
    else if( s .eq. 2 ) then
      fftpack_new_fft%normalization_factor_backward = 1.0_f64/(real(n,kind=f64))
      fftpack_new_fft%normalization_factor_forward = 1.0_f64
    else if( s .eq. 3 ) then
      fftpack_new_fft%normalization_factor_forward = 1.0_f64/(real(n,kind=f64))
      fftpack_new_fft%normalization_factor_backward = 1.0_f64/(real(n,kind=f64))
    endif

    if(data_type .eq. FFT_COMPLEX) then
      SLL_ALLOCATE(fftpack_new_fft%wsave(4*n + 15),ierr)
      call cffti(n, fftpack_new_fft%wsave)
    else
      SLL_ALLOCATE(fftpack_new_fft%dwsave(2*n + 15),ierr)
      call dffti(n, fftpack_new_fft%dwsave)
    endif
  end function fftpack_new_fft

!  function fftpack_apply_fft_real(fft,data,direction)
!    type(sll_fft_plan), pointer             :: fftpack_apply_fft_real
!    type(sll_fft_plan), pointer, intent(in) :: fft
!    sll_int32, intent(in)                   :: direction
!    sll_real32, dimension(1:fft%N), intent(inout) :: data
!    
!    if( direction .eq. FFT_FORWARD) then  
!      call rfftf(fft%N,data,fft%wsave)
!    else
!      call rfftb(fft%N,data,fft%wsave)
!    endif
!    fftpack_apply_fft_real => fft
!  end function fftpack_apply_fft_real
  
  function fftpack_apply_fft_real64(fft,data,direction)
    type(sll_fft_plan), pointer             :: fftpack_apply_fft_real64
    type(sll_fft_plan), pointer, intent(in) :: fft
    sll_int32, intent(in)                   :: direction
    sll_real64, dimension(1:fft%N), intent(inout) :: data
    
    if( direction .eq. FFT_FORWARD) then
      data = data*fft%normalization_factor_forward
      call dfftf(fft%N,data,fft%dwsave)
    else if( direction .eq. FFT_INVERSE) then
      data = data*fft%normalization_factor_backward
      call dfftb ( fft%N, data, fft%dwsave )
    endif
    fftpack_apply_fft_real64 => fft
  end function fftpack_apply_fft_real64

  function fftpack_apply_fft_complex32(fft,data,direction)
    type(sll_fft_plan), pointer             :: fftpack_apply_fft_complex32
    type(sll_fft_plan), pointer, intent(in) :: fft
    sll_int32, intent(in)                   :: direction
    sll_comp32, dimension(:), intent(inout) :: data
    
    if( direction .eq. FFT_FORWARD ) then  
      call cfftf(fft%N,data,fft%wsave)
    else if( direction .eq. FFT_INVERSE ) then
      call cfftb(fft%N,data,fft%wsave)
    else
      stop 'ERROR : direction ==fftpack_apply_fft_complex32=='
    endif
    fftpack_apply_fft_complex32 => fft
  end function fftpack_apply_fft_complex32

#else
#endif
! FFTPACK FFTPACK FFTPACK FFTPACK FFTPACK FFTPACK FFTPACK FFTPACK FFTPACK

  recursive function compute_mode(k,n) result(i)
    sll_int32, intent(in)  :: k, n
    sll_int32              :: i 
    
    if( k .eq. 0 ) then
      i = 0 
    else if( mod(k,2) .eq. 1 ) then
      i = n/2 + compute_mode(k-1,n)
    else
      i = compute_mode(k/2,n/2)
    endif
  end function

  function compute_index(n) result(ind)
    sll_int32, intent(in)                    :: n
    sll_int32, dimension(0:n-1)              :: ind
    sll_int32, dimension(0:n-1)              :: a
    sll_int32                                :: i
  
    !SLL_ALLOCATE(a(0:n-1),ierr)
    !SLL_ALLOCATE(index(0:n-1),ierr)
    do i=0,n-1
      a(i) = i
    enddo
    call bit_reverse(n,a)

    do i=0,n-1 
      ind(a(i)) = i
    enddo 

    !SLL_DEALLOCATE_ARRAY(a,ierr)
  end function compute_index

  function sll_delete_fft(fft)
      type(sll_fft_plan), pointer                :: sll_delete_fft
      type(sll_fft_plan), pointer, intent(inout) :: fft
      sll_int32                                  :: ierr

      if( fft%mod .eq. FFTW_MOD ) then
        !print *, 'nothing in delete fftw'
        !return
      endif

      if( fft%mod .eq. FFTPACK_MOD ) then
       if( fft%data_type .eq. FFT_COMPLEX ) then
        if( associated(fft%wsave) ) &
          SLL_DEALLOCATE(fft%wsave,ierr)
       else
        if( associated(fft%dwsave) ) &
          SLL_DEALLOCATE(fft%dwsave,ierr)
       endif
      endif

      if( fft%mod .eq. SLLFFT_MOD ) then
       if( fft%data_type .eq. FFT_COMPLEX ) then
        if( associated(fft%mode) ) &
          SLL_DEALLOCATE(fft%mode,ierr)
        if( associated(fft%t) ) &
          SLL_DEALLOCATE(fft%t,ierr)
       else
        if( associated(fft%twiddles) ) &
          SLL_DEALLOCATE(fft%twiddles,ierr)
        if( associated(fft%twiddles_n) ) &
          SLL_DEALLOCATE(fft%twiddles_n,ierr)
       endif
        !if( associated(fft%index) ) &
        !  SLL_DEALLOCATE(fft%index,ierr)      
      endif
      

      SLL_DEALLOCATE(fft,ierr)
      fft => NULL()
      sll_delete_fft => NULL()
  end function sll_delete_fft

  function sll_apply_fft_complex64(fft,data,direction)
    type(sll_fft_plan), pointer             :: sll_apply_fft_complex64
    type(sll_fft_plan), pointer, intent(in) :: fft
    sll_int32, intent(in)                   :: direction
    sll_comp64, dimension(0:fft%N-1), intent(inout) :: data
  
    if( fft%mod .eq. FFTW_MOD ) then
#ifndef _NOFFTW
      sll_apply_fft_complex64 => fftw_apply_fft_complex64(fft,data,direction)
      return
#else
       stop 'FFTW NOT INSTALLED'
#endif
    else if( fft%mod .eq. FFTPACK_MOD ) then
      stop 'no 64 bit version in FFTPACK'
    endif

    call fft_dit_nr(data,fft%t,direction)
    sll_apply_fft_complex64 => fft
    !call bit_reverse_complex(fft%N,data)
  end function sll_apply_fft_complex64

  function sll_apply_fft_complex32(fft,data,direction)
    type(sll_fft_plan), pointer             :: sll_apply_fft_complex32
    type(sll_fft_plan), pointer, intent(in) :: fft
    sll_int32, intent(in)                   :: direction
    sll_comp32, dimension(0:fft%N-1), intent(inout) :: data
   
    if( fft%mod .eq. FFTPACK_MOD ) then
#ifndef _NOFFTPACK
      sll_apply_fft_complex32 => fftpack_apply_fft_complex32(fft,data,direction)
      return
#else
      stop 'FFTPACK NOT INSTALLED'
#endif
    else if( fft%mod .eq. FFTW_MOD ) then
#ifndef _NOFFTW
      stop 'no 32bit version in FFTW'
#else
      stop 'FFTW NOT INSTALLED'
#endif
    else if( fft%mod .eq. SLLFFT_MOD ) then
      stop 'no 32bit version in SELALIB'
    else
      stop 'ERROR IN =SLL_APPLY_FFT_COMPLEX32= invalid fft%mod'
    endif
  end function sll_apply_fft_complex32



  function sll_apply_fft_real64(fft,data,direction)
    type(sll_fft_plan), pointer             :: sll_apply_fft_real64
    type(sll_fft_plan), pointer, intent(in) :: fft
    sll_int32, intent(in)                   :: direction
    sll_real64, dimension(0:fft%N-1), intent(inout) :: data
    
    if( fft%mod .eq. FFTW_MOD) then
#ifndef _NOFFTW
      sll_apply_fft_real64 => fftw_apply_fft_real64(fft,data,direction)
      return
#else
      stop 'FFTW NOT INSTALLED'
#endif
    else if( fft%mod .eq. FFTPACK_MOD) then
#ifndef _NOFFTPACK
      sll_apply_fft_real64 => fftpack_apply_fft_real64(fft,data,direction)
      return
#else
      stop 'FFTPACK NOT INSTALLED'
#endif
    endif

    if( direction .eq. FFT_INVERSE ) then
      data = data*fft%normalization_factor_backward
    else if( direction .eq. FFT_FORWARD ) then 
      data = data*fft%normalization_factor_forward
    endif

    call real_data_fft_dit( data, fft%N, fft%twiddles, fft%twiddles_n, direction )
    sll_apply_fft_real64 => fft
  end function sll_apply_fft_real64


!  function new_fft_1D_plan( n_points, data_type, direction )
!    type(sll_fft_plan), pointer :: new_fft_1D_plan
!    ! n_points - 1 must be a power of two in this implementation.
!    sll_int32                :: n_points
! !   sll_int32                :: style
!    sll_int32, intent(in)    :: direction
!    sll_int32, intent(in)    :: data_type
!    sll_int32                :: ierr
!    SLL_ALLOCATE( new_fft_1D_plan, ierr )
!    new_fft_1D_plan%N         = n_points
!  !  new_fft_1D_plan%style     = style
!    new_fft_1D_plan%data_type = data_type
!    SLL_ALLOCATE( new_fft_1D_plan%index(n_points), ierr )
!    SLL_ALLOCATE( new_fft_1D_plan%mode(n_points),  ierr )
!    ! compute twiddles also allocates the memory...
!    call compute_twiddles( n_points-1, new_fft_1D_plan%t )
!  end function new_fft_1D_plan


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
       theta   = 2.0_f64*sll_pi/real(n,kind=f64)      ! basic angular interval
       ! By whatever means we use to compute the twiddles, some sanity
       ! checks are in order: 
       ! t(1)     = (1,0)
       ! t(n/8+1) = (sqrt(2)/2, sqrt(2)/2)
       ! t(n/4+1) = (0,1)
       ! t(n/2+1) = (0,-1) ... but this one is not stored
       do k = 0,n/2-1
          t(k+1) = exp((0.0_f64,1.0_f64)*theta*real(k,kind=f64))
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
    SLL_ASSERT(is_power_of_two(int(n,i64))) 
    theta   = 2.0_f64*sll_pi/real(n,kind=f64)      ! basic angular interval
    ! By whatever means we use to compute the twiddles, some sanity
    ! checks are in order: 
    ! t(0)   = 1; t(1)      = 0
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
    SLL_ASSERT(is_power_of_two(int(n,i64)));                    \
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
  MAKE_BIT_REVERSE_FUNCTION( bit_reverse_integer32, sll_int32 )
  MAKE_BIT_REVERSE_FUNCTION( bit_reverse_integer64, sll_int64 )
  MAKE_BIT_REVERSE_FUNCTION( bit_reverse_real, sll_real64 )

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
    SLL_ASSERT(is_power_of_two(int(num_pairs,i64)))
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
  ! *************************************************************************
  ! Decimation in time FFT, natural order input, bit-reversed output (=NR). 
  ! Size of the data must be a power of 2. Twiddle array must be in 
  ! bit-reversed order. This implementation is 'cache-oblivious'. This is 
  ! only a placeholder for an FFT really; it is not very efficient since it:
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
  ! the omega_2^1 twiddle; the size-4 problems use omega_4^1 and omega_4^2;
  ! the size-8 problems use omega_8^1, omega_8^2, omega_8^3 and omega_8^4. This
  ! is the expected progression of the twiddle indices as we move deeper into
  ! the recursions.
  !    
  ! *************************************************************************
  
  subroutine fft_dit_nr(data, twiddles, sign)
    sll_comp64, dimension(:), intent(inout) :: data
    sll_int32, intent(in)                   :: sign
    sll_comp64, dimension(:), intent(in)    :: twiddles
    sll_int32                               :: n
    !sll_int32                               :: ierr

    n = size(data)
    SLL_ASSERT(is_power_of_two(int(n,i64)))
    !SLL_ALLOCATE(twiddles(n/2),ierr)
    !call compute_twiddles(n,twiddles) 
    !call bit_reverse(n/2,twiddles)
    call fft_dit_nr_aux(data, n, twiddles, 0, sign)
  end subroutine fft_dit_nr
  
  ! Decimation-in-time, natural-order input, bit-reversed order output:
  recursive subroutine fft_dit_nr_aux( dat, size, twiddles, &
                                       twiddle_index, sign )
    intrinsic ishft, conjg
    sll_comp64, dimension(0:size-1), intent(inout)  :: dat
    integer, intent(in)                             :: size
    ! It is more convenient when the twiddles are 0-indexed
    sll_comp64, dimension(0:), intent(in)           :: twiddles 
    integer, intent(in)                             :: twiddle_index
    integer, intent(in)                             :: sign
    integer                                         :: half
    integer                                         :: j
    integer                                         :: t_index_new
    sll_comp64                                      :: omega
    sll_comp64                                      :: tmp

    half = ishft(size,-1) ! any evidence this is faster?
    if ( sign == FFT_INVERSE ) then
       omega = twiddles(twiddle_index)
    else if ( sign == FFT_FORWARD ) then
       omega = conjg(twiddles(twiddle_index))
    end if
    ! Do the butterflies for this stage
    do j=0,half-1
       tmp          = dat(j+half)*omega
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

  subroutine fft_dit_rn( data, sign )
    sll_comp64, dimension(:), intent(inout) :: data
    sll_int32, intent(in)                   :: sign
    sll_comp64, dimension(:), pointer       :: twiddles
    integer                                 :: n
    sll_int32                               :: ierr
    n = size(data) ! bad
    SLL_ASSERT(is_power_of_two(int(n,i64)))
    SLL_ALLOCATE(twiddles(n/2),ierr)
    call compute_twiddles(n,twiddles) 
    ! This algorithm uses the twiddles in natural order. The '1' 
    ! argument is because fft_dit_rn_aux internally 1-indexes its
    ! arrays, so we are just indicating the first twiddle factor.
    call fft_dit_rn_aux(data, n, twiddles, 1, sign)
  end subroutine fft_dit_rn
  
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
       if ( sign == FFT_FORWARD ) then
          omega = conjg(twiddles(jtwiddle))
       else if ( sign == FFT_INVERSE ) then
          omega = twiddles(jtwiddle)
       end if
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
    sll_real64, dimension(0:num_complex-1), intent(in)     :: twiddles 
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
    if ( sign == FFT_FORWARD ) then
       omega_re =  CREAL0(twiddles, twiddle_index)
       omega_im = -CIMAG0(twiddles, twiddle_index)
    else if ( sign == FFT_INVERSE ) then
       omega_re =  CREAL0(twiddles, twiddle_index)
       omega_im =  CIMAG0(twiddles, twiddle_index)
    end if
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
       if( sign == FFT_FORWARD ) then
          omega_re =  CREAL1(twiddles, jtwiddle)
          omega_im = -CIMAG1(twiddles, jtwiddle)
       else if ( sign == FFT_INVERSE ) then
          omega_re =  CREAL1(twiddles, jtwiddle)
          omega_im =  CIMAG1(twiddles, jtwiddle)
       end if
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
  !   special logic will be required to deliver a requested Fourier mode; 
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
    SLL_ASSERT( is_power_of_two(int(n,i64)) )
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
    if( sign .eq. FFT_FORWARD ) then
       ! we use the following as the 'switch' to flip signs between
       ! FORWARD and INVERSE transforms.
       s = -1.0_f64 
       ! The following call has to change to refer to its natural wrapper...
       call fft_dit_nr_real_array_aux( data(0:n-1), &
                                       n_2,         &
                                       twiddles,    &
                                       0,           &
                                       FFT_FORWARD )
       ! but our _nr_ algorithm bit reverses the result, so, until we have
       ! some way to index the data correctly we have to do this:
       call bit_reverse_in_pairs( n_2, data(0:n-1) )
    else if (sign .eq. FFT_INVERSE) then
       s =  1.0_f64
    else
      stop 'ERROR IN =REAL_DATA_FFT_DIT= invalid argument sign'
    end if
    do i=1,n_2/2 ! the index 'i' corresponds to indexing H
       ! FFT_FORWARD case: We intend to mix the odd/even components that we 
       ! have computed into complex numbers H_n. These Complex numbers will 
       ! give the modes as in:
       !
       ! F_i = 1/2*(H_i + H_(N/2-i)^*) - i/2*(H_i-H_(N/2-i)^*)*exp(-j*2*pi*i/N)
       !
       ! which is the answer we are after.
       !
       ! FFT_INVERSE case: The process of decomposing the H_n's into the
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
       ! Compute tmp2 = i/2*(H_n - H_(N/2-n)^*); the sign depends on the
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
    if ( sign .eq. FFT_FORWARD ) then
       ! Set the first and N/2 values independently and pack them in the
       ! memory space provided by the original data.
       tmp_re = data(0)
       data(0) = tmp_re + data(1)   ! mode 0   is real
       data(1) = tmp_re - data(1)   ! mode N/2 is real
    else if ( sign .eq. FFT_INVERSE ) then
       ! Unpack the modes.
       tmp_re  = data(0)
       data(0) = 0.5_f64*(tmp_re + data(1))
       data(1) = 0.5_f64*(tmp_re - data(1))
       ! The following call has to change to refere to its natural wrapper...
       call fft_dit_nr_real_array_aux( data(0:n-1), &
                                       n_2,         &
                                       twiddles,    &
                                       0,           &
                                       FFT_INVERSE )
       ! but our _nr_ algorithm bit reverses the result, so, until we have
       ! some way to index the data correctly we have to do this:
       call bit_reverse_in_pairs( n_2, data(0:n-1) )
    end if
  end subroutine real_data_fft_dit
  
#undef CREAL0
#undef CIMAG0
#undef CREAL1
#undef CIMAG1


end module sll_fft
