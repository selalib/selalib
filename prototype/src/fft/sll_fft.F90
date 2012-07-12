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
!> 1. Declare a fft plan
!> \code type(sll_fft_plan), pointer :: p \endcode
!> 2. Initialize the plan
!> \code p => fft_new_plan(size,in,out,direction,flags) \endcode
!> The arrays in and out can be real and/or complex, 1d or 2d. The size is only a power of two (radix-2).
!> \warning For complex to real and real to complex transform, there is no direction flag.
!>          \code p => fft_new_plan(size,in,out,flags) \endcode
!>
!> \a direction can take two values : FFT_FORWARD and FFT_INVERSE
!>
!> \a flags optional argument that can be : FFT_NORMALIZE
!>                                FFT_ONLY_FIRST_DIRECTION, FFT_ONLY_SECOND_DIRECTION (2d case only)
!> You can combine flags with '+'.
!>
!> 3. Execute the plan
!> \code call fft_apply_plan(p,in,out) \endcode
!> 4. Delete the plan
!> \code call fft_delete_plan(p) \endcode
!>
! \warning the output of sll_apply_fft is only in-place way and it is scrambled.
!          Thus, if you want the mode k (i.e. X_k) you must call sll_get_mode(k).
!          If you want know which mode is in position i in array data call
!          sll_get_index(i)
!>
!> \section sum Summary:
!>
!> 1D
!> <table border="1">
!> <tr>
!> <th> size's problem </th>
!> <th> type of in </th>
!> <th> type of out </th>
!> <th> size of in </th>
!> <th> size of out </th>
!> <th> direction </th>
!> <th> flags </th>
!> </tr>
!> <tr>
!> <td> n </td>
!> <td> real </td>
!> <td> real </td>
!> <td> n </td>
!> <td> n </td>
!> <td> FFT_FORWARD <br /> FFT_INVERSE </td>
!> <td> FFT_NORMALIZE </td>
!> </tr>
!> <tr>
!> <td> n </td>
!> <td> complex </td>
!> <td> complex </td>
!> <td> n </td>
!> <td> n </td>
!> <td> FFT_FORWARD <br /> FFT_INVERSE </td>
!> <td> FFT_NORMALIZE </td>
!> </tr>
!> <tr>
!> <td> n </td>
!> <td> real </td>
!> <td> complex </td>
!> <td> n </td>
!> <td> n/2 </td>
!> <td> ----- </td>
!> <td> FFT_NORMALIZE </td>
!> </tr>
!> <tr>
!> <td> n </td>
!> <td> complex </td>
!> <td> real </td>
!> <td> n/2 </td>
!> <td> n </td>
!> <td> ----- </td>
!> <td> FFT_NORMALIZE </td>
!> </tr>
!> </table>
!>
!>
!> 2D
!> <table border="1">
!> <tr>
!> <th> size's problem </th>
!> <th> type of in </th>
!> <th> type of out </th>
!> <th> size of in </th>
!> <th> size of out </th>
!> <th> direction present </th>
!> <th> flags </th>
!> </tr>
!> <tr>
!> <td> n,m </td>
!> <td> real </td>
!> <td> real </td>
!> <td> n,m </td>
!> <td> n,m </td>
!> <td> FFT_FORWARD <br /> FFT_INVERSE </td>
!> <td> FFT_NORMALIZE </td>
!> </tr>
!> <tr>
!> <td> n,m </td>
!> <td> complex </td>
!> <td> complex </td>
!> <td> n,m </td>
!> <td> n,m </td>
!> <td> FFT_FORWARD <br /> FFT_INVERSE </td>
!> <td> FFT_NORMALIZE <br /> FFT_ONLY_FIRST_DIRECTION <br /> FFT_ONLY_SECOND_DIRECTION </td>
!> </tr>
!> <tr>
!> <td> n,m </td>
!> <td> real </td>
!> <td> complex </td>
!> <td> n,m </td>
!> <td> n/2,m </td>
!> <td> ----- </td>
!> <td> FFT_NORMALIZE </td>
!> </tr>
!> <tr>
!> <td> n,m </td>
!> <td> complex </td>
!> <td> real </td>
!> <td> n/2,m </td>
!> <td> n,m </td>
!> <td> ----- </td>
!> <td> FFT_NORMALIZE </td>
!> </tr>
!> </table>
!>
!>
!>
!> \section example Examples:
!>
!> In-place transform
!> \code
!> sll_int32, parameter :: n = 2**5
!> sll_comp64, dimension(0,n-1) :: in
!> type(sll_fft_plan), pointer  :: p
!>
!> !** INIT DATA **
!>
!> p => fft_new_plan(n,in,in,FFT_FORWARD,FFT_NORMALIZE)
!> call fft_apply_plan(p,in,in)
!> call fft_delete_plan(p)
!> \endcode
!>
!> Two-dimensional transform  
!> \code
!> sll_int32, parameter :: n = 2**5
!> sll_int32, parameter :: m = 2**3
!> sll_comp64, dimension(n/2,m) :: in
!> sll_real64, dimension(n,m) :: out
!> type(sll_fft_plan), pointer  :: p
!>
!> !** INIT DATA **
!>
!> p => fft_new_plan(n,in,out,FFT_INVERSE)
!> call fft_apply_plan(p,in,out)
!> call fft_delete_plan(p)
!> \endcode
!>
!> Transform in one direction
!> \code
!> sll_int32, parameter :: n = 2**5
!> sll_int32, parameter :: m = 2**3
!> sll_comp64, dimension(n,m) :: in
!> sll_comp64, dimension(n,m) :: out
!> type(sll_fft_plan), pointer  :: p
!>
!> !** INIT DATA **
!>
!> p => fft_new_plan(n,in,out,FFT_FORWARD,FFT_ONLY_FIRST_DIRECTION)
!> call fft_apply_plan(p,in,out)
!> call fft_delete_plan(p)
!> \endcode
!>
!>
! \warning let p = sll_get_index(i), if p is even data(i) is the real part of X_p, else if p is odd data(i) is the imaginary part of X_p
!>
!>
!> \section acc Access the mode
!> 
!> To get the value of a mode call the function fft_get_mode(plan,data,num_mode)
!>
!> \code
!> sll_int32, parameter :: n = 2**5
!> sll_comp64, dimension(0,n-1) :: in
!> type(sll_fft_plan), pointer  :: p
!> sll_comp64 :: mode
!> sll_int32 :: k = 3
!>
!> !** INIT DATA **
!>
!> p => fft_new_plan(n,in,in,FFT_FORWARD,FFT_NORMALIZE)
!> call fft_apply_plan(p,in,in)
!> mode = fft_get_mode(plan,in,k)
!> call fft_delete_plan(p)
!> \endcode
!>
!> \code
!> sll_int32, parameter :: n = 2**5
!> sll_int32, parameter :: m = 2**3
!> sll_comp64, dimension(n/2,m) :: in
!> sll_real64, dimension(n,m) :: out
!> type(sll_fft_plan), pointer  :: p
!> sll_comp64 :: mode
!> sll_int32 :: k=4, l=2
!>
!> !** INIT DATA **
!>
!> p => fft_new_plan(n,in,out,FFT_INVERSE)
!> call fft_apply_plan(p,in,out)
!> mode = fft_get_mode(plan,out,k,l)
!> call fft_delete_plan(p)
!> \endcode
!>
!> \code
!> sll_int32, parameter :: n = 2**5
!> sll_real64, dimension(0,n-1) :: in
!> sll_real64, dimension(0,n-1) :: out
!> type(sll_fft_plan), pointer  :: p
!> sll_comp64 :: mode
!> sll_int32 :: k = 3
!>
!> !** INIT DATA **
!>
!> p => fft_new_plan(n,in,out,FFT_FORWARD,FFT_NORMALIZE)
!> call fft_apply_plan(p,in,out)
!> mode = fft_get_mode(plan,out,k)
!> call fft_delete_plan(p)
!> \endcode
!>
!>
!>
!>
!> To set a mode call the subroutine fft_set_mode(plan,data,new_value,num_mode)
!>
!> \code
!> sll_int32, parameter :: n = 2**5
!> sll_comp64, dimension(0,n-1) :: in
!> type(sll_fft_plan), pointer  :: p
!> sll_comp64 :: new_value
!> sll_int32 :: k = 3
!>
!> !** INIT DATA **
!>
!> p => fft_new_plan(n,in,in,FFT_FORWARD,FFT_NORMALIZE)
!> call fft_apply_plan(p,in,in)
!> new_value = complex(5.0_f64,3.2_f64)
!> call fft_set_mode(plan,in,new_value,k)
!> call fft_delete_plan(p)
!> \endcode
!>
!> \code
!> sll_int32, parameter :: n = 2**5
!> sll_int32, parameter :: m = 2**3
!> sll_comp64, dimension(n/2,m) :: in
!> sll_real64, dimension(n,m) :: out
!> type(sll_fft_plan), pointer  :: p
!> sll_comp64 :: new_value
!> sll_int32 :: k=4, l=2
!>
!> !** INIT DATA **
!>
!> p => fft_new_plan(n,in,out,FFT_INVERSE)
!> call fft_apply_plan(p,in,out)
!> new_value = complex(5.0_f64,3.2_f64)
!> call fft_set_mode(plan,out,new_value,k,l)
!> call fft_delete_plan(p)
!> \endcode
!>
!> \code
!> sll_int32, parameter :: n = 2**5
!> sll_real64, dimension(0,n-1) :: in
!> sll_real64, dimension(0,n-1) :: out
!> type(sll_fft_plan), pointer  :: p
!> sll_comp64 :: new_value
!> sll_int32 :: k = 0
!>
!> !** INIT DATA **
!>
!> p => fft_new_plan(n,in,out,FFT_FORWARD,FFT_NORMALIZE)
!> call fft_apply_plan(p,in,out)
!> new_value = complex(5.0_f64,3.2_f64)
!> call fft_set_mode(plan,out,new_value,k)
!> call fft_delete_plan(p)
!> \endcode
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
! For the real transform, we have
! \f$ (x_0,x_1,\dots,x_{n-1}) \rightarrow
!     (r_0,r_{n/2},r_1,i_1,\dots,r_{n/2-1},i_{n/2-1})\f$
! which must be interpreted as the complex array
! \f[ \begin{pmatrix} r_0 &,& 0
!                     \\ r_1 &,& i_1
!                     \\ \vdots  & & \vdots 
!                     \\ r_{n/2-1} &,& i_{n/2-1}
!                     \\ r_{n/2} &,& 0
!                     \\ r_{n/2-1} &,& -i_{n/2-1}
!                     \\ \vdots    & & \vdots
!                     \\ r_1 &,& -i_1 
! \end{pmatrix}\f] 
! \warning Note that ffw use \f$(r_0,r_1,\dots,r_{n/2-1},r_{n/2},i_{n/2-1},\dots,i_1)\f$
!          convention whereas fftpack use \f$(r_0,r_1,i_1,\dots,r_{n/2-1},i_{n/2-1},r_{n/2})\f$
! 
!
! By example if, the input data is \f$(x_0,x_1,x_2,x_3)\f$ the output is \f$(X_0,X_2,X_1,X_3)\f$.
! Thus, sll_get_index(1) returns 2 (cause data[1]=X_2) and sll_get_mode(1) returns X_1.
!------------------------------------------------------------------------------

module sll_fft 
#include "sll_working_precision.h"
#include "misc_utils.h"
#include "sll_assert.h"
#include "sll_memory.h"
  use numeric_constants
  use sll_timer
#ifndef _NOFFTW
  use, intrinsic :: iso_c_binding
#endif
  implicit none
#ifndef _NOFFTW
  include 'fftw3.f03'
#endif

  ! Basic 1D FFT descriptor:
  ! - N: number of samples
  ! - style : transform style
  ! - index : given a wavenumber on 1:N, return the corresponding array 
  !           index (also on 1:N)
  ! - mode  : given an array index, return the corresponding wavenumber
  !type sll_fft_plan
  !  sll_int32 :: N,N2,N3
  !  sll_int32 :: style !, flags
  !  !sll_int32 :: data_type
  !  sll_int32,  dimension(:), pointer :: index => null()
  !  sll_int32,  dimension(:), pointer :: mode => null()
  !  sll_comp64, dimension(:), pointer :: t => null()          ! twiddle factors complex case
  !  sll_real64, dimension(:), pointer :: twiddles => null()  ! twiddles factors real case 
  !  sll_real64, dimension(:), pointer :: twiddles_n => null() ! twiddles factors real case 
  !  sll_real64, dimension(:), pointer :: dwsave => null() ! for use fftpack
  !  !sll_int32                         :: mod !type of library used
  !  !sll_real64                        :: normalization_factor_forward = 1.0_f64
  !  !sll_real64                        :: normalization_factor_backward = 1.0_f64
  !  sll_int32                         :: library
  !  type(C_PTR)                       :: fftw_plan
  !  sll_int32                         :: direction
  !  sll_comp64, dimension(:), pointer  :: in_comp, out_comp 
  !  sll_real64, dimension(:), pointer  :: in_real, out_real
  !  sll_real64 :: fft_time_execution
  !end type sll_fft_plan

!  type fft_plan
!    type(C_PTR)                     :: fftw_plan
!
!    sll_comp64, dimension(:), pointer :: t => null()          ! twiddle factors complex case
!    sll_real64, dimension(:), pointer :: twiddles => null()  ! twiddles factors real case 
!    sll_real64, dimension(:), pointer :: twiddles_n => null() ! twiddles factors real case 
!
!    sll_int32                       :: style
!    sll_int32                       :: library
!    sll_int32                       :: direction
!    sll_int32                       :: problem_rank
!    sll_int32, allocatable          :: problem_shape(:)
!  end type fft_plan

  type sll_fft_plan
#ifndef _NOFFTW
    type(C_PTR)                     :: fftw_plan
#endif
    sll_comp64, dimension(:), pointer :: t => null()          ! twiddle factors complex case
    sll_real64, dimension(:), pointer :: twiddles => null()  ! twiddles factors real case 
    sll_real64, dimension(:), pointer :: twiddles_n => null() ! twiddles factors real case 

    sll_int32                       :: style
    sll_int32                       :: library
    sll_int32                       :: direction
    sll_int32                       :: problem_rank
    sll_int32, pointer              :: problem_shape(:)
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

  sll_int32, parameter :: FFT_FORWARD = -1
  sll_int32, parameter :: FFT_INVERSE = 1
  sll_int32, parameter :: FFT_NORMALIZE_FORWARD = 2**0
  sll_int32, parameter :: FFT_NORMALIZE_INVERSE = 2**0
  sll_int32, parameter :: FFT_NORMALIZE         = 2**0

  sll_int32, parameter :: FFT_ONLY_FIRST_DIRECTION  = 2**2
  sll_int32, parameter :: FFT_ONLY_SECOND_DIRECTION = 2**3
  sll_int32, parameter :: FFT_ONLY_THIRD_DIRECTION  = 2**4

  sll_int32, parameter :: FFTW_MOD = 1000000000
  sll_int32, parameter :: SLLFFT_MOD = 0
  sll_int32, parameter :: FFTPACK_MOD = 100

!need to be change on simple parameter
#define FFTW_MOD 1000000000
#define SLLFFT_MOD 0
#define FFTPACK_MOD 100

  interface bit_reverse
    module procedure bit_reverse_complex, bit_reverse_integer32, &
                     bit_reverse_integer64
  end interface
  
  interface fft_get_mode
     module procedure fft_get_mode_complx_1d, fft_get_mode_complx_2d, &
                      fft_get_mode_complx_3d, fft_get_mode_real_1d
  end interface

  interface fft_set_mode
     module procedure fft_set_mode_complx_1d, fft_set_mode_complx_2d, &
                      fft_set_mode_complx_3d, fft_set_mode_real_1d
  end interface

!----------------------------------------------------------------------------------
! -                                  NEW INTERFACE                                -
!----------------------------------------------------------------------------------

  interface fft_new_plan
    module procedure new_plan_c2c_2d, new_plan_c2r_2d, new_plan_r2c_2d,&
                     fft_new_r2r_1d, new_plan_r2c_1d, new_plan_c2r_1d,&
                     fft_new_c2c_1d_for_1d
  end interface
  interface fft_apply_plan
    module procedure apply_plan_c2c_2d, fft_apply_r2c_2d, fft_apply_c2r_2d,&
                     fft_apply_r2r_1d, fft_apply_r2c_1d, fft_apply_c2r_1d,&
                     fft_apply_c2c_1d
  end interface 

  interface fft_delete_plan
    module procedure fft_delete
  end interface
  
!----------------------------------------------------------------------------------

! -----------------------------------------------
! -                 OLD INTERFACE               -
! -----------------------------------------------
  interface new_plan_c2c_1d
      module procedure fft_new_c2c_1d_for_1d
  end interface  

  interface apply_fft_c2c_1d
      module procedure fft_apply_c2c_1d
  end interface 

  interface delete_fft_plan1d
    module procedure fft_delete 
  end interface
! -----------------------------------------------

  interface fft_is_present_flag
    module procedure fft_is_present_flag_in_plan,fft_is_present_flag_in_integer
  end interface

contains

  ! Return NO if the flag is disabled in the plan or YES if enabled.
  function fft_is_present_flag_in_plan(plan,flag) result(bool)
    type(sll_fft_plan), pointer             :: plan
    sll_int32, intent(in)                   :: flag
    logical                                 :: bool
    sll_int32                               :: m
   
    SLL_ASSERT( is_power_of_two( int(flag,kind=i64) ) )

    m = iand(plan%style,flag)
    if( m .eq. flag ) then
      bool = .true.
    else
      bool = .false.
    endif 
  end function

  ! Return NO if the flag (bit) is disabled in the integer or YES if enabled.
  function fft_is_present_flag_in_integer(intege,bit) result(bool)
    sll_int32, intent(in)                   :: bit, intege
    logical                                 :: bool
    sll_int32                               :: m
   
    SLL_ASSERT( is_power_of_two( int(bit,kind=i64) ) )

    m = iand(intege,bit)
    if( m .eq. bit ) then
      bool = .true.
    else
      bool = .false.
    endif 
  end function

  function fft_get_rank(plan) result(rank)
    type(sll_fft_plan), pointer :: plan
    sll_int32               :: rank
    rank = plan%problem_rank
  end function

  function fft_get_shape(plan) result(etendue)
    type(sll_fft_plan), pointer :: plan
    sll_int32, pointer  :: etendue(:)
    allocate(etendue(fft_get_rank(plan)))
    etendue = plan%problem_shape
  end function

  function fft_get_mode_complx_1d(plan,array,k) result(mode)
    type(sll_fft_plan), pointer :: plan
    sll_comp64, dimension(0:)   :: array
    sll_int32                   :: k
    sll_comp64                  :: mode
    mode = array(k)
  end function

  function fft_get_mode_complx_2d(plan,array,k,l) result(mode)
    type(sll_fft_plan), pointer   :: plan
    sll_comp64, dimension(0:,0:)  :: array
    sll_int32                     :: k, l
    sll_comp64                    :: mode
    mode = array(k,l)
  end function

  function fft_get_mode_complx_3d(plan,array,k,l,m) result(mode)
    type(sll_fft_plan), pointer     :: plan
    sll_comp64, dimension(0:,0:,0:) :: array
    sll_int32                       :: k, l, m
    sll_comp64                      :: mode
    mode = array(k,l,m)
  end function

  function fft_get_mode_real_1d(plan,data,k) result(mode)
    type(sll_fft_plan), pointer :: plan
    sll_real64, dimension(0:)   :: data
    sll_int32                   :: k, n_2, n
    sll_comp64                  :: mode

    if(fft_get_rank(plan) .ne. 1) then
      print*,'Error in sll_fft.F90'
      print*,'      in function fft_get_mode_real_1d'
      print*,'      fata can be only 1d'
      stop
    endif

    n = plan%problem_shape(1)
    n_2 = n/2 !ishft(n,-1)

    if(plan%library .eq. FFTPACK_MOD) then
      if( k .eq. 0 ) then
        mode = cmplx(data(0),0.0_f64,kind=f64)
      else if( k .eq. n_2 ) then
        mode = cmplx(data(n-1),0.0_f64,kind=f64)
      else if( k .gt. n_2 ) then
        mode = cmplx( data(2*(n-k)-1) , -data(2*(n-k)) ,kind=f64)
      else
        mode = cmplx( data(2*k-1) , data(2*k) ,kind=f64)
      endif
    else if(plan%library .eq. SLLFFT_MOD) then
      if( k .eq. 0 ) then
        mode = cmplx(data(0),0.0_f64,kind=f64)
      else if( k .eq. n_2 ) then
        mode = cmplx(data(1),0.0_f64,kind=f64)
      else if( k .gt. n_2 ) then
        mode = cmplx( data(2*(n-k)) , -data(2*(n-k)+1),kind=f64 )
      else
        mode = cmplx( data(2*k) , data(2*k+1) ,kind=f64)
      endif
    else if(plan%library .eq. FFTW_MOD) then
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
    else
      stop 'ERROR IN =fft_get_mode_real_1d= plan%library invalid'
    endif
  end function

  subroutine fft_set_mode_complx_1d(plan,array,new_value,k)
    type(sll_fft_plan), pointer :: plan
    sll_comp64, dimension(0:)   :: array
    sll_int32                   :: k
    sll_comp64                  :: new_value
    array(k) = new_value
  end subroutine


  subroutine fft_set_mode_complx_2d(plan,array,new_value,k,l)
    type(sll_fft_plan), pointer :: plan
    sll_comp64, dimension(0:,0:)   :: array
    sll_int32                   :: k,l
    sll_comp64                  :: new_value
    array(k,l) = new_value
  end subroutine


  subroutine fft_set_mode_complx_3d(plan,array,new_value,k,l,m)
    type(sll_fft_plan), pointer :: plan
    sll_comp64, dimension(0:,0:,0:)   :: array
    sll_int32                   :: k,l,m
    sll_comp64                  :: new_value
    array(k,l,m) = new_value
  end subroutine

  subroutine fft_set_mode_real_1d(plan,data,new_value,k)
    type(sll_fft_plan), pointer :: plan
    sll_real64, dimension(0:)   :: data
    sll_int32                   :: k, n_2, n, index_mode
    sll_comp64                  :: new_value

    if(fft_get_rank(plan) .ne. 1) then
      print*,'Error in sll_fft.F90'
      print*,'      in function fft_set_mode_real_1d'
      print*,'      data can be only 1d'
      stop
    endif

    n = plan%problem_shape(1)
    n_2 = n/2 !ishft(n,-1)

    if(plan%library .eq. FFTPACK_MOD) then
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
    else if(plan%library .eq. SLLFFT_MOD) then
      if( k .eq. 0 ) then
        data(0) = real(new_value,kind=f64)
      else if( k .eq. n_2 ) then
        data(1) = real(new_value,kind=f64)
      else if( k .gt. n_2 ) then
        data(2*(n-k)) = real(new_value,kind=f64)
        data(2*(n-k)+1) = -dimag(new_value)
      else
        data(2*k) = real(new_value,kind=f64)
        data(2*k+1) = dimag(new_value)
      endif
    else if(plan%library .eq. FFTW_MOD) then
      if( k .eq. 0 ) then
        data(0) = real(new_value,kind=f64)
      else if( k .eq. n_2 ) then
        data(n_2) = real(new_value,kind=f64)
      else if( k .gt. n_2 ) then
        data(n-k) = real(new_value,kind=f64)
        data(k) = -dimag(new_value)
      else
        data(k) = real(new_value,kind=f64)
        data(n-k) = dimag(new_value)
      endif
    else
      stop 'ERROR IN =fft_set_mode_real_1d= plan%library invalid'
    endif
  end subroutine 

#define INFO plan%fft_time_execution = time

!#define INFO print *,'-----------------------------';\
!             if(plan%library .eq. SLLFFT_MOD) then;  \
!               print *, 'LIBRARY :: SELALIB';        \
!             endif;                                  \
!             if(plan%library .eq. FFTPACK_MOD) then; \
!               print *, 'LIBRARY :: FFTPACK';        \
!             endif;                                  \
!             if(plan%library .eq. FFTW_MOD) then;    \
!               print *, 'LIBRARY :: FFTW';           \
!             endif;                                  \
!             open(1,file="info.txt",position='append');\
!             print *, 'TIME : ',time;                \
!             plan%fft_time_execution = time;         \
!             write(1,*) time;                        \
!             close(1);                               \
!             print *,'-----------------------------'


#ifndef _DEFAULTFFTLIB
#ifdef _NOFFTW
#define _DEFAULTFFTLIB SLLFFT_MOD
#else
#define _DEFAULTFFTLIB FFTW_MOD
#endif 
#endif


  subroutine print_defaultfftlib()
    print *, '----------------'
    print *,  _DEFAULTFFTLIB
    print *, '----------------'
  end subroutine

  function fft_default_lib_is(lib) result(bool)
    sll_int32 :: lib
    logical   :: bool
    
    if(_DEFAULTFFTLIB .eq. lib ) then
      bool = .true.
    else
      bool = .false.
    endif 
  end function fft_default_lib_is

#ifndef _NOFFTW
#include "fftw_lib.F90"
#endif

#ifndef _NOFFTPACK
#include "fftpack_lib.F90"
#endif

#ifndef _NOFFTSLL
#include "my_lib.F90"
#endif

!#if ( _DEFAULTFFTLIB == FFTW_MOD ) || ( _DEFAULTFFTLIB == FFTPACK_MOD ) 
!  function get_fft_mode(plan,data,k) result(mode)
!    sll_comp64                                :: mode
!    type(sll_fft_plan)                        :: plan
!    sll_int32, intent(in)                     :: k
!    sll_comp64, dimension(0:) , intent(in)    :: data
!    mode = data(k)
!  end function
!#else
!  function get_fft_mode(plan,data,k) result(mode)
!    sll_comp64                                :: mode
!    type(sll_fft_plan)                        :: plan
!    sll_int32, intent(in)                     :: k
!    sll_comp64, dimension(0:) , intent(in)    :: data
!    !mode = data(plan%mode(k))
!    mode = data(k)
!  end function
!#endif
!
!  function sll_get_mode_complex_array(plan,data,k) result(mode)
!    sll_comp64                                :: mode
!    sll_int32, intent(in)                     :: k
!    type(sll_fft_plan)                        :: plan
!    sll_comp64, dimension(0:) , intent(in)    :: data
!    sll_int32                                 :: n_2, n
!  
!    !mode = data(0) 
!    !return
!
!    n = plan%N
!    n_2 = n/2
!
!    SLL_ASSERT( (k .ge. 0) .and. (k .lt. n) )
!    SLL_ASSERT( size(data) .eq. n )
!   
!    if(plan%library .eq. FFTW_MOD) then
!      mode = data(k)
!    else if(plan%library .eq. FFTPACK_MOD) then
!      mode = data(k)
!    else if(plan%library .eq. SLLFFT_MOD) then
!      mode = data(plan%mode(k))
!    else
!      stop 'ERROR IN =SLL_GET_MODE_COMPLEX= plan%mod invalid'
!    endif
!  end function sll_get_mode_complex_array

!  function sll_get_index_real_array(plan,k) result(index)
!    sll_int32                                 :: index
!    sll_int32, intent(in)                     :: k
!    type(sll_fft_plan) , intent(in)           :: plan
!    sll_int32                                 :: n_2, n
!
!    index = plan%index(k)
!    return  
! 
!    n = plan%N
!    n_2 = n/2
!
!    SLL_ASSERT( (k .ge. 0) .and. (k .lt. n) )
!
!    if(plan%library .eq. FFTPACK_MOD) then
!      if( k .eq. 0 ) then
!        index = 0
!      else if( k .eq. n_2 ) then
!        index = n_2
!      !else if( k/2 .eq. (k+1)/2 ) then !k even
!      !  index = k/2
!      else if( k/2 .eq. (k+1)/2 ) then !k even
!        index = n - k/2
!      else
!        index = (k+1)/2 
!      endif
!    else if(plan%library .eq. SLLFFT_MOD) then
!      if( k .eq. 0 ) then
!        index = 0
!      else if( k .eq. 1 ) then
!        index = n_2
!      else if( k/2 .eq. (k+1)/2 ) then !k even
!        index = k/2
!      else
!        index = n - k/2
!      endif
!    else if(plan%library .eq. FFTW_MOD) then
!      if( k .eq. 0 ) then
!        index = 0
!      else if( k .eq. n-1 ) then
!        index = n_2
!      else if( k/2 .eq. (k+1)/2 ) then !k even
!        index = n - k/2
!      else
!        index = k/2
!      endif
!    else
!      stop 'ERROR IN =SLL_GET_MODE_REAL= plan%mod invalid'
!    endif
!  end function sll_get_index_real_array
!
!  function sll_get_mode_real_array(plan,data,k) result(mode)
!    sll_comp64                                :: mode
!    sll_int32, intent(in)                     :: k
!    type(sll_fft_plan) , pointer              :: plan
!    sll_real64, dimension(0:) , intent(in)    :: data
!    sll_int32                                 :: n_2, n
!   
!    n = plan%N
!    n_2 = n/2 !ishft(n,-1)
!
!    if( k .eq. 0 ) then
!      mode = complex(data(0),0.0_f64)
!      return
!    else if( k .eq. n_2 ) then
!      mode = complex(data(1),0.0_f64)
!      return
!    else if( k .gt. n_2 ) then
!      mode = complex( data(2*(n-k)) , -data(2*(n-k)+1) )
!      return
!    else
!      mode = complex( data(2*k) , data(2*k+1) )
!      return
!    endif
!
!    if(plan%library .eq. FFTPACK_MOD) then
!      if( k .eq. 0 ) then
!        mode = complex(data(0),0.0_f64)
!      else if( k .eq. n_2 ) then
!        mode = complex(data(n-1),0.0_f64)
!      else if( k .gt. n_2 ) then
!        mode = cmplx( data(2*(n-k)-1) , -data(2*(n-k)) ,kind=f64)
!      else
!        mode = complex( data(2*k-1) , data(2*k) )
!      endif
!    else if(plan%library .eq. SLLFFT_MOD) then
!      if( k .eq. 0 ) then
!        mode = complex(data(0),0.0_f64)
!      else if( k .eq. n_2 ) then
!        mode = complex(data(1),0.0_f64)
!      else if( k .gt. n_2 ) then
!        mode = complex( data(2*(n-k)) , -data(2*(n-k)+1) )
!      else
!        mode = complex( data(2*k) , data(2*k+1) )
!      endif
!    else if(plan%library .eq. FFTW_MOD) then
!      if( k .eq. 0 ) then
!        mode = complex(data(0),0.0_f64)
!      else if( k .eq. n_2 ) then
!        mode = complex(data(n_2),0.0_f64)
!      else if( k .gt. n_2 ) then
!        !mode = complex( data(k-n_2) , -data(n-k+n_2) )
!        mode = complex( data(n-k) , -data(k) )
!      else
!        mode = complex( data(k) , data(n-k) )
!      endif
!    else
!      stop 'ERROR IN =SLL_GET_MODE_REAL= plan%library invalid'
!    endif
!  end function sll_get_mode_real_array


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


                                  !!!!!!!!!!!!!!!!!!!!!!!!!!
                                  !!  COMPLEX TO COMPLEX  !!
                                  !!!!!!!!!!!!!!!!!!!!!!!!!!

! ---------------------
! - DEFAULT INTERFACE -
! ---------------------
 function fft_new_c2c_1d_for_1d(size_problem,array_in,array_out,direction,flags) result(plan)
    sll_int32, intent(in)                        :: size_problem
    sll_comp64, dimension(:), target, intent(inout) :: array_in, array_out
    sll_int32, intent(in)                        :: direction
    sll_int32, optional, intent(in)              :: flags
    type(sll_fft_plan), pointer                       :: plan

#if defined(_NOFFTW) && _DEFAULTFFTLIB==FFTW_MOD
    stop 'The default library cannot be FFTW because she is not installed'
#endif
    plan => fft_plan_c2c_1d(_DEFAULTFFTLIB,size_problem,array_in,array_out,direction,flags)
  end function

! function fft_new_c2c_1d_for_3d(NX,NY,NZ,array_in,array_out,direction,flags) result(plan)
!    sll_int32, intent(in)                            :: NX,NY,NZ
!    sll_comp64, dimension(:,:,:), target, intent(in) :: array_in, array_out
!    sll_int32, intent(in)                        :: direction
!    sll_int32, optional, intent(in)              :: flags
!    type(sll_fft_plan), pointer                      :: plan
!
!#if defined(_NOFFTW) && _DEFAULTFFTLIB==FFTW_MOD
!    stop 'The default library cannot be FFTW because she is not installed'
!#endif
!    if( _DEFAULTFFTLIB .eq. FFTW_MOD ) then
!      stop 'FFTW not available for DFT in one direction on multidimensional array'
!    endif
!    plan => fft_plan_c2c_1d_for_3d(_DEFAULTFFTLIB,NX,NY,NZ,array_in,array_out,direction,flags)
!  end function

 function new_plan_c2c_2d(NX,NY,array_in,array_out,direction,flags) result(plan)
    sll_int32, intent(in)                          :: NX,NY
    sll_comp64, dimension(0:,0:), target, intent(in) :: array_in, array_out
    sll_int32, intent(in)                          :: direction
    sll_int32, optional, intent(in)                :: flags
    type(sll_fft_plan), pointer                        :: plan
    
#if defined(_NOFFTW) && _DEFAULTFFTLIB==FFTW_MOD
    stop 'The default library cannot be FFTW because she is not installed'
#endif

    if( (size(array_in(:,0)) .ne. NX) .or. (size(array_in(0,:)) .ne. NY) &
          .or. (size(array_out(:,0)) .ne. NX) .or. (size(array_out(0,:)) .ne. NY) ) then
      stop 'Error in new_plan_c2c_2d size problem'
    endif 

    if( present(flags) .and. &
      ( fft_is_present_flag(flags,FFT_ONLY_FIRST_DIRECTION) .or. &
        fft_is_present_flag(flags,FFT_ONLY_SECOND_DIRECTION) ) ) then
      if( _DEFAULTFFTLIB .ne. SLLFFT_MOD ) &
        stop 'option FFT_ONLY_FIRST_DIRECTION available only with Selalib FFT '     
    endif

    plan => fft_plan_c2c_2d(_DEFAULTFFTLIB,NX,NY,array_in,array_out,direction,flags)
  end function

! ---------------------
! -  DEBUG INTERFACE  -
! ---------------------
  function fft_plan_c2c_1d(library,size_problem,array_in,array_out,direction,flags) result(plan)
    sll_int32, intent(in)                        :: library
    sll_int32, intent(in)                        :: size_problem
    sll_comp64, dimension(:), intent(inout)      :: array_in, array_out
    sll_int32, intent(in)                        :: direction
    sll_int32, optional, intent(in)              :: flags
    type(sll_fft_plan), pointer                  :: plan

#ifndef _NOFFTW
    if(library .eq. FFTW_MOD) then
      plan => fftw_new_plan_c2c_1d(library,size_problem,array_in,array_out,direction,flags)
      return
    endif
#endif
#ifndef _NOFFTSLL
    if(library .eq. SLLFFT_MOD) then
      plan => sll_new_plan_c2c_1d(library,size_problem,array_in,array_out,direction,flags)
      return
    endif
#endif
#ifndef _NOFFTPACK
    if(library .eq. FFTPACK_MOD) then
      plan => fftpack_new_plan_c2c_1d(library,size_problem,array_in,array_out,direction,flags)
      return
    endif
#endif
    stop 'ERROR in =fft_plan_c2c_1d= library unknown'
  end function

!  function fft_plan_c2c_1d_for_3d(library,NX,NY,NZ,array_in,array_out,direction,flags) result(plan)
!    sll_int32, intent(in)                        :: library
!    sll_int32, intent(in)                        :: NX,NY,NZ
!    sll_comp64, dimension(:,:,:), target, intent(in) :: array_in, array_out
!    sll_int32, intent(in)                        :: direction
!    sll_int32, optional, intent(in)              :: flags
!    type(sll_fft_plan), pointer                  :: plan
!
!#ifndef _NOFFTW
!    if(library .eq. FFTW_MOD) then
!      plan => fftw_new_plan_c2c_1d(library,size_problem,array_in,array_out,direction,flags)
!      return
!    endif
!#endif
!#ifndef _NOFFTSLL
!    if(library .eq. SLLFFT_MOD) then
!      plan => sll_new_plan_c2c_1d_for_3d(library,NX,NY,NZ,array_in,array_out,direction,flags)
!      return
!    endif
!#endif
!#ifndef _NOFFTPACK
!    if(library .eq. FFTPACK_MOD) then
!      plan => fftpack_new_plan_c2c_1d(library,size_problem,array_in,array_out,direction,flags)
!      return
!    endif
!#endif
!    stop 'ERROR in =fft_plan_c2c_1d= library unknown'
!  end function

  function fft_plan_c2c_2d(library,NX,NY,array_in,array_out,direction,flags) result(plan)
    sll_int32, intent(in)                          :: library
    sll_int32, intent(in)                          :: NX,NY
    sll_comp64, dimension(0:,0:), target           :: array_in, array_out
    sll_int32, intent(in)                          :: direction
    sll_int32, optional, intent(in)                :: flags
    type(sll_fft_plan), pointer                        :: plan

#ifndef _NOFFTW
    if(library .eq. FFTW_MOD) then
      plan => fftw_new_plan_c2c_2d(library,NX,NY,array_in,array_out,direction,flags)
      return
    endif
#endif
#ifndef _NOFFTSLL
    if(library .eq. SLLFFT_MOD) then
      plan => sll_new_plan_c2c_2d(library,NX,NY,array_in,array_out,direction,flags)
      return
    endif
#endif
#ifndef _NOFFTPACK
    if(library .eq. FFTPACK_MOD) then
      print*, 'not available in fftpack'
      stop ''
    endif
#endif
    stop 'ERROR in =fft_plan_c2c_2d= library unknown'
  end function

! ---------------------
! -       APPLY       -
! ---------------------
  subroutine fft_apply_c2c_1d(plan,array_in,array_out)
    type(sll_fft_plan), pointer                     :: plan
    sll_comp64, dimension(0:), intent(inout)        :: array_in, array_out
    type(time_mark), pointer                        :: mark 
    sll_real64                                      :: factor, time
    sll_int32 :: nx

   if( .not. associated(plan) ) then
     print*,'Error in the sll_fft.F90'
     print*,'      in subroutine fft_apply_c2c_1d'
     print*,'      plan not associated'
     stop ''
   endif

    nx = plan%problem_shape(1)
    mark => new_time_mark() 

#ifndef _NOFFTW
    if(plan%library .eq. FFTW_MOD) then
      mark => start_time_mark(mark)
      call fftw_apply_fft_complex(plan, array_in, array_out)
      time = time_elapsed_since(mark)
    endif
#endif
#ifndef _NOFFTSLL
    if(plan%library .eq. SLLFFT_MOD) then
      mark => start_time_mark(mark)
      call sll_apply_fft_complex(plan, array_in, array_out)
      time = time_elapsed_since(mark)
    endif
#endif
#ifndef _NOFFTPACK
    if(plan%library .eq. FFTPACK_MOD) then
      mark => start_time_mark(mark)
      call fftpack_apply_fft_complex(plan, array_in, array_out)
      time = time_elapsed_since(mark)
    endif
#endif
#ifdef _FFTINFO
INFO
#endif
    if( fft_is_present_flag(plan,FFT_NORMALIZE) ) then
      factor = 1.0_f64/real(nx,kind=f64)
      array_out = factor*array_out
    endif
  end subroutine

!  subroutine fft_apply_c2c_1d_3d(plan,array_in,array_out)
!    type(sll_fft_plan), pointer, intent(in)         :: plan
!    sll_comp64, dimension(0:,0:,0:), intent(inout)  :: array_in, array_out
!#ifndef _NOFFTW
!    if(plan%library .eq. FFTW_MOD) &
!      stop 'apply_plan_c2c_1d with 3d data not available with fftw'
!#endif
!#ifndef _NOFFTSLL
!    if(plan%library .eq. SLLFFT_MOD) &
!      call sll_apply_fft_c2c_1d_for_3d(plan, array_in, array_out)
!#endif
!#ifndef _NOFFTPACK
!    if(plan%library .eq. FFTPACK_MOD) &
!      stop 'apply_plan_c2c_1d with 3d data not available with fftpack'
!#endif
!  end subroutine


  subroutine apply_plan_c2c_2d(plan,array_in,array_out)
    type(sll_fft_plan), pointer                     :: plan
    sll_comp64, dimension(0:,0:), intent(inout)     :: array_in, array_out
    sll_real64                                      :: factor

    if( .not. associated(plan) ) then
      print*,'Error in the file sll_fft.F90'
      print*,'      in function apply_plan_c2c_2d'
      print*,'      plan not associated'
      stop ''
    endif

#ifndef _NOFFTW
    if(plan%library .eq. FFTW_MOD) then
      call fftw_apply_plan_c2c_2d(plan, array_in, array_out)
   endif
#endif
#ifndef _NOFFTSLL
    if(plan%library .eq. SLLFFT_MOD) then
        call sll_apply_fft_c2c_2d(plan, array_in, array_out)
    endif
#endif
#ifndef _NOFFTPACK
    if(plan%library .eq. FFTPACK_MOD) &
      stop 'apply_plan_c2c_2d not available with fftpack'
#endif
    if( fft_is_present_flag(plan,FFT_NORMALIZE) ) then
      if( fft_is_present_flag(plan,FFT_ONLY_FIRST_DIRECTION) .and. fft_is_present_flag(plan,FFT_ONLY_SECOND_DIRECTION) ) then
        factor = 1.0_f64/real( plan%problem_shape(1)*plan%problem_shape(2) ,kind=f64)
      else if( .not. fft_is_present_flag(plan,FFT_ONLY_FIRST_DIRECTION) .and. .not. fft_is_present_flag(plan,FFT_ONLY_SECOND_DIRECTION) ) then
        factor = 1.0_f64/real( plan%problem_shape(1)*plan%problem_shape(2) ,kind=f64)
      else if( fft_is_present_flag(plan,FFT_ONLY_FIRST_DIRECTION) ) then
        factor = 1.0_f64/real( plan%problem_shape(1),kind=f64)
      else
        factor = 1.0_f64/real( plan%problem_shape(2),kind=f64)
      endif
      array_out = factor*array_out
    endif
  end subroutine











                 !!!!!!!!!!!!!!!!!!!!!!!
                 !!    REAL TO REAL   !!
                 !!!!!!!!!!!!!!!!!!!!!!!

! ---------------------
! - DEFAULT INTERFACE -
! ---------------------
  function fft_new_r2r_1d(n,array_in,array_out,direction,flags) result(plan)
    !sll_int32, intent(in)                        :: rank
    sll_int32, intent(in)                        :: n
    sll_real64, dimension(:), intent(inout)      :: array_in, array_out
    sll_int32, intent(in)                        :: direction
    sll_int32, optional, intent(in)              :: flags
    type(sll_fft_plan), pointer                  :: plan
    sll_int32                                    :: nx

#if defined(_NOFFTW) && _DEFAULTFFTLIB==FFTW_MOD
    print * ,'Error in file sll_fft.F90'
    print * ,'      function fft_new_r2r_1d'
    stop 'The default library cannot be FFTW because she is not installed'
#endif
  
    nx = n

    if( size(array_in).ne.nx .or. size(array_out).ne.nx) then
      print * ,'Error in file sll_fft.F90'
      print * ,'      function fft_new_r2r_1d'
      stop 'size problem'
    endif
    plan => fft_plan_r2r_1d(_DEFAULTFFTLIB,nx,array_in,array_out,direction,flags)
  end function

! ---------------------
! -  DEBUG INTERFACE  -
! ---------------------
  function fft_plan_r2r_1d(library,size_problem,array_in,array_out,direction,flags) result(plan)
    sll_int32, intent(in)                        :: library
    sll_int32, intent(in)                        :: size_problem
    sll_real64, dimension(:), intent(inout)      :: array_in
    sll_real64, dimension(:), intent(inout)      :: array_out
    sll_int32, intent(in)                        :: direction
    sll_int32, optional, intent(in)              :: flags
    type(sll_fft_plan), pointer                      :: plan

#ifndef _NOFFTW 
    if(library .eq. FFTW_MOD) then
      plan => fftw_new_plan_r2r_1d(library,size_problem,array_in,array_out,direction,flags)
      return
    endif
#endif
#ifndef _NOFFTSLL
    if(library .eq. SLLFFT_MOD) then
      plan => sll_new_plan_r2r_1d(library,size_problem,array_in,array_out,direction,flags)
      return
    endif
#endif
#ifndef _NOFFTPACK
    if(library .eq. FFTPACK_MOD) then
      plan => fftpack_new_plan_r2r_1d(library,size_problem,array_in,array_out,direction,flags)
      return
    endif
#endif

    print * ,'Error in file sll_fft.F90'
    print * ,'      function fft_plan_r2r_1d'
    stop 'library unknown'
  end function

! ---------------------
! -       APPLY       -
! ---------------------
  subroutine fft_apply_r2r_1d(plan,array_in,array_out)
    type(sll_fft_plan), pointer                     :: plan
    sll_real64, dimension(:), intent(inout)         :: array_in
    sll_real64, dimension(:), intent(inout)         :: array_out
    sll_real64                                      :: factor, time
    type(time_mark), pointer                        :: mark
    sll_int32                                       :: nx

    if( .not. associated(plan) ) then
      print*,'Error in the sll_fft.F90'
      print*,'      in subroutine fft_apply_r2r_1d'
      print*,'      plan not associated'
      stop ''
    endif

    nx = plan%problem_shape(1)
    mark => new_time_mark()

    if( size(array_in).ne.nx .or. size(array_out).ne.nx) then
      print * ,'Error in file sll_fft.F90'
      print * ,'      function fft_new_r2r_1d'
      stop 'size problem'
    endif

#ifndef _NOFFTW
    if(plan%library .eq. FFTW_MOD) then
      mark => start_time_mark(mark)
      call fftw_apply_fft_real(plan, array_in, array_out)
      time = time_elapsed_since(mark)
    endif
#endif
#ifndef _NOFFTSLL
    if(plan%library .eq. SLLFFT_MOD) then
      mark => start_time_mark(mark)
      call sll_apply_fft_real(plan, array_in, array_out)
      time = time_elapsed_since(mark)
    endif
#endif
#ifndef _NOFFTPACK
    if(plan%library .eq. FFTPACK_MOD) then
      mark => start_time_mark(mark)
      call fftpack_apply_fft_real(plan, array_in, array_out)
      time = time_elapsed_since(mark)
    endif
#endif

#ifdef _FFTINFO
INFO
#endif

    if( fft_is_present_flag(plan,FFT_NORMALIZE) ) then
      factor = 1.0_f64/real(nx,kind=f64)
      array_out = factor*array_out
    endif
  end subroutine







                                      !!!!!!!!!!!!!!!!!!!!!!!
                                      !!  REAL TO COMPLEX  !!
                                      !!!!!!!!!!!!!!!!!!!!!!!

! ---------------------
! - DEFAULT INTERFACE -
! ---------------------
 function new_plan_r2c_1d(nx,array_in,array_out,flags) result(plan)
    sll_int32, intent(in)                :: nx
    sll_real64, dimension(:)             :: array_in
    sll_comp64, dimension(:)             :: array_out
    sll_int32, optional, intent(in)      :: flags
    type(sll_fft_plan), pointer              :: plan

#if defined(_NOFFTW) && _DEFAULTFFTLIB==FFTW_MOD
    print * ,'Error in file sll_fft.F90'
    print * ,'      function fft_new_r2c_1d_for_1d'
    stop     '      The default library cannot be FFTW because she is not installed'
#endif

    if(_DEFAULTFFTLIB .eq. FFTPACK_MOD) then
      print * ,'Error in file sll_fft.F90'
      print * ,'      function fft_new_r2c_1d_for_1d'
      stop     '      This function doesn''t work with FFTAPCK library'
    endif

    if(size(array_in) .ne. nx) then
      print * ,'Error in file sll_fft.F90'
      print * ,'      function fft_new_r2c_1d_for_1d'
      stop     '      array_in size problem'
    else if( size(array_out) .ne. (nx/2 + 1) ) then
      print * ,'Error in file sll_fft.F90'
      print * ,'      function fft_new_r2c_1d_for_1d'
      stop     '      array_out size problem'
    endif

    plan => fft_plan_r2c_1d(_DEFAULTFFTLIB,nx,array_in(1:nx),array_out(1:nx/2+1),flags)
  end function

 function new_plan_r2c_2d(nx,ny,array_in,array_out,flags) result(plan)
    sll_int32, intent(in)                         :: nx,ny
    sll_real64, dimension(0:,0:), intent(inout)   :: array_in
    sll_comp64, dimension(0:,0:), intent(inout)   :: array_out
    sll_int32, optional                           :: flags
    type(sll_fft_plan), pointer                   :: plan

    if( (size(array_in,dim=1).ne.nx) .and. (size(array_in,dim=2).ne.ny) ) then
      print * ,'Error in file sll_fft.F90'
      print * ,'      in function new_plan_r2c_2d'
      print * ,'      array_in size problem'
      stop ''
    else if( size(array_in,dim=1).ne.nx/2+1 .and. size(array_in,dim=2).ne.ny ) then
      print * ,'Error in file sll_fft.F90'
      print * ,'      in function new_plan_r2c_2d'
      print * ,'      array_out size problem'
      stop ''
    endif

    plan => fft_plan_r2c_2d(_DEFAULTFFTLIB,nx,ny,array_in,array_out,flags)
  end function

! ---------------------
! -  DEBUG INTERFACE  -
! ---------------------
  function fft_plan_r2c_1d(library,nx,array_in,array_out,flags) result(plan)
    sll_int32, intent(in)                            :: library
    sll_int32, intent(in)                            :: nx
    sll_real64, dimension(0:), intent(inout)         :: array_in
    sll_comp64, dimension(0:), intent(out)           :: array_out
    sll_int32, optional, intent(in)                  :: flags
    type(sll_fft_plan), pointer                          :: plan

#ifndef _NOFFTW 
    if(library .eq. FFTW_MOD) then
      plan => fftw_new_plan_r2c_1d(library,nx,array_in,array_out,flags)
      return
    endif
#endif
#ifndef _NOFFTSLL
    if(library .eq. SLLFFT_MOD) then
      plan => sll_new_plan_r2c_1d(library,nx,array_in,array_out,flags)
      return
    endif
#endif
#ifndef _NOFFTPACK
    if(library .eq. FFTPACK_MOD) then
      stop 'not available in fftpack'
    endif
#endif

    stop 'ERROR in =fft_plan_r2c_1d= library unknown'
  end function

  function fft_plan_r2c_2d(library,nx,ny,array_in,array_out,flags) result(plan)
    sll_int32, intent(in)                          :: library
    sll_int32, intent(in)                          :: nx, ny
    sll_real64, dimension(0:,0:)                   :: array_in
    sll_comp64, dimension(0:,0:)                   :: array_out
    sll_int32, optional, intent(in)                :: flags
    type(sll_fft_plan), pointer                        :: plan
    sll_int32                                      :: ierr

#ifndef _NOFFTW 
    if(library .eq. FFTW_MOD) then
      SLL_ALLOCATE(plan,ierr)
      plan%library = FFTW_MOD 
      plan%direction = FFT_FORWARD
      if( present(flags) )then
        plan%style = flags
      else
        plan%style = 0_f32
      endif
      plan%problem_rank = 2
      SLL_ALLOCATE(plan%problem_shape(2),ierr)
      plan%problem_shape = (/ nx , ny /)
      plan%fftw_plan = fftw_plan_dft_r2c_2d(ny,nx,array_in,array_out,FFTW_ESTIMATE + FFTW_UNALIGNED)
      return
    endif
#endif
#ifndef _NOFFTSLL
    if(library .eq. SLLFFT_MOD) then
      plan => sll_new_plan_r2c_2d(library,nx,ny,array_in,array_out,flags)
      return
    endif
#endif
#ifndef _NOFFTPACK
    !if(library .eq. FFTPACK_MOD) then
    !  plan => fftpack_new_plan_r2r_1d(library,size_problem,array_in,array_out,direction,flags)
    !endif
#endif

    stop 'ERROR in =fft_plan_r2c_2d= library unknown'
  end function

! ---------------------
! -       APPLY       -
! ---------------------
  subroutine fft_apply_r2c_1d(plan,array_in,array_out)
    type(sll_fft_plan), pointer                     :: plan
    sll_real64, dimension(:), intent(inout)        :: array_in
    sll_comp64, dimension(:), intent(out)          :: array_out
    sll_int32                                       :: nx
    sll_real64                                      :: time, factor
    type(time_mark), pointer                        :: mark

    if( .not. associated(plan) ) then
      print*,'Error in the sll_fft.F90'
      print*,'      in subroutine fft_apply_r2c_1d'
      print*,'      plan not associated'
      stop ''
    endif

    nx = plan%problem_shape(1)

    if(size(array_in) .ne. nx) then
      print * ,'Error in file sll_fft.F90'
      print * ,'      function fft_new_r2c_1d_for_1d'
      stop     '      array_in size problem'
    else if( size(array_out) .ne. (nx/2 + 1) ) then
      print * ,'Error in file sll_fft.F90'
      print * ,'      function fft_new_r2c_1d_for_1d'
      stop     '      array_out size problem'
    endif

    mark => new_time_mark()

#ifndef _NOFFTW
    if(plan%library .eq. FFTW_MOD) then
      mark => start_time_mark(mark)
      call fftw_apply_fft_r2c_1d(plan, array_in(1:nx), array_out(1:nx/2+1))
      time = time_elapsed_since(mark)
    endif
#endif
#ifndef _NOFFTSLL
    if(plan%library .eq. SLLFFT_MOD) then
      mark => start_time_mark(mark)
      call sll_fft_apply_r2c_1d(plan, array_in(1:nx), array_out(1:nx/2+1))
      time = time_elapsed_since(mark)
    endif
#endif
#ifndef _NOFFTPACK
!    if(plan%library .eq. FFTPACK_MOD) then
!      call fftpack_apply_fft_complex(plan, array_in, array_out)
!    endif
#endif

   if( fft_is_present_flag(plan,FFT_NORMALIZE) ) then
     factor = 1.0_f64/real(nx,kind=f64)
     array_out = array_out*factor
   endif
#ifdef _FFTINFO
INFO
#endif
    mark => delete_time_mark(mark)
  end subroutine

  subroutine fft_apply_r2c_2d(plan,array_in,array_out)
    type(sll_fft_plan), pointer                       :: plan
    sll_real64, dimension(:,:), intent(inout)         :: array_in
    sll_comp64, dimension(:,:), intent(inout)         :: array_out
    sll_real64                                        :: factor
    sll_int32 :: nx, ny

    if( .not. associated(plan) ) then
      print*,'Error in the sll_fft.F90'
      print*,'      in subroutine fft_apply_r2c_2d'
      print*,'      plan not associated'
      stop ''
    endif

    nx = plan%problem_shape(1)
    ny = plan%problem_shape(2)

    if( (size(array_in,dim=1).ne.nx) .and. (size(array_in,dim=2).ne.ny) ) then
      print * ,'Error in file sll_fft.F90'
      print * ,'      in subroutine fft_apply_r2c_2d'
      print * ,'      array_in size problem'
      stop ''
    else if( size(array_in,dim=1).ne.nx/2+1 .and. size(array_in,dim=2).ne.ny ) then
      print * ,'Error in file sll_fft.F90'
      print * ,'      in subroutine fft_apply_r2c_2d'
      print * ,'      array_out size problem'
      stop ''
    endif

#ifndef _NOFFTW
    if(plan%library .eq. FFTW_MOD) then
      call fftw_execute_dft_r2c(plan%fftw_plan, array_in, array_out )
    endif
#endif
#ifndef _NOFFTSLL
    if(plan%library .eq. SLLFFT_MOD) then
      call sll_fft_apply_r2c_2d(plan,array_in,array_out)
    endif
#endif
#ifndef _NOFFTPACK
!    else if(plan%library .eq. FFTPACK_MOD) &
!      call fftpack_apply_fft_complex(plan, array_in, array_out)
#endif
  if( fft_is_present_flag(plan,FFT_NORMALIZE) ) then
    factor = 1.0_f64/real(nx*ny,kind=f64)
    array_out = factor*array_out
  endif
  end subroutine








                                      !!!!!!!!!!!!!!!!!!!!!!!
                                      !!  COMPLEX TO REAL  !!
                                      !!!!!!!!!!!!!!!!!!!!!!!

! ---------------------
! - DEFAULT INTERFACE -
! ---------------------
 function new_plan_c2r_1d(size_problem,array_in,array_out,flags) result(plan)
    sll_int32, intent(in)                    :: size_problem
    sll_comp64, dimension(0:), intent(in)    :: array_in
    sll_real64, dimension(0:), intent(in)    :: array_out
    sll_int32, optional, intent(in)          :: flags
    type(sll_fft_plan), pointer                  :: plan

#if defined(_NOFFTW) && _DEFAULTFFTLIB==FFTW_MOD
    print * ,'Error in file sll_fft.F90'
    print * ,'      function fft_new_c2r_1d_for_1d'
    stop     '      The default library cannot be FFTW because she is not installed'
#endif

    if(_DEFAULTFFTLIB .eq. FFTPACK_MOD) then
      print * ,'Error in file sll_fft.F90'
      print * ,'      function fft_new_c2r_1d_for_1d'
      stop     '      This function doesn''t work with FFTAPCK library'
    endif

    if(size(array_out) .ne. size_problem) then
      print * ,'Error in file sll_fft.F90'
      print * ,'      function fft_new_c2r_1d_for_1d'
      stop     '      array_out size problem'
    else if( size(array_in) .ne. (size_problem/2 + 1) ) then
      print * ,'Error in file sll_fft.F90'
      print * ,'      function fft_new_c2r_1d_for_1d'
      stop     '      array_in size problem'
    endif

    plan => fft_plan_c2r_1d(_DEFAULTFFTLIB,size_problem,array_in(0:size_problem/2),array_out(0:size_problem-1),flags)
  end function

 function new_plan_c2r_2d(nx,ny,array_in,array_out,flags) result(plan)
    sll_int32, intent(in)                         :: nx,ny
    sll_comp64, dimension(0:,0:), intent(inout)      :: array_in
    sll_real64, dimension(0:,0:), intent(inout)      :: array_out
    sll_int32, optional, intent(in)               :: flags
    type(sll_fft_plan), pointer                       :: plan

    if( (size(array_in,dim=1).ne.nx/2+1) .and. (size(array_in,dim=2).ne.ny) ) then
      print * ,'Error in file sll_fft.F90'
      print * ,'      in function new_plan_c2r_2d'
      print * ,'      array_in size problem'
      stop ''
    else if( size(array_in,dim=1).ne.nx .and. size(array_in,dim=2).ne.ny ) then
      print * ,'Error in file sll_fft.F90'
      print * ,'      in function new_plan_c2r_2d'
      print * ,'      array_out size problem'
      stop ''
    endif

    plan => fft_plan_c2r_2d(_DEFAULTFFTLIB,nx,ny,array_in,array_out,flags)
  end function
! ---------------------
! -  DEBUG INTERFACE  -
! ---------------------
  function fft_plan_c2r_1d(library,size_problem,array_in,array_out,flags) result(plan)
    sll_int32, intent(in)                 :: library
    sll_int32, intent(in)                 :: size_problem
    sll_comp64, dimension(:), intent(in)  :: array_in
    sll_real64, dimension(:), intent(in)  :: array_out
    sll_int32, optional, intent(in)       :: flags
    type(sll_fft_plan), pointer               :: plan

#ifndef _NOFFTW 
    if(library .eq. FFTW_MOD) then
      plan => fftw_new_plan_c2r_1d(library,size_problem,array_in,array_out,flags)
      return
    endif
#endif
#ifndef _NOFFTSLL
    if(library .eq. SLLFFT_MOD) then
      plan => sll_new_plan_c2r_1d(library,size_problem,array_in,array_out,flags)
      return
    endif
#endif
#ifndef _NOFFTPACK
!    if(library .eq. FFTPACK_MOD) then
!      stop 'c2r transform not available with fftpack'    
!    endif
#endif
    stop 'ERROR in =fft_plan_c2r_1d= library unknown'
  end function

  function fft_plan_c2r_2d(library,nx,ny,array_in,array_out,flags) result(plan)
    sll_int32, intent(in)                               :: library
    sll_int32, intent(in)                               :: nx, ny
    sll_comp64, dimension(0:,0:), intent(inout)         :: array_in
    sll_real64, dimension(0:,0:), intent(inout)         :: array_out
    sll_int32, optional, intent(in)                     :: flags
    type(sll_fft_plan), pointer                         :: plan 
    sll_int32                                           :: ierr

#ifndef _NOFFTW 
    if(library .eq. FFTW_MOD) then
      SLL_ALLOCATE(plan,ierr)
      plan%library = FFTW_MOD 
      plan%direction = FFT_INVERSE
      if( present(flags) )then
        plan%style = flags
      else
        plan%style = 0_f32
      endif
      plan%problem_rank = 2
      SLL_ALLOCATE(plan%problem_shape(2),ierr)
      plan%problem_shape = (/ nx , ny /)
      plan%fftw_plan = fftw_plan_dft_c2r_2d(ny,nx,array_in,array_out,FFTW_ESTIMATE + FFTW_UNALIGNED)
      return
    endif
#endif
#ifndef _NOFFTSLL
    if(library .eq. SLLFFT_MOD) then
      plan => sll_new_plan_c2r_2d(library,nx,ny,array_in,array_out,flags)
      return
    endif
#endif
#ifndef _NOFFTPACK
    !if(library .eq. FFTPACK_MOD) then
    !  plan => fftpack_new_plan_r2r_1d(library,size_problem,array_in,array_out,direction,flags)
    !endif
#endif

    stop 'ERROR in =fft_plan_c2r_2d= library unknown'
  end function

! ---------------------
! -       APPLY       -
! ---------------------
  subroutine fft_apply_c2r_1d(plan,array_in,array_out)
    type(sll_fft_plan), pointer                 :: plan
    sll_comp64, dimension(0:), intent(inout)    :: array_in
    sll_real64, dimension(0:), intent(inout)    :: array_out
    sll_int32                                   ::  nx
    sll_real64                                  :: factor, time
    type(time_mark), pointer                    :: mark

    nx = plan%problem_shape(1)   
 
    if(size(array_out) .ne. nx) then
      print * ,'Error in file sll_fft.F90'
      print * ,'      function fft_apply_c2r_1d'
      stop     '      array_out size problem'
    else if( size(array_in) .ne. (nx/2 + 1) ) then
      print * ,'Error in file sll_fft.F90'
      print * ,'      function fft_apply_c2r_1d'
      stop     '      array_in size problem'
    endif

    mark => new_time_mark()

#ifndef _NOFFTW
    if(plan%library .eq. FFTW_MOD) then
      mark => start_time_mark(mark)
      call fftw_apply_fft_c2r_1d(plan, array_in(0:nx/2), array_out(0:nx-1))
      time = time_elapsed_since(mark)
    endif
#endif
#ifndef _NOFFTSLL
    if(plan%library .eq. SLLFFT_MOD) then
      mark => start_time_mark(mark)
      call sll_fft_apply_c2r_1d(plan,array_in,array_out)
      time = time_elapsed_since(mark)
    endif
#endif
#ifndef _NOFFTPACK
    if(plan%library .eq. FFTPACK_MOD) then
      stop 'c2r transform don''t work with fftpack'
    endif
#endif

   if( fft_is_present_flag(plan,FFT_NORMALIZE) ) then
     factor = 1.0_f64/real(nx,kind=f64)
     array_out = factor*array_out
   endif

#ifdef _FFTINFO
INFO
#endif
  mark => delete_time_mark(mark)
  end subroutine


  subroutine fft_apply_c2r_2d(plan,array_in,array_out)
    type(sll_fft_plan), pointer                       :: plan
    sll_comp64, dimension(:,:), intent(inout)         :: array_in
    sll_real64, dimension(:,:), intent(inout)         :: array_out
    sll_real64                                        :: factor
    sll_int32                                         :: nx, ny

    if( .not. associated(plan) ) then
      print*,'Error in the sll_fft.F90'
      print*,'      in subroutine fft_apply_c2r_2d'
      print*,'      plan not associated'
      stop ''
    endif

    nx = plan%problem_shape(1)
    ny = plan%problem_shape(2)

#ifndef _NOFFTW
    if(plan%library .eq. FFTW_MOD) then
      call fftw_execute_dft_c2r(plan%fftw_plan, array_in(1:nx/2+1,1:ny), array_out(1:nx,1:ny) )
    endif
#endif
#ifndef _NOFFTSLL
    if(plan%library .eq. SLLFFT_MOD) then
      call sll_fft_apply_c2r_2d(plan,array_in,array_out)
    endif
#endif
!#ifndef _NOFFTPACK
!    else if(plan%library .eq. FFTPACK_MOD) &
!      call fftpack_apply_fft_complex(plan, array_in, array_out)
!#endif
  
  if( fft_is_present_flag(plan,FFT_NORMALIZE) ) then
    factor = 1.0_f64/real(nx*ny,kind=f64)
    array_out = factor*array_out    
  endif
  end subroutine





!  function get_mode(input,x,y) result(mode)
!    !type(sll_fft_plan) :: plan
!    sll_comp64, dimension(0:,0:), target :: input
!    sll_int32 :: x, y
!    sll_comp64, pointer :: mode
!
!    mode => input(x,y)
!  end function

  function which_mode_is_stored_in(x,y) result(k)
    !type(sll_fft_plan) :: plan
    sll_int32 :: x, y
    sll_int32, dimension(2) :: k

    k = (/ x , y /)
  end function

  function where_is_stored_mode(x,y) result(k)
    sll_int32 :: x, y
    sll_int32, dimension(2) :: k

    k = (/ x , y /)
  end function




  !subroutine CFFTF(x)
  ! sll_comp64, dimension(:) :: x
  ! sll_int32 :: n
  ! type(sll_fft_plan), pointer :: plan
  !
  ! n = size(x)
  ! plan => fft_plan_c2c_1d(FFTW_MOD,n,x,x,FFT_FORWARD)
  ! call fft_apply_c2c_1d(plan,x,x)
  ! call fft_delete_plan(plan)
  !end subroutine







            !!!!!!!!!!!!!!!!!!!!!!!
            !! DELETE SUBROUTINE !!
            !!!!!!!!!!!!!!!!!!!!!!!



  subroutine fft_delete(plan)
   type(sll_fft_plan), pointer :: plan
   sll_int32 :: ierr

    if( .not. associated(plan) ) then
      print * , 'Error in file sll_fft.F90'
      print * , '      subroutine fft_delete'
      print * , '      plan is not associated'
      stop 
    endif

#ifndef _NOFFTW
     if(plan%library .eq. FFTW_MOD) then
      call fftw_destroy_plan(plan%fftw_plan)
     endif
#endif
    if(plan%library .eq. SLLFFT_MOD) then
      if(associated(plan%t)) then
        deallocate(plan%t)
        plan%t => null()
      endif
      if(associated(plan%twiddles)) then
        deallocate(plan%twiddles)
        plan%twiddles => null()
      endif
      if(associated(plan%twiddles_n)) then
        deallocate(plan%twiddles_n)
        plan%twiddles_n => null()
      endif
      if(associated(plan%problem_shape)) then
        SLL_DEALLOCATE(plan%problem_shape,ierr)
      endif
    endif
    if(plan%library .eq. FFTPACK_MOD) then
      if(associated(plan%twiddles)) then
        deallocate(plan%twiddles)
        plan%twiddles => null()
      endif
    endif
    plan => null()
  end subroutine










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

  !Conjugate complex array that is represented by an array of reals.
  subroutine conjg_in_pairs(n,array)
    sll_real64, dimension(1:) :: array
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
    integer, intent(in)                             :: size
    sll_comp64, dimension(0:size-1), intent(inout)  :: dat
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
    !if ( sign == FFT_INVERSE ) then
    !   omega = twiddles(twiddle_index)
    !else if ( sign == FFT_FORWARD ) then
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
       !if ( sign == FFT_FORWARD ) then
       !   omega = conjg(twiddles(jtwiddle))
       !else if ( sign == FFT_INVERSE ) then
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
    !if ( sign == FFT_FORWARD ) then
    !   omega_re =  CREAL0(twiddles, twiddle_index)
    !   omega_im = -CIMAG0(twiddles, twiddle_index)
    !else if ( sign == FFT_INVERSE ) then
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
       !if( sign == FFT_FORWARD ) then
       !   omega_re =  CREAL1(twiddles, jtwiddle)
       !   omega_im = -CIMAG1(twiddles, jtwiddle)
       !else if ( sign == FFT_INVERSE ) then
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
       data = data*2.0_f64
    end if
  end subroutine real_data_fft_dit
  
#undef CREAL0
#undef CIMAG0
#undef CREAL1
#undef CIMAG1


end module sll_fft
