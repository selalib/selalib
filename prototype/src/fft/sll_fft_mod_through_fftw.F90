module sll_fft
#include "sll_working_precision.h"
#include "misc_utils.h"
#include "sll_assert.h"
#include "sll_memory.h"
  use fftw_module
  use fft_utils
  implicit none 



  type sll_fft_plan
    type(C_PTR)                     :: fftw_plan

    sll_int32                        :: style
    sll_int32                        :: library
    sll_int32                        :: direction
    sll_int32                        :: problem_rank
    sll_int32, dimension(:), pointer :: problem_shape
  end type sll_fft_plan

  interface fft_new_plan
    module procedure fftw_new_plan_c2c_1d, fftw_new_plan_c2c_2d, &
       fftw_new_plan_r2r_1d, &
       fftw_new_plan_r2c_1d, fftw_new_plan_c2r_1d, &
       fftw_new_plan_r2c_2d, fftw_new_plan_c2r_2d 
  end interface
  interface fft_apply_plan
    module procedure fftw_apply_plan_c2c_1d, fftw_apply_plan_c2c_2d, &
       fftw_apply_plan_r2r_1d, &
       fftw_apply_plan_r2c_1d, fftw_apply_plan_c2r_1d, &
       fftw_apply_plan_r2c_2d, fftw_apply_plan_c2r_2d
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

  integer, parameter :: FFTW_MOD = 1000000000

  interface fft_get_mode
     module procedure fft_get_mode_complx_1d, fft_get_mode_complx_2d, &
                      fft_get_mode_complx_3d, fft_get_mode_real_1d
  end interface

  interface fft_set_mode
     module procedure fft_set_mode_complx_1d, fft_set_mode_complx_2d, &
                      fft_set_mode_complx_3d, fft_set_mode_real_1d
  end interface

contains



  subroutine print_defaultfftlib()
    print *, 'The library used is FFTW'
  end subroutine


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

    n = plan%problem_shape(1)
    n_2 = n/2 !ishft(n,-1)

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
  end subroutine 




! COMPLEX
  function fftw_new_plan_c2c_1d(nx,array_in,array_out,direction,flags) result(plan)
    sll_int32, intent(in)                        :: nx
    sll_comp64, dimension(:), intent(inout)      :: array_in
    sll_comp64, dimension(:), intent(inout)      :: array_out
    sll_int32, intent(in)                        :: direction
    sll_int32, optional, intent(in)              :: flags
    type(sll_fft_plan), pointer                  :: plan
    sll_int32                                    :: ierr

    SLL_ALLOCATE(plan,ierr)
    plan%library = FFTW_MOD 
    plan%direction = direction
    if( present(flags) )then
      plan%style = flags
    else
      plan%style = 0_f32
    endif
    plan%problem_rank = 1
    SLL_ALLOCATE(plan%problem_shape(1),ierr)
    plan%problem_shape = (/ nx  /)

    plan%fftw_plan = fftw_plan_dft_1d(nx,array_in,array_out,direction,FFTW_ESTIMATE + FFTW_UNALIGNED)
  end function

  subroutine fftw_apply_plan_c2c_1d(plan,array_in,array_out)
    type(sll_fft_plan), pointer, intent(in)         :: plan
    sll_comp64, dimension(:), intent(inout)         :: array_in, array_out
    sll_real64 :: factor

    call fftw_execute_dft(plan%fftw_plan, array_in, array_out)

    if( fft_is_present_flag(plan%style,FFT_NORMALIZE) ) then
      factor = 1.0_f64/real(plan%problem_shape(1),kind=f64)
      array_out = factor*array_out
    endif
  end subroutine 
! END COMPLEX

! COMPLEX 2D
  function fftw_new_plan_c2c_2d(NX,NY,array_in,array_out,direction,flags) result(plan)
    sll_int32, intent(in)                            :: NX,NY
    sll_comp64, dimension(0:,0:)                     :: array_in, array_out
    sll_int32, intent(in)                            :: direction
    sll_int32, optional, intent(in)                  :: flags
    type(sll_fft_plan), pointer                      :: plan
    sll_int32                                        :: ierr

    SLL_ALLOCATE(plan,ierr)
    plan%library = FFTW_MOD
    plan%direction = direction
    if( present(flags) )then
      plan%style = flags
    else
      plan%style = 0_f32
    endif
    plan%problem_rank = 2
    SLL_ALLOCATE(plan%problem_shape(2),ierr)
    plan%problem_shape = (/ NX , NY /)
  
    !We must switch the dimension. It's a fftw convention. 
    plan%fftw_plan = fftw_plan_dft_2d(NY,NX,array_in,array_out,direction,FFTW_ESTIMATE + FFTW_UNALIGNED)
  end function

  subroutine fftw_apply_plan_c2c_2d(plan,array_in,array_out)
    type(sll_fft_plan), pointer, intent(in)      :: plan
    sll_comp64, dimension(0:,0:), intent(inout)  :: array_in, array_out
    sll_real64                                   :: factor

    call fftw_execute_dft(plan%fftw_plan, array_in, array_out)

    if( fft_is_present_flag(plan%style,FFT_NORMALIZE) ) then
      factor = 1.0_f64/real(plan%problem_shape(1)*plan%problem_shape(2),kind=f64)
      array_out = factor*array_out
    endif
  end subroutine 
! END COMPLEX 2D

! REAL
  function fftw_new_plan_r2r_1d(nx,array_in,array_out,direction,flags) result(plan)
    sll_int32, intent(in)                        :: nx
    sll_real64, dimension(:), intent(inout)      :: array_in
    sll_real64, dimension(:), intent(inout)      :: array_out
    sll_int32, intent(in)                        :: direction
    sll_int32, optional, intent(in)              :: flags
    type(sll_fft_plan), pointer                  :: plan
    sll_int32 :: ierr

    SLL_ALLOCATE(plan,ierr)
    plan%library = FFTW_MOD
    plan%direction = 0
    if( present(flags) )then
      plan%style = flags
    else
      plan%style = 0_f32
    endif
    plan%problem_rank = 1
    SLL_ALLOCATE(plan%problem_shape(1),ierr)
    plan%problem_shape = (/ nx /)

    if(direction .eq. FFT_FORWARD) then
      plan%fftw_plan = fftw_plan_r2r_1d(nx,array_in,array_out,FFTW_R2HC,FFTW_ESTIMATE + FFTW_UNALIGNED)
    else if(direction .eq. FFT_INVERSE) then
      plan%fftw_plan = fftw_plan_r2r_1d(nx,array_in,array_out,FFTW_HC2R,FFTW_ESTIMATE + FFTW_UNALIGNED)
    endif
  end function

  subroutine fftw_apply_plan_r2r_1d(plan,array_in,array_out)
    type(sll_fft_plan), pointer, intent(in)     :: plan
    sll_real64, dimension(:), intent(inout)     :: array_in, array_out
    sll_real64                                  :: factor

    call fftw_execute_r2r(plan%fftw_plan, array_in, array_out)

    if( fft_is_present_flag(plan%style,FFT_NORMALIZE) ) then
      factor = 1.0_f64/real(plan%problem_shape(1),kind=f64)
      array_out = factor*array_out
    endif
  end subroutine


! R2C
  function fftw_new_plan_r2c_1d(nx,array_in,array_out,flags) result(plan)
    sll_int32, intent(in)                   :: nx
    sll_real64, dimension(:), intent(inout) :: array_in
    sll_comp64, dimension(:), intent(out)   :: array_out
    sll_int32, optional, intent(in)         :: flags
    type(sll_fft_plan), pointer             :: plan
    sll_int32 :: ierr

    SLL_ALLOCATE(plan,ierr)
    plan%library = FFTW_MOD
    plan%direction = 0
    if( present(flags) )then
      plan%style = flags
    else
      plan%style = 0_f32
    endif
    plan%problem_rank = 1
    SLL_ALLOCATE(plan%problem_shape(1),ierr)
    plan%problem_shape = (/ nx /)

    plan%fftw_plan = fftw_plan_dft_r2c_1d(nx,array_in,array_out,FFTW_ESTIMATE + FFTW_UNALIGNED)
  end function

  subroutine fftw_apply_plan_r2c_1d(plan,array_in,array_out)
    type(sll_fft_plan), pointer, intent(in)         :: plan
    sll_real64, dimension(:), intent(inout)         :: array_in
    sll_comp64, dimension(:), intent(out)           :: array_out
    sll_real64                                      :: factor

    call fftw_execute_dft_r2c(plan%fftw_plan, array_in, array_out)

    if( fft_is_present_flag(plan%style,FFT_NORMALIZE) ) then
      factor = 1.0_f64/real(plan%problem_shape(1),kind=f64)
      array_out = factor*array_out
    endif
  end subroutine

  function fftw_new_plan_r2c_2d(nx,ny,array_in,array_out,flags) result(plan)
    sll_int32, intent(in)                     :: nx,ny
    sll_real64, dimension(:,:), intent(inout) :: array_in
    sll_comp64, dimension(:,:), intent(out)   :: array_out
    sll_int32, optional, intent(in)           :: flags
    type(sll_fft_plan), pointer               :: plan
    sll_int32 :: ierr

    SLL_ALLOCATE(plan,ierr)
    plan%library = FFTW_MOD
    plan%direction = 0
    if( present(flags) )then
      plan%style = flags
    else
      plan%style = 0_f32
    endif
      
    plan%problem_rank = 2
    SLL_ALLOCATE(plan%problem_shape(2),ierr)
    plan%problem_shape = (/ nx , ny /)
    plan%fftw_plan = fftw_plan_dft_r2c_2d(ny,nx,array_in,array_out,FFTW_ESTIMATE + FFTW_UNALIGNED)
  end function

  subroutine fftw_apply_plan_r2c_2d(plan,array_in,array_out)
    type(sll_fft_plan), pointer, intent(in)           :: plan
    sll_real64, dimension(:,:), intent(inout)         :: array_in
    sll_comp64, dimension(:,:), intent(out)           :: array_out
    sll_int32     :: nx, ny
    sll_real64    :: factor

    if( .not. associated(plan) ) then
      print*,'Eroor in subroutine fftw_apply_plan_r2c_2d'
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

    call fftw_execute_dft_r2c(plan%fftw_plan, array_in, array_out)

    if( fft_is_present_flag(plan%style,FFT_NORMALIZE) ) then
      factor = 1.0_f64/real(nx*ny,kind=f64)
      array_out = factor*array_out
    endif
  end subroutine

!END R2C

! C2R
  function fftw_new_plan_c2r_1d(nx,array_in,array_out,flags) result(plan)
    sll_int32, intent(in)                  :: nx
    sll_comp64, dimension(:)               :: array_in
    sll_real64, dimension(:)               :: array_out
    sll_int32, optional, intent(in)        :: flags
    type(sll_fft_plan), pointer            :: plan
    sll_int32 :: ierr

    SLL_ALLOCATE(plan,ierr)
    plan%library = FFTW_MOD
    plan%direction = 0
    if( present(flags) )then
      plan%style = flags
    else
      plan%style = 0_f32
    endif
    plan%problem_rank = 1
    SLL_ALLOCATE(plan%problem_shape(1),ierr)
    plan%problem_shape = (/ nx /)

    plan%fftw_plan = fftw_plan_dft_c2r_1d(nx,array_in,array_out,FFTW_ESTIMATE + FFTW_UNALIGNED)
  end function

  subroutine fftw_apply_plan_c2r_1d(plan,array_in,array_out)
    type(sll_fft_plan), pointer :: plan
    sll_comp64, dimension(:)    :: array_in
    sll_real64, dimension(:)    :: array_out
    sll_real64                  :: factor

    call fftw_execute_dft_c2r(plan%fftw_plan, array_in, array_out)

    if( fft_is_present_flag(plan%style,FFT_NORMALIZE) ) then
      factor = 1.0_f64/real(plan%problem_shape(1),kind=f64)
      array_out = factor*array_out
    endif
  end subroutine

  function fftw_new_plan_c2r_2d(nx,ny,array_in,array_out,flags) result(plan)
    sll_int32, intent(in)                     :: nx, ny
    sll_comp64, dimension(:,:), intent(inout) :: array_in
    sll_real64, dimension(:,:), intent(out)   :: array_out
    sll_int32, optional, intent(in)           :: flags
    type(sll_fft_plan), pointer               :: plan
    sll_int32 :: ierr

    SLL_ALLOCATE(plan,ierr)
    plan%library = FFTW_MOD
    plan%direction = 0
    if( present(flags) )then
      plan%style = flags
    else
      plan%style = 0_f32
    endif
      
    plan%problem_rank = 2
    SLL_ALLOCATE(plan%problem_shape(2),ierr)
    plan%problem_shape = (/ nx , ny /)
    plan%fftw_plan = fftw_plan_dft_c2r_2d(ny,nx,array_in,array_out,FFTW_ESTIMATE + FFTW_UNALIGNED)
  end function

  subroutine fftw_apply_plan_c2r_2d(plan,array_in,array_out)
    type(sll_fft_plan), pointer, intent(in)           :: plan
    sll_comp64, dimension(1:,1:), intent(inout)       :: array_in
    sll_real64, dimension(1:,1:), intent(out)         :: array_out
    sll_int32                                         :: nx, ny
    sll_real64                                        :: factor

    if( .not. associated(plan) ) then
      print*,'Error in subroutine fftw_apply_plan_c2r_2d'
      print*,'      plan not associated'
      stop ''
    endif

    nx = plan%problem_shape(1)
    ny = plan%problem_shape(2)
    call fftw_execute_dft_c2r(plan%fftw_plan, array_in(1:nx/2+1,1:ny), array_out(1:nx,1:ny) )

    if( fft_is_present_flag(plan%style,FFT_NORMALIZE) ) then
      factor = 1.0_f64/real(nx*ny,kind=f64)
      array_out = factor*array_out
    endif
  end subroutine
!END C2R
! FFTW FFTW FFTW FFTW FFTW FFTW FFTW FFTW FFTW FFTW FFTW FFTW FFTW FFTW
! FFTW FFTW FFTW FFTW FFTW FFTW FFTW FFTW FFTW FFTW FFTW FFTW FFTW FFTW
! FFTW FFTW FFTW FFTW FFTW FFTW FFTW FFTW FFTW FFTW FFTW FFTW FFTW FFTW


  subroutine fft_delete_plan(plan)
   type(sll_fft_plan), pointer :: plan
   sll_int32 :: ierr

    if( .not. associated(plan) ) then
      print * , '  Error in fft_delete_plan subroutine'
      print * , '  you try to delete a plan not associated'
      stop 
    endif

    call fftw_destroy_plan(plan%fftw_plan)
    if(associated(plan%problem_shape)) then
      SLL_DEALLOCATE(plan%problem_shape,ierr)
    endif
    plan => null()
  end subroutine

end module sll_fft
