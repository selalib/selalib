! FFTPACK FFTPACK FFTPACK FFTPACK FFTPACK FFTPACK FFTPACK FFTPACK FFTPACK
! FFTPACK FFTPACK FFTPACK FFTPACK FFTPACK FFTPACK FFTPACK FFTPACK FFTPACK
! FFTPACK FFTPACK FFTPACK FFTPACK FFTPACK FFTPACK FFTPACK FFTPACK FFTPACK
! COMPLEX
  function fftpack_new_plan_c2c_1d(library,nx,array_in,array_out,direction,flags) result(plan)
    sll_int32, intent(in)                        :: library
    sll_int32, intent(in)                        :: nx
    sll_comp64, dimension(:)                     :: array_in, array_out
    sll_int32, intent(in)                        :: direction
    sll_int32, optional,  intent(in)             :: flags
    type(sll_fft_plan), pointer                      :: plan
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

  subroutine fftpack_apply_fft_complex(plan,array_in,array_out)
    type(sll_fft_plan), pointer, intent(in)             :: plan
    sll_comp64, dimension(:), intent(inout)         :: array_in, array_out
    sll_int32 :: nx

!    if( loc(array_in) .ne. loc(array_out)) then ! out-place transform
!       array_out = array_in
!    endif   
    nx = plan%problem_shape(1)
 
    if( plan%direction .eq. FFT_FORWARD ) then
      call zfftf( nx , array_out ,plan%twiddles )
    else if( plan%direction .eq. FFT_INVERSE ) then
      call zfftb( nx, array_out , plan%twiddles )
    endif
  end subroutine
! END COMPLEX

! REAL
  function fftpack_new_plan_r2r_1d(library,nx,array_in,array_out,direction,flags) result(plan)
    sll_int32, intent(in)                        :: library
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

  subroutine fftpack_apply_fft_real(plan,array_in,array_out)
    type(sll_fft_plan), pointer, intent(in)         :: plan
    sll_real64, dimension(:), intent(inout)         :: array_in, array_out
    sll_int32 :: nx

!    if( loc(array_in) .ne. loc(array_out)) then ! out-place transform
!       array_out = array_in
!    endif   

    nx = plan%problem_shape(1)

    if( plan%direction .eq. FFT_FORWARD ) then
      call dfftf( nx , array_out ,plan%twiddles )
    else if( plan%direction .eq. FFT_INVERSE ) then
      call dfftb( nx, array_out , plan%twiddles )
    endif
  end subroutine
! END REAL
! FFTPACK FFTPACK FFTPACK FFTPACK FFTPACK FFTPACK FFTPACK FFTPACK FFTPACK
! FFTPACK FFTPACK FFTPACK FFTPACK FFTPACK FFTPACK FFTPACK FFTPACK FFTPACK
! FFTPACK FFTPACK FFTPACK FFTPACK FFTPACK FFTPACK FFTPACK FFTPACK FFTPACK
