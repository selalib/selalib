! FFTPACK FFTPACK FFTPACK FFTPACK FFTPACK FFTPACK FFTPACK FFTPACK FFTPACK
! FFTPACK FFTPACK FFTPACK FFTPACK FFTPACK FFTPACK FFTPACK FFTPACK FFTPACK
! FFTPACK FFTPACK FFTPACK FFTPACK FFTPACK FFTPACK FFTPACK FFTPACK FFTPACK
! COMPLEX
  function fftpack_new_plan_c2c_1d(library,size_problem,array_in,array_out,direction,flags) result(plan)
    sll_int32, intent(in)                        :: library
    sll_int32, intent(in)                        :: size_problem
    sll_comp64, dimension(:), target, intent(in) :: array_in, array_out
    sll_int32, intent(in)                        :: direction
    sll_int32, optional,  intent(in)             :: flags
    type(sll_fft_plan), pointer                  :: plan
 
    allocate(plan)
    
    plan%N = size_problem
    if( present(flags)) then
      plan%style = flags
    else
      plan%style = 0
    endif
    plan%library = library
    plan%direction = direction
    plan%in_comp => array_in
    plan%out_comp => array_out
    plan%in_real => null()
    plan%out_real => null()

    allocate(plan%dwsave(4*plan%N + 15))
    call zffti(plan%N,plan%dwsave)
  end function

  subroutine fftpack_apply_fft_complex(plan,array_in,array_out)
    type(sll_fft_plan), pointer, intent(in)         :: plan
    sll_comp64, dimension(:), intent(inout)         :: array_in, array_out

    if( loc(array_in) .ne. loc(array_out)) then ! out-place transform
       array_out = array_in
    endif   

    if( plan%direction .eq. FFT_FORWARD ) then
      call zfftf( plan%N , array_out ,plan%dwsave )
    else if( plan%direction .eq. FFT_INVERSE ) then
      call zfftb( plan%N, array_out , plan%dwsave )
    endif
  end subroutine
! END COMPLEX

! REAL
  function fftpack_new_plan_r2r_1d(library,size_problem,array_in,array_out,direction,flags) result(plan)
    sll_int32, intent(in)                        :: library
    sll_int32, intent(in)                        :: size_problem
    sll_real64, dimension(:), target, intent(in) :: array_in, array_out
    sll_int32, intent(in)                        :: direction
    sll_int32, optional,  intent(in)             :: flags
    type(sll_fft_plan), pointer                  :: plan
 
    allocate(plan)
    
    plan%N = size_problem
    if( present(flags)) then
      plan%style = flags
    else
      plan%style = 0
    endif
    plan%library = library
    plan%direction = direction
    plan%in_real => array_in
    plan%out_real => array_out
    plan%in_comp => null()
    plan%out_comp => null()

    allocate(plan%dwsave(2*plan%N + 15))
    call dffti(plan%N,plan%dwsave)
  end function

  subroutine fftpack_apply_fft_real(plan,array_in,array_out)
    type(sll_fft_plan), pointer, intent(in)         :: plan
    sll_real64, dimension(:), intent(inout)         :: array_in, array_out

    if( loc(array_in) .ne. loc(array_out)) then ! out-place transform
       array_out = array_in
    endif   

    if( plan%direction .eq. FFT_FORWARD ) then
      call dfftf( plan%N , array_out ,plan%dwsave )
    else if( plan%direction .eq. FFT_INVERSE ) then
      call dfftb( plan%N, array_out , plan%dwsave )
    endif
  end subroutine
! END REAL
! FFTPACK FFTPACK FFTPACK FFTPACK FFTPACK FFTPACK FFTPACK FFTPACK FFTPACK
! FFTPACK FFTPACK FFTPACK FFTPACK FFTPACK FFTPACK FFTPACK FFTPACK FFTPACK
! FFTPACK FFTPACK FFTPACK FFTPACK FFTPACK FFTPACK FFTPACK FFTPACK FFTPACK
