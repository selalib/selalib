! FFTW FFTW FFTW FFTW FFTW FFTW FFTW FFTW FFTW FFTW FFTW FFTW FFTW FFTW
! FFTW FFTW FFTW FFTW FFTW FFTW FFTW FFTW FFTW FFTW FFTW FFTW FFTW FFTW
! FFTW FFTW FFTW FFTW FFTW FFTW FFTW FFTW FFTW FFTW FFTW FFTW FFTW FFTW
! COMPLEX
  function fftw_new_plan_c2c_1d(library,size_problem,array_in,array_out,direction,flags) result(plan)
    sll_int32, intent(in)                        :: library
    sll_int32, intent(in)                        :: size_problem
    sll_comp64, dimension(:), target, intent(in) :: array_in, array_out
    sll_int32, intent(in)                        :: direction
    sll_int32, optional, intent(in)              :: flags
    type(sll_fft_plan), pointer                  :: plan

    allocate(plan)
    
    plan%N = size_problem
    if( present(flags)) then
      plan%style = flags + FFTW_ESTIMATE + FFTW_UNALIGNED
    else
      plan%style = FFTW_ESTIMATE + FFTW_UNALIGNED
    endif
    plan%library = library
    plan%direction = direction
    plan%in_comp => array_in
    plan%out_comp => array_out
    plan%in_real => null()
    plan%out_real => null()

    plan%fftw_plan = fftw_plan_dft_1d(plan%N,plan%in_comp,plan%out_comp,plan%direction,plan%style)
  end function

  subroutine fftw_apply_fft_complex(plan,array_in,array_out)
    type(sll_fft_plan), pointer, intent(in)         :: plan
    sll_comp64, dimension(:), intent(inout)         :: array_in, array_out

    call fftw_execute_dft(plan%fftw_plan, array_in, array_out)
  end subroutine fftw_apply_fft_complex

  subroutine fftw_apply_fft_comp(plan)
    type(sll_fft_plan), pointer, intent(in)         :: plan
    call fftw_execute_dft(plan%fftw_plan,plan%in_comp,plan%out_comp)
  end subroutine
! END COMPLEX

! REAL
  function fftw_new_plan_r2r_1d(library,size_problem,array_in,array_out,direction,flags) result(plan)
    sll_int32, intent(in)                        :: library
    sll_int32, intent(in)                        :: size_problem
    sll_real64, dimension(:), target, intent(in) :: array_in, array_out
    sll_int32, intent(in)                        :: direction
    sll_int32, optional, intent(in)              :: flags
    type(sll_fft_plan), pointer                  :: plan

    allocate(plan)
    
    plan%N = size_problem
    if( present(flags)) then
      plan%style = flags + FFTW_ESTIMATE + FFTW_UNALIGNED
    else
      plan%style = FFTW_ESTIMATE + FFTW_UNALIGNED
    endif
    plan%library = library
    plan%direction = direction
    plan%in_real => array_in
    plan%out_real => array_out
    plan%in_comp => null()
    plan%out_comp => null()

    if(plan%direction .eq. FFT_FORWARD) then
      plan%fftw_plan = fftw_plan_r2r_1d(plan%N,plan%in_real,plan%out_real,FFTW_HC2R,plan%style)
    else if(plan%direction .eq. FFT_INVERSE) then
      plan%fftw_plan = fftw_plan_r2r_1d(plan%N,plan%in_real,plan%out_real,FFTW_R2HC,plan%style)
    endif
  end function

  subroutine fftw_apply_fft_real(plan,array_in,array_out)
    type(sll_fft_plan), pointer, intent(in)         :: plan
    sll_real64, dimension(:), intent(inout)         :: array_in, array_out

    call fftw_execute_r2r(plan%fftw_plan, array_in, array_out)
  end subroutine

  subroutine fftw_apply_fft_re(plan)
    type(sll_fft_plan), pointer, intent(in)         :: plan
    call fftw_execute_r2r(plan%fftw_plan,plan%in_real,plan%out_real)
  end subroutine
! END REAL
! FFTW FFTW FFTW FFTW FFTW FFTW FFTW FFTW FFTW FFTW FFTW FFTW FFTW FFTW
! FFTW FFTW FFTW FFTW FFTW FFTW FFTW FFTW FFTW FFTW FFTW FFTW FFTW FFTW
! FFTW FFTW FFTW FFTW FFTW FFTW FFTW FFTW FFTW FFTW FFTW FFTW FFTW FFTW
