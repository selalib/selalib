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
      plan%style =  flags !+ FFTW_ESTIMATE + FFTW_UNALIGNED
    else
      plan%style = 0 !FFTW_ESTIMATE + FFTW_UNALIGNED
    endif
    plan%library = library
    plan%direction = direction
    plan%in_comp => array_in
    plan%out_comp => array_out
    plan%in_real => null()
    plan%out_real => null()

    plan%fftw_plan = fftw_plan_dft_1d(plan%N,plan%in_comp,plan%out_comp,plan%direction,FFTW_ESTIMATE + FFTW_UNALIGNED)
  end function

  subroutine fftw_apply_fft_complex(plan,array_in,array_out)
    type(sll_fft_plan), pointer, intent(in)         :: plan
    sll_comp64, dimension(:), intent(inout)         :: array_in, array_out

    call fftw_execute_dft(plan%fftw_plan, array_in, array_out)
  end subroutine 

  subroutine fftw_apply_fft_comp(plan)
    type(sll_fft_plan), pointer, intent(in)         :: plan
    call fftw_execute_dft(plan%fftw_plan,plan%in_comp,plan%out_comp)
  end subroutine
! END COMPLEX

! COMPLEX 2D
  function fftw_new_plan_c2c_2d(library,NX,NY,array_in,array_out,direction,flags) result(plan)
    sll_int32, intent(in)                            :: library
    sll_int32, intent(in)                            :: NX,NY
    sll_comp64, dimension(0:,0:), target             :: array_in, array_out
    sll_int32, intent(in)                            :: direction
    sll_int32, optional, intent(in)                  :: flags
    type(sll_fft_plan_2d), pointer                   :: plan

    allocate(plan)
    
    allocate(plan%plan_x)
    allocate(plan%plan_y)
    plan%plan_x%N = NX
    plan%plan_y%N = NY
    if( present(flags)) then
      plan%plan_x%style =  flags !+ FFTW_ESTIMATE + FFTW_UNALIGNED
    else
      plan%plan_x%style = 0 !FFTW_ESTIMATE + FFTW_UNALIGNED
    endif
    plan%library = library
    plan%plan_x%direction = direction
 
    plan%plan_x%fftw_plan = fftw_plan_dft_2d(NX,NY,array_in,array_out,direction,FFTW_ESTIMATE + FFTW_UNALIGNED)
  end function

  subroutine fftw_apply_plan_c2c_2d(plan,array_in,array_out)
    type(sll_fft_plan_2d), pointer, intent(in)      :: plan
    sll_comp64, dimension(0:,0:), intent(inout)     :: array_in, array_out

    call fftw_execute_dft(plan%plan_x%fftw_plan, array_in, array_out)
  end subroutine 
! END COMPLEX 2D

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
      plan%style = flags !+ FFTW_ESTIMATE + FFTW_UNALIGNED
    else
      plan%style = 0 !FFTW_ESTIMATE + FFTW_UNALIGNED
    endif
    plan%library = library
    plan%direction = direction
    plan%in_real => array_in
    plan%out_real => array_out
    plan%in_comp => null()
    plan%out_comp => null()

    if(plan%direction .eq. FFT_FORWARD) then
      plan%fftw_plan = fftw_plan_r2r_1d(plan%N,plan%in_real,plan%out_real,FFTW_R2HC,FFTW_ESTIMATE + FFTW_UNALIGNED)
    else if(plan%direction .eq. FFT_INVERSE) then
      plan%fftw_plan = fftw_plan_r2r_1d(plan%N,plan%in_real,plan%out_real,FFTW_HC2R,FFTW_ESTIMATE + FFTW_UNALIGNED)
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


! R2C
  function fftw_new_plan_r2c_1d(library,size_problem,array_in,array_out,flags) result(plan)
    sll_int32, intent(in)                           :: library
    sll_int32, intent(in)                           :: size_problem
    sll_real64, dimension(:), target, intent(inout) :: array_in
    sll_comp64, dimension(:), target, intent(out)   :: array_out
    sll_int32, optional, intent(in)                 :: flags
    type(sll_fft_plan), pointer                     :: plan

    allocate(plan)
    
    plan%N = size_problem
    if( present(flags)) then
      plan%style = flags !+ FFTW_ESTIMATE + FFTW_UNALIGNED
    else
      plan%style = 0 !FFTW_ESTIMATE + FFTW_UNALIGNED
    endif
    plan%library = library
    plan%in_real => array_in
    plan%out_real => null()
    plan%in_comp => null()
    plan%out_comp => array_out

    plan%fftw_plan = fftw_plan_dft_r2c_1d(size_problem,array_in,array_out,FFTW_ESTIMATE + FFTW_UNALIGNED)
  end function

  subroutine fftw_apply_fft_r2c_1d(plan,array_in,array_out)
    type(sll_fft_plan), pointer, intent(in)         :: plan
    sll_real64, dimension(:), intent(inout)         :: array_in
    sll_comp64, dimension(:), intent(out)           :: array_out

    call fftw_execute_dft_r2c(plan%fftw_plan, array_in, array_out)
  end subroutine
!END R2C

! C2R
  function fftw_new_plan_c2r_1d(library,size_problem,array_in,array_out,flags) result(plan)
    sll_int32, intent(in)                           :: library
    sll_int32, intent(in)                           :: size_problem
    sll_comp64, dimension(:), target                :: array_in
    sll_real64, dimension(:), target                :: array_out
    sll_int32, optional, intent(in)                 :: flags
    type(sll_fft_plan), pointer                     :: plan

    allocate(plan)
    
    plan%N = size_problem
    if( present(flags)) then
      plan%style = flags !+ FFTW_ESTIMATE + FFTW_UNALIGNED
    else
      plan%style = 0 !FFTW_ESTIMATE + FFTW_UNALIGNED
    endif
    plan%library = library
    plan%in_real => null()
    plan%out_real => array_out
    plan%in_comp => array_in
    plan%out_comp => null()

    plan%fftw_plan = fftw_plan_dft_c2r_1d(size_problem,array_in,array_out,FFTW_ESTIMATE + FFTW_UNALIGNED)
  end function

  subroutine fftw_apply_fft_c2r_1d(plan,array_in,array_out)
    type(sll_fft_plan), pointer :: plan
    sll_comp64, dimension(:)    :: array_in
    sll_real64, dimension(:)    :: array_out

    call fftw_execute_dft_c2r(plan%fftw_plan, array_in, array_out)
  end subroutine
!END C2R
! FFTW FFTW FFTW FFTW FFTW FFTW FFTW FFTW FFTW FFTW FFTW FFTW FFTW FFTW
! FFTW FFTW FFTW FFTW FFTW FFTW FFTW FFTW FFTW FFTW FFTW FFTW FFTW FFTW
! FFTW FFTW FFTW FFTW FFTW FFTW FFTW FFTW FFTW FFTW FFTW FFTW FFTW FFTW
