! FFTW FFTW FFTW FFTW FFTW FFTW FFTW FFTW FFTW FFTW FFTW FFTW FFTW FFTW
! FFTW FFTW FFTW FFTW FFTW FFTW FFTW FFTW FFTW FFTW FFTW FFTW FFTW FFTW
! FFTW FFTW FFTW FFTW FFTW FFTW FFTW FFTW FFTW FFTW FFTW FFTW FFTW FFTW
! COMPLEX
  function fftw_new_plan_c2c_1d(library,nx,array_in,array_out,direction,flags) result(plan)
    sll_int32, intent(in)                        :: library
    sll_int32, intent(in)                        :: nx
    sll_comp64, dimension(:), intent(inout)         :: array_in
    sll_comp64, dimension(:), intent(inout)        :: array_out
    sll_int32, intent(in)                        :: direction
    sll_int32, optional, intent(in)              :: flags
    type(sll_fft_plan), pointer                      :: plan
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

  subroutine fftw_apply_fft_complex(plan,array_in,array_out)
    type(sll_fft_plan), pointer, intent(in)             :: plan
    sll_comp64, dimension(:), intent(inout)         :: array_in, array_out

    call fftw_execute_dft(plan%fftw_plan, array_in, array_out)
  end subroutine 
! END COMPLEX

! COMPLEX 2D
  function fftw_new_plan_c2c_2d(library,NX,NY,array_in,array_out,direction,flags) result(plan)
    sll_int32, intent(in)                            :: library
    sll_int32, intent(in)                            :: NX,NY
    sll_comp64, dimension(0:,0:)                     :: array_in, array_out
    sll_int32, intent(in)                            :: direction
    sll_int32, optional, intent(in)                  :: flags
    type(sll_fft_plan), pointer                          :: plan
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
  
    !We must swith the dimension. It's a fftw convention. 
    plan%fftw_plan = fftw_plan_dft_2d(NY,NX,array_in,array_out,direction,FFTW_ESTIMATE + FFTW_UNALIGNED)
  end function

  subroutine fftw_apply_plan_c2c_2d(plan,array_in,array_out)
    type(sll_fft_plan), pointer, intent(in)      :: plan
    sll_comp64, dimension(0:,0:), intent(inout)     :: array_in, array_out

    call fftw_execute_dft(plan%fftw_plan, array_in, array_out)
  end subroutine 
! END COMPLEX 2D

! REAL
  function fftw_new_plan_r2r_1d(library,nx,array_in,array_out,direction,flags) result(plan)
    sll_int32, intent(in)                        :: library
    sll_int32, intent(in)                        :: nx
    sll_real64, dimension(:), intent(inout)      :: array_in
    sll_real64, dimension(:), intent(inout)      :: array_out
    sll_int32, intent(in)                        :: direction
    sll_int32, optional, intent(in)              :: flags
    type(sll_fft_plan), pointer                     :: plan
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

  subroutine fftw_apply_fft_real(plan,array_in,array_out)
    type(sll_fft_plan), pointer, intent(in)         :: plan
    sll_real64, dimension(:), intent(inout)     :: array_in, array_out

    call fftw_execute_r2r(plan%fftw_plan, array_in, array_out)
  end subroutine

!  subroutine fftw_apply_fft_re(plan)
!    type(sll_sll_fft_plan), pointer, intent(in)         :: plan
!    call fftw_execute_r2r(plan%fftw_plan,plan%in_real,plan%out_real)
!  end subroutine
! END REAL


! R2C
  function fftw_new_plan_r2c_1d(library,nx,array_in,array_out,flags) result(plan)
    sll_int32, intent(in)                           :: library
    sll_int32, intent(in)                           :: nx
    sll_real64, dimension(:), intent(inout) :: array_in
    sll_comp64, dimension(:), intent(out)   :: array_out
    sll_int32, optional, intent(in)                 :: flags
    type(sll_fft_plan), pointer                     :: plan
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

  subroutine fftw_apply_fft_r2c_1d(plan,array_in,array_out)
    type(sll_fft_plan), pointer, intent(in)             :: plan
    sll_real64, dimension(:), intent(inout)         :: array_in
    sll_comp64, dimension(:), intent(out)           :: array_out

    call fftw_execute_dft_r2c(plan%fftw_plan, array_in, array_out)
  end subroutine
!END R2C

! C2R
  function fftw_new_plan_c2r_1d(library,nx,array_in,array_out,flags) result(plan)
    sll_int32, intent(in)                           :: library
    sll_int32, intent(in)                           :: nx
    sll_comp64, dimension(:)               :: array_in
    sll_real64, dimension(:)                :: array_out
    sll_int32, optional, intent(in)                 :: flags
    type(sll_fft_plan), pointer                     :: plan
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
