! SELALIB SELALIB SELALIB SELALIB SELALIB SELALIB SELALIB SELALIB SELALIB SELALIB
! SELALIB SELALIB SELALIB SELALIB SELALIB SELALIB SELALIB SELALIB SELALIB SELALIB
! SELALIB SELALIB SELALIB SELALIB SELALIB SELALIB SELALIB SELALIB SELALIB SELALIB
! COMPLEX
  function sll_new_plan_c2c_1d(library,size_problem,array_in,array_out,direction,flags) result(plan)
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

    allocate(plan%t(1:plan%N/2))
    call compute_twiddles(plan%N,plan%t)
    if ( plan%direction == FFT_FORWARD ) then
       plan%t = conjg(plan%t)
    end if
    call bit_reverse(plan%N/2,plan%t)
  end function

  subroutine sll_apply_fft_complex(plan,array_in,array_out)
    type(sll_fft_plan), pointer, intent(in)         :: plan
    sll_comp64, dimension(:), intent(inout)         :: array_in, array_out

    if( loc(array_in) .ne. loc(array_out)) then ! out-place transform
       array_out = array_in
    endif   
    
    call fft_dit_nr(array_out,plan%t,plan%direction)
    call bit_reverse_complex(plan%N,array_out)
  end subroutine

! END COMPLEX

! REAL
  function sll_new_plan_r2r_1d(library,size_problem,array_in,array_out,direction,flags) result(plan)
    sll_int32, intent(in)                        :: library
    sll_int32, intent(in)                        :: size_problem
    sll_real64, dimension(:), target, intent(in) :: array_in, array_out
    sll_int32, intent(in)                        :: direction
    sll_int32, optional,  intent(in)             :: flags
    type(sll_fft_plan), pointer                  :: plan
    sll_int32 :: i
 
    allocate(plan)
    
    plan%N = size_problem
    if( present(flags)) then
      plan%style = flags
    else
      plan%style = 0
    endif

    allocate(plan%index(0:plan%N-1))
    plan%index(0) = 0;
    plan%index(1) = plan%N/2
    do i=1,plan%N/2-1
      plan%index(2*i) = i
      plan%index(2*i + 1) = plan%N - i
    enddo


    plan%library = library
    plan%direction = direction
    plan%in_real => array_in
    plan%out_real => array_out
    plan%in_comp => null()
    plan%out_comp => null()
    allocate(plan%twiddles(0:plan%N/2-1))
    allocate(plan%twiddles_n(0:plan%N-1))
    call compute_twiddles_real_array( plan%N, plan%twiddles_n )
    call compute_twiddles_real_array( plan%N/2, plan%twiddles(0:plan%N/2-1) )
    if( plan%direction .eq. FFT_FORWARD ) then
      call conjg_in_pairs(plan%N/2,plan%twiddles)
    endif
    call bit_reverse_in_pairs( plan%N/4, plan%twiddles(0:plan%N/2-1))
  end function

  subroutine sll_apply_fft_real(plan,array_in,array_out)
    type(sll_fft_plan), pointer, intent(in)         :: plan
    sll_real64, dimension(:), intent(inout)         :: array_in, array_out

    if( loc(array_in) .ne. loc(array_out)) then ! out-place transform
       array_out = array_in
    endif

    call real_data_fft_dit( array_out, plan%N, plan%twiddles, plan%twiddles_n, plan%direction )
  end subroutine
! END REAL
! SELALIB SELALIB SELALIB SELALIB SELALIB SELALIB SELALIB SELALIB SELALIB SELALIB
! SELALIB SELALIB SELALIB SELALIB SELALIB SELALIB SELALIB SELALIB SELALIB SELALIB
! SELALIB SELALIB SELALIB SELALIB SELALIB SELALIB SELALIB SELALIB SELALIB SELALIB
