! SELALIB SELALIB SELALIB SELALIB SELALIB SELALIB SELALIB SELALIB SELALIB SELALIB
! SELALIB SELALIB SELALIB SELALIB SELALIB SELALIB SELALIB SELALIB SELALIB SELALIB
! SELALIB SELALIB SELALIB SELALIB SELALIB SELALIB SELALIB SELALIB SELALIB SELALIB
! COMPLEX

! ------
! - 1D -
! ------
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

! ------
! - 3D -
! ------
  function sll_new_plan_c2c_1d_for_3d(library,N,N2,N3,array_in,array_out,direction,flags) result(plan)
    sll_int32, intent(in)                            :: library
    sll_int32, intent(in)                            :: N,N2,N3
    sll_comp64, dimension(0:,0:,0:), target, intent(in) :: array_in, array_out
    sll_int32, intent(in)                            :: direction
    sll_int32, optional,  intent(in)                 :: flags
    type(sll_fft_plan), pointer                      :: plan
    sll_int32 :: size_problem
 
    allocate(plan)
    
    plan%N = N
    plan%N2 = N2
    plan%N3 = N3
    if( present(flags)) then
      plan%style = flags
    else
      plan%style = 0
    endif
    plan%library = library
    plan%direction = direction
    !plan%in_comp => array_in
    !plan%out_comp => array_out
    !plan%in_real => null()
    !plan%out_real => null()

    size_problem = 0
    if( (direction==FFT_FORWARD_X) .or. (direction==FFT_INVERSE_X) ) then
      size_problem = N
    else if( (direction==FFT_FORWARD_Y) .or. (direction==FFT_INVERSE_Y) ) then
      size_problem = N2
    else if( (direction==FFT_FORWARD_Z) .or. (direction==FFT_INVERSE_Z) ) then
      size_problem = N3
    endif

    allocate(plan%t(1:size_problem/2))
    call compute_twiddles(size_problem,plan%t)
    if ( plan%direction >= FFT_FORWARD ) then
       plan%t = conjg(plan%t)
    end if
    call bit_reverse(size_problem/2,plan%t)
  end function

  subroutine sll_apply_fft_c2c_1d_for_3d(plan,array_in,array_out)
    type(sll_fft_plan), pointer, intent(in)         :: plan
    sll_comp64, dimension(0:,0:,0:), intent(inout)  :: array_in, array_out
    sll_int32                                       :: i, j, direction

    if( loc(array_in) .ne. loc(array_out)) then ! out-place transform
       array_out = array_in
    endif   
     
    if( plan%direction > 0 ) then
      direction = FFT_FORWARD
    else
      direction = FFT_INVERSE
    endif    

    if( (plan%direction .eq. FFT_FORWARD_X) .or. (plan%direction .eq. FFT_INVERSE_X) ) then
      do i=0,plan%N2-1
       do j=0,plan%N3-1
        call fft_dit_nr(array_out(:,i,j),plan%t,direction)
        call bit_reverse_complex(plan%N,array_out(:,i,j))
       enddo
      enddo
    else if( (plan%direction==FFT_FORWARD_Y) .or. (plan%direction==FFT_INVERSE_Y) ) then
      do i=0,plan%N-1
       do j=0,plan%N3-1
        call fft_dit_nr(array_out(i,:,j),plan%t,direction)
        call bit_reverse_complex(plan%N,array_out(i,:,j))
       enddo
      enddo
    else if( (plan%direction==FFT_FORWARD_Z) .or. (plan%direction==FFT_INVERSE_Z) ) then
      do i=0,plan%N-1
       do j=0,plan%N2-1
        call fft_dit_nr(array_out(i,j,:),plan%t,direction)
        call bit_reverse_complex(plan%N,array_out(i,j,:))
       enddo
      enddo
    endif
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
