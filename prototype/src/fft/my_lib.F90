! SELALIB SELALIB SELALIB SELALIB SELALIB SELALIB SELALIB SELALIB SELALIB SELALIB
! SELALIB SELALIB SELALIB SELALIB SELALIB SELALIB SELALIB SELALIB SELALIB SELALIB
! SELALIB SELALIB SELALIB SELALIB SELALIB SELALIB SELALIB SELALIB SELALIB SELALIB
! COMPLEX

! ------
! - 1D -
! ------
  function sll_new_plan_c2c_1d(library,nx,array_in,array_out,direction,flags) result(plan)
    sll_int32, intent(in)                        :: library
    sll_int32, intent(in)                        :: nx
    sll_comp64, dimension(:), intent(in)         :: array_in, array_out
    sll_int32, intent(in)                        :: direction
    sll_int32, optional,  intent(in)             :: flags
    type(sll_fft_plan), pointer                      :: plan
    sll_int32                                    :: ierr

    SLL_ASSERT(size(array_in).eq.nx)
    SLL_ASSERT(size(array_out).eq.nx)
    SLL_ALLOCATE(plan,ierr)
    plan%library = library 
    plan%direction = direction
    if( present(flags) )then
      plan%style = flags
    else
      plan%style = 0_f32
    endif
    plan%problem_rank = 1
    SLL_ALLOCATE(plan%problem_shape(1),ierr)
    plan%problem_shape = (/ nx  /)

    SLL_ALLOCATE(plan%t(1:nx/2),ierr)
    call compute_twiddles(nx,plan%t)
    if ( direction == FFT_FORWARD ) then
       plan%t = conjg(plan%t)
    end if
    call bit_reverse(nx/2,plan%t)
  end function

  subroutine sll_apply_fft_complex(plan,array_in,array_out)
    type(sll_fft_plan), pointer                 :: plan
    sll_comp64, dimension(:), intent(inout)     :: array_in, array_out

    if( loc(array_in) .ne. loc(array_out)) then ! out-place transform
       array_out = array_in
    endif   
    
    call fft_dit_nr(array_out,plan%t,plan%direction)
    call bit_reverse_complex(plan%problem_shape(1),array_out)
  end subroutine


#ifdef CACA 
! --------------------
! - 2D One direction -
! --------------------
  function sll_new_plan_c2c_1d_for_2d(library,N,N2,array_in,array_out,direction,flags) result(plan)
    sll_int32, intent(in)                            :: library
    sll_int32, intent(in)                            :: N,N2
    sll_comp64, dimension(0:,0:), target, intent(in) :: array_in, array_out
    sll_int32, intent(in)                            :: direction
    sll_int32, optional,  intent(in)                 :: flags
    type(sll_fft_plan), pointer                      :: plan
    sll_int32 :: size_problem
 
    allocate(plan)
    
    plan%N = N
    plan%N2 = N2
    plan%N3 = -1
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
      stop 'direction error in sll_new_plan_c2c_1d_for_2d'
    endif
 
    allocate(plan%t(1:size_problem/2))
    call compute_twiddles(size_problem,plan%t)
    if ( plan%direction >= FFT_FORWARD ) then
       plan%t = conjg(plan%t)
    end if
    call bit_reverse(size_problem/2,plan%t)
  end function

  subroutine sll_apply_fft_c2c_1d_for_2d(plan,array_in,array_out)
    type(sll_fft_plan), pointer                     :: plan
    sll_comp64, dimension(0:,0:), intent(inout)     :: array_in, array_out
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
        call fft_dit_nr(array_out(:,i),plan%t,direction)
        call bit_reverse_complex(plan%N,array_out(:,i))
      enddo
    else if( (plan%direction==FFT_FORWARD_Y) .or. (plan%direction==FFT_INVERSE_Y) ) then
      do i=0,plan%N-1
        call fft_dit_nr(array_out(i,:),plan%t,direction)
        call bit_reverse_complex(plan%N2,array_out(i,:))
      enddo
    endif
  end subroutine

#endif

! --------------------
! - 2D one direction -
! --------------------
  function sll_new_plan_c2c_1d_for_2d(library,NX,NY,array_in,array_out,direction,flags) result(plan)
    sll_int32, intent(in)                            :: library
    sll_int32, intent(in)                            :: NX,NY
    sll_comp64, dimension(0:,0:), target, intent(in) :: array_in, array_out
    sll_int32, intent(in)                            :: direction
    sll_int32, optional,  intent(in)                 :: flags
    type(sll_fft_plan), pointer                          :: plan
    sll_int32                                        :: ierr    

    SLL_ASSERT(size(array_in,dim=1).eq.NX)
    SLL_ASSERT(size(array_in,dim=2).eq.NY)
    SLL_ASSERT(size(array_out,dim=1).eq.NX)
    SLL_ASSERT(size(array_out,dim=2).eq.NY)
    SLL_ALLOCATE(plan,ierr)
    plan%library = library 
    plan%direction = direction
    if( present(flags) )then
      plan%style = flags
    else
      plan%style = 0_f32
    endif
    plan%problem_rank = 2
    SLL_ALLOCATE(plan%problem_shape(2),ierr)
    plan%problem_shape = (/ NX , NY /)

    if( fft_is_present_flag(plan,FFT_ONLY_FIRST_DIRECTION) ) then
      SLL_ALLOCATE(plan%t(1:NX/2),ierr)
      call compute_twiddles(NX,plan%t(1:NX/2))
      if ( direction == FFT_FORWARD ) then
         plan%t(1:NX/2) = conjg(plan%t(1:NX/2))
      end if
      call bit_reverse(NX/2,plan%t(1:NX/2))
    else if( fft_is_present_flag(plan,FFT_ONLY_SECOND_DIRECTION) ) then
      SLL_ALLOCATE(plan%t(1:NY/2),ierr)
      call compute_twiddles(NY,plan%t(1:NY/2))
      if ( direction == FFT_FORWARD ) then
         plan%t(1:NY/2) = conjg(plan%t(1:NY/2))
      end if
      call bit_reverse(NY/2,plan%t(1:NY/2))
    endif
  end function

  subroutine sll_apply_fft_c2c_1d_for_2d(plan,array_in,array_out)
    type(sll_fft_plan), pointer                     :: plan
    sll_comp64, dimension(0:,0:), intent(inout)     :: array_in, array_out
    sll_int32                                       :: i, nx, ny
    sll_int32 :: fft_shape(2)

    if( loc(array_in) .ne. loc(array_out)) then ! out-place transform
       array_out = array_in
    endif   
      
    fft_shape = fft_get_shape(plan)
    nx = fft_shape(1)
    ny = fft_shape(2)

    if( fft_is_present_flag(plan,FFT_ONLY_FIRST_DIRECTION) ) then    
      do i=0,ny-1
        call fft_dit_nr(array_out(0:nx-1,i),plan%t(1:nx/2),plan%direction)
        call bit_reverse_complex(nx,array_out(0:nx-1,i))
      enddo
    else if( fft_is_present_flag(plan,FFT_ONLY_SECOND_DIRECTION) ) then
      do i=0,nx-1
        call fft_dit_nr(array_out(i,0:ny-1),plan%t(1:ny/2),plan%direction)
        call bit_reverse_complex(ny,array_out(i,0:ny-1))
      enddo
    endif
  end subroutine

! --------------------
! - 2D Two direction -
! --------------------
  function sll_new_plan_c2c_2d(library,NX,NY,array_in,array_out,direction,flags) result(plan)
    sll_int32, intent(in)                            :: library
    sll_int32, intent(in)                            :: NX,NY
    sll_comp64, dimension(0:,0:), target, intent(in) :: array_in, array_out
    sll_int32, intent(in)                            :: direction
    sll_int32, optional,  intent(in)                 :: flags
    type(sll_fft_plan), pointer                          :: plan
    sll_int32                                        :: ierr    
    logical                                          :: two_direction !true if dft in the two directions, false otherwise.

    SLL_ASSERT(size(array_in,dim=1).eq.NX)
    SLL_ASSERT(size(array_in,dim=2).eq.NY)
    SLL_ASSERT(size(array_out,dim=1).eq.NX)
    SLL_ASSERT(size(array_out,dim=2).eq.NY)
    SLL_ALLOCATE(plan,ierr)
    plan%library = library
    plan%direction = direction
    if( present(flags) )then
      plan%style = flags
    else
      plan%style = 0_f32
    endif
    plan%problem_rank = 2
    SLL_ALLOCATE(plan%problem_shape(2),ierr)
    plan%problem_shape = (/ NX , NY /)

    two_direction = .false.
    if( fft_is_present_flag(plan%style,FFT_ONLY_FIRST_DIRECTION) .and. fft_is_present_flag(plan%style,FFT_ONLY_SECOND_DIRECTION) ) then
      SLL_ALLOCATE(plan%t(1:NX/2 + NY/2),ierr)
    else if( fft_is_present_flag(plan%style,FFT_ONLY_FIRST_DIRECTION) ) then
      SLL_ALLOCATE(plan%t(1:NX/2),ierr)
    else if( fft_is_present_flag(plan%style,FFT_ONLY_SECOND_DIRECTION) ) then
      SLL_ALLOCATE(plan%t(NX/2+1:NX/2+NY/2),ierr)
    else
      !If we are here, there is no FFT_ONLY_XXXXX_DIRECTION flags. So we want a 2D FFT in all direction.
      SLL_ALLOCATE(plan%t(1:NX/2 + NY/2),ierr)
      two_direction = .true.
    endif
    
    if( fft_is_present_flag(plan%style,FFT_ONLY_FIRST_DIRECTION) .or. two_direction ) then
      call compute_twiddles(NX,plan%t(1:NX/2))
      if ( direction == FFT_FORWARD ) then
        plan%t(1:NX/2) = conjg(plan%t(1:NX/2))
      end if
      call bit_reverse(NX/2,plan%t(1:NX/2))
    endif

    if( fft_is_present_flag(plan%style,FFT_ONLY_SECOND_DIRECTION) .or. two_direction ) then
      call compute_twiddles(NY,plan%t(NX/2+1:NX/2 + NY/2))
      if ( direction == FFT_FORWARD ) then
        plan%t(NX/2+1:NX/2 + NY/2) = conjg(plan%t(NX/2+1:NX/2 + NY/2))
      end if
      call bit_reverse(NY/2,plan%t(NX/2+1:NX/2 + NY/2))
    endif
  end function

  subroutine sll_apply_fft_c2c_2d(plan,array_in,array_out)
    type(sll_fft_plan), pointer                     :: plan
    sll_comp64, dimension(0:,0:), intent(inout)     :: array_in, array_out
    sll_int32                                       :: i, nx, ny
    sll_int32                                       :: fft_shape(2)
    logical                                         :: two_direction

    if( loc(array_in) .ne. loc(array_out)) then ! out-place transform
       array_out = array_in
    endif   
      
    fft_shape = fft_get_shape(plan)
    nx = fft_shape(1)
    ny = fft_shape(2)
    
    two_direction = .true.
    if( fft_is_present_flag(plan%style,FFT_ONLY_FIRST_DIRECTION) .or. fft_is_present_flag(plan%style,FFT_ONLY_SECOND_DIRECTION) ) then
      two_direction = .false.
    endif

    
    if( fft_is_present_flag(plan%style,FFT_ONLY_FIRST_DIRECTION) .or. two_direction ) then
      do i=0,ny-1
        call fft_dit_nr(array_out(0:nx-1,i),plan%t(1:nx/2),plan%direction)
        call bit_reverse_complex(nx,array_out(0:nx-1,i))
      enddo
    endif
    if( fft_is_present_flag(plan%style,FFT_ONLY_SECOND_DIRECTION) .or. two_direction ) then
      do i=0,nx-1
        call fft_dit_nr(array_out(i,0:ny-1),plan%t(nx/2+1:nx/2+ny/2),plan%direction)
        call bit_reverse_complex(ny,array_out(i,0:ny-1))
      enddo
    endif
  end subroutine

! ------
! - 3D -
! ------
!  function sll_new_plan_c2c_1d_for_3d(library,N,N2,N3,array_in,array_out,direction,flags) result(plan)
!    sll_int32, intent(in)                            :: library
!    sll_int32, intent(in)                            :: N,N2,N3
!    sll_comp64, dimension(0:,0:,0:), target, intent(in) :: array_in, array_out
!    sll_int32, intent(in)                            :: direction
!    sll_int32, optional,  intent(in)                 :: flags
!    type(sll_fft_plan), pointer                      :: plan
!    sll_int32 :: size_problem
! 
!    allocate(plan)
!    
!    plan%N = N
!    plan%N2 = N2
!    plan%N3 = N3
!    if( present(flags)) then
!      plan%style = flags
!    else
!      plan%style = 0
!    endif
!    plan%library = library
!    plan%direction = direction
!    !plan%in_comp => array_in
!    !plan%out_comp => array_out
!    !plan%in_real => null()
!    !plan%out_real => null()
!
!    size_problem = 0
!    if( (direction==FFT_FORWARD_X) .or. (direction==FFT_INVERSE_X) ) then
!      size_problem = N
!    else if( (direction==FFT_FORWARD_Y) .or. (direction==FFT_INVERSE_Y) ) then
!      size_problem = N2
!    else if( (direction==FFT_FORWARD_Z) .or. (direction==FFT_INVERSE_Z) ) then
!      size_problem = N3
!    endif
!
!    allocate(plan%t(1:size_problem/2))
!    call compute_twiddles(size_problem,plan%t)
!    if ( plan%direction >= FFT_FORWARD ) then
!       plan%t = conjg(plan%t)
!    end if
!    call bit_reverse(size_problem/2,plan%t)
!  end function
!
!  subroutine sll_apply_fft_c2c_1d_for_3d(plan,array_in,array_out)
!    type(sll_fft_plan), pointer, intent(in)         :: plan
!    sll_comp64, dimension(0:,0:,0:), intent(inout)  :: array_in, array_out
!    sll_int32                                       :: i, j, direction
!
!    if( loc(array_in) .ne. loc(array_out)) then ! out-place transform
!       array_out = array_in
!    endif   
!     
!    if( plan%direction > 0 ) then
!      direction = FFT_FORWARD
!    else
!      direction = FFT_INVERSE
!    endif    
!
!    if( (plan%direction .eq. FFT_FORWARD_X) .or. (plan%direction .eq. FFT_INVERSE_X) ) then
!      do i=0,plan%N2-1
!       do j=0,plan%N3-1
!        call fft_dit_nr(array_out(:,i,j),plan%t,direction)
!        call bit_reverse_complex(plan%N,array_out(:,i,j))
!       enddo
!      enddo
!    else if( (plan%direction==FFT_FORWARD_Y) .or. (plan%direction==FFT_INVERSE_Y) ) then
!      do i=0,plan%N-1
!       do j=0,plan%N3-1
!        call fft_dit_nr(array_out(i,:,j),plan%t,direction)
!        call bit_reverse_complex(plan%N,array_out(i,:,j))
!       enddo
!      enddo
!    else if( (plan%direction==FFT_FORWARD_Z) .or. (plan%direction==FFT_INVERSE_Z) ) then
!      do i=0,plan%N-1
!       do j=0,plan%N2-1
!        call fft_dit_nr(array_out(i,j,:),plan%t,direction)
!        call bit_reverse_complex(plan%N,array_out(i,j,:))
!       enddo
!      enddo
!    endif
!  end subroutine
!
! END COMPLEX






! REAL
  function sll_new_plan_r2r_1d(library,nx,array_in,array_out,direction,flags) result(plan)
    sll_int32, intent(in)                        :: library
    sll_int32, intent(in)                        :: nx, direction
    sll_real64, dimension(:), intent(in)         :: array_in
    sll_real64, dimension(:), intent(in)         :: array_out
    sll_int32, optional,  intent(in)             :: flags
    type(sll_fft_plan), pointer                      :: plan
    sll_int32                                    :: ierr

    SLL_ASSERT(size(array_in).eq.nx)
    SLL_ASSERT(size(array_out).eq.nx) 
    SLL_ALLOCATE(plan,ierr)
    plan%library = library
    plan%direction = direction 
    if( present(flags) )then
      plan%style = flags
    else
      plan%style = 0_f32
    endif
    plan%problem_rank = 1
    SLL_ALLOCATE(plan%problem_shape(1),ierr)
    plan%problem_shape = (/ nx /)

    SLL_ALLOCATE(plan%twiddles(0:nx/2-1),ierr)
    SLL_ALLOCATE(plan%twiddles_n(0:nx-1),ierr)
    call compute_twiddles_real_array( nx, plan%twiddles_n )
    call compute_twiddles_real_array( nx/2, plan%twiddles(0:nx/2-1) )
    if( direction .eq. FFT_FORWARD ) then
      call conjg_in_pairs(nx/2,plan%twiddles)
    endif
    call bit_reverse_in_pairs( nx/4, plan%twiddles(0:nx/2-1))
  end function

  subroutine sll_apply_fft_real(plan,array_in,array_out)
    implicit none
    type(sll_fft_plan), pointer                     :: plan
    sll_real64, dimension(:), intent(inout)         :: array_in, array_out
    sll_int32 :: nx

    nx = plan%problem_shape(1)

    if( loc(array_in) .ne. loc(array_out)) then ! out-place transform
       array_out = array_in
    endif
    call real_data_fft_dit( array_out, nx, plan%twiddles, plan%twiddles_n, plan%direction )
  end subroutine
! END REAL







! REAL TO COMPLEX
  function sll_new_plan_r2c_1d(library,nx,array_in,array_out,flags) result(plan)
    sll_int32, intent(in)                        :: library
    sll_int32, intent(in)                        :: nx
    sll_real64, dimension(:), intent(in)         :: array_in
    sll_comp64, dimension(:), intent(in)         :: array_out
    sll_int32, optional,  intent(in)             :: flags
    type(sll_fft_plan), pointer                      :: plan
    sll_int32                                    :: ierr

    SLL_ASSERT(size(array_in).eq.nx)
    SLL_ASSERT(size(array_out).eq.nx/2+1) 
    SLL_ALLOCATE(plan,ierr)
    plan%library = library
    plan%direction = FFT_FORWARD
    if( present(flags) )then
      plan%style = flags
    else
      plan%style = 0_f32
    endif
    plan%problem_rank = 1
    SLL_ALLOCATE(plan%problem_shape(1),ierr)
    plan%problem_shape = (/ nx /)

    SLL_ALLOCATE(plan%twiddles(0:nx/2-1),ierr)
    SLL_ALLOCATE(plan%twiddles_n(0:nx-1),ierr)
    call compute_twiddles_real_array( nx, plan%twiddles_n )
    call compute_twiddles_real_array( nx/2, plan%twiddles(0:nx/2-1) )
    call conjg_in_pairs(nx/2,plan%twiddles)
    call bit_reverse_in_pairs( nx/4, plan%twiddles(0:nx/2-1))
  end function

  subroutine sll_fft_apply_r2c_1d(plan,array_in,array_out)
    implicit none
    type(sll_fft_plan), pointer                     :: plan
    sll_real64, dimension(0:), intent(inout)        :: array_in
    sll_comp64, dimension(0:), intent(out)          :: array_out
    sll_int32                                       :: nx, i

    nx = plan%problem_shape(1)
 
    call real_data_fft_dit( array_in, nx , plan%twiddles, plan%twiddles_n, plan%direction )
    !mode k=0
    array_out(0) = cmplx(array_in(0),0.0_f64,kind=f64)
    !mode k=n/2
    array_out(nx/2) = cmplx(array_in(1),0.0_f64,kind=f64)
    !mode k=1 to k= n-2
    do i=1,nx/2-1
        array_out(i) = cmplx(array_in(2*i),array_in(2*i+1),kind=f64)
        !array_out(plan%N-i) = cmplx(array_in(2*i),-array_in(2*i+1),kind=f64)
    enddo
  end subroutine
! END REAL TO COMPLEX

! COMPLEX TO REAL
  function sll_new_plan_c2r_1d(library,nx,array_in,array_out,flags) result(plan)
    sll_int32, intent(in)                        :: library
    sll_int32, intent(in)                        :: nx
    sll_comp64, dimension(:), intent(in)         :: array_in
    sll_real64, dimension(:), intent(in)         :: array_out
    sll_int32, optional,  intent(in)             :: flags
    type(sll_fft_plan), pointer                      :: plan
    sll_int32                                    :: ierr

    SLL_ASSERT(size(array_in).eq.nx/2+1)
    SLL_ASSERT(size(array_out).eq.nx)
    SLL_ALLOCATE(plan,ierr)
    plan%library = library
    plan%direction = FFT_INVERSE
    if( present(flags) )then
      plan%style = flags
    else
      plan%style = 0_f32
    endif
    plan%problem_rank = 1
    SLL_ALLOCATE(plan%problem_shape(1),ierr)
    plan%problem_shape = (/ nx /)

    SLL_ALLOCATE(plan%twiddles(0:nx/2-1),ierr)
    SLL_ALLOCATE(plan%twiddles_n(0:nx-1),ierr)
    call compute_twiddles_real_array( nx, plan%twiddles_n )
    call compute_twiddles_real_array( nx/2, plan%twiddles(0:nx/2-1) )
    call bit_reverse_in_pairs( nx/4, plan%twiddles(0:nx/2-1))
  end function

  subroutine sll_fft_apply_c2r_1d(plan,array_in,array_out)
    implicit none
    type(sll_fft_plan), pointer                   :: plan
    sll_comp64, dimension(0:), intent(inout)      :: array_in
    sll_real64, dimension(0:), intent(out)        :: array_out
    sll_int32                                     :: nx, i

    nx = plan%problem_shape(1)
 
    !mode k=0
    array_out(0) = real(array_in(0),kind=f64)
    !mode k=n/2
    array_out(1) = real(array_in(nx/2),kind=f64)
    !mode k=1 to k= n-2
    do i=1,nx/2-1
      array_out(2*i) = real(array_in(i),kind=f64)
      array_out(2*i+1) = dimag(array_in(i))
    enddo
    call real_data_fft_dit( array_out, nx , plan%twiddles, plan%twiddles_n, plan%direction )
  end subroutine
! END COMPLEX TO REAL

! REAL TO COMPLEX 2D
  function sll_new_plan_r2c_2d(library,nx,ny,array_in,array_out,flags) result(plan)
    sll_int32, intent(in)                        :: library
    sll_int32, intent(in)                        :: nx,ny
    sll_real64, dimension(:,:), intent(in)       :: array_in
    sll_comp64, dimension(:,:), intent(in)       :: array_out
    sll_int32, optional,  intent(in)             :: flags
    type(sll_fft_plan), pointer                      :: plan
    sll_int32                                    :: ierr

    SLL_ASSERT(size(array_in,dim=1).eq.nx)
    SLL_ASSERT(size(array_in,dim=2).eq.ny)
    SLL_ASSERT(size(array_out,dim=1).eq.nx/2+1)
    SLL_ASSERT(size(array_out,dim=2).eq.ny) 
    SLL_ALLOCATE(plan,ierr)
    plan%library = library
    plan%direction = FFT_FORWARD
    if( present(flags) )then
      plan%style = flags
    else
      plan%style = 0_f32
    endif
    plan%problem_rank = 2
    SLL_ALLOCATE(plan%problem_shape(2),ierr)
    plan%problem_shape = (/ nx , ny /)

    SLL_ALLOCATE(plan%t(1:ny/2),ierr)
    call compute_twiddles(ny,plan%t)
    plan%t = conjg(plan%t)
    call bit_reverse(ny/2,plan%t)

    SLL_ALLOCATE(plan%twiddles(0:nx/2-1),ierr)
    SLL_ALLOCATE(plan%twiddles_n(0:nx-1),ierr)
    call compute_twiddles_real_array( nx, plan%twiddles_n )
    call compute_twiddles_real_array( nx/2, plan%twiddles(0:nx/2-1) )
    call conjg_in_pairs(nx/2,plan%twiddles)
    call bit_reverse_in_pairs( nx/4, plan%twiddles(0:nx/2-1))
  end function

  subroutine sll_fft_apply_r2c_2d(plan,array_in,array_out)
    implicit none
    type(sll_fft_plan), pointer                     :: plan
    sll_real64, dimension(0:,0:), intent(inout)     :: array_in
    sll_comp64, dimension(0:,0:), intent(out)       :: array_out
    sll_int32                                       :: nx, i, ny, k

    nx = plan%problem_shape(1)
    ny = plan%problem_shape(2)
 
    do k=0,ny-1
      call real_data_fft_dit( array_in(0:nx-1,k), nx , plan%twiddles, plan%twiddles_n, plan%direction )
      array_out(0,k) = cmplx(array_in(0,k),0.0_f64,kind=f64)
      array_out(nx/2,k) = cmplx(array_in(1,k),0.0_f64,kind=f64)
      do i=1,nx/2-1
          array_out(i,k) = cmplx(array_in(2*i,k),array_in(2*i+1,k),kind=f64)
      enddo
    enddo
    do k=0,nx/2
      call fft_dit_nr(array_out(k,0:ny-1),plan%t,plan%direction)
      call bit_reverse_complex(ny,array_out(k,0:ny-1))
    enddo

  end subroutine
! END REAL TO COMPLEX 2D

! COMPLEX TO REAL 2D
  function sll_new_plan_c2r_2d(library,nx,ny,array_in,array_out,flags) result(plan)
    sll_int32, intent(in)                        :: library
    sll_int32, intent(in)                        :: nx,ny
    sll_comp64, dimension(:,:), intent(in)       :: array_in
    sll_real64, dimension(:,:), intent(in)       :: array_out
    sll_int32, optional,  intent(in)             :: flags
    type(sll_fft_plan), pointer                      :: plan
    sll_int32                                    :: ierr
 
    SLL_ASSERT(size(array_in,dim=1).eq.nx/2+1)
    SLL_ASSERT(size(array_in,dim=2).eq.ny)
    SLL_ASSERT(size(array_out,dim=1).eq.nx)
    SLL_ASSERT(size(array_out,dim=2).eq.ny)

    SLL_ALLOCATE(plan,ierr)
    plan%library = library
    plan%direction = FFT_INVERSE
    if( present(flags) )then
      plan%style = flags
    else
      plan%style = 0_f32
    endif
    plan%problem_rank = 2
    SLL_ALLOCATE(plan%problem_shape(2),ierr)
    plan%problem_shape = (/ nx , ny /)

    SLL_ALLOCATE(plan%t(1:ny/2),ierr)
    call compute_twiddles(ny,plan%t)
    call bit_reverse(ny/2,plan%t)

    SLL_ALLOCATE(plan%twiddles(0:nx/2-1),ierr)
    SLL_ALLOCATE(plan%twiddles_n(0:nx-1),ierr)
    call compute_twiddles_real_array( nx, plan%twiddles_n )
    call compute_twiddles_real_array( nx/2, plan%twiddles(0:nx/2-1) )
    call bit_reverse_in_pairs( nx/4, plan%twiddles(0:nx/2-1))
  end function

  subroutine sll_fft_apply_c2r_2d(plan,array_in,array_out)
    implicit none
    type(sll_fft_plan), pointer                     :: plan
    sll_comp64, dimension(0:,0:), intent(inout)     :: array_in
    sll_real64, dimension(0:,0:), intent(out)       :: array_out
    sll_int32                                       :: nx, i, ny, k, j

    nx = plan%problem_shape(1)
    ny = plan%problem_shape(2)

    do j=0,nx/2
      call fft_dit_nr(array_in(j,0:ny-1),plan%t,plan%direction)
      call bit_reverse_complex(ny,array_in(j,0:ny-1))
    enddo
    do i=0,ny-1
      array_out(0,i) = real(array_in(0,i),kind=f64)
      array_out(1,i) = real(array_in(nx/2,i),kind=f64)
      do k=1,nx/2-1
        array_out(2*k,i) = real(array_in(k,i),kind=f64)
        array_out(2*k+1,i) = dimag(array_in(k,i))
      enddo
      call real_data_fft_dit( array_out(0:nx-1,i), nx , plan%twiddles, plan%twiddles_n, plan%direction )
    enddo
  end subroutine
! END COMPLEX TO REAL 2D

! SELALIB SELALIB SELALIB SELALIB SELALIB SELALIB SELALIB SELALIB SELALIB SELALIB
! SELALIB SELALIB SELALIB SELALIB SELALIB SELALIB SELALIB SELALIB SELALIB SELALIB
! SELALIB SELALIB SELALIB SELALIB SELALIB SELALIB SELALIB SELALIB SELALIB SELALIB
