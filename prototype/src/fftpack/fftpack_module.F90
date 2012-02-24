module fftpack

 contains

  subroutine fftpack_cffti(n,wsave)
    integer , intent(in):: n
    real, dimension(4*n+15), intent(out) :: wsave
    call cffti(n,wsave)
  end subroutine fftpack_cffti

  subroutine fftpack_cfftb(n,c,wsave)
    integer, intent(in)                  :: n
    complex, dimension(n), intent(inout) :: c 
    real, dimension(4*n+15), intent(out) :: wsave
    call cfftb(n,c,wsave)
  end subroutine fftpack_cfftb

  subroutine fftpack_cfftf(n,c,wsave)
    integer, intent(in)                  :: n
    complex, dimension(n), intent(inout) :: c 
    real, dimension(4*n+15), intent(out) :: wsave
    call cfftf(n,c,wsave)
  end subroutine fftpack_cfftf

  subroutine fftpack_rffti(n,wsave)
    integer , intent(in):: n
    real, dimension(2*n+15), intent(out) :: wsave
    call rffti(n,wsave)
  end subroutine fftpack_rffti

  subroutine fftpack_rfftb(n,c,wsave)
    integer, intent(in)                  :: n
    real, dimension(n), intent(inout)    :: c 
    real, dimension(2*n+15), intent(out) :: wsave
    call rfftb(n,c,wsave)
  end subroutine fftpack_rfftb

  subroutine fftpack_rfftf(n,c,wsave)
    integer, intent(in)                  :: n
    real, dimension(n), intent(inout)    :: c 
    real, dimension(2*n+15), intent(out) :: wsave
    call rfftf(n,c,wsave)
  end subroutine fftpack_rfftf
end module fftpack

