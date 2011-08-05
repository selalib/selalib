module fft1d_module
  implicit none
  private
  type, public :: fft1dclass
     real, dimension(:), pointer :: coefs
     double precision, dimension(:), pointer :: coefd
     integer  :: n  ! number of samples in each sequence
  end type  fft1dclass
!  interface initfft
!     module procedure initsimpfft, initdoubfft
!  end interface
  interface fft
     module procedure simpfft, doubfft
  end interface
  interface fftinv
     module procedure simpfftinv, doubfftinv
  end interface

  public :: initdoubfft, initsimpfft, fft, fftinv
  contains

    subroutine initsimpfft(this,l)
      type( fft1dclass) :: this
      integer :: l  

      this%n = l
      allocate(this%coefs(2*this%n+15))
      call rffti(this%n,this%coefs)

    end subroutine initsimpfft

    subroutine initdoubfft(this,l)
      type( fft1dclass) :: this
      integer :: l 

      this%n = l
      allocate(this%coefd(2*this%n+15))
      call dffti(this%n,this%coefd)

    end subroutine initdoubfft
  
    subroutine simpfft(this,array)
      type( fft1dclass) :: this
      real, dimension(:) :: array

      call rfftf( this%n, array, this%coefs )

      array = array /this%n      ! normalize FFT
    end subroutine simpfft

    subroutine doubfft(this,array)
      type( fft1dclass) :: this
      double precision, dimension(:) :: array

    call dfftf( this%n, array, this%coefd)

      array = array /this%n      ! normalize FFT
    end subroutine doubfft

    subroutine simpfftinv(this,array)
      type( fft1dclass) :: this
      real, dimension(:) :: array

      call rfftb( this%n, array, this%coefs )

    end subroutine simpfftinv

    subroutine doubfftinv(this,array)
      type( fft1dclass) :: this
      double precision, dimension(:) :: array

      call dfftb( this%n, array,  this%coefd )

    end subroutine doubfftinv
  end module fft1d_module
