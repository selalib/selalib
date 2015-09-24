!> Module to compute FFT using FFTPACK library
module fft_module
  implicit none
  !> fft plan for FFTPACK library
  type, public :: fftclass
     real(4), dimension(:), pointer :: coefc  !< simple precision
     real(4), dimension(:), pointer :: work   !< simple precision
     real(4), dimension(:), pointer :: workc  !< simple precision
     real(8), dimension(:), pointer :: coefd  !< double precision
     real(8), dimension(:), pointer :: workd  !< double precision
     real(8), dimension(:), pointer :: coefcd !< double precision
     integer                        :: n      !< number of samples in each sequence
  end type fftclass

  !> Initialize the fftpack plan
  interface initfft
     module procedure initdoubfft,  initdoubcfft
  end interface
  !> Forward fft with fftpack
  interface fft
     module procedure doubfft, doubcfft
  end interface
  !> Inverse fft with fftpack
  interface fftinv
     module procedure doubfftinv,  doubcfftinv
  end interface

  contains

    !> fftpack initialization
    subroutine initdoubfft(this,f,l)
      type(fftclass) :: this !< fft plan
      double precision, dimension(:,:) :: f !< data array
      integer :: l !< array size
      this%n = l 
      allocate(this%coefd(2*this%n+15))
      call dffti(this%n,this%coefd)
    end subroutine initdoubfft

    !> fftpack initialization
    subroutine initdoubcfft(this,f,l)
      type(fftclass) :: this !< fft plan
      double complex, dimension(:,:) :: f !< data array
      integer :: l !< array size
      this%n = l
      allocate(this%coefcd(4*this%n+15))
      call zffti(this%n,this%coefcd)
    end subroutine initdoubcfft

    !> forward fft
    subroutine doubfft(this,array)
      type(fftclass) :: this !< fft plan
      integer :: i,p
      double precision, dimension(:,:) :: array !< data array
      p = size(array,2)   ! number of 1d transforms

      do i=1,p
         call dfftf( this%n, array(:,i), this%coefd)
      end do

      array = array /this%n      ! normalize FFT
    end subroutine doubfft

    !> forward fft
    subroutine doubcfft(this,array)
      type(fftclass) :: this !< fft plan
      integer :: i, p
      double complex, dimension(:,:) :: array !< data array

      p = size(array,2)   ! number of 1d transforms

      do i=1,p
         call zfftf( this%n, array(:,i), this%coefcd)
      end do
      array = array /this%n      ! normalize FFT
    end subroutine doubcfft

    !> inverse fft
    subroutine doubfftinv(this,array)
      type(fftclass) :: this !< fft plan
      integer :: i, p
      double precision, dimension(:,:) :: array !< data array

      p = size(array,2)   ! number of 1d transforms

      do i=1,p
         call dfftb( this%n, array(:,i),  this%coefd )
      end do
    end subroutine doubfftinv

    !> inverse fft
    subroutine doubcfftinv(this,array)
      type(fftclass) :: this !< fft plan
      integer :: i, p
      double complex, dimension(:,:) :: array !< data array

      p = size(array,2)   ! number of 1d transforms

      do i=1,p
         call zfftb( this%n, array(:,i),  this%coefcd )
      end do
    end subroutine doubcfftinv

  end module fft_module
