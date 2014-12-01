module fft_module
  implicit none

  type, public :: fftclass
     real, dimension(:), pointer ::  coefc, work, workc
     double precision, dimension(:), pointer :: coefd, workd, coefcd
     integer  :: n  ! number of samples in each sequence
  end type fftclass

  interface initfft
     module procedure initdoubfft,  initdoubcfft
  end interface
  interface fft
     module procedure doubfft, doubcfft
  end interface
  interface fftinv
     module procedure doubfftinv,  doubcfftinv
  end interface

  contains

    subroutine initdoubfft(this,f,l)
      type(fftclass) :: this
      double precision, dimension(:,:) :: f
      integer :: l 
      this%n = l 
      allocate(this%coefd(2*this%n+15))
      call dffti(this%n,this%coefd)
    end subroutine initdoubfft


    subroutine initdoubcfft(this,f,l)
      type(fftclass) :: this
      double complex, dimension(:,:) :: f
      integer :: l 
      this%n = l
      allocate(this%coefcd(4*this%n+15))
      call zffti(this%n,this%coefcd)
    end subroutine initdoubcfft

    subroutine doubfft(this,array)
      type(fftclass) :: this
      integer :: i,p
      double precision, dimension(:,:) :: array
      p = size(array,2)   ! number of 1d transforms

      do i=1,p
         call dfftf( this%n, array(:,i), this%coefd)
      end do

      array = array /this%n      ! normalize FFT
    end subroutine doubfft

    subroutine doubcfft(this,array)
      type(fftclass) :: this
      integer :: i, p
      double complex, dimension(:,:) :: array

      p = size(array,2)   ! number of 1d transforms

      do i=1,p
         call zfftf( this%n, array(:,i), this%coefcd)
      end do
      array = array /this%n      ! normalize FFT
    end subroutine doubcfft

    subroutine doubfftinv(this,array)
      type(fftclass) :: this
      integer :: i, p
      double precision, dimension(:,:) :: array

      p = size(array,2)   ! number of 1d transforms

      do i=1,p
         call dfftb( this%n, array(:,i),  this%coefd )
      end do
    end subroutine doubfftinv

    subroutine doubcfftinv(this,array)
      type(fftclass) :: this
      integer :: i, p
      double complex, dimension(:,:) :: array

      p = size(array,2)   ! number of 1d transforms

      do i=1,p
         call zfftb( this%n, array(:,i),  this%coefcd )
      end do
    end subroutine doubcfftinv

  end module fft_module
