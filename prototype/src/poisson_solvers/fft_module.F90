module fft_module
  implicit none
  public :: initdfft, initcfft, fft, fftinv
  type, public :: fftclass
     real, dimension(:), pointer ::  coefc, work, workc
     double precision, dimension(:), pointer :: coefd, workd, coefcd
     integer  :: n  ! number of samples in each sequence
  end type fftclass
  interface initdfft
     module procedure initdoubfft
  end interface
  interface initcfft
     module procedure initdoubcfft
  end interface
  interface fft
     module procedure doubfft, doubcfft
  end interface
  interface fftinv
     module procedure doubfftinv,  doubcfftinv
  end interface

  contains

    subroutine initdoubfft(this,l)
      type(fftclass) :: this
      integer :: l 
      this%n = l 
      allocate(this%coefd(2*this%n+15))
      call dffti(this%n,this%coefd)
    end subroutine initdoubfft

    subroutine initdoubcfft(this,l)
      type(fftclass) :: this
      integer :: l 
      this%n = l
      allocate(this%coefcd(4*this%n+15))
      call zffti(this%n,this%coefcd)
    end subroutine initdoubcfft

    subroutine doubfft(this,array)
      type(fftclass) :: this
      integer :: i,j,p, inc, lda
      double precision, dimension(:,:) :: array

      p = size(array,2)   ! number of 1d transforms
      inc = 1         ! all data are samples
      lda = size(array,1) ! leading dimension of array

      do i=1,p
         call dfftf( this%n, array(:,i), this%coefd)
      end do
      do j= 1,p
         do i=1,lda
            array(i,j) = array(i,j) /this%n      ! normalize FFT
         end do
      end do     
    end subroutine doubfft

    subroutine doubcfft(this,array)
      type(fftclass) :: this
      integer :: i, p, inc, lda
      complex(8), dimension(:,:) :: array

      p = size(array,2)   ! number of 1d transforms
      inc = 1             ! all data are samples
      lda = size(array,1) ! leading dimension of array

      do i=1,p
         call zfftf( this%n, array(:,i), this%coefcd)
      end do
      array = array /this%n      ! normalize FFT
    end subroutine doubcfft

    subroutine doubfftinv(this,array)
      type(fftclass) :: this
      integer :: i, p
      double precision, dimension(:,:) :: array

      do i=1,p
         call dfftb( this%n, array(:,i),  this%coefd )
      end do
    end subroutine doubfftinv

    subroutine doubcfftinv(this,array)
      type(fftclass) :: this
      integer :: i, p, inc, lda
      complex(8), dimension(:,:) :: array

      p = size(array,2)   ! number of 1d transforms
      inc = 1             ! all data are samples
      lda = size(array,1) ! leading dimension of array

      do i=1,p
         call zfftb( this%n, array(:,i),  this%coefcd )
      end do
    end subroutine doubcfftinv

  end module fft_module
