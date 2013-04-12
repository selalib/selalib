module fft_module
#include "sll_working_precision.h"
  implicit none
  private
  type, public :: fftclass
     real, dimension(:), pointer ::  coefc, work, workc
     sll_real64, dimension(:), pointer :: coefd, workd, coefcd
     sll_int32  :: n  ! number of samples in each sequence
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

  public :: initfft, fft, fftinv
  contains

    subroutine initdoubfft(this,f,l)
      type(fftclass) :: this
      sll_real64, dimension(:,:) :: f
      sll_int32 :: l 
      this%n = l 
      allocate(this%coefd(2*this%n+15))
      call dffti(this%n,this%coefd)
    end subroutine initdoubfft


    subroutine initdoubcfft(this,f,l)
      type(fftclass) :: this
      sll_comp64, dimension(:,:) :: f
      sll_int32 :: l 
      this%n = l
      allocate(this%coefcd(4*this%n+15))
      call zffti(this%n,this%coefcd)
    end subroutine initdoubcfft

    subroutine doubfft(this,array)
      type(fftclass) :: this
      sll_int32, parameter :: sign = -1   ! we choose this for direct transform
      sll_int32 :: i,j,p, inc, lda
      sll_real64, dimension(:,:) :: array
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
      sll_int32, parameter :: sign = -1   ! we choose this for direct transform
      sll_int32 :: i,j, p, inc, lda
      sll_comp64, dimension(:,:) :: array

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
      sll_int32, parameter :: sign = 1   ! we choose this for inverse transform
      sll_int32 :: i,j, p, inc, lda
      sll_real64, dimension(:,:) :: array

      sll_real64, dimension(size(array,2),size(array,1)) :: DX

      p = size(array,2)   ! number of 1d transforms
      inc = 1             ! all data are samples
      lda = size(array,1) ! leading dimension of array

      do i=1,p
         call dfftb( this%n, array(:,i),  this%coefd )
      end do
    end subroutine doubfftinv

    subroutine doubcfftinv(this,array)
      type(fftclass) :: this
      sll_int32, parameter :: sign = 1   ! we choose this for inverse transform
      sll_int32 :: i, p, inc, lda
      sll_comp64, dimension(:,:) :: array

      p = size(array,2)   ! number of 1d transforms
      inc = 1             ! all data are samples
      lda = size(array,1) ! leading dimension of array

      do i=1,p
         call zfftb( this%n, array(:,i),  this%coefcd )
      end do
    end subroutine doubcfftinv

  end module fft_module
