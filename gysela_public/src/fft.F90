!**************************************************************
!  Copyright Euratom-CEA
!  Authors : 
!     Virginie Grandgirard (virginie.grandgirard@cea.fr)
!     Chantal Passeron (chantal.passeron@cea.fr)
!     Guillaume Latu (guillaume.latu@cea.fr)
!     Xavier Garbet (xavier.garbet@cea.fr)
!     Philippe Ghendrih (philippe.ghendrih@cea.fr)
!     Yanick Sarazin (yanick.sarazin@cea.fr)
!  
!  This code GYSELA (for GYrokinetic SEmi-LAgrangian) 
!  is a 5D gyrokinetic global full-f code for simulating 
!  the plasma turbulence in a tokamak.
!  
!  This software is governed by the CeCILL-B license 
!  under French law and abiding by the rules of distribution 
!  of free software.  You can  use, modify and redistribute 
!  the software under the terms of the CeCILL-B license as 
!  circulated by CEA, CNRS and INRIA at the following URL
!  "http://www.cecill.info". 
!**************************************************************
      
!-------------------------------------------------------
! file : fft.f90 
! date : 30/04/2001
!  - fast fourier transform using :
!   1) Numerical Recipes Version for fortran 90
!   2) fft4f2d
!
! ATTENTION : to have a correct Fourier Transform, you 
!  must not take into account the last periodic point
! ex : if the function f(1:nx) is periodic, you have to
!       use the (nx-1) first points, and (nx-1) must
!       be a power of 2
!-------------------------------------------------------
module fft_module
  use prec_const
  use fft_NRF90_module
      
  implicit none
  public :: NR_fft1_2D, NR_fft1_3D, NR_fft2_3D
      
  !******************************
  contains
  !******************************
      
  !****************************************************
  !  FFT using Numerical Recipes
  !****************************************************
  !-------------------------------------------------------------
  ! Computation of the wave vector of the Fourier space 
  !  associated to the storage of the FFT basis vector
  !  in Numerical Recipes. 
  !  example for N=8 :
  !     ---------------------------------
  !     indx  | 1 | 2  3  4  5 |  6  7  8 
  !     ---------------------------------
  !     modes | 0 | 1  2  3  4 | -3 -2 -1
  !     ---------------------------------
  ! 
  ! Arguments:
  ! ==========
  !  . Nx   : (input) INTEGER 
  !             Number of Fourier modes
  !  . dx   : (input) REAL(RKIND)
  !             discretisation step
  !  . kfft : (output) REAL(RKIND) array, dimension(N)
  !             wave vector of the Fourier space
  !-------------------------------------------------------------
  subroutine NR_fft_vector(Nx,dx,kfft)
    integer                     , intent(in)  :: Nx
    real(RKIND)                 , intent(in)  :: dx
    real(RKIND), dimension(1:Nx), intent(out) :: kfft
      
    integer     :: ix
    real(RKIND) :: wkx
    
    !-> wave vector
    wkx = TWOPI/(real(Nx)*dx)
    do ix = 1,(Nx/2)+1
      kfft(ix) = (ix-1) * wkx
    end do
    do ix = (Nx/2)+2,Nx
      kfft(ix) = -1._RKIND * (Nx-(ix-1)) * wkx
    end do
  end subroutine NR_fft_vector
      
  !-------------------------------------------------------------
  ! Computation of a FFT in 1D (with complex function)
  !  input/output array = 2D array.
  ! the FFT 1D is done on the second indexes
  !  with the FFT routines of Numerical Recipes
  !  (by using fourrow)
  ! -> isign = 1 for a FFT and isign = -1 for an inverse FFT
  !   . n2  = size(data,2)
  !-------------------------------------------------------------
  subroutine NR_fft1_2D(data,n2,isign)
    use prec_const
    use mem_alloc_module
    complex(CKIND), dimension(:,:), intent(inout) :: data
    integer                       , intent(in)    :: n2
    integer                       , intent(in)    :: isign
      
    integer                         :: i, n1
    real(RKIND)                     :: rem
    complex(CKIND), dimension(1:n2) :: data_i
      
    n1  = size(data,1)
    rem = modulo(log(real(n2)),real(log(2._RKIND)))
    if (rem.ne.0._RKIND) then
      write(6,*) 'Warning in NR_fft1_2D : n2 = ', &
        n2, ' is not a power of 2'
      stop
    end if
 
    do i = 1,n1
      data_i = data(i,:)
      call four1(data_i,isign)
      data(i,:) = data_i
    end do
    if (isign.eq.-1) data = data/n2
  end subroutine NR_fft1_2D
      
  !-------------------------------------------------------------
  ! Computation of a FFT in 1D (with complex function)
  !  input/output array = 3D array.
  ! the FFT 1D is done on the second indexes
  !  with the FFT routines of Numerical Recipes 
  !  (by using fourrow)
  ! -> isign = 1 for a FFT and isign = -1 for an inverse FFT
  !   . n2  = size(data,2)
  !-------------------------------------------------------------
  subroutine NR_fft1_3D(data,n2,isign)
    use prec_const
    use mem_alloc_module
    complex(CKIND), dimension(:,:,:), intent(inout) :: data
    integer                         , intent(in)    :: n2
    integer                         , intent(in)    :: isign
      
    integer                         :: i, k, n1, n3
    real(RKIND)                     :: rem
    complex(CKIND), dimension(1:n2) :: data_i
      
    n1  = size(data,1)
    n3  = size(data,3)
    rem = modulo(log(real(n2)),real(log(2._RKIND)))
    if (rem.ne.0._RKIND) then
      write(6,*) 'Warning in NR_fft1_3D : n2 = ', &
        n2, ' is not a power of 2'
      stop
    end if
 
    do k = 1,n3
      do i = 1,n1
        data_i = data(i,:,k)
        call four1(data_i,isign)
        data(i,:,k) = data_i
      end do
    end do
    if (isign.eq.-1) data = data/n2
  end subroutine NR_fft1_3D
      
  !-------------------------------------------------------------
  ! Computation of a FFT in 2D (with complex function)
  !  input/output array = 3D array.
  ! the FFT 2D is done on the second and third indexes
  !  with the FFT routines of Numerical Recipes 
  !  (by using fourrow)
  ! -> isign = 1 for a FFT and isign = -1 for an inverse FFT
  !   . n2 = size(data,2)
  !   . n3 = size(data,3)
  !-------------------------------------------------------------
 subroutine NR_fft2_3D(data,n2,n3,isign)
    use prec_const
    use mem_alloc_module
    complex(CKIND), dimension(:,:,:), intent(inout) :: data
    integer                         , intent(in)    :: n2, n3
    integer                         , intent(in)    :: isign
      
    complex(CKIND), dimension(1:n2,1:n3) :: data_i
    complex(CKIND), dimension(1:n3,1:n2) :: temp
    integer :: i, n1
      
    n1 = size(data,1)
 
    do i = 1,n1 
      data_i = data(i,:,:)
      call fourrow(data_i,isign)
      temp = transpose(data_i)
      call fourrow(temp,isign)
      data_i = transpose(temp)
      data(i,:,:) = data_i
    end do
    if (isign.eq.-1) data = data/(n2*n3)
  end subroutine NR_fft2_3D
end module fft_module
