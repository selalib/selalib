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
! file : fft_NRF90.f90 
! date : 30/04/2001
!  - fast fourier transform
! ( Numerical Recipes Version for fortran 90 )
!
! The basis vector of the Fourier space are stored as
! follow :
!  example for N=8 :
!     ---------------------------------
!     indx  | 1 | 2  3  4  5 |  6  7  8 
!     ---------------------------------
!     modes | 0 | 1  2  3  4 | -3 -2 -1
!     ---------------------------------
!
! ATTENTION : to have a correct Fourier Transform, you 
!  must not take into account the last periodic point
! ex : if the function f(1:nx) is periodic, you have to
!       use the (nx-1) first points, and (nx-1) must
!       be a power of 2
!-------------------------------------------------------
module fft_NRF90_module
  implicit none
      
  integer, parameter :: npar_arth  = 16
  integer, parameter :: npar2_arth = 8
      
  interface assert
     module procedure assert1
  end interface
      
  interface assert_eq
     module procedure assert_eq2
  end interface
      
  interface arth
     module procedure arth_d, arth_i
  end interface
      
  interface swap
     module procedure swap_zv
  end interface
      
  !******************************
  contains
  !******************************
      
  !------------------------------------------------------
  ! numerical recipes : functions extracted of nrutil.f90
  !------------------------------------------------------
  subroutine assert1(n1,string)
    character(len=*), intent(in) :: string
    logical         , intent(in) :: n1
      
    if (.not. n1) then
      write (*,*) 'nrerror: an assertion failed with this tag:', &
        string
      stop 'program terminated by assert1'
    end if
  end subroutine assert1
      
  !------------------------------------------------------
  function assert_eq2(n1,n2,string)
    character(len=*), intent(in) :: string
    integer         , intent(in) :: n1,n2
      
    integer :: assert_eq2
      
    if (n1 == n2) then
      assert_eq2 = n1
    else
      write (*,*) 'nrerror: an assert_eq failed with this tag:', &
        string
      stop 'program terminated by assert_eq2'
    end if
  end function assert_eq2
      
  !------------------------------------------------------
  function arth_d(first,increment,n)
    use prec_const
    real(RKIND), intent(in) :: first, increment
    integer    , intent(in) :: n
      
    real(RKIND), dimension(n) :: arth_d
    integer                   :: k, k2
    real(RKIND)               :: temp
      
    if (n > 0) arth_d(1) = first
    if (n <= npar_arth) then
      do k = 2,n
        arth_d(k) = arth_d(k-1)+increment
      end do
    else
      do k = 2,npar2_arth
        arth_d(k) = arth_d(k-1)+increment
      end do
      temp = increment*npar2_arth
      k    = npar2_arth
      do
        if (k >= n) exit
        k2                    = k+k
        arth_d(k+1:min(k2,n)) = temp+arth_d(1:min(k,n-k))
        temp                  = temp+temp
        k                     = k2
      end do
    end if
  end function arth_d
      
  !------------------------------------------------------
  function arth_i(first,increment,n)
    use prec_const
    integer, intent(in) :: first, increment, n
      
    integer, dimension(n) :: arth_i
    integer               :: k, k2, temp
      
    if (n > 0) arth_i(1) = first
    if (n <= npar_arth) then
      do k = 2,n
        arth_i(k) = arth_i(k-1)+increment
      end do
    else
      do k = 2,npar2_arth
        arth_i(k) = arth_i(k-1)+increment
      end do
      temp = increment*npar2_arth
      k    = npar2_arth
      do
        if (k >= n) exit
        k2                    = k+k
        arth_i(k+1:min(k2,n)) = temp+arth_i(1:min(k,n-k))
        temp                  = temp+temp
        k                     = k2
      end do
    end if
  end function arth_i
      
  !------------------------------------------------------
  subroutine swap_zv(a,b)
    use prec_const
    complex(CKIND), dimension(:), intent(inout) :: a,b
      
    complex(CKIND), dimension(size(a)) :: dum
      
    dum = a
    a   = b
    b   = dum
  end subroutine swap_zv
      
  !------------------------------------------------------
  ! numerical recipes : functions extracted of fourrow.f90
  !------------------------------------------------------
  subroutine fourrow(data,isign)
    use prec_const
    complex(CKIND), dimension(:,:), intent(inout) :: data
    integer                       , intent(in)    :: isign
      
    integer        :: n, i, istep, j, m, mmax, n2
    real(RKIND)    :: theta
    complex(CKIND) :: w, wp
    complex(CKIND) :: ws
    complex(CKIND), dimension(size(data,1)) :: temp
      
    n = size(data,2)
    call assert(iand(n,n-1)==0, &
      'n must be a power of 2 in fourrow_dp')
    n2 = n/2
    j  = n2
    do i = 1,n-2
      if (j > i) call swap(data(:,j+1),data(:,i+1))
      m = n2
      do
        if (m < 2 .or. j < m) exit
        j = j-m
        m = m/2
      end do
      j = j+m
    end do
    mmax = 1
    do
      if (n <= mmax) exit
      istep = 2*mmax
      theta = PI/(isign*mmax)
      wp    = cmplx(-2.0_RKIND * &
        sin(0.5_RKIND*theta)**2,sin(theta),kind=CKIND)
      w     = cmplx(1.0_RKIND,0.0_RKIND,kind=CKIND)
      do m = 1,mmax
        ws = w
        do i = m,n,istep
          j         = i+mmax
          temp      = ws*data(:,j)
          data(:,j) = data(:,i)-temp
          data(:,i) = data(:,i)+temp
        end do
        w = w*wp+w
      end do
      mmax = istep
    end do
  end subroutine fourrow
      
  !--------------------------------------------------------
  ! Numerical Recipes : functions extracted of fourcol.f90
  !--------------------------------------------------------
  subroutine fourcol(data,isign)
    use prec_const
    implicit none
    complex(CKIND), dimension(:,:), intent(inout) :: data
    integer                       , intent(in)    :: isign
      
    integer        :: n, i, istep, j, m, mmax, n2
    real(RKIND)    :: theta
    complex(CKIND) :: w, wp
    complex(CKIND) :: ws
    complex(CKIND), dimension(size(data,2)) :: temp
      
    n = size(data,1)
    call assert(iand(n,n-1)==0, 'n must be a power of 2 in fourcol')
    n2 = n/2
    j  = n2
    do i = 1,n-2
      if (j > i) call swap(data(j+1,:),data(i+1,:))
      m = n2
      do
        if (m < 2 .or. j < m) exit
        j = j-m
        m = m/2
      end do
      j = j+m
    end do
    mmax = 1
    do
      if (n <= mmax) exit
      istep = 2*mmax
      theta = PI/(isign*mmax)
      wp    = cmplx(-2.0_RKIND * &
        sin(0.5_RKIND*theta)**2,sin(theta),kind=CKIND)
      w     = cmplx(1.0_RKIND,0.0_RKIND,kind=CKIND)
      do m = 1,mmax
        ws = w
        do i = m,n,istep
          j         = i+mmax
          temp      = ws*data(j,:)
          data(j,:) = data(i,:)-temp
          data(i,:) = data(i,:)+temp
        end do
        w = w*wp+w
      end do
      mmax = istep
    end do
  end subroutine fourcol
      
  !------------------------------------------------------
  ! numerical recipes : functions extracted of four1.f90
  !------------------------------------------------------
  subroutine four1(data,isign)
    use prec_const
    use mem_alloc_module
    complex(CKIND), dimension(:), intent(inout) :: data
    integer                     , intent(in)    :: isign
      
    complex(CKIND), dimension(:,:), pointer :: dat, temp
    complex(CKIND), dimension(:)  , pointer :: w, wp
    real(RKIND)   , dimension(:)  , pointer :: theta
    integer                                 :: n, m1, m2, j
      
    n = size(data)
    call assert(iand(n,n-1)==0, &
      'n must be a power of 2 in four1_dp')
    m1 = 2**ceiling(0.5_RKIND*log(real(n,RKIND))/0.693147_RKIND)
    m2 = n/m1
      
    call temp_allocate(dat,1,m1,1,m2,'dat')
    call temp_allocate(theta,1,m1,'theta')
    call temp_allocate(w,1,m1,'w')
    call temp_allocate(wp,1,m1,'wp')
    call temp_allocate(temp,1,m2,1,m1,'temp')
      
    dat = reshape(data,shape(dat))
    call fourrow(dat,isign)
    theta = arth(0,isign,m1)*TWOPI/n
    wp    = cmplx(-2.0_RKIND * &
      sin(0.5_RKIND*theta)**2,sin(theta),kind=CKIND)
    w     = cmplx(1.0_RKIND,0.0_RKIND,kind=CKIND)
    do j = 2,m2
      w        = w*wp+w
      dat(:,j) = dat(:,j)*w
    end do
    temp = transpose(dat)
    call fourrow(temp,isign)
    data = reshape(temp,shape(data))
    call temp_deallocate(dat)
    call temp_deallocate(w)
    call temp_deallocate(wp)
    call temp_deallocate(theta)
    call temp_deallocate(temp)
  end subroutine four1
end module fft_NRF90_module
