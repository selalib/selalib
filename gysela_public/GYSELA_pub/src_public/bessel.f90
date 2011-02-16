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
      
module bessel_module
  use prec_const
  use mem_alloc_module
      
  implicit none
      
  integer, parameter      :: I4B = SELECTED_INT_KIND(9)
  integer(I4B), parameter :: NPAR_POLY=8
      
  interface assert
     module procedure assert1 
  end interface
      
  interface poly
     module procedure poly_dd
  end interface
      
  interface bessj0
     module procedure bessj0_s
  end interface
      
  interface bessj1
     module procedure bessj1_s
  end interface
      
  interface bessj
     module procedure bessj_s
  end interface
      
  !******************************
  contains
  !******************************
      
  !---------------------------------------------
  !  function in nrutil.f90 (Numerical Recipes)
  !---------------------------------------------
  subroutine nrerror(string)
    character(len=*), intent(in) :: string
    write (*,*) 'nrerror: ',string
    stop 'program terminated by nrerror'
  end subroutine nrerror
      
  !-----------------------------------------------
  !  subroutine in nrutil.f90 (Numerical Recipes)
  !-----------------------------------------------
  subroutine assert1(n1,string)
    character(len=*), intent(in) :: string
    logical, intent(in) :: n1
    if (.not. n1) then
      write (*,*) 'nrerror: an assertion failed with this tag:', &
        string
      stop 'program terminated by assert1'
    end if
  end subroutine assert1
    
      
  !---------------------------------------------
  !  function in nrutil.f90 (Numerical Recipes)
  !---------------------------------------------
  function poly_dd(x,coeffs)
    real(RKIND)              , intent(in) :: x
    real(RKIND), dimension(:), intent(in) :: coeffs
      
    real(RKIND)                        :: poly_dd
    real(RKIND)                        :: pow
    real(RKIND), dimension(:), pointer :: vec
    integer(I4B)                       :: i,n,nn
      
    n=size(coeffs)
    if (n <= 0) then
      poly_dd = 0.0_RKIND
    else if (n < NPAR_POLY) then
      poly_dd = coeffs(n)
      do i = n-1,1,-1
        poly_dd = x*poly_dd+coeffs(i)
      end do
    else
      call temp_allocate(vec,1,n+1,'vec')
      pow      = x
      vec(1:n) = coeffs
      do
        vec(n+1)  = 0.0_RKIND
        nn        = ishft(n+1,-1)
        vec(1:nn) = vec(1:n:2)+pow*vec(2:n+1:2)
        if (nn == 1) exit
        pow = pow*pow
        n   = nn
      end do
      poly_dd = vec(1)
      call temp_deallocate(vec)
    end if
  end function poly_dd
      
  !--------------------------------------------
  !  function in bessj0.f90 (Numerical Recipes) 
  !   for the computation of j0(x) 
  !--------------------------------------------
  function bessj0_s(x)
    implicit none
    real(RKIND), intent(in) :: x
    real(RKIND)             :: bessj0_s
      
    real(RKIND) :: ax,xx,z
    real(RKIND) :: y
    real(RKIND), dimension(5) :: &
      p = (/1.0_RKIND,-0.1098628627e-2_RKIND, &
      0.2734510407e-4_RKIND,-0.2073370639e-5_RKIND, &
      0.2093887211e-6_RKIND/)
    real(RKIND), dimension(5) :: &
      q = (/-0.1562499995e-1_RKIND, &
      0.1430488765e-3_RKIND,-0.6911147651e-5_RKIND, &
      0.7621095161e-6_RKIND,&
      -0.934945152e-7_RKIND/)
    real(RKIND), dimension(6) :: &
      r = (/57568490574.0_RKIND,-13362590354.0_RKIND, &
      651619640.7_RKIND,-11214424.18_RKIND,77392.33017_RKIND, &
      -184.9052456_RKIND/)
    real(RKIND), dimension(6) :: &
      s = (/57568490411.0_RKIND,1029532985.0_RKIND, &
      9494680.718_RKIND,59272.64853_RKIND, &
      267.8532712_RKIND,1.0_RKIND/)
      
    if (abs(x) < 8.0) then
      y        = x**2
      bessj0_s = poly(y,r)/poly(y,s)
    else
      ax        = abs(x)
      z         = 8.0_RKIND/ax
      y         = z**2
      xx        = ax-0.785398164_RKIND
      bessj0_s  = sqrt(0.636619772_RKIND/ax)*(cos(xx)* &
        poly(y,p)-z*sin(xx)*poly(y,q))
    end if
  end function bessj0_s
      
  !-------------------------------------------
  !  function in bessj1.f90 (Numerical Recipes) 
  !    for the computation of j1(x) 
  !-------------------------------------------
  function bessj1_s(x)
    implicit none
    real(RKIND), intent(in) :: x
    real(RKIND)             :: bessj1_s
      
    real(RKIND) :: ax,xx,z
    real(RKIND) :: y
    real(RKIND), dimension(6) :: &
      r = (/72362614232.0_RKIND, &
      -7895059235.0_RKIND,242396853.1_RKIND,-2972611.439_RKIND, &
      15704.48260_RKIND,-30.16036606_RKIND/)
    real(RKIND), dimension(6) :: &
      s = (/144725228442.0_RKIND,2300535178.0_RKIND, &
      18583304.74_RKIND,99447.43394_RKIND, &
      376.9991397_RKIND,1.0_RKIND/)
    real(RKIND), dimension(5) :: &
      p = (/1.0_RKIND,0.183105e-2_RKIND, &
      -0.3516396496e-4_RKIND,0.2457520174e-5_RKIND, &
      -0.240337019e-6_RKIND/)
    real(RKIND), dimension(5) :: &
      q = (/0.04687499995_RKIND, &
      -0.2002690873e-3_RKIND,0.8449199096e-5_RKIND, &
      -0.88228987e-6_RKIND,0.105787412e-6_RKIND/)
      
    if (abs(x) < 8.0) then
      y        = x**2
      bessj1_s = x*(poly(y,r)/poly(y,s))
    else
      ax       = abs(x)
      z        = 8.0_RKIND/ax
      y        = z**2
      xx       = ax-2.356194491_RKIND
      bessj1_s = sqrt(0.636619772_RKIND/ax)*(cos(xx)* &
        poly(y,p)-z*sin(xx)*poly(y,q))*sign(1.0_RKIND,x)
    end if
  end function bessj1_s
  
      
  !----------------------------------------------------
  !  function in bessj0.f90 (Numerical Recipes) 
  !   for the computation of jn(x) with n>=2
  !---------------------------------------------------
  function bessj_s(n,x)
    implicit none
    integer(I4B), intent(in) :: n
    real(RKIND) , intent(in) :: x
    real(RKIND)              :: bessj_s
      
    integer(I4B), parameter :: iacc = 40
    integer(I4B), parameter :: iexp = maxexponent(x)/2
    integer(I4B)            :: j, jsum, m
    real(RKIND)             :: ax, bj, bjm, bjp, summ, tox
      
    call assert(n >= 2, 'bessj_s args')
    ax = abs(x)
    if (ax*ax <= 8.0_RKIND*tiny(x)) then
      bessj_s = 0.0_RKIND
    else if (ax > real(n,RKIND)) then
      tox = 2.0_RKIND/ax
      bjm = bessj0(ax)
      bj  = bessj1(ax)
      do j = 1,n-1
        bjp = j*tox*bj-bjm
        bjm = bj
        bj  = bjp
      end do
      bessj_s = bj
    else
      tox     = 2.0_RKIND/ax
      m       = 2*((n+int(sqrt(real(iacc*n,RKIND))))/2)
      bessj_s = 0.0_RKIND
      jsum    = 0
      summ    = 0.0_RKIND
      bjp     = 0.0_RKIND
      bj      = 1.0_RKIND
      do j = m,1,-1
        bjm = j*tox*bj-bjp
        bjp = bj
        bj  = bjm
        if (exponent(bj) > iexp) then
          bj      = scale(bj,-iexp)
          bjp     = scale(bjp,-iexp)
          bessj_s = scale(bessj_s,-iexp)
          summ    = scale(summ,-iexp)
        end if
        if (jsum /= 0) summ = summ+bj
        jsum = 1-jsum
        if (j == n) bessj_s=bjp
      end do
      summ    = 2.0_RKIND*summ-bj
      bessj_s = bessj_s/summ
    end if
    if (x < 0.0 .and. mod(n,2) == 1) bessj_s =- bessj_s
  end function bessj_s
      
  !----------------------------------------------------
  !  function in bessel.f90 
  !   for the computation of jm(x) with m>=2
  !---------------------------------------------------
  function bessel(m,x)
    implicit none
    integer(I4B), intent(in) :: m
    real(RKIND) , intent(in) :: x
    real(RKIND)              :: bessel
    
    select case(m)
      case(0) 
        bessel = bessj0_s(x)
      case(1)
        bessel = bessj1_s(x)
      case default
        bessel = bessj_s(m,x)
    end select
  end function bessel
      
  !--------------------------------------------------------
  ! The same kind of function than the function 'zbrent'
  !  (from page 253 of Numerical Recipes) which
  !  finds the root known to lie between two numbers 
  !  x1 and x2, i.e f(x1)<0<f(x2)
  !--------------------------------------------------------
  function zbrent_Jmx(m,x1,x2,tol)
    use prec_const
    implicit none
    integer    , intent(in) :: m
    real(RKIND), intent(in) :: x1, x2, tol
    real(RKIND)             :: zbrent_Jmx
    
    integer    , parameter :: itmax=100
    real(RKIND), parameter :: eps=epsilon(x1)
    integer                :: iter
    real(RKIND)            :: a, b, c, d, e
    real(RKIND)            :: fa, fb, fc
    real(RKIND)            :: p, q, r, s, tol1, xm
    
    a  = x1
    b  = x2
    fa = bessel(m,a)
    fb = bessel(m,b)
    if ((fa > 0.0 .and. fb > 0.0) .or. (fa < 0.0 .and. fb < 0.0)) &
      call nrerror('root must be bracketed for zbrent')
    c  = b
    fc = fb
    do iter = 1,itmax
      if ((fb > 0.0 .and. fc > 0.0) .or. &
        (fb < 0.0 .and. fc < 0.0)) then
        c  = a
        fc = fa
        d  = b-a
        e  = d
      end if
      if (abs(fc) < abs(fb)) then
        a  = b
        b  = c
        c  = a
        fa = fb
        fb = fc
        fc = fa
      end if
      tol1 = 2.0_RKIND*eps*abs(b)+0.5_RKIND*tol
      xm   = 0.5_RKIND*(c-b)
      if (abs(xm) <= tol1 .or. fb == 0.0) then
        zbrent_Jmx = b
        return
      end if
      if (abs(e) >= tol1 .and. abs(fa) > abs(fb)) then
        s = fb/fa
        if (a == c) then
          p = 2.0_RKIND*xm*s
          q = 1.0_RKIND-s
        else
          q = fa/fc
          r = fb/fc
          p = s*(2.0_RKIND*xm*q*(q-r)-(b-a)*(r-1.0_RKIND))
          q = (q-1.0_RKIND)*(r-1.0_RKIND)*(s-1.0_RKIND)
        end if
        if (p > 0.0) q = -q
        p = abs(p)
        if (2.0_RKIND*p<min(3.0_RKIND*xm*q-abs(tol1*q),abs(e*q))) &
          then
          e = d
          d = p/q
        else
          d = xm
          e = d
        end if
      else
        d = xm
        e = d
      end if
      a  = b
      fa = fb
      b  = b + merge(d,sign(tol1,xm), abs(d) > tol1 )
      fb = bessel(m,b)
    end do
    call nrerror('zbrent: exceeded maximum iterations')
    zbrent_Jmx = b  
  end function zbrent_Jmx
      
  !------------------------------------------
  ! computes the zeros of Jm(x)
  !------------------------------------------
  subroutine zeros_bessel_Jmx(m,nb_zeros,Jmx_zeros)
    use prec_const
    implicit none
    integer                           , intent(in)  :: m
    integer                           , intent(in)  :: nb_zeros
    real(RKIND), dimension(1:nb_zeros), intent(out) :: Jmx_zeros
    
    integer     :: izero
    real(RKIND) :: del, tol
    real(RKIND) :: x, x1, x2, f1, f2, root
    
    !*** Stepping increment ***
    del = 0.1
    !*** Tolerance for the zero finder ***
    tol = 1e-9
    x   = del
    do izero = 1,nb_zeros
      x1 = x
      x2 = x+del
      f1 = bessel(m,x1)
      f2 = bessel(m,x2)
      do while ((f1*f2).gt.0._RKIND)
        x1 = x
        x2 = x + del
        f1 = bessel(m,x1)
        f2 = bessel(m,x2)
        x  = x2
      end do
      root             = zbrent_Jmx(m,x1,x2,tol)
      Jmx_zeros(izero) = root
      x                = x2
    end do
  end subroutine zeros_bessel_Jmx
      
  !------------------------------------------
  ! computes the maximum of abs(Jm(x))
  !------------------------------------------
  function max_Jmx(m,xgrid)
    use prec_const
    implicit none
    integer                   , intent(in) :: m
    real(RKIND), dimension(0:), intent(in) :: xgrid
      
    integer     :: Nx, ix
    real(RKIND) :: max_tmp, max_Jmx
      
    Nx      = size(xgrid)
    max_tmp = 0._RKIND
    do ix = 0,Nx-1
      max_tmp = max(max_tmp,bessel(m,xgrid(ix)))
    end do
    max_Jmx = max_tmp
  end function max_Jmx
      
  !*******************************************************
  !            SPECIAL FUNCTIONS
  !*******************************************************
  !***********************************************************
  !   Adapted from a Fortran 77 routine
  !   Computes the modified Bessel functions I0(z) and K0(z)
  !     for a complex argument
  !
  !   Input :  z --- Complex argument
  !   Output:  cbess_I0 --- I0(z)
  !            cbess_K0 --- K0(z)
  !**********************************************************
  subroutine compute_modified_bessel(Z,cbess_I0,cbess_K0)
    implicit none
      
    complex(CKIND), intent(in)   :: Z
    complex(CKIND), intent(out)  :: cbess_I0, cbess_K0
    complex(CKIND)               :: Z1, Z2, CI, CA, CS, CT, CR, CW
    complex(CKIND)               :: tmp, ZR, ZR2
    real(RKIND), dimension(1:12) :: A, B
    real(RKIND), dimension(1:10) :: A1
    real(RKIND)                  :: A0, W0
    integer                      :: k, k0
      
    A0 = abs(Z)
    Z2 = Z*Z
    Z1 = Z
    CI = (0._RKIND,1._RKIND)
    if ( A0.eq.0._RKIND ) then
      cbess_I0 = (1._RKIND,0._RKIND)
      cbess_K0 = (1.e300_RKIND,0._RKIND)
    end if
    if ( real(Z) < 0._RKIND ) then
      Z1 = -Z
    end if
    if ( A0 <= 18._RKIND ) then
      cbess_I0 = (1._RKIND,0._RKIND)
      CR       = (1._RKIND,0._RKIND)
      do k =1,50
        CR       = 0.25_RKIND*CR*Z2/(k*k)
        cbess_I0 = cbess_I0+CR
      end do
    else
      A(1:12) = (/0.125_RKIND,7.03125e-2_RKIND, &
        7.32421875e-2_RKIND,1.1215209960938e-1_RKIND, &
        2.2710800170898e-1_RKIND,5.7250142097473e-1_RKIND, &
        1.7277275025845_RKIND,6.0740420012735_RKIND, &
        2.4380529699556e1_RKIND,1.1001714026925e2_RKIND, &
        5.5133589612202e2_RKIND,3.0380905109224e3_RKIND/)
      B(1:12) = (/-0.375_RKIND,-1.171875e-1_RKIND, &
        -1.025390625e-1_RKIND,-1.4419555664063e-1_RKIND, &
        -2.7757644653320e-1_RKIND,-6.7659258842468e-1_RKIND, &
        -1.9935317337513_RKIND,-6.8839142681099_RKIND, &
        -2.7248827311269e1_RKIND,-1.2159789187654e2_RKIND, &
        -6.0384407670507e2_RKIND,-3.3022722944809e3_RKIND/)
      k0 = 12
      if ( A0 >= 35._RKIND ) k0 = 9
      if ( A0 >= 50._RKIND ) k0 = 7
      CA       = exp(Z1)/sqrt(2._RKIND*PI*Z1)
      cbess_I0 = (1._RKIND,0._RKIND)
      ZR       = 1._RKIND/Z1
      do k = 1,k0
        cbess_I0 = cbess_I0+A(k)*ZR**k
      end do
        cbess_I0 = CA*cbess_I0
    end if
    if ( A0 <= 9._RKIND ) then
      CS = (0._RKIND,0._RKIND)
      CT = -log(0.5_RKIND*Z1)-0.5772156649015329_RKIND
      W0 = 0._RKIND
      CR = (1._RKIND,0._RKIND)
      CW = (0._RKIND,0._RKIND)
      do k = 1,50
        W0 = W0+1._RKIND/k
        CR = 0.25_RKIND*CR/(k*k)*Z2
        CS = CS+CR*(W0+CT)
        if ( abs((CS-CW)/CS) < 1.e-15_RKIND ) then
          cbess_K0 = CT+CS
        else        
          CW = CS
        end if
      end do
    else
      A1(1:10) = (/0.125_RKIND,0.2109375_RKIND, &
        1.0986328125_RKIND,1.1775970458984e1_RKIND, &
        2.1461706161499e2_RKIND,5.9511522710323e3_RKIND, &
        2.3347645606175e5_RKIND,1.2312234987631e7_RKIND, &
        8.401390346421e8_RKIND,7.2031420482627e10_RKIND/)
      tmp      = 0.5_RKIND/Z1
      ZR2      = 1._RKIND/Z2
      cbess_K0 =(1._RKIND,0._RKIND)
      do k = 1,10
        cbess_K0 = cbess_K0+A1(k)*ZR2**k
      end do
      cbess_K0 = tmp*cbess_K0/cbess_I0
    end if
    if ( real(Z) < 0._RKIND ) then
      if ( aimag(Z) < 0._RKIND ) cbess_K0 = cbess_K0+CI*PI*cbess_I0
      if ( aimag(Z) > 0._RKIND ) cbess_K0 = cbess_K0-CI*PI*cbess_I0
    end if
  end subroutine compute_modified_bessel
      
  !**************************************************************
  !   Adapted from C++ routines for FORTRAN 90
  !   Purpose: Compute modified Bessel function I0(x), called by
  !            BesselI[0, x] in Mathematica
  !      Input :  x   ---  Argument REAL
  !      Output:  fx  ---  I0(x)
  !**************************************************************
  subroutine compute_bessel_I0_Real(x,fx)
    implicit none
      
    real(RKIND)   , intent(in)  :: x
    complex(CKIND), intent(out) :: fx
    real(RKIND)                 :: y
      
    if ( x.eq.0._RKIND ) then
      fx = 1._RKIND
    else if ( x < 3.75_RKIND ) then
      y  = x / 3.75_RKIND
      fx = 1._RKIND + &
        y*(3.5156229_RKIND+y*(3.0899424_RKIND+y*(1.2067492_RKIND+ &
        y*(0.2659732_RKIND+y*(0.0360768_RKIND+y*0.0045813_RKIND)))))
    else
      y  = 3.75_RKIND / x
      fx = (exp(x)/sqrt(x)) * &
        (0.39894228_RKIND+y*(0.01328592_RKIND+y*(0.00225319_RKIND+ &
        y*(-0.00157565_RKIND+y*(0.00916281_RKIND) + &
        y*(-0.02057706_RKIND + &
        y*(0.02635537_RKIND + &
        y*(-0.01647633_RKIND+y*0.00392377_RKIND)))))))
    end if
  end subroutine compute_bessel_I0_Real
      
  !**************************************************************
  !   Adapted from C++ routines for FORTRAN 90
  !   Purpose: Compute modified Bessel function K0(x), called by
  !            BesselK[0, x] in Mathematica
  !      Input :  x   ---  Argument REAL AND POSITIVE
  !      Output:  fx  ---  K0(x)
  !**************************************************************
  subroutine compute_bessel_K0_Real(x,fx)
    implicit none
      
    real(RKIND)   , intent(in)  :: x
    complex(CKIND), intent(out) :: fx
    real(RKIND)                 :: y, tempI0r
    complex(CKIND)              :: tempI0c
      
    if ( x.eq.0._RKIND ) then
      fx = (1.e300_RKIND,0._RKIND)
    else if ( x <= 2._RKIND ) then
      call compute_bessel_I0_Real(x,tempI0c)
      tempI0r = dble(tempI0c)
      y       = x*x*0.25_RKIND
      fx      = (-log(x*0.5_RKIND)*tempI0r) + &
        (-0.5721566_RKIND+y*(0.42278420_RKIND+ &
        y*(0.23069756+y*(0.03488590_RKIND+y*(0.00262698_RKIND+ &
        y*(0.00010750_RKIND+y*0.74e-5_RKIND))))))
    else
      y  = 2._RKIND / x
      fx = (exp(-x)/sqrt(x)) * &
        (1.25331414_RKIND+y*(-0.07832358_RKIND+ &
        y*(0.02189568_RKIND + &
        y*(-0.01062446_RKIND+y*(0.00587872_RKIND+ &
        y*(-0.00251540_RKIND+y*0.53208e-3_RKIND))))))
    end if
  end subroutine compute_bessel_K0_Real
      
end module bessel_module
