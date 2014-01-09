!------------------------------------------------------------------
! computes the Fried and Conte plasma dispersion function
! Z(z) = i * sqrt(pi) * w(z) = i * sqrt(pi) * exp(-z**2)*erfc(-i*z)
!      with erfc(z)= 1 - erf(z)
! Eric Sonnendrucker 2008/12/13
!------------------------------------------------------------------
subroutine FriedConte(z,zeta,dzeta)
  double complex, intent(in)  :: z
  double complex, intent(out) :: zeta
  double complex, intent(out) :: dzeta

  double precision :: sqrtpi   ! square root of pi
  double precision :: xi,yi,u,v
  integer          :: iflag

  sqrtpi = 1.7724538509055159

  xi = real(z)
  yi = imag(z)

  call wofz(xi,yi,u,v,iflag)

  zeta  = dcmplx(-sqrtpi*v,sqrtpi*u)
  dzeta = -2.0d0*(1+z*zeta)

  return
end subroutine FriedConte
