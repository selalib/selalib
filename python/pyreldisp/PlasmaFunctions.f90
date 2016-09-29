module plasma_functions

  use precision_module, only: dp

  interface
    subroutine wofz(xi,yi,u,v,flag)
      use precision_module, only: dp
      real(dp), intent(in   ) :: xi,yi
      real(dp), intent(  out) :: u,v
      logical , intent(  out) :: flag
    end subroutine wofz
  end interface

contains

  !-----------------------------------------------------------------------------
  ! Computes the Fried and Conte plasma dispersion function
  !
  ! Z(z)=i*sqrt(pi)*exp(-z**2)*(1-erf(-i*z))
  !
  ! Eric Sonnendrucker 2008/12/13
  ! Modified by Yaman Gü¢lü and Edoardo Zoni: 2016/09
  !-----------------------------------------------------------------------------
  subroutine FriedConte(z,zeta,dzeta)

    complex(dp), intent(in   ) :: z
    complex(dp), intent(  out) :: zeta
    complex(dp), intent(  out) :: dzeta

    real(dp), parameter :: sqrtpi = 1.7724538509055159_dp
    real(dp) :: u,v
    logical  :: flag

    call wofz(xi=real(z),yi=imag(z),u=u,v=v,flag=flag)

    zeta  = cmplx(-sqrtpi*v,sqrtpi*u,kind=dp)
    dzeta = -2.0_dp*(1.0_dp+z*zeta)

  end subroutine FriedConte

end module plasma_functions
