module sll_common_array_initializers_module
#include "sll_assert.h" 
#include "sll_working_precision.h"
#include "sll_constants.h"
#include "sll_assert.h"

  implicit none

  ! The functions specified here are meant to have the specific signature
  ! described in the sll_parallel_array_initializer_module. Else, they could
  ! not be used with this module.

contains


  ! -------------------------------------------------------------------------
  !
  !             Landau damping 4d initialization function
  !
  ! -------------------------------------------------------------------------
  !
  ! The params array is declared optional to conform with the expected 
  ! function signature of the initializer subroutines, but in the particular
  ! case of the landau initializer, the params array must be passed.



  function sll_landau_initializer_2d( x, vx, params ) 
    sll_real64 :: sll_landau_initializer_2d
    sll_real64, intent(in) :: x
    sll_real64, intent(in) :: vx
 
    sll_real64, dimension(:), intent(in), optional :: params
    sll_real64 :: eps
    sll_real64 :: kx
    sll_real64 :: factor1

    if( .not. present(params) ) then
       print *, 'sll_landau_initializer_4d, error: the params array must ', &
            'be passed. params(1) = epsilon, params(2) = kx, params(3) = ky.'
       stop
    end if
    SLL_ASSERT(size(params)>=2)
    kx = params(1)
    eps = params(2)

    !Normalization
    !sagemath command
    !sage : var('u v epsilon a b c d x y')
    !sage : f(a,b,c,d,epsilon) =integral(integral(integral(integral((1+epsilon*cos(2*pi/(b-a)*x))*exp(-(u*u+v*v)/2),u,-oo,oo),v,-oo,oo),x,a,b),y,c,d)
    
!!$    factor1 =  1./( (eta2_min - eta2_max) &
!!$               *(((eta1_min - eta1_max)* &
!!$               sin(2*sll_pi*eta1_min/(eta1_min - eta1_max)) &
!!$                - (eta1_min - eta1_max)* &
!!$               sin(2*sll_pi*eta1_max/(eta1_min - eta1_max)))*eps  &
!!$               + 2*sll_pi*eta1_min - 2*sll_pi*eta1_max))
    factor1 = 1.0_f64/sqrt(2.0_f64*sll_pi)
!!$    sll_landau_initializer_4d = factor1 * &
!!$         (1.0_f64/((eta2_max-eta2_min)*(eta1_max-eta1_min))+eps*cos(kx*x))*exp(-0.5_f64*(vx**2+vy**2))
    sll_landau_initializer_2d = factor1 * &
         (1.0_f64+eps*cos(kx*x))*exp(-0.5_f64*vx**2)
  end function sll_landau_initializer_2d






  ! This is a simplistic initializer aimed at a 4d cartesian distribution
  ! function, periodic in x and y, and compact-ish in vx and vy.
  !
  ! Basically:
  !                          
  ! f(x,y,vx,vy) = alpha*exp(-0.5*((x -xc )^2+(y - yc)^2)) + 
  !                beta* exp(-0.5*((vx-vxc)^2+(vy-vyc)^2))
  !                          
  !
  ! It is meant to be used in the intervals:
  ! x:  [ 0,2*pi/kx]
  ! y:  [ 0,2*pi/ky]
  ! vx: [-6,6]
  ! vy: [-6,6]

  ! convention for the params array:
  ! params(1) = eta1_min
  ! params(2) = eta1_max
  ! params(3) = eta2_min
  ! params(4) = eta2_max
  ! params(5) = epsilon

  function sll_gaussian_initializer_4d( x, y, vx, vy, params ) 
    sll_real64 :: sll_gaussian_initializer_4d
    sll_real64, intent(in) :: x
    sll_real64, intent(in) :: y
    sll_real64, intent(in) :: vx
    sll_real64, intent(in) :: vy

    sll_real64, dimension(:), intent(in), optional :: params
    sll_real64 :: xc
    sll_real64 :: yc
    sll_real64 :: vxc
    sll_real64 :: vyc
    sll_real64 :: alpha
    sll_real64 :: beta 

    if( .not. present(params) ) then
       print *, 'sll_gaussian_initializer_4d, error: the params array must ', &
            'be passed: ', &
            'params(1) = xc, params(2) = yc, params(3) = vxc, params(4) = vyc',&
            'params(5) = alpha, params(6) = beta'
       stop
    end if

    xc    = params(1)
    yc    = params(2)
    vxc   = params(3)
    vyc   = params(4)
    alpha = params(5)
    beta  = params(6)

    sll_gaussian_initializer_4d = alpha*exp(-0.5_f64*((x-xc)**2+(y-yc)**2)) + &
                                  beta *exp(-0.5_f64*((vx-vxc)**2+(vy-vyc)**2))

  end function sll_gaussian_initializer_4d


  ! -------------------------------------------------------------------------
  !
  !             Landau damping 4d initialization function
  !
  ! -------------------------------------------------------------------------
  !
  ! The params array is declared optional to conform with the expected 
  ! function signature of the initializer subroutines, but in the particular
  ! case of the landau initializer, the params array must be passed.



  function sll_landau_initializer_4d( x, y, vx, vy, params ) 
    sll_real64 :: sll_landau_initializer_4d
    sll_real64, intent(in) :: x
    sll_real64, intent(in) :: y
    sll_real64, intent(in) :: vx
    sll_real64, intent(in) :: vy

    sll_real64, dimension(:), intent(in), optional :: params
    sll_real64 :: eta1_min
    sll_real64 :: eta1_max
    sll_real64 :: eta2_min
    sll_real64 :: eta2_max

    sll_real64 :: eps
    sll_real64 :: kx
    sll_real64 :: factor1

    if( .not. present(params) ) then
       print *, 'sll_landau_initializer_4d, error: the params array must ', &
            'be passed. params(1) = eta1_min, params(2) = eta1_max, ', &
            'params(3) = eta2_min, params(4) = eta2_max, params(5) = epsilon.'
       stop
    end if

    SLL_ASSERT( size(params) >= 5 )

    eta1_min = params(1)
    eta1_max = params(2)
    eta2_min = params(3)
    eta2_max = params(4)
    eps      = params(5)
    kx  =  2. * sll_pi / (eta1_max - eta1_min)

    !Normalization
    !sagemath command
    !sage : var('u v epsilon a b c d x y')
    !sage : f(a,b,c,d,epsilon) =integral(integral(integral(integral((1+epsilon*cos(2*pi/(b-a)*x))*exp(-(u*u+v*v)/2),u,-oo,oo),v,-oo,oo),x,a,b),y,c,d)
    
!!$    factor1 =  1./( (eta2_min - eta2_max) &
!!$               *(((eta1_min - eta1_max)* &
!!$               sin(2*sll_pi*eta1_min/(eta1_min - eta1_max)) &
!!$                - (eta1_min - eta1_max)* &
!!$               sin(2*sll_pi*eta1_max/(eta1_min - eta1_max)))*eps  &
!!$               + 2*sll_pi*eta1_min - 2*sll_pi*eta1_max))
    factor1 = 1.0_f64/(2.0*sll_pi)
!!$    sll_landau_initializer_4d = factor1 * &
!!$         (1.0_f64/((eta2_max-eta2_min)*(eta1_max-eta1_min))+eps*cos(kx*x))*exp(-0.5_f64*(vx**2+vy**2))
    sll_landau_initializer_4d = factor1 * &
         (1.0_f64+eps*cos(kx*x))*exp(-0.5_f64*(vx**2+vy**2))
  end function sll_landau_initializer_4d

  function sll_test_x_transport_initializer_v1v2x1x2( vx, vy, x, y, params ) 
    sll_real64 :: sll_test_x_transport_initializer_v1v2x1x2
    sll_real64, intent(in) :: x
    sll_real64, intent(in) :: y
    sll_real64, intent(in) :: vx
    sll_real64, intent(in) :: vy
    sll_real64 :: t

    sll_real64, dimension(:), intent(in), optional :: params
    sll_real64 :: eta1_min
    sll_real64 :: eta1_max
    sll_real64 :: kx

    if( .not. present(params) ) then
       print *, ' sll_test_x_transport_initializer, error: the params array', & 
            'must be passed. params(1) = epsilon, params(2)=kx, params(3) = ky.'
       stop
    end if

    eta1_min = params(1)
    eta1_max = params(2)
    t=params(11)
    kx  =  2.0_f64 * sll_pi / (eta1_max - eta1_min)


    !sll_test_x_transport_initializer_v1v2x1x2 =sin(kx*(x-t))
    sll_test_x_transport_initializer_v1v2x1x2 =sin(kx*(x-vx*t))
    !sll_test_x_transport_initializer_v1v2x1x2 = exp(-4*x**2)

  end function sll_test_x_transport_initializer_v1v2x1x2

  function sll_test_y_transport_initializer_v1v2x1x2( vx, vy, x, y, params ) 
    sll_real64 :: sll_test_y_transport_initializer_v1v2x1x2
    sll_real64, intent(in) :: x
    sll_real64, intent(in) :: y
    sll_real64, intent(in) :: vx
    sll_real64, intent(in) :: vy
    sll_real64 :: t

    sll_real64, dimension(:), intent(in), optional :: params
    sll_real64 :: eta2_min
    sll_real64 :: eta2_max
    sll_real64 :: kx

    eta2_min = params(3)
    eta2_max = params(4)
    t=params(11)
    kx  =  2.0_f64 * sll_pi / (eta2_max - eta2_min)
    if( .not. present(params) ) then
       print *, ' sll_test_y_transport_initializer, error: the params array', & 
            'must be passed. params(1) = epsilon, params(2)=kx, params(3) = ky.'
       stop
    end if

    sll_test_y_transport_initializer_v1v2x1x2 =sin(kx*(y-vy*t))
    !sll_test_y_transport_initializer_v1v2x1x2 =exp(-4*y**2)
    !sll_test_x_transport_initializer_v1v2x1x2 = 2_f64

  end function sll_test_y_transport_initializer_v1v2x1x2


  function sll_test_vx_transport_initializer_v1v2x1x2( vx, vy, x, y, params ) 
    sll_real64 :: sll_test_vx_transport_initializer_v1v2x1x2
    sll_real64, intent(in) :: x
    sll_real64, intent(in) :: y
    sll_real64, intent(in) :: vx
    sll_real64, intent(in) :: vy
    sll_real64 :: t

    sll_real64, dimension(:), intent(in), optional :: params
    sll_real64 :: eta1_min
    sll_real64 :: eta1_max
    sll_real64 :: kx,v,hh
    if( .not. present(params) ) then
       print *, ' sll_test_vx_transport_initializer, error: the params array', &
            ' mustbe passed. params(1)=epsilon,params(2)=kx, params(3) = ky.'
       stop
    end if

    eta1_min = params(1)
    eta1_max = params(2)
    t=params(11)
    kx  =  2. * sll_pi / (eta1_max - eta1_min)

    !sll_test_vx_transport_initializer_v1v2x1x2 =exp(-4*(vx-t)**2)
    v=vx-t
    hh=1
    if ((v.gt.1).or.(v.lt.-1)) hh=0
    sll_test_vx_transport_initializer_v1v2x1x2=(1+(-3+(3-v**2)*v**2)*v**2)*hh
    !write(*,*) 'tini=',t

!!$    sll_landau_initializer_v1v2x1x2 = 1
!!$    if (x < 0) sll_landau_initializer_v1v2x1x2 = 0

  end function sll_test_vx_transport_initializer_v1v2x1x2
  function sll_test_vy_transport_initializer_v1v2x1x2( vx, vy, x, y, params ) 
    sll_real64 :: sll_test_vy_transport_initializer_v1v2x1x2
    sll_real64, intent(in) :: x
    sll_real64, intent(in) :: y
    sll_real64, intent(in) :: vx
    sll_real64, intent(in) :: vy
    sll_real64 :: t

    sll_real64, dimension(:), intent(in), optional :: params
    sll_real64 :: eta1_min
    sll_real64 :: eta1_max
    sll_real64 :: eta2_min
    sll_real64 :: eta2_max

    sll_real64 :: eps
    sll_real64 :: kx
    sll_real64 :: factor1
    sll_real64 :: v,hh
    t=params(11)

    if( .not. present(params) ) then
       print *, ' sll_test_vy_transport_initializer, error: the params array', & 
            'must be passed. params(1) = epsilon, params(2)=kx, params(3) = ky.'
       stop
    end if
    v=vy-t
    hh=1
    if ((v.gt.1).or.(v.lt.-1)) hh=0
    sll_test_vy_transport_initializer_v1v2x1x2 =(1+(-3+(3-v**2)*v**2)*v**2)*hh
    !sll_test_vy_transport_initializer_v1v2x1x2 =exp(-4*vy**2)
    !sll_test_vy_transport_initializer_v1v2x1x2 = 2_f64

  end function sll_test_vy_transport_initializer_v1v2x1x2

  function sll_test_xvx_transport_initializer_v1v2x1x2( vx, vy, x, y, params ) 
    sll_real64 :: sll_test_xvx_transport_initializer_v1v2x1x2
    sll_real64, intent(in) :: x
    sll_real64, intent(in) :: y
    sll_real64, intent(in) :: vx
    sll_real64, intent(in) :: vy
    sll_real64 :: t

    sll_real64, dimension(:), intent(in), optional :: params
    sll_real64 :: eta1_min
    sll_real64 :: eta1_max
    sll_real64 :: kx,v,hh,xx
    if( .not. present(params) ) then
       print *, ' sll_test_xvx_transport_initializer, error: the params array', &
            ' mustbe passed. params(1)=epsilon,params(2)=kx, params(3) = ky.'
       stop
    end if

    eta1_min = params(1)
    eta1_max = params(2)
    t=params(11)
    kx  =  2. * sll_pi / (eta1_max - eta1_min)

    !sll_test_vx_transport_initializer_v1v2x1x2 =exp(-4*(vx-t)**2)
    v=vx-t
    xx=x-t*vx+t*t/2
    hh=1
    if ((v.gt.1).or.(v.lt.-1)) hh=0
    sll_test_xvx_transport_initializer_v1v2x1x2=(1+(-3+(3-v**2)*v**2)*v**2)*hh*sin(kx*xx)

  end function sll_test_xvx_transport_initializer_v1v2x1x2

function sll_test_yvy_transport_initializer_v1v2x1x2( vx, vy, x, y, params ) 
    sll_real64 :: sll_test_yvy_transport_initializer_v1v2x1x2
    sll_real64, intent(in) :: x
    sll_real64, intent(in) :: y
    sll_real64, intent(in) :: vx
    sll_real64, intent(in) :: vy
    sll_real64 :: t

    sll_real64, dimension(:), intent(in), optional :: params
    sll_real64 :: eta2_min
    sll_real64 :: eta2_max
    sll_real64 :: ky,v,hh,yy
    if( .not. present(params) ) then
       print *, ' sll_test_yvy_transport_initializer, error: the params array', &
            ' mustbe passed. params(1)=epsilon,params(2)=kx, params(3) = ky.'
       stop
    end if

    eta2_min = params(3)
    eta2_max = params(4)
    t=params(11)
    ky  =  2. * sll_pi / (eta2_max - eta2_min)

    !sll_test_vx_transport_initializer_v1v2x1x2 =exp(-4*(vx-t)**2)
    v=vy-t
    yy=y-t*vy+t*t/2
    hh=1
    if ((v.gt.1).or.(v.lt.-1)) hh=0
    sll_test_yvy_transport_initializer_v1v2x1x2=(1+(-3+(3-v**2)*v**2)*v**2)*hh*sin(ky*yy)

  end function sll_test_yvy_transport_initializer_v1v2x1x2


!!$function sll_test_transport_2d_initializer_v1v2x1x2( vx, vy, x, y, params ) 
!!$    sll_real64 :: sll_test_transport_2d_initializer_v1v2x1x2
!!$    sll_real64, intent(in) :: x
!!$    sll_real64, intent(in) :: y
!!$    sll_real64, intent(in) :: vx
!!$    sll_real64, intent(in) :: vy
!!$    sll_real64 :: t
!!$
!!$    sll_real64, dimension(:), intent(in), optional :: params
!!$    sll_real64 :: eta1_min
!!$    sll_real64 :: eta1_max
!!$    sll_real64 :: kx,v,hh,xx
!!$    if( .not. present(params) ) then
!!$       print *, ' sll_test_xvx_transport_initializer, error: the params array', &
!!$            ' mustbe passed. params(1)=epsilon,params(2)=kx, params(3) = ky.'
!!$       stop
!!$    end if
!!$
!!$    eta1_min = params(1)
!!$    eta1_max = params(2)
!!$    eta2_min = params(3)
!!$    eta2_max = params(4)
!!$    t=params(11)
!!$    kx  =  2. * sll_pi / (eta1_max - eta1_min)
!!$    ky  =  2. * sll_pi / (eta2_max - eta2_min)
!!$
!!$    !sll_test_vx_transport_initializer_v1v2x1x2 =exp(-4*(vx-t)**2)
!!$    vvx=vx-t
!!$    xx=x-t*vx+t*t/2
!!$    vvy=vx-t
!!$    yy=x-t*vy+t*t/2
!!$    hh=1
!!$    if ((v.gt.1).or.(v.lt.-1)) hh=0
!!$    sll_test_transport_2d_initializer_v1v2x1x2=(1+(-3+(3-v**2)*v**2)*v**2)*hh*sin(kx*(xx-v*t))
!!$
!!$  end function sll_test_transport_2d_initializer_v1v2x1x2


  function sll_landau_1d_xvx_initializer_v1v2x1x2( vx, vy, x, y, params ) 
    sll_real64 :: sll_landau_1d_xvx_initializer_v1v2x1x2
    sll_real64, intent(in) :: x
    sll_real64, intent(in) :: y
    sll_real64, intent(in) :: vx
    sll_real64, intent(in) :: vy
    sll_real64  :: t

    sll_real64, dimension(:), intent(in), optional :: params
    sll_real64 :: eta1_min
    sll_real64 :: eta1_max
    sll_real64 :: eta2_min
    sll_real64 :: eta2_max

    sll_real64 :: eps
    sll_real64 :: kx
    sll_real64 :: factor1

    if( .not. present(params) ) then
       print *, 'sll_landau_1d_initializer_v1v2x1x2, error: the params array', &
            'must be passed params(1)= epsilon, params(2) = kx, params(3) = ky.'
       stop
    end if

    eta1_min = params(1)
    eta1_max = params(2)
    eta2_min = params(3)
    eta2_max = params(4)

    eps = params(5)
    !kx  =  2. * sll_pi / (eta1_max - eta1_min)
    kx=0.2_f64
    factor1 = 1.0_f64/sqrt((2.0*sll_pi))

    sll_landau_1d_xvx_initializer_v1v2x1x2 = factor1 * &
         (1.0_f64+eps*cos(kx*x))*exp(-0.5_f64*(vx**2))
  end function sll_landau_1d_xvx_initializer_v1v2x1x2

  function sll_landau_1d_yvy_initializer_v1v2x1x2( vx, vy, x, y, params ) 
    sll_real64 :: sll_landau_1d_yvy_initializer_v1v2x1x2
    sll_real64, intent(in) :: x
    sll_real64, intent(in) :: y
    sll_real64, intent(in) :: vx
    sll_real64, intent(in) :: vy
    sll_real64  :: t

    sll_real64, dimension(:), intent(in), optional :: params
    sll_real64 :: eta1_min
    sll_real64 :: eta1_max
    sll_real64 :: eta2_min
    sll_real64 :: eta2_max

    sll_real64 :: eps
    sll_real64 :: ky
    sll_real64 :: factor1

    if( .not. present(params) ) then
       print *, 'sll_landau_1d_initializer_v1v2x1x2, error: the params array', &
            'must be passed params(1)= epsilon, params(2) = kx, params(3) = ky.'
       stop
    end if

    eta1_min = params(1)
    eta1_max = params(2)
    eta2_min = params(3)
    eta2_max = params(4)

    eps = params(5)
    ky  =  2. * sll_pi / (eta2_max - eta2_min)
    factor1 = 1.0_f64/sqrt((2.0*sll_pi))

    sll_landau_1d_yvy_initializer_v1v2x1x2 = factor1 * &
         (1.0_f64+eps*cos(ky*y))*exp(-0.5_f64*(vy**2))
  end function sll_landau_1d_yvy_initializer_v1v2x1x2


  function sll_landau_2d_initializer_v1v2x1x2( vx, vy, x, y, params ) 
    sll_real64 :: sll_landau_2d_initializer_v1v2x1x2
    sll_real64, intent(in) :: x
    sll_real64, intent(in) :: y
    sll_real64, intent(in) :: vx
    sll_real64, intent(in) :: vy

    sll_real64, dimension(:), intent(in), optional :: params
    sll_real64 :: eta1_min
    sll_real64 :: eta1_max
    sll_real64 :: eta2_min
    sll_real64 :: eta2_max

    sll_real64 :: eps
    sll_real64 :: kx
    sll_real64 :: ky
    sll_real64 :: factor1

    if( .not. present(params) ) then
       print *, 'sll_landau_initializer_4d, error: the params array must ', &
            'be passed. params(1) = epsilon, params(2) = kx, params(3) = ky.'
       stop
    end if

    eta1_min = params(1)
    eta1_max = params(2)
    eta2_min = params(3)
    eta2_max = params(4)

    eps = params(5)
    kx  =  2. * sll_pi / (eta1_max - eta1_min)
    ky = kx
    factor1 = 1.0_f64/(2.0*sll_pi)
    sll_landau_2d_initializer_v1v2x1x2 = factor1 * &
         (1.0_f64+eps*cos(kx*x)*cos(ky*y))*exp(-0.5_f64*(vx**2+vy**2))

  end function sll_landau_2d_initializer_v1v2x1x2


  ! this function is a 1D landau initializer used for debugging
  ! 4D drift kinetic simulations in variables x1,x2,x3 ,v1
  ! the function is constant with respect to x2 and x3

  function sll_landau_initializer_dk_test_4d(v1,x1,x2,x3,params ) 
    sll_real64 :: sll_landau_initializer_dk_test_4d
    sll_real64, intent(in) :: x1
    sll_real64, intent(in) :: x2
    sll_real64, intent(in) :: x3
    sll_real64, intent(in) :: v1
    sll_real64, dimension(:), intent(in), optional :: params

    sll_real64 :: epsilon
    sll_real64 :: kx
    sll_real64 :: factor1

    if( .not. present(params) ) then
       print *, 'sll_landau_initializer_dk_test_4d, error: the params array must ', &
            'be passed. params(1) = epsilon, params(2) = kx'
       stop
    end if

    epsilon = params(1)
    kx      = params(2)
    factor1 = 0.5_f64/sll_pi

    !write(*,*) 'x1 v1',x1,v1

    sll_landau_initializer_dk_test_4d = factor1*&
         (1.0_f64+epsilon*cos(kx*x1))*exp(-0.5_f64*(v1**2))
  end function sll_landau_initializer_dk_test_4d

  !---------------------------------------------------------------------------
  !
  !                         Periodic Maxwellian
  !
  ! 4D distribution in [0,1]X[0,1][-6,6]X[-6,6]  with the property of being 
  ! periodic in the spatial directions (x1,x2) and gaussian in velocity space.
  !
  ! f(x,y,vx,vy) = sin(kx*x)*sin(ky*y) + 
  !                beta* exp(-0.5*((vx-vxc)^2+(vy-vyc)^2))
  !                          
  !
  ! It is meant to be used in the intervals:
  ! x:  [ 0,2*pi/kx]
  ! y:  [ 0,2*pi/ky]
  ! vx: [-6,6]
  ! vy: [-6,6]

  ! convention for the params array:
  ! params(1) = eta1_min
  ! params(2) = eta1_max
  ! params(3) = eta2_min
  ! params(4) = eta2_max
  ! params(5) = epsilon
  !
  !---------------------------------------------------------------------------

  function sll_periodic_gaussian_initializer_4d( x, y, vx, vy, params ) &
    result(val)

    sll_real64 :: val !sll_gaussian_initializer_4d
    sll_real64, intent(in) :: x
    sll_real64, intent(in) :: y
    sll_real64, intent(in) :: vx
    sll_real64, intent(in) :: vy

    sll_real64, dimension(:), intent(in), optional :: params
    sll_real64 :: xc
    sll_real64 :: yc
    sll_real64 :: vxc
    sll_real64 :: vyc
    sll_real64 :: alpha
    sll_real64 :: beta 

    if( .not. present(params) ) then
       print *, 'sll_gaussian_initializer_4d, error: the params array must ', &
            'be passed: ', &
            'params(1) = xc, params(2) = yc, params(3) = vxc, params(4) = vyc',&
            'params(5) = alpha, params(6) = beta'
       stop
    end if

    xc    = params(1)
    yc    = params(2)
    vxc   = params(3)
    vyc   = params(4)
    alpha = params(5)
    beta  = params(6)

    val = alpha*exp(-0.5_f64*((x-xc)**2+(y-yc)**2)) + &
          beta *exp(-0.5_f64*((vx-vxc)**2+(vy-vyc)**2))

  end function sll_periodic_gaussian_initializer_4d




 !---------------------------------------------------------------------------
  !
  !                         Periodic-Periodic initializer
  !
  ! 4D distribution in [0,2*pi/kx]X[0,2*pi/ky][-6,6]X[-6,6]  with the property of being 
  ! periodic in the spatial directions (x1,x2) and gaussian in velocity space.
  !
  ! f(x,y,vx,vy) = (1 + alpha * cos(kx*x)*cos(ky*y) )
  !                * exp(-0.5*((vx-vxc)^2+(vy-vyc)^2))/ 2*pi
  !                          
  !  This function is described in the article of Crouseilles and al. 2009
  ! ' A parallel Vlasov solver based on local cubic spline interpolation on patches'
  !
  ! It is meant to be used in the intervals:
  ! x:  [ 0,2*pi/kx]
  ! y:  [ 0,2*pi/ky]
  ! vx: [-6,6]
  ! vy: [-6,6]
  
  ! convention for the params array:
  ! params(1) = eta1_min
  ! params(2) = eta1_max
  ! params(3) = eta2_min
  ! params(4) = eta2_max
  ! params(5) = vxc
  ! params(6) = vyc
  ! params(7) = alpha
  ! params(8) = beta
  !---------------------------------------------------------------------------

  function sll_periodic_periodic_gaussian2009_initializer_4d( x, y, vx, vy, params ) &
    result(val)
    
    sll_real64 :: val !sll_gaussian_initializer_4d
    sll_real64, intent(in) :: x
    sll_real64, intent(in) :: y
    sll_real64, intent(in) :: vx
    sll_real64, intent(in) :: vy
    
    sll_real64, dimension(:), intent(in), optional :: params
    sll_real64 :: eta1_min,eta1_max
    sll_real64 :: eta2_min,eta2_max
    sll_real64 :: vxc
    sll_real64 :: vyc
    sll_real64 :: alpha
    sll_real64 :: beta
    sll_real64 :: kx,ky
    
    if( .not. present(params) ) then
       print *, 'sll_periodic_periodic_gaussian2009_initializer_4d, error: the params array must ', &
            'be passed: ', &
            'params(1) = eta1_min, params(2) = eta1_max, params(3) = eta2_min, params(4) = eta2_max',&
            'params(5) = vxc, params(6) = vyx, params(7) = alpha, params(8) = beta'
       stop
    end if
    
    eta1_min = params(1)
    eta1_max = params(2)
    eta2_min = params(3)
    eta2_max = params(4)
    vxc   = params(5)
    vyc   = params(6)
    alpha      = params(7)
    beta      = params(8)
    kx  =  2. * sll_pi / (eta1_max - eta1_min)
    ky  =  2. * sll_pi / (eta2_max - eta2_min)
    


    val = beta*(1 + alpha * cos(kx*x)*cos(ky*y) )&
         * exp(-0.5*((vx-vxc)**2+(vy-vyc)**2)) &
         / (2*sll_pi)
    
  end function sll_periodic_periodic_gaussian2009_initializer_4d
  

   !---------------------------------------------------------------------------
  !
  !                         Periodic-Periodic initializer another
  !
  ! 4D distribution in [0,2*pi/kx]X[0,2*pi/ky][-6,6]X[-6,6]  with the property of being 
  ! periodic in the spatial directions (x1,x2) and gaussian in velocity space.
  !
  ! f(x,y,vx,vy) = (1 + alpha * (cos(kx*x) + cos(ky*y)) )
  !                * exp(-0.5*((vx-vxc)^2+(vy-vyc)^2))/ (2*sll_pi)
  !                          
  !  This function is described in the article of Filbet and Sonnendrucker 2002
  ! 'Comparison of Eulerian Vlasov solvers'
  !
  ! It is meant to be used in the intervals:
  ! x:  [ 0,2*pi/kx]
  ! y:  [ 0,2*pi/ky]
  ! vx: [-6,6]
  ! vy: [-6,6]
  
  ! convention for the params array:
  ! params(1) = eta1_min
  ! params(2) = eta1_max
  ! params(3) = eta2_min
  ! params(4) = eta2_max
  ! params(5) = vxc
  ! params(6) = vyc
  ! params(7) = alpha
  ! params(8) = beta
  !---------------------------------------------------------------------------

  function sll_periodic_periodic_gaussian2002_initializer_4d( x, y, vx, vy, params ) &
    result(val)
    
    sll_real64 :: val !sll_gaussian_initializer_4d
    sll_real64, intent(in) :: x
    sll_real64, intent(in) :: y
    sll_real64, intent(in) :: vx
    sll_real64, intent(in) :: vy
    
    sll_real64, dimension(:), intent(in), optional :: params
    sll_real64 :: eta1_min,eta1_max
    sll_real64 :: eta2_min,eta2_max
    sll_real64 :: vxc
    sll_real64 :: vyc
    sll_real64 :: alpha
    sll_real64 :: beta
    sll_real64 :: kx,ky
    
    if( .not. present(params) ) then
       print *, 'sll_periodic_periodic_gaussian2002_initializer_4d, error: the params array must ', &
            'be passed: ', &
            'params(1) = eta1_min, params(2) = eta1_max, params(3) = eta2_min, params(4) = eta2_max',&
            'params(5) = vxc, params(6) = vyx, params(7) = alpha, params(8) = beta'
       stop
    end if
    
    eta1_min = params(1)
    eta1_max = params(2)
    eta2_min = params(3)
    eta2_max = params(4)
    vxc   = params(5)
    vyc   = params(6)
    alpha      = params(7)
    beta      = params(8)
    kx  =  2. * sll_pi / (eta1_max - eta1_min)
    ky  =  2. * sll_pi / (eta2_max - eta2_min)
    


    val = beta*(1 + alpha * (cos(kx*x)+cos(ky*y)) )&
         * exp(-0.5*((vx-vxc)**2+(vy-vyc)**2)) &
         / (2*sll_pi)
    
  end function sll_periodic_periodic_gaussian2002_initializer_4d


end module sll_common_array_initializers_module
