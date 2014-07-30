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



  function sll_gaussian_initializer_2d( x_1, x_2, params ) 
    sll_real64 :: sll_gaussian_initializer_2d
    sll_real64, intent(in) :: x_1
    sll_real64, intent(in) :: x_2 
    sll_real64, dimension(:), intent(in), optional :: params
    sll_real64 :: xc_1
    sll_real64 :: xc_2
    sll_real64 :: sigma_1
    sll_real64 :: sigma_2
    !sll_real64 :: kx
    !sll_real64 :: factor1

    if( .not. present(params) ) then
       print *, 'sll_gaussian_initializer_2d, error: ', &
         'the params array must be passed.', &
         ' params(1) = xc_1', &
         ' params(2) = xc_2', &
         ' params(3) = sigma_1', &
         ' params(4) = sigma_2'
    end if
    SLL_ASSERT(size(params)>=4)
    xc_1 = params(1)
    xc_2 = params(2)
    sigma_1 = params(3)
    sigma_2 = params(4)
    sll_gaussian_initializer_2d =  &
         exp(-0.5_f64*(x_1-xc_1)**2/sigma_1**2 &
         -0.5_f64*(x_2-xc_2)**2/sigma_2**2)
  end function sll_gaussian_initializer_2d


  function sll_cos_bell_initializer_2d( x_1, x_2, params ) result(res)
    sll_real64 :: res
    sll_real64, intent(in) :: x_1
    sll_real64, intent(in) :: x_2 
    sll_real64, dimension(:), intent(in), optional :: params
    sll_real64 :: xc_1
    sll_real64 :: xc_2
    sll_real64 :: r

    if( .not. present(params) ) then
       print *, 'sll_cos_bell_initializer_2d, error: ', &
         'the params array must be passed.', &
         ' params(1) = xc_1', &
         ' params(2) = xc_2'
    end if
    SLL_ASSERT(size(params)>=2)
    xc_1 = params(1)
    xc_2 = params(2)
    
    r = sqrt((x_1-xc_1)**2+(x_2-xc_2)**2)
    
    if(r<0.5_f64*sll_pi)then
      res = cos(r)**6
    else
      res = 0._f64
    endif
    
  end function sll_cos_bell_initializer_2d




  !swirling deformation flow

  function sll_SDF_A1_initializer_2d( x_1, x_2, params ) result(res)
    sll_real64 :: res
    sll_real64, intent(in) :: x_1
    sll_real64, intent(in) :: x_2 
    sll_real64, dimension(:), intent(in), optional :: params
    if(present(params))then
      if(size(params)>=100)then
        print *,'#params needs not to have size >=100'
      endif
    endif
    res = -cos(0.5_f64*x_1)**2*sin(x_2)
  end function sll_SDF_A1_initializer_2d

  function sll_SDF_A2_initializer_2d( x_1, x_2, params ) result(res)
    sll_real64 :: res
    sll_real64, intent(in) :: x_1
    sll_real64, intent(in) :: x_2 
    sll_real64, dimension(:), intent(in), optional :: params
    if(present(params))then
      if(size(params)>=100)then
        print *,'#params needs not to have size >=100'
      endif
    endif    
    res = cos(0.5_f64*x_2)**2*sin(x_1)
  end function sll_SDF_A2_initializer_2d


  function sll_SDF_time_initializer_1d( t, params ) result(res)
    sll_real64 :: res
    sll_real64, intent(in) :: t
    sll_real64, dimension(:), intent(in), optional :: params
    sll_real64 :: time_period

    if( .not. present(params) ) then
      print *, 'sll_SDF_time_initializer_1d, error: ', &
        'the params array must be passed.', &
        ' params(1) = time_period'
      stop
    end if
    SLL_ASSERT(size(params)>=1)
    time_period = params(1)
    res = sll_pi*cos(sll_pi*t/time_period)
  end function sll_SDF_time_initializer_1d




  


  ! -------------------------------------------------------------------------
  !
  !             Landau damping 2d initialization function
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
       print *, 'sll_landau_initializer_2d, error: the params array must ', &
            'be passed. params(1) = epsilon, params(2) = kx.'
       stop
    end if
    SLL_ASSERT(size(params)>=2)
    kx = params(1)
    eps = params(2)

    factor1 = 1.0_f64/sqrt(2.0_f64*sll_pi)
    sll_landau_initializer_2d = factor1 * &
         (1.0_f64+eps*cos(kx*x))*exp(-0.5_f64*vx**2)
  end function sll_landau_initializer_2d


  function sll_bump_on_tail_initializer_2d( x, vx, params ) result(res)
    sll_real64 :: res
    sll_real64, intent(in) :: x
    sll_real64, intent(in) :: vx
 
    sll_real64, dimension(:), intent(in), optional :: params
    sll_real64 :: eps
    sll_real64 :: kx
    sll_real64 :: factor1

    if( .not. present(params) ) then
       print *, 'sll_bump_on_tail_initializer_2d, error: the params array must ', &
            'be passed. params(1) = epsilon, params(2) = kx.'
       stop
    end if
    SLL_ASSERT(size(params)>=2)
    kx = params(1)
    eps = params(2)
    factor1 = 1.0_f64/sqrt(2.0_f64*sll_pi)
    res = factor1 * (1._f64+eps*cos(kx*x))&
      *(0.9_f64*exp(-0.5_f64*vx**2)+0.2_f64*exp(-0.5_f64*(vx-4.5_f64)**2/0.5**2))
    !nbox = 3     
  end function sll_bump_on_tail_initializer_2d

  function sll_two_stream_instability_initializer_2d( x, vx, params ) result(res)
    sll_real64 :: res
    sll_real64, intent(in) :: x
    sll_real64, intent(in) :: vx
 
    sll_real64, dimension(:), intent(in), optional :: params
    sll_real64 :: eps
    sll_real64 :: kx
    sll_real64 :: factor1

    if( .not. present(params) ) then
       print *, 'sll_two_stream_instability_initializer_2d, error: the params array must ', &
            'be passed. params(1) = epsilon, params(2) = kx.'
       stop
    end if
    SLL_ASSERT(size(params)>=2)
    kx = params(1)
    eps = params(2)
    factor1 = 1.0_f64/sqrt(2.0_f64*sll_pi)
    res = factor1 * (1._f64+eps*cos(kx*x))*vx**2*exp(-0.5_f64*vx**2)
  end function sll_two_stream_instability_initializer_2d


  
  
  

  function sll_diocotron_initializer_2d( r, theta, params ) result(res)
    sll_real64 :: res
    sll_real64, intent(in) :: r
    sll_real64, intent(in) :: theta
 
    sll_real64, dimension(:), intent(in), optional :: params
    sll_real64 :: eps
    sll_real64 :: k_mode
    sll_real64 :: r_minus
    sll_real64 :: r_plus

    if( .not. present(params) ) then
       print *, '#sll_diocotron_initializer_2d, error: the params array must ', &
            'be passed. params(1) = r_minus'
       print *,'#params(2)= r_plus params(3)=epsilon param(4)=k_mode'     
       stop
    end if
    SLL_ASSERT(size(params)>=4)
    r_minus = params(1) 
    r_plus =  params(2)
    eps = params(3) 
    k_mode = params(4) 
    
    
    if((r>=r_minus).and.(r<=r_plus))then
      res = (1.0_f64+eps*cos(k_mode*theta))
    else
      res = 0._f64  
    endif 
       
  end function sll_diocotron_initializer_2d
  
  
  function sll_diocotron_initializer_2d2( x, y, params ) result(res)
    sll_real64 :: res
    sll_real64, intent(in) :: x
    sll_real64, intent(in) :: y
 
    sll_real64, dimension(:), intent(in), optional :: params
    sll_real64 :: eps
    sll_real64 :: k_mode
    sll_real64 :: r_minus
    sll_real64 :: r_plus
    sll_real64 :: r
    sll_real64 :: theta

    if( .not. present(params) ) then
       print *, '#sll_diocotron_initializer_2d, error: the params array must ', &
            'be passed. params(1) = r_minus'
       print *,'#params(2)= r_plus params(3)=epsilon param(4)=k_mode'     
       stop
    end if
    SLL_ASSERT(size(params)>=4)
    r_minus = params(1) 
    r_plus =  params(2)
    eps = params(3) 
    k_mode = params(4) 
    
    r= sqrt(x**2+y**2)
    
    if (y>=0) then
      theta = acos(x/r)
    else
      theta = 2._f64*sll_pi-acos(x/r)
    endif
    if((r>=r_minus).and.(r<=r_plus))then
      res = (1.0_f64+eps*cos(k_mode*theta))
    else
      res = 0._f64  
    endif 
  
  end function sll_diocotron_initializer_2d2
 
  function sll_beam_initializer_2d( x, vx, params ) result(res)
    sll_real64 :: res
    sll_real64, intent(in) :: x
    sll_real64, intent(in) :: vx
 
    sll_real64, dimension(:), intent(in), optional :: params
    sll_real64 :: alpha
    sll_real64 :: x_plus
    sll_real64 :: x_minus
    

    if( .not. present(params) ) then
       print *, '#sll_beam_initializer_2d, error: the params array must ', &
            'be passed. params(1) = alpha'
       stop
    end if
    SLL_ASSERT(size(params)>=1)
    alpha = params(1) 
    
    x_plus = (x+1.2_f64)/0.3_f64
    x_minus = (x-1.2_f64)/0.3_f64
    res = 4._f64/sqrt(2._f64*sll_pi*alpha)
    res = res*(0.5_f64*erf(x_plus)-0.5_f64*erf(x_minus))
    res = res*exp(-vx**2/(2*alpha))    

  end function sll_beam_initializer_2d


  
  
  
  


  function sll_KHP1_2d( x, y, params ) result(res)
   sll_real64 :: res
   sll_real64, intent(in) :: x
   sll_real64, intent(in) :: y

   sll_real64, dimension(:), intent(in), optional :: params
   sll_real64 :: eps
   sll_real64 :: k_mode_x
   sll_real64 :: k_mode_y

   if( .not. present(params) ) then
      print *, '#sll_KHP1_2d, error: the params array must ', &
           'be passed.'
      print *,'#params(1)= eps  param(2)=k_mode'
      stop
   end if
   SLL_ASSERT(size(params)>=3)
   eps = params(1)
   k_mode_x = params(2)
   k_mode_y = params(3)

   res = sin(k_mode_y*y)+eps*cos(k_mode_x*x)

  end function sll_KHP1_2d


  function sll_DSG_2d( eta1, eta2, params ) result(res)
    sll_real64  :: res
    sll_real64, intent(in)   :: eta1
    sll_real64, intent(in)   :: eta2
    sll_real64, dimension(:), intent(in), optional :: params    
    
    sll_real64  :: eta1n
    sll_real64  :: eta2n
    sll_real64  :: eta1_min
    sll_real64  :: eta2_min
    sll_real64  :: eta1_max
    sll_real64  :: eta2_max
    if( .not. present(params) ) then
       print *, '#sll_sll_D_sharped_Geo_2d, error: the params array must ', &
            'be passed. params(1) = eta1_min, params(2) = eta2_min', &
            'be passed. params(3) = eta1_max, params(4) = eta2_max'
       stop
    end if
    SLL_ASSERT(size(params)>=4)
    eta1_min =params(1)
    eta2_min =params(2)
    eta1_max =params(3)
    eta2_max =params(4)
    eta1n = (eta1 - eta1_min)/(eta1_max - eta1_min)
    eta2n = (eta2 - eta2_min)/(eta2_max - eta2_min)
    res =  4._f64*eta1n*(1._f64 - eta1n)* (1._f64 + 0.1_f64*sin(8.*sll_pi*eta2n))
  end function sll_DSG_2d


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
       print *,'does not depend on y',y
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

  function sll_landau_mode_initializer_4d( x, y, vx, vy, params ) result(res)
    sll_real64 :: res
    sll_real64, intent(in) :: x
    sll_real64, intent(in) :: y
    sll_real64, intent(in) :: vx
    sll_real64, intent(in) :: vy
    sll_real64, dimension(:), intent(in), optional :: params
    sll_real64 :: eps
    sll_real64 :: kx
    sll_real64 :: ky
    sll_real64 :: factor1

    if( .not. present(params) ) then
       print *, 'sll_landau_initializer_4d, error: the params array must ', &
            'be passed. params(1) = kx, params(2) = ky, ', &
            'params(3) = eps.'
       stop
    end if

    SLL_ASSERT( size(params) >= 3 )

    kx = params(1)
    ky = params(2)
    eps      = params(3)

    factor1 = 1.0_f64/(2.0*sll_pi)
    res = factor1 * &
         (1.0_f64+eps*cos(kx*x)*cos(ky*y))*exp(-0.5_f64*(vx**2+vy**2))
  end function sll_landau_mode_initializer_4d


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
       print *,'does not depend on y,vy',y,vy     
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
       print *,'does not depend on x vx',x,vx     
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
       print *,'does not depend on x,y,vy',x,y,vy
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
    !sll_real64 :: eta1_min
    !sll_real64 :: eta1_max
    !sll_real64 :: eta2_min
    !sll_real64 :: eta2_max

    !sll_real64 :: eps
    !sll_real64 :: kx
    !sll_real64 :: factor1
    sll_real64 :: v,hh
    t=params(11)

    if( .not. present(params) ) then
       print *, ' sll_test_vy_transport_initializer, error: the params array', & 
            'must be passed. params(1) = epsilon, params(2)=kx, params(3) = ky.'
       print *, 'does not depend on x,y,vx',x,y,vx     
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
       print *, 'does not depend on y and vy',y,vy
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
       print *,'does not depend on x and vx',x,vx
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
    !sll_real64  :: t

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
       print *,'does not depend on y and vy',y,vy
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

function sll_twostream_1d_xvx_initializer_v1v2x1x2( vx, vy, x, y, params ) 
    sll_real64 :: sll_twostream_1d_xvx_initializer_v1v2x1x2
    sll_real64, intent(in) :: x
    sll_real64, intent(in) :: y
    sll_real64, intent(in) :: vx
    sll_real64, intent(in) :: vy
    sll_real64  :: v0

    sll_real64, dimension(:), intent(in), optional :: params
    sll_real64 :: eta1_min
    sll_real64 :: eta1_max
    sll_real64 :: eta2_min
    sll_real64 :: eta2_max

    sll_real64 :: eps
    sll_real64 :: kx
    sll_real64 :: factor1

    if( .not. present(params) ) then
       print *, 'sll_twostream_1d_initializer_v1v2x1x2, error: the params array', &
            'must be passed params(1)= epsilon, params(2) = kx, params(3) = ky.'
       stop
    end if

    eta1_min = params(1)
    eta1_max = params(2)
    eta2_min = params(3)
    eta2_max = params(4)

    eps = params(5)
    v0=3.0_f64
    !kx  =  2. * sll_pi / (eta1_max - eta1_min)
    kx=0.2_f64
    factor1 = 1.0_f64/sqrt((2.0*sll_pi))

    sll_twostream_1d_xvx_initializer_v1v2x1x2 = factor1/2 * &
         (1.0_f64+eps*cos(kx*x))*(exp(-0.5_f64*((vx-v0)**2))+ &
         exp(-0.5_f64*((vx+v0)**2)))
  end function sll_twostream_1d_xvx_initializer_v1v2x1x2

 function sll_galaxy_1d_xvx_initializer_v1v2x1x2( vx, vy, x, y, params ) 
    sll_real64 :: sll_galaxy_1d_xvx_initializer_v1v2x1x2
    sll_real64, intent(in) :: x
    sll_real64, intent(in) :: y
    sll_real64, intent(in) :: vx
    sll_real64, intent(in) :: vy
    !sll_real64  :: t

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
       print *,'does not depend on y and vy',y,vy
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
    sll_galaxy_1d_xvx_initializer_v1v2x1x2=0.0_f64
    if ((x.lt.1.0_f64).and.(x.gt.-1.0_f64)) then
       sll_galaxy_1d_xvx_initializer_v1v2x1x2 = factor1*exp(-0.5_f64*(vx**2))
    end if
  end function sll_galaxy_1d_xvx_initializer_v1v2x1x2

 function sll_galaxy_2d_initializer_v1v2x1x2( vx, vy, x, y, params ) 
    sll_real64 :: sll_galaxy_2d_initializer_v1v2x1x2
    sll_real64, intent(in) :: x
    sll_real64, intent(in) :: y
    sll_real64, intent(in) :: vx
    sll_real64, intent(in) :: vy
    sll_real64 :: v0x
    sll_real64 :: v0y
    !sll_real64  :: t

    sll_real64, dimension(:), intent(in), optional :: params
    sll_real64 :: eta1_min
    sll_real64 :: eta1_max
    sll_real64 :: eta2_min
    sll_real64 :: eta2_max
    sll_real64 :: a1_xmin,a1_xmax,a1_ymin,a1_ymax,x1mil
    sll_real64 :: a2_xmin,a2_xmax,a2_ymin,a2_ymax,x2mil
    sll_real64 :: eps
    sll_real64 :: kx
    sll_real64 :: factor1

    if( .not. present(params) ) then
       print *, 'sll_landau_1d_initializer_v1v2x1x2, error: the params array', &
            'must be passed params(1)= epsilon, params(2) = kx, params(3) = ky.'
       print *,'does not depend on y and vy',y,vy
       stop
    end if
    !lmin=0.0_f64
    eta1_min = params(1)
    eta1_max = params(2)
    eta2_min = params(3)
    eta2_max = params(4)
    !amas 1
    a1_xmin = -2.0_f64
    a1_xmax = 2.0_f64
    x1mil = 0.0_f64
    a1_ymin = -a1_xmax
    a1_ymax = -a1_xmin
    !amas 2
    a2_xmin = -6.0_f64
    a2_xmax = -2.0_f64
    x2mil=-4.0_f64
    a2_ymin = -a2_xmax
    a2_ymax = -a2_xmin
    
    eps = params(5)
    !kx  =  2. * sll_pi / (eta1_max - eta1_min)
    kx=0.2_f64
    factor1 = 1.0_f64 /sqrt((2.0*sll_pi))**2
    sll_galaxy_2d_initializer_v1v2x1x2=0.0_f64
    !un bloc de particles
!!$    if ((x.le.1.0_f64).and.(x.ge.-1.0_f64).and.(y.le.1.0_f64).and.(y.ge.-1.0_f64)) then
!!$       sll_galaxy_2d_initializer_v1v2x1x2 = factor1*exp(-0.5_f64*(vx**2+vy**2))
!!$    end if
    !2 blocs de particles
    !amas 1
    if ((x.le.a1_xmax).and.(x.ge.a1_xmin).and.(y.le.a1_ymax).and.(y.ge.a1_ymin)) then
       v0x=0.0_f64
       v0y=0.0_f64
    sll_galaxy_2d_initializer_v1v2x1x2 = factor1*exp(-0.5_f64*((vx-v0x)**2+(vy-v0y)**2))*exp(-0.5_f64*((x-x1mil)**2+(y+x1mil)**2)*16)
    end if
    !amas 2
     if ((x.le.a2_xmax).and.(x.ge.a2_xmin).and.(y.le.a2_ymax).and.(y.ge.a2_ymin)) then
       v0x=0 ! 1.0_f64 !0.3109793990   !0.2433894805
       v0y=-0 !.5
        sll_galaxy_2d_initializer_v1v2x1x2 =factor1*exp(-0.5_f64*((vx-v0x)**2+(vy-v0y)**2))*exp(-0.5_f64*((x-x2mil)**2+(y+x2mil)**2)*16)
    end if
!!$ if ((x.le.lmax).and.(x.ge.lmin).and.(y.le.-lmin).and.(y.ge.-lmax)) then
!!$       v0x=-1.0
!!$       v0y=0
!!$       sll_galaxy_2d_initializer_v1v2x1x2 = factor1*exp(-0.5_f64*((vx-v0x)**2+(vy-v0y)**2))*exp(-((x-lmil)**2+(y+lmil)**2))
!!$    end if
!!$    if ((y.le.lmax).and.(y.ge.lmin).and.(x.le.-lmin).and.(x.ge.-lmax)) then
!!$       v0x=1.0
!!$       v0y=0
!!$        sll_galaxy_2d_initializer_v1v2x1x2 = factor1*exp(-0.5_f64*((vx-v0x)**2+(vy-v0y)**2))*exp(-((x+lmil)**2+(y-lmil)**2))
!!$     end if
  !  end if
!!$    if ((x.le.lmax).and.(x.ge.lmin).and.(y.le.-lmin).and.(y.ge.-lmax)) then
!!$       v0x=-1.0
!!$       v0y=0
!!$       sll_galaxy_2d_initializer_v1v2x1x2 = factor1*exp(-0.5_f64*((vx-v0x)**2+(vy-v0y)**2))*exp(-((x-lmil)**2+(y+lmil)**2))
!!$    end if
!!$ if ((y.le.lmax).and.(y.ge.lmin).and.(x.le.-lmin).and.(x.ge.-lmax)) then
!!$       v0x=1.0
!!$       v0y=0
!!$        sll_galaxy_2d_initializer_v1v2x1x2 = factor1*exp(-0.5_f64*((vx-v0x)**2+(vy-v0y)**2))*exp((x+lmil)**2+(y-lmil)**2)
!!$    end if
!!$  
!!$    if ((x.le.3.0_f64).and.(x.ge.1.5_f64).and.(y.le.-1.5_f64).and.(y.ge.-3.0_f64)) then
!!$       v0x=-1.0_f64
!!$       v0y=0.0
!!$       sll_galaxy_2d_initializer_v1v2x1x2 = factor1*exp(-0.5_f64*((vx-v0x)**2+(vy-v0y)**2))
!!$    end if
!!$ if ((y.le.3.0_f64).and.(y.ge.1.5_f64).and.(x.le.-1.5_f64).and.(x.ge.-3.0_f64)) then
!!$       v0x=1.0_f64
!!$       v0y=-0.0
!!$        sll_galaxy_2d_initializer_v1v2x1x2 = factor1*exp(-0.5_f64*((vx-v0x)**2+(vy-v0y)**2))
!!$    end if

!!$    if ((x.le.2.5_f64).and.(x.ge.0.5_f64).and.(y.le.-0.5_f64).and.(y.ge.-2.5_f64)) then
!!$       sll_galaxy_2d_initializer_v1v2x1x2 = factor1*exp(-0.5_f64*(vx**2+vy**2))
!!$    end if
!!$ if ((y.le.2.5_f64).and.(y.ge.0.5_f64).and.(x.le.-0.5_f64).and.(x.ge.-2.5_f64)) then
!!$       sll_galaxy_2d_initializer_v1v2x1x2 = factor1*exp(-0.5_f64*(vx**2+vy**2))
!!$    end if
  end function sll_galaxy_2d_initializer_v1v2x1x2


  function sll_landau_1d_yvy_initializer_v1v2x1x2( vx, vy, x, y, params ) 
    sll_real64 :: sll_landau_1d_yvy_initializer_v1v2x1x2
    sll_real64, intent(in) :: x
    sll_real64, intent(in) :: y
    sll_real64, intent(in) :: vx
    sll_real64, intent(in) :: vy
    !sll_real64  :: t

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
       print *,'does not depend on vx and x',vx,x
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
       print *,'does not depend on x2,x3',x2,x3     
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
    vxc      = params(5)
    vyc      = params(6)
    alpha    = params(7)
    beta     = params(8)
    kx  =  2. * sll_pi / (eta1_max - eta1_min)
    ky  =  2. * sll_pi / (eta2_max - eta2_min)
    


    val = beta*(1 + alpha * (cos(kx*x)+cos(ky*y)) )&
         * exp(-0.5*((vx-vxc)**2+(vy-vyc)**2)) &
         / (2*sll_pi)
    
  end function sll_periodic_periodic_gaussian2002_initializer_4d



 !---------------------------------------------------------------------------
  !
  !                         Gaussian beam 4d initializer
  !  
  ! convention for the params array:
  ! params(1) = vth
  ! params(2) = xth
  ! params(3) = vxc
  ! params(4) = vyc
  ! params(5) = xc
  ! params(6) = yc
  ! params(7) = n0
  ! params(8) = radius
  !---------------------------------------------------------------------------
  
  function sll_gaussian_beam_initializer_4d( x, y, vx, vy, params ) &
       result(val)
    
    sll_real64 :: val
    sll_real64, intent(in) :: x
    sll_real64, intent(in) :: y
    sll_real64, intent(in) :: vx
    sll_real64, intent(in) :: vy
    
    sll_real64, dimension(:), intent(in), optional :: params
    sll_real64 :: vxc
    sll_real64 :: vyc
    sll_real64 :: xc,yc,radius
    sll_real64 :: vt,xt,n0
    
    if( .not. present(params) ) then
       print *, 'sll_gaussian_initializer_4d, error: the params array must ', &
            'be passed: ', &
            'params(1) = vt, params(2) = xt, params(3) = sigma_x, params(4) = sigma_v',&
            ' params(5) = vxc, params(6) = vyc, params(7) = xc, params(8) = yc, params(9) = n0'
       stop
    end if
    
    vt      = params(1)
    xt      = params(2)
    vxc     = params(3)
    vyc     = params(4)
    xc      = params(5)
    yc      = params(6)
    n0      = params(7)
    radius   = params(8)
   
    
    val = n0 *exp(-0.5*(radius**2*(x-xc)**2  + radius**2*(y-yc)**2)/(xt*xt)  ) &
         / (2*sll_pi*xt**2) &
         *exp(-0.5*((vx-vxc)**2+(vy-vyc)**2)/(vt*vt))/ (2*sll_pi*vt**2)
    
  end function sll_gaussian_beam_initializer_4d

end module sll_common_array_initializers_module
