module m_init_functions

#include "sll_working_precision.h"

  use sll_m_constants, only: sll_pi

contains

  function func_one( eta1, eta2, params ) result(res)
    real(8), intent(in) :: eta1
    real(8), intent(in) :: eta2
    real(8), intent(in) :: params(:)
    real(8) :: res
#ifdef DEBUG
    real(8) :: dummy
    if(size(params)>0) dummy = eta1+eta2
#endif
    res = 1.0_8
  end function func_one

  function func_zero( eta1, eta2, params ) result(res)
    real(8), intent(in) :: eta1
    real(8), intent(in) :: eta2
    real(8), intent(in) :: params(:)
    real(8) :: res
#ifdef DEBUG
    real(8) :: dummy
    if(size(params)>0) dummy = eta1+eta2
#endif
    res = 0.0_8
  end function func_zero

  function func_four( eta1, eta2, params ) result(res)
    real(8), intent(in) :: eta1
    real(8), intent(in) :: eta2
    real(8), dimension(:), intent(in) :: params
    real(8) :: res
#ifdef DEBUG
    real(8) :: dummy
    if (size(params)>0) dummy = eta1+eta2
#endif
    res = 4.0_8
  end function func_four

  function func_epsi( eta1, eta2, params ) result(res)
    real(8), intent(in) :: eta1
    real(8), intent(in) :: eta2
    real(8), intent(in) :: params(:)
    real(8) :: res
#ifdef DEBUG
    real(8) :: dummy
    if (size(params)>0) dummy = eta1+eta2
#endif

    res = 0.0_8
  end function func_epsi

  !----------------------------------------------------------
  !  Solution for a identity change of coordinates 
  !   and periodic-periodic conditions
  !   the matrix A is equal to identity 
  !   the scalar c is equal to zero 
  !-------------------------------------------------------------
  function source_term_perper( eta1, eta2, params ) result(res)
    real(8), intent(in) :: eta1
    real(8), intent(in) :: eta2
    real(8), intent(in) :: params(:)
    real(8) :: res

    res =  cos(2*sll_pi*eta1)
    !!-2*(2.0*sll_pi)**2*cos(2.0*sll_pi*eta1)*cos(2.0*sll_pi*eta2)! 0.001*cos(2*sll_pi*eta1)!
  end function source_term_perper

  real(8) function sol_exacte_perper( eta1, eta2 )
    real(8), intent(in) :: eta1
    real(8), intent(in) :: eta2
    
    sol_exacte_perper = -1.0_8/((2*sll_pi)**2)*cos(2*sll_pi*eta1)!cos(2.0*sll_pi*eta1)*cos(2.0*sll_pi*eta2)!-0.001/((2*sll_pi)**2)*cos(2*sll_pi*eta1)
  end function sol_exacte_perper

  real(8) function sol_exacte_perper_der1( eta1, eta2 )
    real(8), intent(in) :: eta1
    real(8), intent(in) :: eta2
    
    sol_exacte_perper_der1 = 1.0_8/(2*sll_pi)*sin(2*sll_pi*eta1) !-2.0*sll_pi*sin(2.0*sll_pi*eta1)*cos(2.0*sll_pi*eta2)
  end function sol_exacte_perper_der1

  real(8) function sol_exacte_perper_der2( eta1, eta2 )
    real(8), intent(in) :: eta1
    real(8), intent(in) :: eta2

    sol_exacte_perper_der2 = 0.0_8!-2.0*sll_pi*cos(2.0*sll_pi*eta1)*sin(2.0*sll_pi*eta2)
  end function sol_exacte_perper_der2

  !----------------------------------------------------------
  !  Solution for a identity change of coordinates 
  !   and periodic-dirichlet conditions
  !   and also dirichlet-dirichlet conditons
  !   the matrix A is equal to identity 
  !   the scalar c is equal to zero 
  !-------------------------------------------------------------
  real(8) function source_term_perdir( eta1, eta2, params ) ! in the path
    real(8), intent(in) :: eta1
    real(8), intent(in) :: eta2
    real(8), intent(in) :: params(:)
!#ifdef DEBUG
!    real(8) :: dummy
!    if (present(params)) dummy = eta1+eta2
!#endif

    source_term_perdir = -2*(2*sll_pi)**2* sin(2*sll_pi*eta1)*sin(2*sll_pi*eta2)
        ! -(16.0*sll_pi**2*eta2**4 &
        ! - 16.0*sll_pi**2*eta2**2 &
        ! - 12.0*eta2**2 + 2.0)*cos(2*sll_pi*eta1)*sin(2*sll_pi*eta1)
  end function source_term_perdir


  real(8) function sol_exacte_perdir( eta1, eta2 )
    real(8), intent(in) :: eta1
    real(8), intent(in) :: eta2

    sol_exacte_perdir = sin(2.0*sll_pi*eta1)*sin(2.0*sll_pi*eta2)!eta2 ** 2 * (eta2**2-1)&
        ! * cos(2.0*sll_pi*eta1)*sin(2.0*sll_pi*eta1)
  end function sol_exacte_perdir


  real(8) function sol_exacte_perdir_der1( eta1, eta2 )
    real(8), intent(in) :: eta1
    real(8), intent(in) :: eta2

    sol_exacte_perdir_der1 = 2.0*sll_pi*cos(2.0*sll_pi*eta1)*sin(2.0*sll_pi*eta2)
  end function sol_exacte_perdir_der1


  real(8) function sol_exacte_perdir_der2( eta1, eta2 )
    real(8), intent(in) :: eta1
    real(8), intent(in) :: eta2

    sol_exacte_perdir_der2 = 2.0*sll_pi*sin(2.0*sll_pi*eta1)*cos(2.0*sll_pi*eta2)
  end function sol_exacte_perdir_der2

  !-------------------------------------------------------------
  !  Solution for a identity change of coordinates 
  !   and also dirichlet-periodicconditons
  !   the matrix A is equal to identity 
  !   the scalar c is equal to zero 
  !-------------------------------------------------------------
  real(8) function source_term_dirper( eta1, eta2, params ) ! in the path
    real(8), intent(in) :: eta1
    real(8), intent(in) :: eta2
    real(8), intent(in) :: params(:)
!#ifdef DEBUG
!    real(8) :: dummy
!    if (present(params)) dummy = params(1)
!#endif

    source_term_dirper = -2*(2*sll_pi)**2* sin(2*sll_pi*eta1)*cos(2*sll_pi*eta2)
       ! -(16.0*sll_pi**2*eta1**4 &
       ! - 16.0*sll_pi**2*eta1**2 &
       ! - 12.0*eta1**2 + 2.0)*sin(2*sll_pi*eta2)*cos(2*sll_pi*eta2)
  end function source_term_dirper


  real(8) function sol_exacte_dirper( eta1, eta2 )
    real(8), intent(in) :: eta1
    real(8), intent(in) :: eta2

    sol_exacte_dirper = sin(2.0*sll_pi*eta1)*cos(2.0*sll_pi*eta2)
         !eta1 ** 2 * (eta1**2-1)* cos(2*sll_pi*eta2)*sin(2*sll_pi*eta2)
  end function sol_exacte_dirper

  real(8) function sol_exacte_dirper_der1( eta1, eta2 )
    real(8), intent(in) :: eta1
    real(8), intent(in) :: eta2

    sol_exacte_dirper_der1 = 2*sll_pi*cos(2.0*sll_pi*eta1)*cos(2.0*sll_pi*eta2)
         !eta1 ** 2 * (eta1**2-1)* cos(2*sll_pi*eta2)*sin(2*sll_pi*eta2)
  end function sol_exacte_dirper_der1

  real(8) function sol_exacte_dirper_der2( eta1, eta2 )
    real(8), intent(in) :: eta1
    real(8), intent(in) :: eta2

    sol_exacte_dirper_der2 = -2.0*sll_pi*sin(2.0*sll_pi*eta1)*sin(2.0*sll_pi*eta2)
         !eta1 ** 2 * (eta1**2-1)* cos(2*sll_pi*eta2)*sin(2*sll_pi*eta2)
  end function sol_exacte_dirper_der2

  !----------------------------------------------------------
  !  Solution for a r theta change of coordinates 
  !   and periodic-dirichlet conditions
  !   and also dirivhlet-dirichlet conditons
  !   the matrix A is equal to identity 
  !   the scalar c is equal to zero 
  !-------------------------------------------------------------
  real(8) function rho_rtheta( eta1, eta2, params ) ! in the path
    real(8), intent(in) :: eta1
    real(8), intent(in) :: eta2
    real(8), intent(in) :: params(:)
    real(8) :: x, y
    
!#ifdef DEBUG
!    if (present(params)) print*, params
!#endif

    x = eta2*cos(2*sll_pi*eta1)
    y = eta2*sin(2*sll_pi*eta1)
    
    rho_rtheta = x*y*(-32.0*x**2 - 32.0*y**2 + 15.0)  
    
  end function rho_rtheta


  real(8) function sol_exacte_rtheta( eta1, eta2, params ) ! in the path
    real(8), intent(in) :: eta1
    real(8), intent(in) :: eta2
    real(8), intent(in) :: params(:)
    
!#ifdef DEBUG
!    if (present(params)) print*, params
!#endif
    
    sol_exacte_rtheta = ( eta2**2-1)*(eta2**2-0.5**2)*eta2**2&
         *cos(2*sll_pi*eta1)*sin(2*sll_pi*eta1)
    
  end function sol_exacte_rtheta

  !----------------------------------------------------------
  !  Solution for a colella change of coordinates 
  !   and periodic-periodic conditons
  !   the matrix A is equal to identity 
  !   the scalar c is equal to zero 
  !-------------------------------------------------------------
  real(8) function source_term_chgt_perper( eta1, eta2, params ) ! in the path
    real(8), intent(in) :: eta1
    real(8), intent(in) :: eta2
    real(8), intent(in) :: params(:)

    real(8) :: x, y
    
    x =   eta1 + 0.1_8*sin(2*sll_pi*eta1) * sin(2*sll_pi*eta2)
    y =   eta2 + 0.1_8*sin(2*sll_pi*eta1) * sin(2*sll_pi*eta2)
    
    source_term_chgt_perper = -8.0*sll_pi**2*cos(2*sll_pi*x)*cos(2*sll_pi*y) 
    
  end function source_term_chgt_perper

  real(8) function sol_exacte_chgt_perper( eta1, eta2 )
    real(8), intent(in) :: eta1
    real(8), intent(in) :: eta2

    real(8) :: x,y
    
    x =   eta1 + 0.1*sin(2* sll_pi*eta1) * sin(2*sll_pi*eta2)
    y =   eta2 + 0.1*sin(2* sll_pi*eta1) * sin(2*sll_pi*eta2)
    
    sol_exacte_chgt_perper = cos(2*sll_pi*x)*cos(2*sll_pi*y)

  end function sol_exacte_chgt_perper

  real(8) function sol_exacte_chgt_perper_der1( eta1, eta2 )
    real(8), intent(in) :: eta1
    real(8), intent(in) :: eta2

    real(8) :: x,y
    
    x =   eta1 + 0.1*sin(2* sll_pi*eta1) * sin(2*sll_pi*eta2)
    y =   eta2 + 0.1*sin(2* sll_pi*eta1) * sin(2*sll_pi*eta2)
    
    sol_exacte_chgt_perper_der1 = -2*sll_pi*sin(2*sll_pi*x)*cos(2*sll_pi*y)&
         * ( 1.0_8 + 0.1*2*sll_pi*cos(2* sll_pi*eta1) * sin(2*sll_pi*eta2) )&
         -2*sll_pi*cos(2*sll_pi*x)*sin(2*sll_pi*y)&
         * ( 0.1*2*sll_pi*cos(2* sll_pi*eta1) * sin(2*sll_pi*eta2) )
 
  end function sol_exacte_chgt_perper_der1

  real(8) function sol_exacte_chgt_perper_der2( eta1, eta2 )
    real(8), intent(in) :: eta1
    real(8), intent(in) :: eta2

    real(8) :: x,y
    
    x =   eta1 + 0.1*sin(2* sll_pi*eta1) * sin(2*sll_pi*eta2)
    y =   eta2 + 0.1*sin(2* sll_pi*eta1) * sin(2*sll_pi*eta2)
    
    !sol_exacte_chgt_perper_der2 = -2*sll_pi*cos(2*sll_pi*x)*sin(2*sll_pi*y)
    
    sol_exacte_chgt_perper_der2 = -2*sll_pi*sin(2*sll_pi*x)*cos(2*sll_pi*y)&
         * ( 0.1*2*sll_pi*sin(2* sll_pi*eta1) * cos(2*sll_pi*eta2) )&
         -2*sll_pi*cos(2*sll_pi*x)*sin(2*sll_pi*y)&
         * ( 1.0_8 + 0.1*2*sll_pi*sin(2* sll_pi*eta1)*cos(2*sll_pi*eta2) )
  end function sol_exacte_chgt_perper_der2

  !----------------------------------------------------------
  !  Solution for a colella change of coordinates 
  !   and periodic-dirichlet conditons
  !   and dircihlet-diichlet conditions
  !   the matrix A is equal to identity 
  !   the scalar c is equal to zero 
  !-------------------------------------------------------------
  real(8) function source_term_chgt_perdir( eta1, eta2, params ) ! in the path
    real(8), intent(in) :: eta1
    real(8), intent(in) :: eta2
    real(8), intent(in) :: params(:)

    real(8) :: x, y

    x =   eta1 + 0.1*sin(2* sll_pi*eta1) * sin(2*sll_pi*eta2)
    y =   eta2 + 0.1*sin(2* sll_pi*eta1) * sin(2*sll_pi*eta2)
    
    source_term_chgt_perdir= -2*(2*sll_pi)**2 * sin(2*sll_pi*y)*cos(2*sll_pi*x)
  end function source_term_chgt_perdir

  real(8) function sol_exacte_chgt_perdir( eta1, eta2 )
    real(8), intent(in) :: eta1
    real(8), intent(in) :: eta2

    real(8) :: x,y

    x =   eta1 + 0.1*sin(2* sll_pi*eta1) * sin(2*sll_pi*eta2)
    y =   eta2 + 0.1*sin(2* sll_pi*eta1) * sin(2*sll_pi*eta2)
    
    sol_exacte_chgt_perdir = cos(2*sll_pi*x)*sin(2*sll_pi*y)
  end function sol_exacte_chgt_perdir

  real(8) function sol_exacte_chgt_perdir_der1( eta1, eta2 )
    real(8), intent(in) :: eta1
    real(8), intent(in) :: eta2

    real(8) :: x,y

    x =   eta1 + 0.1*sin(2* sll_pi*eta1) * sin(2*sll_pi*eta2)
    y =   eta2 + 0.1*sin(2* sll_pi*eta1) * sin(2*sll_pi*eta2)
    
    sol_exacte_chgt_perdir_der1 = -2*sll_pi*sin(2*sll_pi*x)*sin(2*sll_pi*y)&
         * ( 1.0_8 + 0.1*2*sll_pi*cos(2*sll_pi*eta1) * sin(2*sll_pi*eta2) )&
         + 2*sll_pi*cos(2*sll_pi*x)*cos(2*sll_pi*y)&
         * ( 2*sll_pi*0.1*cos(2* sll_pi*eta1) * sin(2*sll_pi*eta2) ) 
    
  end function sol_exacte_chgt_perdir_der1

  real(8) function sol_exacte_chgt_perdir_der2( eta1, eta2 )
    real(8), intent(in) :: eta1
    real(8), intent(in) :: eta2

    real(8) :: x,y

    x =   eta1 + 0.1*sin(2* sll_pi*eta1) * sin(2*sll_pi*eta2)
    y =   eta2 + 0.1*sin(2* sll_pi*eta1) * sin(2*sll_pi*eta2)
    
    sol_exacte_chgt_perdir_der2 = -2*sll_pi*sin(2*sll_pi*x)*sin(2*sll_pi*y)&
         * ( 0.1*2*sll_pi*sin(2*sll_pi*eta1) * cos(2*sll_pi*eta2) ) &
         + 2*sll_pi*cos(2*sll_pi*x)*cos(2*sll_pi*y)&
         * ( 1.0_8 + 2*sll_pi*0.1*sin(2* sll_pi*eta1) *cos(2*sll_pi*eta2) ) 
    
  end function sol_exacte_chgt_perdir_der2

  !----------------------------------------------------------
  !  Solution for a colella change of coordinates 
  !   and dirchlet-periodic conditions
  !   the matrix A is equal to identity 
  !   the scalar c is equal to zero 
  !-------------------------------------------------------------
  real(8) function source_term_chgt_dirdir( eta1, eta2, params ) ! in the path
    real(8), intent(in) :: eta1
    real(8), intent(in) :: eta2
    real(8), intent(in) :: params(:)

    real(8) :: x, y
    ! -------------------------------------------------
    ! In the case without change of coordinates
    ! -------------------------------------------------
    x =   eta1 + 0.1_8*sin(2* sll_pi*eta1) * sin(2*sll_pi*eta2)
    y =   eta2 + 0.1_8*sin(2* sll_pi*eta1) * sin(2*sll_pi*eta2)
    
    source_term_chgt_dirdir = &
         -2*(2.0*sll_pi)**2*sin(2*sll_pi*x)*sin(2*sll_pi*y)
    
  end function source_term_chgt_dirdir

  real(8) function sol_exacte_chgt_dirdir( eta1, eta2 )
    real(8), intent(in) :: eta1
    real(8), intent(in) :: eta2

    real(8) :: x,y

    x =   eta1 + 0.1*sin(2* sll_pi*eta1) * sin(2*sll_pi*eta2)
    y =   eta2 + 0.1*sin(2* sll_pi*eta1) * sin(2*sll_pi*eta2)

    sol_exacte_chgt_dirdir = sin(2* sll_pi*y)*sin(2* sll_pi*x)

  end function sol_exacte_chgt_dirdir

  real(8) function sol_exacte_chgt_dirdir_der1( eta1, eta2 )
    real(8), intent(in) :: eta1
    real(8), intent(in) :: eta2

    real(8) :: x,y

    x =   eta1 + 0.1*sin(2* sll_pi*eta1) * sin(2*sll_pi*eta2)
    y =   eta2 + 0.1*sin(2* sll_pi*eta1) * sin(2*sll_pi*eta2)

    sol_exacte_chgt_dirdir_der1 = 2*sll_pi*cos(2* sll_pi*x)*sin(2* sll_pi*y)&
         * ( 1.0_8 + 0.1*2*sll_pi*cos(2*sll_pi*eta1) * sin(2*sll_pi*eta2) )&
         + 2*sll_pi*sin(2* sll_pi*x)*cos(2* sll_pi*y) &
         * ( 2*sll_pi*0.1*cos(2* sll_pi*eta1) * sin(2*sll_pi*eta2) )
  end function sol_exacte_chgt_dirdir_der1

  real(8) function sol_exacte_chgt_dirdir_der2( eta1, eta2 )
    real(8), intent(in) :: eta1
    real(8), intent(in) :: eta2
 
    real(8) :: x,y
    
    x =   eta1 + 0.1*sin(2* sll_pi*eta1) * sin(2*sll_pi*eta2)
    y =   eta2 + 0.1*sin(2* sll_pi*eta1) * sin(2*sll_pi*eta2)
    
    sol_exacte_chgt_dirdir_der2 =  2*sll_pi*cos(2* sll_pi*x)*sin(2* sll_pi*y)&
         * ( 0.1*2*sll_pi*sin(2*sll_pi*eta1) * cos(2*sll_pi*eta2)  )&
         + 2*sll_pi*sin(2* sll_pi*x)*cos(2* sll_pi*y) &
         * ( 1.0_8 + 2*sll_pi*0.1*sin(2* sll_pi*eta1) *cos(2*sll_pi*eta2) )
    
  end function sol_exacte_chgt_dirdir_der2

  !----------------------------------------------------------
  !  Solution for a colella change of coordinates 
  !   and dirchlet-periodic conditions
  !   the matrix A is equal to identity 
  !   the scalar c is equal to zero 
  !-------------------------------------------------------------
  real(8) function source_term_chgt_dirper( eta1, eta2, params ) ! in the path
    real(8), intent(in) :: eta1
    real(8), intent(in) :: eta2
    real(8), intent(in) :: params(:)

    real(8) :: x, y
    ! -------------------------------------------------
    ! In the case without change of coordinates
    ! -------------------------------------------------
    x =   eta1 + 0.1_8*sin(2* sll_pi*eta1) * sin(2*sll_pi*eta2)
    y =   eta2 + 0.1_8*sin(2* sll_pi*eta1) * sin(2*sll_pi*eta2)
    
    source_term_chgt_dirper = -2*(2*sll_pi)**2*sin(2*sll_pi*x)*cos(2*sll_pi*y)
    
  end function source_term_chgt_dirper


  real(8) function sol_exacte_chgt_dirper( eta1, eta2 )
    real(8), intent(in) :: eta1
    real(8), intent(in) :: eta2

    real(8) :: x,y

    x =   eta1 + 0.1*sin(2* sll_pi*eta1) * sin(2*sll_pi*eta2)
    y =   eta2 + 0.1*sin(2* sll_pi*eta1) * sin(2*sll_pi*eta2)
    
    sol_exacte_chgt_dirper = sin(2* sll_pi*x)*cos(2* sll_pi*y)
    
  end function sol_exacte_chgt_dirper


  real(8) function sol_exacte_chgt_dirper_der1( eta1, eta2 )
    real(8), intent(in) :: eta1
    real(8), intent(in) :: eta2

    real(8) :: x,y
    
    x =   eta1 + 0.1*sin(2* sll_pi*eta1) * sin(2*sll_pi*eta2)
    y =   eta2 + 0.1*sin(2* sll_pi*eta1) * sin(2*sll_pi*eta2)
    
    sol_exacte_chgt_dirper_der1 = 2*sll_pi*cos(2* sll_pi*x)*cos(2* sll_pi*y) &
         * ( 1.0_8 + 0.1*2*sll_pi*cos(2*sll_pi*eta1) * sin(2*sll_pi*eta2) )&
         - 2*sll_pi*sin(2* sll_pi*x)*sin(2* sll_pi*y)&
         * ( 2*sll_pi*0.1*cos(2* sll_pi*eta1) * sin(2*sll_pi*eta2) ) 
  end function sol_exacte_chgt_dirper_der1


  real(8) function sol_exacte_chgt_dirper_der2( eta1, eta2 )
    real(8), intent(in) :: eta1
    real(8), intent(in) :: eta2

    real(8) :: x,y
    
    x =   eta1 + 0.1*sin(2* sll_pi*eta1) * sin(2*sll_pi*eta2)
    y =   eta2 + 0.1*sin(2* sll_pi*eta1) * sin(2*sll_pi*eta2)
    
    sol_exacte_chgt_dirper_der2 = 2*sll_pi*cos(2* sll_pi*x)*cos(2* sll_pi*y) &
         * ( 0.1*2*sll_pi*sin(2*sll_pi*eta1) * cos(2*sll_pi*eta2)  )&
         - 2*sll_pi*sin(2* sll_pi*x)*sin(2* sll_pi*y)&
         * (1.0_8 + 2*sll_pi*0.1*sin(2* sll_pi*eta1) *cos(2*sll_pi*eta2) ) 
    
  end function sol_exacte_chgt_dirper_der2


  !!!!!! test case with F(theta,phi) = (2pi theta , 2pi phi)

  real(8) function adimension_chgt_x( eta1, eta2, params )
    real(8), intent(in) :: eta1
    real(8), intent(in) :: eta2
    real(8), intent(in) :: params(:)

    print*, eta1, eta2
    adimension_chgt_x = 2*sll_pi*eta1 !+ eta2)
  end function adimension_chgt_x

  real(8) function adimension_chgt_y( eta1, eta2, params )
    real(8), intent(in) :: eta1
    real(8), intent(in) :: eta2
    real(8), intent(in) :: params(:)

    print*, eta1, eta2
    adimension_chgt_y = 2*sll_pi*eta2
  end function adimension_chgt_y

  real(8) function jac11_adimension_chgt( eta1, eta2, params )
    real(8), intent(in) :: eta1
    real(8), intent(in) :: eta2
    real(8), intent(in) :: params(:)

    print*, eta1, eta2
    jac11_adimension_chgt = 2*sll_pi
  end function jac11_adimension_chgt

  real(8) function jac12_adimension_chgt( eta1, eta2, params )
    real(8), intent(in) :: eta1
    real(8), intent(in) :: eta2
    real(8), intent(in) :: params(:)

    print*, eta1, eta2
    jac12_adimension_chgt = 0.0_8!sll_pi
  end function jac12_adimension_chgt

  real(8) function jac21_adimension_chgt( eta1, eta2, params )
    real(8), intent(in) :: eta1
    real(8), intent(in) :: eta2
    real(8), intent(in) :: params(:)

    print*, eta1, eta2
    jac21_adimension_chgt = 0.0_8!2*sll_pi!0.0
  end function jac21_adimension_chgt

  real(8) function jac22_adimension_chgt( eta1, eta2, params )
    real(8), intent(in) :: eta1
    real(8), intent(in) :: eta2
    real(8), intent(in) :: params(:)

    print*, eta1, eta2
    jac22_adimension_chgt = 2*sll_pi
  end function jac22_adimension_chgt

  real(8) function sol_exacte_chgt_adim( eta1, eta2, params )
    real(8), intent(in) :: eta1
    real(8), intent(in) :: eta2
    real(8), intent(in) :: params(:)

    real(8) :: x,y
    
    x =   2*sll_pi*eta1!+eta2)
    y =   2* sll_pi*eta2
    
    sol_exacte_chgt_adim = cos(x)*cos(y)
    
  end function sol_exacte_chgt_adim

  real(8) function source_term_chgt_adim( eta1, eta2, params )
    real(8), intent(in) :: eta1
    real(8), intent(in) :: eta2
    real(8), intent(in) :: params(:)

    real(8) :: x,y
    
    x =   2.0_8*sll_pi*eta1 !+eta2)
    y =   2.0_8*sll_pi*eta2
    
    source_term_chgt_adim = -2.0_8*cos(x)*cos(y)
    
  end function source_term_chgt_adim

#define R_MIN 1.0_8
#define R_MAX 2.0_8
#define N_MOD 4

  real(8) function f_cos( eta1, eta2, params )
    real(8), intent(in) :: eta1
    real(8), intent(in) :: eta2
    real(8), intent(in) :: params(:)

    !sage: assume(r>=1)
    !sage: assume(r<=2)
    !sage: phi = (r-r_min)*(r-r_max)*r*cos(n*theta)
    !sage: diff(r*diff(phi,r),r)/r + diff(phi,theta,theta)/(r*r)

    real(8) :: r, theta
    integer :: n = N_MOD 

    r = eta1 * (R_MAX-R_MIN) + R_MIN
    theta = eta2 * 2.0_8 * sll_pi

    f_cos = -(r-R_MAX)*(r-R_MIN)*n*n*cos(n*theta)/r &
            + ((r-R_MAX)*(r-R_MIN)*cos(n*theta)  &
            + (r-R_MAX)*r*cos(n*theta) + (r-R_MIN)*r*cos(n*theta) &
            + 2.0_8*((r-R_MAX)*cos(n*theta) + (r-R_MIN)*cos(n*theta) &
            + r*cos(n*theta))*r)/r

  end function f_cos

  real(8) function f_sin( eta1, eta2, params )
    real(8), intent(in) :: eta1
    real(8), intent(in) :: eta2
    real(8), intent(in) :: params(:)

    !sage: assume(r>=1)
    !sage: assume(r<=2)
    !sage: phi = (r-r_min)*(r-r_max)*r*sin(n*theta)
    !sage: diff(r*diff(phi,r),r)/r + diff(phi,theta,theta)/(r*r)

    real(8) :: r, theta
    integer :: n = N_MOD 

    r = eta1 * (R_MAX-R_MIN) + R_MIN
    theta = eta2 * 2 * sll_pi

    f_sin = -(r-R_MAX)*(r-R_MIN)*n*n*sin(n*theta)/r &
          + ((r-R_MAX)*(r-R_MIN)*sin(n*theta) &
          + (r-R_MAX)*r*sin(n*theta) + (r-R_MIN)*r*sin(n*theta) &
          + 2.0_8*((r-R_MAX)*sin(n*theta) + (r-R_MIN)*sin(n*theta)  &
          + r*sin(n*theta))*r)/r

  end function f_sin
          
  real(8) function u_sin_der1( eta1, eta2 )
    real(8), intent(in) :: eta1
    real(8), intent(in) :: eta2

    !sage: assume(r>=1)
    !sage: assume(r<=2)
    !sage: phi = (r-r_min)*(r-r_max)*r*sin(n*theta)
    !sage: diff(r*diff(phi,r),r)/r + diff(phi,theta,theta)/(r*r)

    real(8) :: r, theta
    integer :: n = N_MOD 
    
    r = eta1 * (R_MAX-R_MIN) + R_MIN
    theta = eta2 * 2 * sll_pi

    u_sin_der1 = (r - R_MAX)*(r - R_MIN)*sin(n*theta) &
               + (r - R_MAX)*r*sin(n*theta) &
               + (r - R_MIN)*r*sin(n*theta)

  end function u_sin_der1

  real(8) function u_sin_der2( eta1, eta2 )
    real(8), intent(in) :: eta1
    real(8), intent(in) :: eta2

    !sage: assume(r>=1)
    !sage: assume(r<=2)
    !sage: phi = (r-r_min)*(r-r_max)*r*sin(n*theta)
    !sage: diff(r*diff(phi,r),r)/r + diff(phi,theta,theta)/(r*r)

    real(8) :: r, theta
    integer :: n = N_MOD 
    
    r = eta1 + R_MIN
    theta = eta2 * 2.0_8 * sll_pi

    u_sin_der2 = n*(r - R_MAX)*(r - R_MIN)*r*cos(n*theta)

  end function u_sin_der2

  real(8) function u_sin( eta1, eta2 )
    real(8), intent(in) :: eta1
    real(8), intent(in) :: eta2

    !sage: assume(r>=1)
    !sage: assume(r<=2)
    !sage: phi = (r-r_min)*(r-r_max)*r*sin(n*theta)
    !sage: diff(r*diff(phi,r),r)/r + diff(phi,theta,theta)/(r*r)

    real(8) :: r, theta
    integer :: n = N_MOD 
    
    r = eta1 * (R_MAX-R_MIN) + R_MIN
    theta = eta2 * 2 * sll_pi

    u_sin = (r-R_MIN)*(r-R_MAX)*sin(n*theta)*r

  end function u_sin

  real(8) function u_cos( eta1, eta2 )
    real(8), intent(in) :: eta1
    real(8), intent(in) :: eta2

    !sage: assume(r>=1)
    !sage: assume(r<=2)
    !sage: phi = (r-r_min)*(r-r_max)*r*cos(n*theta)
    !sage: diff(r*diff(phi,r),r)/r + diff(phi,theta,theta)/(r*r)

    real(8) :: r, theta
    integer :: n = N_MOD
    
    r = eta1 * (R_MAX-R_MIN) + R_MIN
    theta = eta2 * 2 * sll_pi

    u_cos = (r-R_MIN)*(r-R_MAX)*cos(n*theta)*r

  end function u_cos

  real(8) function u_cos_der1( eta1, eta2 )
    real(8), intent(in) :: eta1
    real(8), intent(in) :: eta2

    !sage: assume(r>=1)
    !sage: assume(r<=2)
    !sage: phi = (r-r_min)*(r-r_max)*r*sin(n*theta)
    !sage: diff(r*diff(phi,r),r)/r + diff(phi,theta,theta)/(r*r)

    real(8) :: r, theta
    integer :: n = N_MOD
    
    r = eta1 * (R_MAX-R_MIN) + R_MIN
    theta = eta2 * 2 * sll_pi
    u_cos_der1 = (r - R_MAX)*(r - R_MIN)*cos(n*theta) &
               + (r - R_MAX)*r*cos(n*theta) &
               + (r - R_MIN)*r*cos(n*theta)

  end function u_cos_der1

  real(8) function u_cos_der2( eta1, eta2 )
    real(8), intent(in) :: eta1
    real(8), intent(in) :: eta2

    !sage: assume(r>=1)
    !sage: assume(r<=2)
    !sage: phi = (r-r_min)*(r-r_max)*r*sin(n*theta)
    !sage: diff(r*diff(phi,r),r)/r + diff(phi,theta,theta)/(r*r)

    real(8) :: r, theta
    integer :: n = N_MOD
    
    r = eta1 * (R_MAX-R_MIN) + R_MIN
    theta = eta2 * 2 * sll_pi
    u_cos_der2 = - n*(r - R_MAX)*(r - R_MIN)*r*sin(n*theta)

  end function u_cos_der2

end module m_init_functions
