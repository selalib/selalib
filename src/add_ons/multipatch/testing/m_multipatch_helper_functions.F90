module m_multipatch_helper_functions

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  use sll_m_constants, only: &
    sll_p_pi

  implicit none

  public :: &
    func_one, &
    func_zero

  private
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

contains

! External functions used as parameters in the above unit test:

function func_one( eta1, eta2) result(res)
  real(8), intent(in) :: eta1
  real(8), intent(in) :: eta2
  real(8) :: res
  real(8) :: dummy
  dummy = eta1+eta2
  res = 1.0_8
end function func_one

function func_zero( eta1, eta2) result(res)
  real(8), intent(in) :: eta1
  real(8), intent(in) :: eta2
  real(8) :: res
  real(8) :: dummy
  dummy = eta1+eta2
  res = 0.0_8
end function func_zero

function func_epsi( eta1, eta2, params ) result(res)
  real(8), intent(in) :: eta1
  real(8), intent(in) :: eta2
  real(8), dimension(:), intent(in) :: params
  real(8) :: res
  real(8) :: dummy
  dummy = eta1+eta2+params(1)

  res = 0.0_8
end function func_epsi

!----------------------------------------------------------
!  Solution for a identity change of coordinates 
!   and periodic-periodic conditions
!   the matrix A is equal to identity 
!   the scalar c is equal to zero 
!-------------------------------------------------------------
function source_term_perper( eta1, eta2) result(res)
  real(8), intent(in) :: eta1
  real(8), intent(in) :: eta2
  ! real(8), dimension(:), intent(in), optional :: params
  real(8) :: res
  real(8) :: dummy
  dummy = eta1+eta2

  res =  0.001*cos(2*sll_p_pi*eta1)
  !!-2*(2.0*sll_p_pi)**2*cos(2.0*sll_p_pi*eta1)*cos(2.0*sll_p_pi*eta2)! 0.001*cos(2*sll_p_pi*eta1)!
end function source_term_perper

real(8) function sol_exacte_perper(eta1,eta2)
  real(8), intent(in) :: eta1
  real(8), intent(in) :: eta2
  real(8) :: dummy
  dummy = eta1+eta2
  
  !real(8), dimension(:), intent(in), optional :: params
  sol_exacte_perper = -0.001/((2*sll_p_pi)**2)*cos(2*sll_p_pi*eta1)!cos(2.0*sll_p_pi*eta1)*cos(2.0*sll_p_pi*eta2)!-0.001/((2*sll_p_pi)**2)*cos(2*sll_p_pi*eta1)
end function sol_exacte_perper

real(8) function sol_exacte_perper_der1(eta1,eta2)
  real(8), intent(in) :: eta1
  real(8), intent(in) :: eta2
  real(8) :: dummy
  dummy = eta1+eta2
  
  !real(8), dimension(:), intent(in), optional :: params
  sol_exacte_perper_der1 = 0.001/(2*sll_p_pi)*sin(2*sll_p_pi*eta1) !-2.0*sll_p_pi*sin(2.0*sll_p_pi*eta1)*cos(2.0*sll_p_pi*eta2)
end function sol_exacte_perper_der1

real(8) function sol_exacte_perper_der2(eta1,eta2)
  real(8), intent(in) :: eta1
  real(8), intent(in) :: eta2
  real(8) :: dummy
  dummy = eta1+eta2
  
  !real(8), dimension(:), intent(in), optional :: params
  sol_exacte_perper_der2 = 0.0_8!-2.0*sll_p_pi*cos(2.0*sll_p_pi*eta1)*sin(2.0*sll_p_pi*eta2)
end function sol_exacte_perper_der2

!----------------------------------------------------------
!  Solution for a identity change of coordinates 
!   and periodic-dirichlet conditions
!   and also dirichlet-dirichlet conditons
!   the matrix A is equal to identity 
!   the scalar c is equal to zero 
!-------------------------------------------------------------

real(8) function source_term_perdir(eta1,eta2) ! in the path
  real(8), intent(in) :: eta1
  real(8), intent(in) :: eta2
  
  source_term_perdir = -2*(0.5*sll_p_pi)**2* sin(0.5*sll_p_pi*eta1)*sin(0.5*sll_p_pi*eta2)
      ! -(16.0*sll_p_pi**2*eta2**4 &
      ! - 16.0*sll_p_pi**2*eta2**2 &
      ! - 12.0*eta2**2 + 2.0)*cos(2*sll_p_pi*eta1)*sin(2*sll_p_pi*eta1)
  
end function source_term_perdir


real(8) function sol_exacte_perdir(eta1,eta2)
  real(8) :: eta1,eta2
  !real(8), dimension(:), intent(in), optional :: params
  sol_exacte_perdir = sin(0.5*sll_p_pi*eta1)*sin(0.5*sll_p_pi*eta2)!eta2 ** 2 * (eta2**2-1)&
      ! * cos(2.0*sll_p_pi*eta1)*sin(2.0*sll_p_pi*eta1)
  
  !print*, 'heho'
end function sol_exacte_perdir


real(8) function sol_exacte_perdir_der1(eta1,eta2)
  real(8) :: eta1,eta2
  !real(8), dimension(:), intent(in), optional :: params
  sol_exacte_perdir_der1 = 2.0*sll_p_pi*cos(2.0*sll_p_pi*eta1)*sin(2.0*sll_p_pi*eta2)
end function sol_exacte_perdir_der1


real(8) function sol_exacte_perdir_der2(eta1,eta2)
  real(8) :: eta1,eta2
  !real(8), dimension(:), intent(in), optional :: params
  sol_exacte_perdir_der2 = 2.0*sll_p_pi*sin(2.0*sll_p_pi*eta1)*cos(2.0*sll_p_pi*eta2)
end function sol_exacte_perdir_der2

!  Solution for a identity change of coordinates 
!   and also dirichlet-periodicconditons
!   the matrix A is equal to identity 
!   the scalar c is equal to zero 
!-------------------------------------------------------------

real(8) function source_term_dirper(eta1,eta2,params) ! in the path
  real(8),intent(in) :: eta1,eta2
  real(8), dimension(:), intent(in), optional :: params
  real(8) :: dummy
  dummy = params(1)

  source_term_dirper = -2*(2*sll_p_pi)**2* sin(2*sll_p_pi*eta1)*cos(2*sll_p_pi*eta2)
     ! -(16.0*sll_p_pi**2*eta1**4 &
     ! - 16.0*sll_p_pi**2*eta1**2 &
     ! - 12.0*eta1**2 + 2.0)*sin(2*sll_p_pi*eta2)*cos(2*sll_p_pi*eta2)
end function source_term_dirper


real(8) function sol_exacte_dirper(eta1,eta2)
  use sll_m_constants
  real(8) :: eta1,eta2
  !real(8), dimension(:), intent(in), optional :: params

  sol_exacte_dirper = sin(2.0*sll_p_pi*eta1)*cos(2.0*sll_p_pi*eta2)
       !eta1 ** 2 * (eta1**2-1)* cos(2*sll_p_pi*eta2)*sin(2*sll_p_pi*eta2)
end function sol_exacte_dirper

real(8) function sol_exacte_dirper_der1(eta1,eta2)
  real(8) :: eta1,eta2
  !real(8), dimension(:), intent(in), optional :: params

  sol_exacte_dirper_der1 = 2*sll_p_pi*cos(2.0*sll_p_pi*eta1)*cos(2.0*sll_p_pi*eta2)
       !eta1 ** 2 * (eta1**2-1)* cos(2*sll_p_pi*eta2)*sin(2*sll_p_pi*eta2)
end function sol_exacte_dirper_der1

real(8) function sol_exacte_dirper_der2(eta1,eta2)
  real(8) :: eta1,eta2
  !real(8), dimension(:), intent(in), optional :: params

  sol_exacte_dirper_der2 = -2.0*sll_p_pi*sin(2.0*sll_p_pi*eta1)*sin(2.0*sll_p_pi*eta2)
       !eta1 ** 2 * (eta1**2-1)* cos(2*sll_p_pi*eta2)*sin(2*sll_p_pi*eta2)
end function sol_exacte_dirper_der2

!----------------------------------------------------------
!  Solution for a r theta change of coordinates 
!   and periodic-dirichlet conditions
!   and also dirivhlet-dirichlet conditons
!   the matrix A is equal to identity 
!   the scalar c is equal to zero 
!-------------------------------------------------------------
real(8) function rho_rtheta(eta1,eta2,params) ! in the path
  real(8),intent(in) :: eta1,eta2
  real(8) :: x, y
  real(8), dimension(:), intent(in), optional :: params
  
  if (present(params)) print*, params

  x = eta2*cos(2*sll_p_pi*eta1)
  y = eta2*sin(2*sll_p_pi*eta1)
  
  rho_rtheta = x*y*(-32.0*x**2 - 32.0*y**2 + 15.0)  
  
end function rho_rtheta


real(8) function sol_exacte_rtheta(eta1,eta2,params) ! in the path
  real(8),intent(in) :: eta1,eta2
  real(8), dimension(:), intent(in), optional :: params
  
  if (present(params)) print*, params
  
  sol_exacte_rtheta = ( eta2**2-1)*(eta2**2-0.5**2)*eta2**2&
       *cos(2*sll_p_pi*eta1)*sin(2*sll_p_pi*eta1)
  
end function sol_exacte_rtheta

!----------------------------------------------------------
!  Solution for a colella change of coordinates 
!   and periodic-periodic conditons
!   the matrix A is equal to identity 
!   the scalar c is equal to zero 
!-------------------------------------------------------------
real(8) function source_term_chgt_perper(eta1,eta2) ! in the path
  real(8):: eta1,eta2
  real(8) :: x, y
  ! real(8), dimension(:), intent(in), optional :: params
  
  x =   eta1 + 0.1_8*sin(2*sll_p_pi*eta1) * sin(2*sll_p_pi*eta2)
  y =   eta2 + 0.1_8*sin(2*sll_p_pi*eta1) * sin(2*sll_p_pi*eta2)
  
  source_term_chgt_perper = -8.0*sll_p_pi**2*cos(2*sll_p_pi*x)*cos(2*sll_p_pi*y) 
  
end function source_term_chgt_perper

real(8) function sol_exacte_chgt_perper(eta1,eta2)
  real(8) :: eta1,eta2
  real(8) :: x,y
  
  x =   eta1 + 0.1*sin(2* sll_p_pi*eta1) * sin(2*sll_p_pi*eta2)
  y =   eta2 + 0.1*sin(2* sll_p_pi*eta1) * sin(2*sll_p_pi*eta2)
  
  sol_exacte_chgt_perper = cos(2*sll_p_pi*x)*cos(2*sll_p_pi*y)
  
  
end function sol_exacte_chgt_perper

real(8) function sol_exacte_chgt_perper_der1(eta1,eta2)
  real(8) :: eta1,eta2
  real(8) :: x,y
  
  x =   eta1 + 0.1*sin(2* sll_p_pi*eta1) * sin(2*sll_p_pi*eta2)
  y =   eta2 + 0.1*sin(2* sll_p_pi*eta1) * sin(2*sll_p_pi*eta2)
  
  sol_exacte_chgt_perper_der1 = -2*sll_p_pi*sin(2*sll_p_pi*x)*cos(2*sll_p_pi*y)&
       * ( 1.0_8 + 0.1*2*sll_p_pi*cos(2* sll_p_pi*eta1) * sin(2*sll_p_pi*eta2) )&
       -2*sll_p_pi*cos(2*sll_p_pi*x)*sin(2*sll_p_pi*y)&
       * ( 0.1*2*sll_p_pi*cos(2* sll_p_pi*eta1) * sin(2*sll_p_pi*eta2) )
  
end function sol_exacte_chgt_perper_der1

real(8) function sol_exacte_chgt_perper_der2(eta1,eta2)
  real(8) :: eta1,eta2
  real(8) :: x,y
  
  x =   eta1 + 0.1*sin(2* sll_p_pi*eta1) * sin(2*sll_p_pi*eta2)
  y =   eta2 + 0.1*sin(2* sll_p_pi*eta1) * sin(2*sll_p_pi*eta2)
  
  !sol_exacte_chgt_perper_der2 = -2*sll_p_pi*cos(2*sll_p_pi*x)*sin(2*sll_p_pi*y)
  
  sol_exacte_chgt_perper_der2 = -2*sll_p_pi*sin(2*sll_p_pi*x)*cos(2*sll_p_pi*y)&
       * ( 0.1*2*sll_p_pi*sin(2* sll_p_pi*eta1) * cos(2*sll_p_pi*eta2) )&
       -2*sll_p_pi*cos(2*sll_p_pi*x)*sin(2*sll_p_pi*y)&
       * ( 1.0_8 + 0.1*2*sll_p_pi*sin(2* sll_p_pi*eta1)*cos(2*sll_p_pi*eta2) )
end function sol_exacte_chgt_perper_der2


!----------------------------------------------------------
!  Solution for a colella change of coordinates 
!   and periodic-dirichlet conditons
!   and dircihlet-diichlet conditions
!   the matrix A is equal to identity 
!   the scalar c is equal to zero 
!-------------------------------------------------------------
real(8) function source_term_chgt_perdir(eta1,eta2) ! in the path
  real(8),intent(in) :: eta1,eta2
  real(8) :: x, y
    
  x =   eta1 + 0.1*sin(2* sll_p_pi*eta1) * sin(2*sll_p_pi*eta2)
  y =   eta2 + 0.1*sin(2* sll_p_pi*eta1) * sin(2*sll_p_pi*eta2)
  
  source_term_chgt_perdir= -2*(2*sll_p_pi)**2 * sin(2*sll_p_pi*y)*cos(2*sll_p_pi*x)
  
end function source_term_chgt_perdir

real(8) function sol_exacte_chgt_perdir(eta1,eta2)
  real(8) :: eta1,eta2
  real(8) :: x,y

  x =   eta1 + 0.1*sin(2* sll_p_pi*eta1) * sin(2*sll_p_pi*eta2)
  y =   eta2 + 0.1*sin(2* sll_p_pi*eta1) * sin(2*sll_p_pi*eta2)
  
  sol_exacte_chgt_perdir = cos(2*sll_p_pi*x)*sin(2*sll_p_pi*y)
  
end function sol_exacte_chgt_perdir

real(8) function sol_exacte_chgt_perdir_der1(eta1,eta2)
  real(8) :: eta1,eta2
  real(8) :: x,y

  x =   eta1 + 0.1*sin(2* sll_p_pi*eta1) * sin(2*sll_p_pi*eta2)
  y =   eta2 + 0.1*sin(2* sll_p_pi*eta1) * sin(2*sll_p_pi*eta2)
  
  sol_exacte_chgt_perdir_der1 = -2*sll_p_pi*sin(2*sll_p_pi*x)*sin(2*sll_p_pi*y)&
       * ( 1.0_8 + 0.1*2*sll_p_pi*cos(2*sll_p_pi*eta1) * sin(2*sll_p_pi*eta2) )&
       + 2*sll_p_pi*cos(2*sll_p_pi*x)*cos(2*sll_p_pi*y)&
       * ( 2*sll_p_pi*0.1*cos(2* sll_p_pi*eta1) * sin(2*sll_p_pi*eta2) ) 
  
end function sol_exacte_chgt_perdir_der1

real(8) function sol_exacte_chgt_perdir_der2(eta1,eta2)
  real(8) :: eta1,eta2
  real(8) :: x,y

  x =   eta1 + 0.1*sin(2* sll_p_pi*eta1) * sin(2*sll_p_pi*eta2)
  y =   eta2 + 0.1*sin(2* sll_p_pi*eta1) * sin(2*sll_p_pi*eta2)
  
  sol_exacte_chgt_perdir_der2 = -2*sll_p_pi*sin(2*sll_p_pi*x)*sin(2*sll_p_pi*y)&
       * ( 0.1*2*sll_p_pi*sin(2*sll_p_pi*eta1) * cos(2*sll_p_pi*eta2) ) &
       + 2*sll_p_pi*cos(2*sll_p_pi*x)*cos(2*sll_p_pi*y)&
       * ( 1.0_8 + 2*sll_p_pi*0.1*sin(2* sll_p_pi*eta1) *cos(2*sll_p_pi*eta2) ) 
  
end function sol_exacte_chgt_perdir_der2

!----------------------------------------------------------
!  Solution for a colella change of coordinates 
!   and dirchlet-periodic conditions
!   the matrix A is equal to identity 
!   the scalar c is equal to zero 
!-------------------------------------------------------------
real(8) function source_term_chgt_dirdir(eta1,eta2) ! in the path
  real(8),intent(in) :: eta1,eta2
  real(8) :: x, y

  ! -------------------------------------------------
  ! In the case without change of coordinates
  ! -------------------------------------------------
  x =   eta1 + 0.1_8*sin(2* sll_p_pi*eta1) * sin(2*sll_p_pi*eta2)
  y =   eta2 + 0.1_8*sin(2* sll_p_pi*eta1) * sin(2*sll_p_pi*eta2)
  
  
  source_term_chgt_dirdir = &
       -2*(2.0*sll_p_pi)**2*sin(2*sll_p_pi*x)*sin(2*sll_p_pi*y)
  
end function source_term_chgt_dirdir

real(8) function sol_exacte_chgt_dirdir(eta1,eta2)
  real(8) :: eta1,eta2
  real(8) :: x,y
  
  x =   eta1 + 0.1*sin(2* sll_p_pi*eta1) * sin(2*sll_p_pi*eta2)
  y =   eta2 + 0.1*sin(2* sll_p_pi*eta1) * sin(2*sll_p_pi*eta2)
  
  sol_exacte_chgt_dirdir = sin(2* sll_p_pi*y)*sin(2* sll_p_pi*x)
  
end function sol_exacte_chgt_dirdir


real(8) function sol_exacte_chgt_dirdir_der1(eta1,eta2)
  real(8) :: eta1,eta2
  real(8) :: x,y
  
  x =   eta1 + 0.1*sin(2* sll_p_pi*eta1) * sin(2*sll_p_pi*eta2)
  y =   eta2 + 0.1*sin(2* sll_p_pi*eta1) * sin(2*sll_p_pi*eta2)
  
  
  sol_exacte_chgt_dirdir_der1 = 2*sll_p_pi*cos(2* sll_p_pi*x)*sin(2* sll_p_pi*y)&
       * ( 1.0_8 + 0.1*2*sll_p_pi*cos(2*sll_p_pi*eta1) * sin(2*sll_p_pi*eta2) )&
       + 2*sll_p_pi*sin(2* sll_p_pi*x)*cos(2* sll_p_pi*y) &
       * ( 2*sll_p_pi*0.1*cos(2* sll_p_pi*eta1) * sin(2*sll_p_pi*eta2) )
end function sol_exacte_chgt_dirdir_der1


real(8) function sol_exacte_chgt_dirdir_der2(eta1,eta2)
  real(8) :: eta1,eta2
  real(8) :: x,y
  
  x =   eta1 + 0.1*sin(2* sll_p_pi*eta1) * sin(2*sll_p_pi*eta2)
  y =   eta2 + 0.1*sin(2* sll_p_pi*eta1) * sin(2*sll_p_pi*eta2)
  
  sol_exacte_chgt_dirdir_der2 =  2*sll_p_pi*cos(2* sll_p_pi*x)*sin(2* sll_p_pi*y)&
       * ( 0.1*2*sll_p_pi*sin(2*sll_p_pi*eta1) * cos(2*sll_p_pi*eta2)  )&
       + 2*sll_p_pi*sin(2* sll_p_pi*x)*cos(2* sll_p_pi*y) &
       * ( 1.0_8 + 2*sll_p_pi*0.1*sin(2* sll_p_pi*eta1) *cos(2*sll_p_pi*eta2) )
  
end function sol_exacte_chgt_dirdir_der2

!----------------------------------------------------------
!  Solution for a colella change of coordinates 
!   and dirchlet-periodic conditions
!   the matrix A is equal to identity 
!   the scalar c is equal to zero 
!-------------------------------------------------------------

real(8) function source_term_chgt_dirper(eta1,eta2) ! in the path
  real(8),intent(in) :: eta1,eta2
  real(8) :: x, y
  !real(8), dimension(:), intent(in), optional :: params
  ! -------------------------------------------------
  ! In the case without change of coordinates
  ! -------------------------------------------------
  x =   eta1 + 0.1_8*sin(2* sll_p_pi*eta1) * sin(2*sll_p_pi*eta2)
  y =   eta2 + 0.1_8*sin(2* sll_p_pi*eta1) * sin(2*sll_p_pi*eta2)

  source_term_chgt_dirper = -2*(2*sll_p_pi)**2*sin(2*sll_p_pi*x)*cos(2*sll_p_pi*y)
  
end function source_term_chgt_dirper

real(8) function sol_exacte_chgt_dirper(eta1,eta2)
  real(8) :: eta1,eta2
  real(8) :: x,y

  x =   eta1 + 0.1*sin(2* sll_p_pi*eta1) * sin(2*sll_p_pi*eta2)
  y =   eta2 + 0.1*sin(2* sll_p_pi*eta1) * sin(2*sll_p_pi*eta2)
  
  sol_exacte_chgt_dirper = sin(2* sll_p_pi*x)*cos(2* sll_p_pi*y)
  
end function sol_exacte_chgt_dirper

real(8) function sol_exacte_chgt_dirper_der1(eta1,eta2)
  real(8) :: eta1,eta2
  real(8) :: x,y
  
  x =   eta1 + 0.1*sin(2* sll_p_pi*eta1) * sin(2*sll_p_pi*eta2)
  y =   eta2 + 0.1*sin(2* sll_p_pi*eta1) * sin(2*sll_p_pi*eta2)
  
  sol_exacte_chgt_dirper_der1 = 2*sll_p_pi*cos(2* sll_p_pi*x)*cos(2* sll_p_pi*y) &
       * ( 1.0_8 + 0.1*2*sll_p_pi*cos(2*sll_p_pi*eta1) * sin(2*sll_p_pi*eta2) )&
       - 2*sll_p_pi*sin(2* sll_p_pi*x)*sin(2* sll_p_pi*y)&
       * ( 2*sll_p_pi*0.1*cos(2* sll_p_pi*eta1) * sin(2*sll_p_pi*eta2) ) 
end function sol_exacte_chgt_dirper_der1

real(8) function sol_exacte_chgt_dirper_der2(eta1,eta2)
  real(8) :: eta1,eta2
  real(8) :: x,y
  
  x =   eta1 + 0.1*sin(2* sll_p_pi*eta1) * sin(2*sll_p_pi*eta2)
  y =   eta2 + 0.1*sin(2* sll_p_pi*eta1) * sin(2*sll_p_pi*eta2)
  
  sol_exacte_chgt_dirper_der2 = 2*sll_p_pi*cos(2* sll_p_pi*x)*cos(2* sll_p_pi*y) &
       * ( 0.1*2*sll_p_pi*sin(2*sll_p_pi*eta1) * cos(2*sll_p_pi*eta2)  )&
       - 2*sll_p_pi*sin(2* sll_p_pi*x)*sin(2* sll_p_pi*y)&
       * (1.0_8 + 2*sll_p_pi*0.1*sin(2* sll_p_pi*eta1) *cos(2*sll_p_pi*eta2) ) 
  
end function sol_exacte_chgt_dirper_der2


!!!!!! test case with F(theta,phi) = (2pi theta , 2pi phi)

real(8) function adimension_chgt_x(eta1,eta2)
  real(8) :: eta1,eta2
  print*, eta1, eta2
  adimension_chgt_x = 2_8*sll_p_pi*eta1 !+ eta2)
end function adimension_chgt_x

real(8) function adimension_chgt_y(eta1,eta2)
  real(8) :: eta1,eta2
  print*, eta1, eta2
  adimension_chgt_y = 2_8*sll_p_pi*eta2
end function adimension_chgt_y

real(8) function jac11_adimension_chgt(eta1,eta2)
  real(8) :: eta1,eta2
  print*, eta1, eta2
  jac11_adimension_chgt = 2_8*sll_p_pi
end function jac11_adimension_chgt

real(8) function jac12_adimension_chgt(eta1,eta2)
  real(8) :: eta1,eta2
  print*, eta1, eta2
  jac12_adimension_chgt = 0.0_8!sll_p_pi
end function jac12_adimension_chgt

real(8) function jac21_adimension_chgt(eta1,eta2)
  real(8) :: eta1,eta2
  print*, eta1, eta2
  jac21_adimension_chgt = 0.0_8!2*sll_p_pi!0.0
end function jac21_adimension_chgt

real(8) function jac22_adimension_chgt(eta1,eta2)
  real(8) :: eta1,eta2
  print*, eta1, eta2
  jac22_adimension_chgt = 2_8*sll_p_pi
end function jac22_adimension_chgt

real(8) function sol_exacte_chgt_adim(eta1,eta2)
  real(8) :: eta1,eta2
  real(8) :: x,y
  
  x =   2_8*sll_p_pi*eta1!+eta2)
  y =   2_8* sll_p_pi*eta2
  
  sol_exacte_chgt_adim = cos(x)*cos(y)
  
end function sol_exacte_chgt_adim

real(8) function source_term_chgt_adim(eta1,eta2)
  real(8) :: eta1,eta2
  real(8) :: x,y
  
  x =   2*sll_p_pi*eta1 !+eta2)
  y =   2* sll_p_pi*eta2
  
  source_term_chgt_adim = -2.0_8*cos(x)*cos(y)
  
end function source_term_chgt_adim

end module m_multipatch_helper_functions
