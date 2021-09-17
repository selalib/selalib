!! eta1, eta2 are the logical coordinates which are transformed
!! to the physical coordinates x1, x2
!! The Jacobian matrix is the derivate matrix of the transformation function
!! the determinant of the Jacobian is needed in an integral if the transformation is used
!! the logical paramters are in the interval [0,1]
module sll_m_mapping_2d
  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#include "sll_assert.h"
#include "sll_memory.h"
#include "sll_working_precision.h"
  implicit none

  public :: &
       sll_t_mapping_2d
     
  
  private

  abstract interface
     function sll_i_eval_function( eta1, eta2, params ) result(res)
       use sll_m_working_precision
       sll_real64, intent(in) :: eta1
       sll_real64, intent(in) :: eta2
       sll_real, dimension(:), intent(in) :: params
       sll_real64             :: res
     end function sll_i_eval_function
  end interface
  
  type matrix_element
     procedure(sll_i_eval_function), pointer, nopass :: f
  end type matrix_element
  
  type :: sll_t_mapping_2d
     
     type(matrix_element), dimension(:,:), pointer :: j_matrix
     procedure(sll_i_eval_function), pointer, nopass :: x1_func
     procedure(sll_i_eval_function), pointer, nopass :: x2_func
     sll_real64, dimension(:), pointer :: params
     
   contains
     procedure :: get_x1

     procedure :: get_x2
     
     procedure :: get_x

     procedure :: jacobian

     procedure :: jacobian_matrix

     procedure :: jacobian_matrix_inverse

     procedure :: jacobian_matrix_inverse_transposed

     procedure :: metric

     procedure :: metric_inverse

     procedure :: init

     procedure :: free

  end type sll_t_mapping_2d
  
contains
   function get_x1( self, eta1, eta2 ) result(x1)
    class(sll_t_mapping_2d), intent(in) :: self
    sll_real64, intent(in)   :: eta1
    sll_real64, intent(in)   :: eta2
    sll_real64               :: x1

    x1=self%x1_func(eta1,eta2,self%params)
 
  end function get_x1

  function get_x2( self, eta1, eta2 ) result(x2)
    class(sll_t_mapping_2d), intent(in) :: self
    sll_real64, intent(in)   :: eta1
    sll_real64, intent(in)   :: eta2
    sll_real64               :: x2
    x2=self%x2_func(eta1,eta2,self%params)

  end function get_x2 
  
  function get_x( self, eta1, eta2 ) result(x)
    class(sll_t_mapping_2d), intent(in) :: self
    sll_real64, intent(in)   :: eta1
    sll_real64, intent(in)   :: eta2
    sll_real64               :: x(2)
    x(1)=self%x1_func(eta1,eta2,self%params)
    x(2)=self%x2_func(eta1,eta2,self%params)
  end function get_x
  

  function jacobian(self, eta1, eta2)result(x)
    class(sll_t_mapping_2d), intent(in) :: self
    sll_real64               :: x
    sll_real64, intent(in)   :: eta1
    sll_real64, intent(in)   :: eta2
    
    x=self%j_matrix(1,1)%f(eta1,eta2,self%params)*&
         self%j_matrix(2,2)%f(eta1,eta2,self%params)-&
         self%j_matrix(1,2)%f(eta1,eta2,self%params)*&
         self%j_matrix(2,1)%f(eta1,eta2,self%params)
  end function jacobian



  function jacobian_matrix(self, eta1, eta2)result(y)
    class(sll_t_mapping_2d), intent(in) :: self
    sll_real64, intent(in)   :: eta1
    sll_real64, intent(in)   :: eta2
    sll_real64               :: y(2,2)

    y(1,1)=self%j_matrix(1,1)%f(eta1,eta2,self%params)
    y(1,2)=self%j_matrix(1,2)%f(eta1,eta2,self%params)
    y(2,1)=self%j_matrix(2,1)%f(eta1,eta2,self%params)
    y(2,2)=self%j_matrix(2,2)%f(eta1,eta2,self%params)
  end function jacobian_matrix

  function jacobian_matrix_inverse(self, eta1, eta2)result(y)
    class(sll_t_mapping_2d), intent(in) :: self
    sll_real64, intent(in)   :: eta1
    sll_real64, intent(in)   :: eta2
    sll_real64               :: y(2,2)

    y(1,1)=self%j_matrix(2,2)%f(eta1,eta2,self%params)
    y(1,2)=-self%j_matrix(1,2)%f(eta1,eta2,self%params)
    y(2,1)=-self%j_matrix(2,1)%f(eta1,eta2,self%params)
    y(2,2)=self%j_matrix(1,1)%f(eta1,eta2,self%params)
    y=y/self%jacobian(eta1, eta2)
  end function jacobian_matrix_inverse

  function jacobian_matrix_inverse_transposed(self, eta1, eta2)result(y)
    class(sll_t_mapping_2d), intent(in) :: self
    sll_real64, intent(in)   :: eta1
    sll_real64, intent(in)   :: eta2
    sll_real64               :: y(2,2)

    y(1,1)=self%j_matrix(2,2)%f(eta1,eta2,self%params)
    y(1,2)=-self%j_matrix(2,1)%f(eta1,eta2,self%params)
    y(2,1)=-self%j_matrix(1,2)%f(eta1,eta2,self%params)
    y(2,2)=self%j_matrix(1,1)%f(eta1,eta2,self%params)
    y=y/self%jacobian(eta1, eta2)
  end function jacobian_matrix_inverse_transposed

  function metric(self, eta1, eta2)result(g)
    class(sll_t_mapping_2d), intent(in) :: self
    sll_real64, intent(in)   :: eta1
    sll_real64, intent(in)   :: eta2
    sll_real64               :: y(2,2)
    sll_real64               :: g(2,2)
    y=self%jacobian_matrix(eta1, eta2)
    g(1,1)=y(1,1)**2+y(2,1)**2
    g(1,2)=y(1,1)*y(1,2)+y(2,1)*y(2,2)
    g(2,1)=g(1,2)
    g(2,2)=y(1,2)**2+y(2,2)**2
  end function metric


  function metric_inverse(self, eta1, eta2)result(y)
    class(sll_t_mapping_2d), intent(in) :: self
    sll_real64, intent(in)   :: eta1
    sll_real64, intent(in)   :: eta2
    sll_real64               :: y(2,2)
    sll_real64               :: g(2,2)
    y=self%jacobian_matrix_inverse_transposed(eta1, eta2)
    g(1,1)=y(1,1)**2+y(2,1)**2
    g(1,2)=y(1,1)*y(1,2)+y(2,1)*y(2,2)
    g(2,1)=g(1,2)
    g(2,2)=y(1,2)**2+y(2,2)**2
  end function metric_inverse
  
  
  subroutine init(self,j11,j12,j21,j22,x1_func,x2_func,params)
    class(sll_t_mapping_2d), intent(out) :: self
    procedure(sll_i_eval_function) :: j11
    procedure(sll_i_eval_function) :: j12
    procedure(sll_i_eval_function) :: j21
    procedure(sll_i_eval_function) :: j22
    procedure(sll_i_eval_function) :: x1_func
    procedure(sll_i_eval_function) :: x2_func
    sll_real64, dimension(:), intent(in) :: params 
    !local variables
    sll_int32 :: ierr
    
    SLL_ALLOCATE(self%j_matrix(2,2), ierr)
    SLL_ALLOCATE(self%params(size(params)),ierr)
    self%x1_func => x1_func
    self%x2_func => x2_func
    self%j_matrix(1,1)%f=>j11
    self%j_matrix(1,2)%f=>j12
    self%j_matrix(2,1)%f=>j21
    self%j_matrix(2,2)%f=>j22
    self%params=params
  end subroutine init

  subroutine free(self)
    class(sll_t_mapping_2d), intent(inout) :: self
    DEALLOCATE(self%j_matrix)
    DEALLOCATE(self%params)
  end subroutine free
  

end module sll_m_mapping_2d
