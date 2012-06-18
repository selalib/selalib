module base_class
  implicit none

  type, abstract :: shape
   contains
     procedure(no_argument_real_method), deferred, pass :: area
     procedure(no_argument_subroutine), deferred, pass  :: print
  end type shape

  type :: two_real
     real :: a, b
  end type two_real

  abstract interface 
     function no_argument_real_method(this) result(res)
       import shape
       real :: res
       class(shape), intent(in) :: this
     end function no_argument_real_method
     function one_argument_real_method(this,a) result(res)
       import shape
       real :: res
       class(shape), intent(in) :: this
       real, intent(in) :: a
     end function one_argument_real_method
     subroutine no_argument_subroutine(this)
       import shape
       class(shape), intent(in) :: this
     end subroutine no_argument_subroutine
  end interface
end module base_class

module class_Circle
  use base_class
  implicit none
  real :: pi = 3.1415926535897931d0 ! Class-wide private constant
  type, extends(shape) :: Circle
     real :: radius
   contains
     procedure, pass :: area => circle_area
     procedure, pass :: print => circle_print
  end type Circle
contains
  function circle_area(this) result(area)
    class(Circle), intent(in) :: this
    real :: area
    area = pi * this%radius**2
  end function circle_area

subroutine circle_print(this)
 class(Circle), intent(in) :: this
 real :: area
 area = this%area()  ! Call the type-bound function
 print *, 'Circle: r = ', this%radius, ' area = ', area
end subroutine circle_print
end module class_Circle

module class_square
  use base_class
  implicit none
  type, extends(shape) :: square
     real :: side
   contains
     procedure, pass :: area => square_area
     procedure, pass :: print => square_print
     procedure, pass :: set_side
     procedure, pass :: square_sum
     generic :: operator(+) => square_sum
  end type square
contains
  function square_sum(this,b) result(ssum)
    type(square), pointer :: ssum
    class(square),  intent(in) :: this,b
    allocate(ssum)
    ssum%side = this%side + b%side
  end function square_sum
  function square_area(this) result(area)
    class(square), intent(in) :: this
    real :: area
    area = this%side**2
  end function square_area

  subroutine square_print(this)
    class(square), intent(in) :: this
    real :: area
    area = this%area()  ! Call the type-bound function
    print *, 'Square: a = ', this%side, ' area = ', area
  end subroutine square_print

  subroutine set_side(this, side)
    class(square), intent(inout) :: this
    real, intent(in) :: side
    
    this%side = side
  end subroutine set_side

  elemental function sum_squares(a,b) result(c)
    real, intent(in) :: a, b
    real :: c
    c = a*a + b*b
  end function sum_squares
end module class_square

program circle_test
  use base_class
  use class_Circle
  use class_square
  implicit none

  type(Circle), target :: c
  type(square), target :: s
  class(shape), pointer :: sh

  real, dimension(3) :: a =[real(8):: 1, 2, 1./3]
  integer :: i,j
  real, dimension(6), target  :: X = [ (I, I = 2, 17, 3) ]
  real, dimension(:,:), pointer :: b
  real, dimension(3,2) :: tab1=reshape([1,2,3,1,2,3],[3,2]), tab2=reshape([4,5,6,4,5,6],[3,2])
  type(two_real) :: t
  
     
  print*, 'type constructor'
  t = two_real(2.3,1)
  print*, t%a, t%b
  
  !X = [1, 2, 3]
  print*, 'x ',x

  !b(1:3,1:2) => x(:)
  !print*, b
  
  c = circle(radius=1.5)
  c%radius = 1.5
  call c%print         
 
  s%side = 2.0
  call s%print

  call s%set_side(3.0)
  call s%print
  
  print*, 'shape is a circle'
  sh => c
  call sh%print
  print*, 'now shape is a square'
  sh => s
  call sh%print

  print*, 'test elemental '
  print*, tab1
  print*, tab2
  print*, sum_squares(tab1,tab2)

end program circle_test
