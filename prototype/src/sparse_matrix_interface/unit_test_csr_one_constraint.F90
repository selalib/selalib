!gfortran -O3 sll_csr_one_constraint.F90

! we add one constraint for a matrix in CSR format
! form a matrix A
! we get a matrix
!   A   b
!   b^T 0
! where b in the vector of constraint (constraint_vec)
! so that all the component of b are nonzero
! b^T stands for the transposition of b
! we suppose that A is square
!
! example of CSR format see 
! http://netlib.org/linalg/html_templates/node91.html#SECTION00931100000000000000
!A= 10  0  0  0 -2  0
!    3  9  0  0  0  3
!    0  7  8  7  0  0
!    3  0  8  7  5  0
!    0  8  0  9  9 13 
!
!val     = 10 -2  3  9  3  7 ...
!col_ind =  1  5  1  2  6  2 ...
!
!row_ptr =  1 3 6 9 13 17 20
!
! here ia stands for row_prt
! ja for col_ind
! and a for val


program csr_one_constraint
  use csr_one_constraint_module
  implicit none
  !local variables
  integer :: num_rows
  integer :: num_rows2
  integer :: num_nz
  integer :: num_nz2
  real(8), dimension(:), allocatable :: a
  real(8), dimension(:), allocatable :: a2
  real(8), dimension(:), allocatable :: constraint_vec
  integer, dimension(:), allocatable :: ia
  integer, dimension(:), allocatable :: ia2
  integer, dimension(:), allocatable :: ja
  integer, dimension(:), allocatable :: ja2
  integer, parameter :: TEST_ADNANE = 0 
  integer, parameter :: TEST_NETLIB = 1 
  integer :: test_matrix
  integer :: i
  integer :: k
  
  test_matrix = TEST_ADNANE
  !test_matrix = TEST_NETLIB
  
  select case  (test_matrix)
    case (TEST_ADNANE)
      num_rows = 4
      num_nz = 9
      allocate(ia(num_rows+1))
      allocate(ja(num_nz))
      allocate(a(num_nz))
      allocate(constraint_vec(num_rows))
      ia(1:num_rows+1) = (/1,2,5,8,10/) !change with previous (/1,2,5,8,9/)
      ja(1:num_nz) = (/1,1,2,4,1,3,4,3,4/)
      a(1:num_nz) = (/1._8,3._8,4._8,5._8,6._8,7._8,8._8,10._8,11._8/)
      constraint_vec(1:num_rows) = 1._8      
    case (TEST_NETLIB)
      num_rows = 6
      num_nz = 19
      allocate(ia(num_rows+1))
      allocate(ja(num_nz))
      allocate(a(num_nz))
      allocate(constraint_vec(num_rows))
      ia(1:num_rows+1) = (/1,3,6,9,13,17,20/)
      ja(1:num_nz) = (/&
                        1,5, &
                        1,2,6, &
                        2,3,4, &
                        1,3,4,5, &
                        2,4,5,6, &
                        2,5,6 &
                      /)
      a(1:num_nz) = (/&
                     10._8,-2._8, &
                     3._8,9._8,3._8, &
                     7._8,8._8,7._8, &
                     3._8,8._8,7._8,5._8, &
                     8._8,9._8,9._8,13._8, &
                     4._8,2._8,-1._8 &
                     /)
      constraint_vec(1:num_rows) = 1._8            
    case default
  end select
  

  num_rows2 = num_rows+1
  num_nz2 = num_nz+2*num_rows

  allocate(ia2(num_rows2+1))
  allocate(ja2(num_nz2))
  allocate(a2(num_nz2))
  
  call csr_add_one_constraint(ia,ja,a,num_rows,num_nz,constraint_vec,ia2,ja2,a2)
  
  print*,'#input'
  
  do i=1,num_rows
    do k=ia(i),ia(i+1)-1
      print *,'A(',i,',',ja(k),')=',a(k)
    enddo
  enddo

  print*,'#output'
  
  do i=1,num_rows2
    do k=ia2(i),ia2(i+1)-1
      print *,'A(',i,',',ja2(k),')=',a2(k)
    enddo
  enddo

  print *,'ia=',ia
  print *,'ia2=',ia2

  print *,'ja=',ja
  print *,'ja2=',ja2

  print *,'a=',floor(a)
  print *,'a2=',floor(a2)



  
  

end program csr_one_constraint







