program unit_test_sll_quadtree
 #include "sll_working_precision.h"
 #include "sll_memory.h"

  use quadtree_mod

  sll_real64 :: xmin, xmax, ymin, ymax, sigma
  integer, parameter :: max_num_el = 10
  integer, parameter :: npart = 5000
  type(quadtree_box), pointer :: root
  sll_real64, dimension(npart,3) :: part
  integer :: i

  ! Define outer bounding box
  xmin = -1.
  xmax = 1
  ymin = -1
  ymax = 1.
  ! Define smoothing distance
  sigma = 0.1
  ! initialise root pointer
  call initialise_quadtree(root, NULL(root), xmin, xmax, ymin, ymax, max_num_el)

  ! generate the part array with a uniform random distribution on [0,1]
  CALL RANDOM_NUMBER(part)
  ! scale to the box
  do i=1,npart
     part(i,1) = xmin + part(i,1)*(xmax-xmin)
     part(i,2) = ymin + part(i,2)*(ymax-ymin)
     print*, i, part(i,1), part(i,2)
  enddo
    print*, '--------'

  ! insert part elements into quadtree starting from roo box
  do i = 1, npart
     call quadtree_insert (root, part, i)
  enddo

  ! traverse quadtree
  call quadtree_traversal(root, part, sigma)

  ! deallocate quadtree
  call deallocate_quadtree(root)

  print*, 'test_quadtree has exited normally'

end program unit_test_sll_quadtree
