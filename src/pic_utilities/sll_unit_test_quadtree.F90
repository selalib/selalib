program unit_test_sll_quadtree
#include "sll_working_precision.h"
#include "sll_memory.h"
#include "sll_assert.h"

  use sll_quadtree


 implicit none

  sll_real64 :: xmin, xmax, ymin, ymax, sigma
  sll_int32, parameter :: max_num_el = 10
  sll_int32, parameter :: npart = 5000
  type(sll_quadtree_box), pointer :: root
  sll_real64, dimension(npart,3) :: part
  type(sll_particle_1d1v_group) :: particlespecies
  sll_int32 :: i,ierr

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

  !-----------------------------------------------------------------------------
  !Generate random weights for second test
  call random_number(part(:,3))

  allocate(particlespecies%particle(1:npart))
  particlespecies%particle%dx=part(:,1)
  particlespecies%particle%vx=part(:,2)
  particlespecies%particle%w=part(:,3)

  call quadtree_smoothing_1d1v( particlespecies, xmin, xmax, ymin, ymax,sigma,max_num_el)



  print*, 'test_quadtree has exited normally'

end program unit_test_sll_quadtree
