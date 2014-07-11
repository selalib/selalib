module sll_quadtree
#include "sll_working_precision.h"
#include "sll_memory.h"
#include "sll_assert.h"

  ! Written by Eric Sonnendrucker 26 March 2014

  !Replace this with pz as particle%dx, pv as particle%pv,
  !use global, only : POS1 => PZ,  POS2 => PV, WEIGHT => PF
  use sll_particle_representations


  implicit none

  sll_int32,parameter:: POS1 = 1, POS2 = 2, WEIGHT = 3


  type sll_quadtree_box
     ! bounding box of quadtree
     sll_real64 :: xmin, xmax, ymin, ymax
     ! maximum number of elements stored in box
     sll_int32 :: max_num_el
     ! actual number of elements stored in box
     sll_int32 :: num_el
     ! array containing the elements
     sll_int32, dimension(:), pointer :: element_list
     ! pointers linking the quadtree box to parent and four children
     type (sll_quadtree_box), pointer :: parent
     type (sll_quadtree_box), pointer :: north_west
     type (sll_quadtree_box), pointer :: north_east
     type (sll_quadtree_box), pointer :: south_west
     type (sll_quadtree_box), pointer :: south_east
  end type  sll_quadtree_box

  contains

!> @brief Initialize sll_quadtree_box
!> @param[in] xmin Minimum x coordinate of box
!> @param[in] xmax Maximum x coordinate of box
!> @param[in] ymin Minimum y coordinate of box
!> @param[in] ymax Maximum y coordinate of box
!> @param[in] max_num_el Maximum number of elements in the box
!> @param[out] this_quadtree Root quadtree box with no elements inside
subroutine initialise_quadtree(this_quadtree, parent, xmin, xmax, ymin, ymax, max_num_el)

  type(sll_quadtree_box), pointer :: this_quadtree
  type(sll_quadtree_box), pointer :: parent
  sll_real64, intent(in) :: xmin, xmax, ymin, ymax
  sll_int32 :: max_num_el
  ! local variables
  sll_int32 :: ierr  ! error flag

  ! allocate the sll_quadtree_box pointer
  allocate(this_quadtree)
  ! Point to parent
  this_quadtree%parent => parent
  ! initialise sll_quadtree_box values
  this_quadtree%xmin =  xmin
  this_quadtree%xmax =  xmax
  this_quadtree%ymin =  ymin
  this_quadtree%ymax =  ymax
  this_quadtree%max_num_el = max_num_el
  ! This box contains no elements for the moment
  this_quadtree%num_el = 0
  ! allocate element_list
  allocate(this_quadtree%element_list(max_num_el),stat=ierr)
  ! Nullify childen pointers
  nullify(this_quadtree%north_west)
  nullify(this_quadtree%north_east)
  nullify(this_quadtree%south_west)
  nullify(this_quadtree%south_east)

end subroutine initialise_quadtree

!> @brief Deallocates an item in the quadtree recursively
!> @param[in]  this_quadtree  Quadtree  which is to be deallocated
recursive subroutine deallocate_quadtree(this_quadtree)
    type (sll_quadtree_box), pointer :: this_quadtree
    ! local variables
    sll_int32 :: ierr

    if (.not.(associated(this_quadtree))) return

    ! Deallocate childen pointers
    call deallocate_quadtree(this_quadtree%north_west)
    call deallocate_quadtree(this_quadtree%north_east)
    call deallocate_quadtree(this_quadtree%south_west)
    call deallocate_quadtree(this_quadtree%south_east)

    ! deallocate element_list array
    deallocate(this_quadtree%element_list,stat=ierr)

    ! deallocate pointer to element
    deallocate(this_quadtree)

end subroutine deallocate_quadtree

!> @brief Inserts a new item in the quadtree.
!> @param[in]  this_quadtree  Quadtree in which item is to be insterted
!> @param[in]  array  Array containing elements to be insterted
!> @param[in]  index    sll_int32 giving access to the array describing element to be inserted
!>
recursive subroutine quadtree_insert (this_quadtree, array, index)

  type (sll_quadtree_box), pointer :: this_quadtree
  sll_real64, dimension(:,:) :: array
  sll_int32, intent(in) :: index
  ! local variables
  sll_real64 :: xmid, ymid, xmin, xmax, ymin, ymax
  sll_int32 :: i, ind

  ! compute box limits and midpoints
  xmin = this_quadtree%xmin
  xmax = this_quadtree%xmax
  ymin = this_quadtree%ymin
  ymax = this_quadtree%ymax
  xmid = (xmin+xmax)/2
  ymid = (ymin+ymax)/2

  ! Ignore item if it is not in the box
  if ((array(index,POS1) < this_quadtree%xmin) .or. (array(index,POS1) > this_quadtree%xmax) &
       .or. (array(index,POS2) < this_quadtree%ymin) .or. (array(index,POS2) > this_quadtree%ymax)) then
     print*, 'Quadtree: particle rejected', index, this_quadtree%xmin,this_quadtree%xmax, &
          this_quadtree%ymin,this_quadtree%ymax
     return
  end if

  ! If there is space in this quadtree box insert item here
  if (this_quadtree%num_el < this_quadtree%max_num_el) then
     this_quadtree%num_el = this_quadtree%num_el + 1
     this_quadtree%element_list(this_quadtree%num_el) = index
!     print*, index, this_quadtree%xmin,this_quadtree%xmax,this_quadtree%ymin,this_quadtree%ymax
  else
     ! try and insert the element in a child
     ! check if the box has already be subdivided)
     if (associated(this_quadtree%north_west)) then
        ! insert present element in appropriate child
        if (array(index,POS1) < xmid) then
           if (array(index,POS2) < ymid) then
              call quadtree_insert(this_quadtree%south_west, array, index)
           else
              call quadtree_insert(this_quadtree%north_west, array, index)
           endif
        else
           if (array(index,POS2) < ymid) then
              call quadtree_insert(this_quadtree%south_east, array, index)
           else
              call quadtree_insert(this_quadtree%north_east, array, index)
           end if
        end if
     else
        ! subdivide the box
        call initialise_quadtree(this_quadtree%north_west,this_quadtree, xmin, xmid, &
             ymid, ymax, this_quadtree%max_num_el)
        call initialise_quadtree(this_quadtree%north_east,this_quadtree, xmid, xmax, &
             ymid, ymax, this_quadtree%max_num_el)
        call initialise_quadtree(this_quadtree%south_west,this_quadtree, xmin, xmid, &
             ymin, ymid, this_quadtree%max_num_el)
        call initialise_quadtree(this_quadtree%south_east,this_quadtree, xmid, xmax, &
             ymin, ymid, this_quadtree%max_num_el)
        ! insert all elements in this box in children boxes
        do i = 1, this_quadtree%num_el
           ind = this_quadtree%element_list(i)
           if (array(ind,POS1) < xmid) then
              if (array(ind,POS2) < ymid) then
                 call quadtree_insert(this_quadtree%south_west, array, ind)
              else
                 call quadtree_insert(this_quadtree%north_west, array, ind)
              endif
           else
              if (array(ind,POS2) < ymid) then
                 call quadtree_insert(this_quadtree%south_east, array, ind)
              else
                 call quadtree_insert(this_quadtree%north_east, array, ind)
              end if
           end if
        end do
        ! insert present element in appropriate child
        if (array(index,POS1) < xmid) then
           if (array(index,POS2) < ymid) then
              call quadtree_insert(this_quadtree%south_west, array, index)
           else
              call quadtree_insert(this_quadtree%north_west, array, index)
           endif
        else
           if (array(index,POS2) < ymid) then
              call quadtree_insert(this_quadtree%south_east, array, index)
           else
              call quadtree_insert(this_quadtree%north_east, array, index)
           end if
        end if
     end if
  end if
end subroutine quadtree_insert

!> @brief Apply function on all elements if at lowest level. Else go to the four children
!> @param[in]  this_quadtree  sll_quadtree_box which is being traversed
!> @param[in]  array Array which will be operated on during quadtree traversal.
!> @param[in]  sigma Decay factor for smoothing of the order of a velocity cell
recursive subroutine quadtree_traversal(this_quadtree, array, sigma)
  type(sll_quadtree_box), pointer :: this_quadtree
  sll_real64, dimension(:,:), intent(inout) :: array
  sll_real64, intent(in) :: sigma
  ! local variables
  sll_int32 :: i, num, num_next
  sll_real64 :: dist, fact, avg
  logical, parameter :: print_quadtree = .false.

  if (.not.(associated(this_quadtree%north_west))) then
     ! current box has no children. So apply function to element list
     ! print quadtree elements if desired
     if (print_quadtree) then
        print*, 'box with bounds:',  this_quadtree%xmin, this_quadtree%xmax, &
             this_quadtree%ymin, this_quadtree%ymax
        do i = 1, this_quadtree%num_el
           num = this_quadtree%element_list(i)
           print*, num, array(num,POS1), array(num,POS2)
        end do
     end if
     ! do close neigbour smoothing
     do i = 1, 2*int(this_quadtree%num_el/2), 2
        num = this_quadtree%element_list(i)
        num_next = this_quadtree%element_list(i+1)
        ! Compute average of two weights
        avg = 0.5*(array(num, WEIGHT)+array(num_next, WEIGHT))
        ! compute weighting factor related to inter-particle distance
        dist = (array(num,POS1)-array(num_next,POS1))**2 + (array(num,POS2)-array(num_next,POS2))**2
        fact = exp(-dist/(sigma**2))
        ! Smooth weights
        array(num, WEIGHT) =  (1-fact)*array(num, WEIGHT) + fact*avg
        array(num_next, WEIGHT) =  (1-fact)*array(num_next, WEIGHT) + fact*avg
     end do
  else
     ! Visit the four children of current quadtree box
     call quadtree_traversal(this_quadtree%north_west, array, sigma)
     call quadtree_traversal(this_quadtree%north_east, array, sigma)
     call quadtree_traversal(this_quadtree%south_west, array, sigma)
     call quadtree_traversal(this_quadtree%south_east, array, sigma)
  end if
end subroutine quadtree_traversal

!==================================================
! Quadtree smoothing algorithm
! Creates a quadtree in phase space
! to pair particles that are close to each other in phase space
! and averages the weights of two particles sharing the same quadtree box
! ES 27 March 2014
   subroutine quadtree_smoothing(part, npart, pos1min,pos1max,pos2min,pos2max,hv,maxppc)
     sll_real64, dimension(:,:) :: part ! array containing the data to insert in the quadtree
     sll_int32 :: npart  ! number of particles (first dimension of part array)
     sll_real64 :: pos1min, pos1max, pos2min, pos2max ! outer box for quadtree
     sll_real64 :: hv    ! parameter representing inter particle distance for smoothing
     sll_int32 :: maxppc  ! maximum number of particles per quadtree cell
     type(sll_quadtree_box), pointer :: root
     ! local variables
     sll_int32 :: i

     ! initialise root quadtree corresponding to box pos1min, pos1max, pos2min, pos2max
     call initialise_quadtree(root, NULL(root), pos1min, pos1max, pos2min, pos2max, maxppc)

     ! insert particles into quadtree
     do i = 1, npart
        call quadtree_insert (root, part, i)
     enddo
     ! traverse quadtree to smooth particle weights in same box
     call quadtree_traversal(root, part, hv)

     ! deallocate quadtree
     call deallocate_quadtree(root)

   end subroutine quadtree_smoothing


!> @brief Quadtree smoothing for one particle species
   subroutine quadtree_smoothing_1d1v( species, pos1min,pos1max,pos2min,pos2max,hv,maxppc)
     type(sll_particle_1d1v_group), intent(inout) :: species  !Group of particles
     sll_real64, intent(in) :: pos1min, pos1max, pos2min, pos2max ! outer box for quadtree
     sll_real64, intent(in) :: hv    ! parameter representing inter particle distance for smoothing
     sll_int32, intent(in) :: maxppc  ! maximum number of particles per quadtree cell
     sll_real64, dimension(size(species%particle),3) :: part !Temporary array
     type(sll_quadtree_box), pointer :: root

     ! local variables
     sll_int32 :: i
     sll_int32 :: npart  ! number of particles (first dimension of part array)


     !Set up working array for the quadtree
     !There maybe a faster solution including pointers
     npart=size(species%particle)
     part(:,POS1)=species%particle%dx
     part(:,POS2)=species%particle%vx
     part(:,WEIGHT)=species%particle%w

     ! initialise root quadtree corresponding to box pos1min, pos1max, pos2min, pos2max
     call initialise_quadtree(root, NULL(root), pos1min, pos1max, pos2min, pos2max, maxppc)

     ! insert particles into quadtree
     do i = 1, npart
        call quadtree_insert (root, part, i)
     enddo

    !   traverse quadtree to smooth particle weights in same box
     call quadtree_traversal(root, part, hv)

     ! deallocate quadtree
     call deallocate_quadtree(root)
   endsubroutine


end module sll_quadtree
