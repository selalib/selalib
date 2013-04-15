!------------------------------------------------------------------------------
! SELALIB
!------------------------------------------------------------------------------
! MODULE: sll_nu_cart_mesh
!
! DESCRIPTION:
!> @file mesh.F90
!> @namespace sll_nu_cart_mesh
!> @author Madaule Eric
!> @brief module for non uniform cartesian mesh
!> @details In this module we define a non uniform cartesian mesh type.
!>          Here also are the constructors and destructor of this type
!>          Note that this type is not made to be a pointer and is build with subroutines
!------------------------------------------------------------------------------
module sll_nu_cart_mesh
#include "sll_working_precision.h"
#include "sll_memory.h"

  implicit none

  type :: non_unif_cart_mesh
     !---------------------------------------------------------------------------
     !> @brief non uniform 2D mesh
     !> @details non uniform 2D mesh
     !>          contains number of cells in each direction, each cell size in each direction
     !>                   and position of each nodes
     sll_int32 :: n_etat1,n_etat2
     sll_real64,dimension(:),pointer :: d_etat1,d_etat2
     sll_real64,dimension(:),pointer :: etat1,etat2
     ! position of node j in direction i, can be usefull for computation, but I'm nor sure
     ! see if we keep this in some more developed version
     logical :: unif1, unif2 ! .true. if the meshing is uniform, just for optimization, .fale. else
  end type non_unif_cart_mesh

contains

  subroutine init_nu_cart_mesh(n_etat1,n_etat2,d_etat1,d_etat2,mesh)
    !---------------------------------------------------------------------------
    !> @brief creation of type non_unif_cart_mesh
    !> @details creation of type non_unif_cart_mesh
    !>          the mesh can be uniform in one direction etat1 or etat2 (both is possible)
    !>          we allocate nodes and d_etati
    !>          if the size of cells in one direction is not constant, just don't enter any value
    !>          the vector will be allocated and you can fill it by yourself
    !>          if the size of the cell is constant in direction i, then size(d_etati) = 1
    !> @param[IN] n_etat1 number of cell in direction etat1
    !> @param[IN] n_etat2 number of cell in direction etat2
    !> @param[IN,OPTIONAL] d_etat1 size of step in direction etat1 if the size is constant
    !> @param[IN,OPTIONAL] d_etat2 size of step in direction etat2 if the size is constant
    !> @param[OUT] mesh the mesh to build
    
    sll_int32,intent(in) :: n_etat1,n_etat2
    sll_real64,intent(in),optional :: d_etat1,d_etat2
    type(non_unif_cart_mesh), intent(out) :: mesh

    integer :: err

    mesh%n_etat1=n_etat1
    mesh%n_etat2=n_etat2
    SLL_ALLOCATE(mesh%etat1(n_etat1+1),err)
    SLL_ALLOCATE(mesh%etat2(n_etat2+1),err)

    if (present(d_etat1)) then
       SLL_ALLOCATE(mesh%d_etat1(1),err)
       mesh%d_etat1=d_etat1
       mesh%unif1=.true.
       if (.not. present(d_etat2)) then
          SLL_ALLOCATE(mesh%d_etat2(n_etat2),err)
          mesh%unif2=.false.
       end if
    end if
    if (present(d_etat2)) then
       SLL_ALLOCATE(mesh%d_etat2(1),err)
       mesh%d_etat2=d_etat2
       mesh%unif2=.true.
       if (.not. present(d_etat1)) then
          SLL_ALLOCATE(mesh%d_etat1(n_etat1),err)
          mesh%unif1=.false.
       end if
    end if
    if ((.not. present(d_etat1)) .and. (.not. present(d_etat2))) then
       print*,"The mesh is non-uniform in both direction"
       SLL_ALLOCATE(mesh%d_etat1(n_etat1),err)
       SLL_ALLOCATE(mesh%d_etat2(n_etat2),err)
       mesh%unif1=.false.
       mesh%unif2=.false.
    end if

  end subroutine init_nu_cart_mesh

  subroutine fill_node_nuc_mesh(etat1_min,etat2_min,mesh)
    !---------------------------------------------------------------------------
    !> @brief filling the mesh%node with nodes coordinates
    !> @details filling the mesh%node with nodes coordinates
    !>          the non_unif_cart_mesh object must already be build and the d_etat1 and d_etat2 be filled
    !> @param[IN] etat1_min minimum in direction etat1
    !> @param[IN] etat2_min minimum in direction etat2
    !> @param[INOUT] mesh the mesh to build

    sll_real64,intent(in) :: etat1_min,etat2_min
    type(non_unif_cart_mesh), intent(inout) :: mesh

    integer :: i,j

    mesh%etat1(1)=etat1_min
    mesh%etat2(1)=etat2_min

    if (mesh%unif1) then
       do i=2,mesh%n_etat1+1
          mesh%etat1(i)=mesh%etat1(i-1)+mesh%d_etat1(1)
       end do
    else
       do i=2,mesh%n_etat1+1
          mesh%etat1(i)=mesh%etat1(i-1)+mesh%d_etat1(i-1)
       end do
    end if

    if (mesh%unif2) then
       do i=2,mesh%n_etat2+1
          mesh%etat2(i)=mesh%etat2(i-1)+mesh%d_etat2(1)
       end do
    else
       do i=2,mesh%n_etat2+1
          mesh%etat2(i)=mesh%etat2(i-1)+mesh%d_etat2(i-1)
       end do
    end if

    ! we test if some node have zero for coordinate in one direction
    ! this can be requiered for some cases (as Vlasov-Poisson with DG)
    j=0
    do i=1,mesh%n_etat1+1
       if (mesh%etat1(i)==epsilon( max(abs(mesh%etat1(1)),abs(mesh%etat1(mesh%n_etat1+1))) )) then
          j=1
       end if
    end do
    if (j==0) then
       print*,'0 does not belong to nodes in direction etat1'
    end if
    j=0
    do i=1,mesh%n_etat2+1
       if (mesh%etat2(i)==epsilon( max(abs(mesh%etat2(1)),abs(mesh%etat2(mesh%n_etat2+1))) )) then
          j=1
       end if
    end do
    if (j==0) then
       print*,'0 does not belong to nodes in direction etat2'
    end if

  end subroutine fill_node_nuc_mesh

  subroutine delete_nu_cart_mesh(mesh)
    !---------------------------------------------------------------------------
    !> @brief delete an object of type non_unif_cart_mesh
    !> @details delete all array in an object of type non_unif_cart_mesh
    !> @param[INOUT] mesh the mesh to delete

    type(non_unif_cart_mesh),intent(inout) :: mesh

    DEALLOCATE(mesh%etat1)
    DEALLOCATE(mesh%etat2)
    DEALLOCATE(mesh%d_etat1)
    DEALLOCATE(mesh%d_etat2)

  end subroutine delete_nu_cart_mesh

end module sll_nu_cart_mesh
