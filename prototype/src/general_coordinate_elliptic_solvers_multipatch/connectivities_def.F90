!     
! File:   connectivities_def.F90
! Author: root
!
! Created on December 20, 2011, 12:10 PM
!

module connectivities_def
    implicit none

    private

    type, public :: CONNECTIVITY
        integer :: oi_npatch    
        integer :: oi_nelt
        integer :: oi_N
        !> this is the local connectivity array
        !> opi_IEN (id_p, b, elt)
        !> id_p : the id of the current patch
        !> b : the local num of the basis
        !> elt : the num of the element in the current patch
        integer, dimension(:,:,:), pointer :: opi_IEN
        integer, dimension(:), pointer :: opi_ID
        !> LM = ID ( IEN )
        integer, dimension(:,:,:), pointer :: opi_LM
        !> nen array : for each patch and element we give the number of element basis
        integer, dimension(:,:), pointer :: opi_nen
        !> for each patch we give the correspondance
        !> between the real elt id and the loop-id
        !> this is to handle boundary grids, but can also be used
        !> in the case of local refinement, or to parallelize the assembling process
        !> when a local mesh has a lot of GL points (compared to others elts)
        integer, dimension(:,:), pointer :: opi_real_elts
    end type CONNECTIVITY

    type, public :: CONNECTIVITIES
        integer :: oi_n
        type(CONNECTIVITY), dimension(:), pointer :: opo_con
    end type CONNECTIVITIES

end module connectivities_def