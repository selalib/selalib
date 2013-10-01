!     
! File:   geometries_def.F90
! Author: ratnani
!
! Created on December 5, 2011, 2:13 PM
!

module sparse_matrix_def
    use used_precision
    implicit none

    private

    !******************************************************************************************
    !		SparseMatrix
    !******************************************************************************************
    type, public :: CSR_MATRIX
        integer :: oi_nR !NUMBER OF ROWS
        integer :: oi_nC !NUMBER OF COLUMNS
        integer :: oi_nel !NUMBER OF NON ZERO ELTS
        integer, dimension(:), pointer :: opi_ia
        integer, dimension(:), pointer :: opi_ja
        real(wp), dimension(:), pointer :: opr_a
        !................
        logical :: ol_use_mm_format
        integer, dimension(:), pointer :: opi_i
        !................
    end type CSR_MATRIX
    
    type, public :: CSR_MATRICES
        integer :: oi_n
        type(CSR_MATRIX), dimension(:), pointer :: opo_M
    end type CSR_MATRICES
    !******************************************************************************************

end module sparse_matrix_def

