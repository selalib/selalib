!> @ingroup maxwell_solvers
!> @brief
!> Helper for spline finite elements utilites
!> 
!> @author
!> Benedikt Perse

module sll_m_spline_fem_utilities_3d_helper
  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#include "sll_assert.h"
#include "sll_errors.h"
#include "sll_memory.h"
#include "sll_working_precision.h"

  use sll_m_matrix_csr, only: &
       sll_t_matrix_csr

  implicit none

  public :: &
       sll_s_spline_fem_sparsity_mass3d, &
       sll_s_spline_fem_sparsity_mass3d_clamped, &
       sll_s_spline_fem_sparsity_mixedmass3d, &
       sll_s_spline_fem_sparsity_mixedmass3d_clamped, &
       assemble_mass3d, &
       assemble_mass3d_clamped, &
       assemble_mass3d_clamped_boundary, &
       assemble_mixedmass3d, &
       assemble_mixedmass3d_clamped_boundary, &
       assemble_mixedmass3d_clamped

  private
  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

contains

  !---------------------------------------------------------------------------!
  !> Helper function to create sparsity pattern of the 3d mass matrix
  subroutine sll_s_spline_fem_sparsity_mass3d( deg, n_cells, spmat )
    sll_int32,              intent( in    ) :: deg(3) !< spline degree for every direction
    sll_int32,              intent( in    ) :: n_cells(3) !< number of cells (and grid points)
    type(sll_t_matrix_csr), intent(   out ) :: spmat     !< sparse matrix

    !local variables
    sll_int32 :: n_nnz
    sll_int32 :: n_nnz_per_row
    sll_int32 :: ind, shift,row,i,j,k,k2

    ! Compute number of non-zero elements
    n_nnz_per_row = (2*deg(1)+1)*(2*deg(2)+1)*(2*deg(3)+1)
    n_nnz = n_nnz_per_row*n_cells(1)*n_cells(2)*n_cells(3)

    ! Create the csr matrix
    call spmat%create( n_rows=n_cells(1)*n_cells(2)*n_cells(3), n_cols=n_cells(1)*n_cells(2)*n_cells(3), n_nnz=n_nnz )

    ! Set the row and column vector for the elements (row vector compressed)
    ! There are n_nnz_per_row elements in each of the rows
    spmat%arr_ia(1) = 1
    do row = 2, n_cells(1)*n_cells(2)*n_cells(3)+1
       spmat%arr_ia(row) = spmat%arr_ia(row-1) + n_nnz_per_row
    end do

    SLL_ASSERT( spmat%arr_ia(n_cells(1)*n_cells(2)*n_cells(3)+1) == n_nnz+1 )


    ind = 1
    !initialise mass matrix 
    !loop over number of rows
    do k = 1, n_cells(3)
       if(k <= deg(3))then
          do j = 1, n_cells(2)
             do i = 1, n_cells(1)
                do k2 = 1, k+deg(3)
                   shift = (k2-1)*n_cells(1)*n_cells(2)
                   call loop2( i, j, shift, deg, n_cells, ind, spmat )
                end do
                do k2 = n_cells(3)+k-deg(3), n_cells(3)
                   shift = (k2-1)*n_cells(1)*n_cells(2)
                   call loop2( i, j, shift, deg, n_cells, ind, spmat )
                end do
                SLL_ASSERT( ind == (i+(j-1)*n_cells(1)+(k-1)*n_cells(1)*n_cells(2))*n_nnz_per_row+1 )
             end do
          end do
       elseif( k > deg(3) .and. k<= n_cells(3)-deg(3))then
          do j = 1, n_cells(2)
             do i = 1, n_cells(1)
                do k2 = k-deg(3), k+deg(3)
                   shift = (k2-1)*n_cells(1)*n_cells(2)
                   call loop2( i, j, shift, deg, n_cells, ind, spmat )
                end do
                SLL_ASSERT( ind == (i+(j-1)*n_cells(1)+(k-1)*n_cells(1)*n_cells(2))*n_nnz_per_row+1 )
             end do
          end do
       else
          do j =1, n_cells(2)
             do i =1, n_cells(1)
                do k2 = 1, k+deg(3)-n_cells(3)
                   shift = (k2-1)*n_cells(1)*n_cells(2)
                   call loop2( i, j, shift, deg, n_cells, ind, spmat )
                end do
                do k2 = k-deg(3), n_cells(3)
                   shift = (k2-1)*n_cells(1)*n_cells(2)
                   call loop2( i, j, shift, deg, n_cells, ind, spmat )
                end do
                SLL_ASSERT( ind == (i+(j-1)*n_cells(1)+(k-1)*n_cells(1)*n_cells(2))*n_nnz_per_row+1 )
             end do
          end do
       end if
    end do
    spmat%arr_a = 0.0_f64

    SLL_ASSERT( maxval(spmat%arr_ja) == spmat%n_cols )

    SLL_ASSERT( ind == spmat%n_nnz+1 )

  end subroutine sll_s_spline_fem_sparsity_mass3d

  !initialise sparse matrix in second direction
  subroutine loop2( i, j, shift, deg, n_cells, ind, spmat )
    sll_int32,              intent( in    ) :: i !< row in first direction
    sll_int32,              intent( in    ) :: j !< row in second direction
    sll_int32,              intent( in    ) :: shift !< shift in the third direction
    sll_int32,              intent( in    ) :: deg(3) !< spline degree for every direction
    sll_int32,              intent( in    ) :: n_cells(3) !< number of cells (and grid points)
    sll_int32,              intent( inout ) :: ind   !< index 
    type(sll_t_matrix_csr), intent( inout ) :: spmat !< sparse matrix
    !local variables
    sll_int32 :: j2, shift1

    !in the first deg(2) rows the entries are splitted in the second direction
    if(j<=deg(2))then
       do j2 = 1, j+deg(2)
          shift1 = (j2-1)*n_cells(1)+shift
          call loop1( i, shift1, deg(1), n_cells(1), ind, spmat )
       end do
       do j2 = n_cells(2)+j-deg(2), n_cells(2)
          shift1 = (j2-1)*n_cells(1)+shift
          call loop1( i, shift1, deg(1), n_cells(1), ind, spmat )
       end do

       !in rows: deg(2)+1-n_cells(2)-deg(2) we have consecutive entries in the second direction  
    elseif( j > deg(2) .and. j<= n_cells(2)-deg(2))then
       do j2 = j-deg(2), j+deg(2)
          shift1 = (j2-1)*n_cells(1)+shift
          call loop1( i, shift1, deg(1), n_cells(1), ind, spmat )
       end do

       !in the last deg(2) rows the entries are splitted in the second direction
    else
       do j2 = 1, j+deg(2)-n_cells(2)
          shift1 = (j2-1)*n_cells(1)+shift
          call loop1( i, shift1, deg(1), n_cells(1), ind, spmat )
       end do
       do j2 = j-deg(2), n_cells(2)
          shift1 = (j2-1)*n_cells(1)+shift
          call loop1(i,shift1, deg(1),n_cells(1),ind,spmat)
       end do
    end if

  end subroutine loop2

  !initialise sparse matrix in first direction
  subroutine loop1( i, shift, deg, n_cells, ind, spmat )
    sll_int32,              intent( in    ) :: i !< row in first direction
    sll_int32,              intent( in    ) :: shift !< shift in the second direction
    sll_int32,              intent( in    ) :: deg !< spline degree in first direction
    sll_int32,              intent( in    ) :: n_cells !< number of cells in the first direction 
    sll_int32,              intent( inout ) :: ind   !< index    
    type(sll_t_matrix_csr), intent( inout ) :: spmat !< sparse matrix
    !local variables
    sll_int32 :: i2

    !in the first deg(1) rows the entries are splitted in the first direction
    if(i<= deg)then
       do i2 = 1, i+deg
          spmat%arr_ja(ind) = i2+shift
          ind = ind+1
       end do
       do i2 = n_cells+i-deg, n_cells 
          spmat%arr_ja(ind) = i2+shift
          ind = ind+1
       end do

       !in rows: deg(1)+1-n_cells(1)-deg(1) we have consecutive entries in the first direction  
    elseif( i > deg .and. i<= n_cells-deg) then
       do i2 = i-deg, i+deg
          spmat%arr_ja(ind) = i2+shift
          ind = ind+1
       end do

       !in the last deg(1) rows the entries are splitted in the first direction
    else
       do i2 = 1, i+deg-n_cells
          spmat%arr_ja(ind) = i2+shift
          ind = ind+1
       end do
       do i2 = i-deg, n_cells
          spmat%arr_ja(ind) = i2+shift
          ind = ind+1
       end do
    end if
  end subroutine loop1

  !---------------------------------------------------------------------------!
  !> Helper function to create sparsity pattern of the 3d clamped mass matrix
  subroutine sll_s_spline_fem_sparsity_mass3d_clamped( deg, n_cells, spmat )
    sll_int32,              intent( in    ) :: deg(3) !< spline deg for every direction
    sll_int32,              intent( in    ) :: n_cells(3) !< number of cells (and grid points)
    type(sll_t_matrix_csr), intent(   out ) :: spmat     !< sparse matrix

    !local variables
    sll_int32 :: n_nnz
    sll_int32 :: n_nnz_per_row
    sll_int32 :: ind, shift,row,i,j,k,k2

    ! Compute number of non-zero elements
    n_nnz_per_row = (2*deg(1)+1)*(2*deg(2)+1)*(2*deg(3)+1)
    n_nnz = (n_nnz_per_row*(n_cells(1)-deg(1))+(3*deg(1)**2+deg(1))*(2*deg(2)+1)*(2*deg(3)+1))*n_cells(2)*n_cells(3)

    ! Create the csr matrix
    call spmat%create( n_rows=(n_cells(1)+deg(1))*n_cells(2)*n_cells(3), n_cols=(n_cells(1)+deg(1))*n_cells(2)*n_cells(3), n_nnz=n_nnz )

    ! Set the row and column vector for the elements (row vector compressed)
    ! There are n_nnz_per_row elements in each of the rows
    spmat%arr_ia(1) = 1
    row=1
    do k = 1, n_cells(3)
       do j = 1, n_cells(2)
          do i=2, deg(1)+1
             row=row+1
             spmat%arr_ia(row) = spmat%arr_ia(row-1) + (deg(1) + i - 1)*(2*deg(2)+1)*(2*deg(3)+1)
          end do
          do i = deg(1)+2, n_cells(1)+1
             row=row+1
             spmat%arr_ia(row) = spmat%arr_ia(row-1) + n_nnz_per_row
          end do
          do i=n_cells(1)+2, n_cells(1)+deg(1)+1
             row=row+1
             spmat%arr_ia(row) = spmat%arr_ia(row-1) + (2*deg(1) + n_cells(1)+2 - i)*(2*deg(2)+1)*(2*deg(3)+1)
          end do
       end do
    end do

    SLL_ASSERT( row == (n_cells(1)+deg(1))*n_cells(2)*n_cells(3)+1)
    SLL_ASSERT( spmat%arr_ia(row) == n_nnz+1 )


    ind = 1
    !initialise mass matrix 
    !loop over number of rows
    do k=1, deg(3)
       do j=1,n_cells(2)
          do i=1,n_cells(1)+deg(1)
             do k2 = 1, k+deg(3)
                shift=(k2-1)*(n_cells(1)+deg(1))*n_cells(2)
                call loop2_clamped( i, j, shift, deg, n_cells, ind, spmat )
             end do

             do k2 = n_cells(3)+k-deg(3), n_cells(3)
                shift=(k2-1)*(n_cells(1)+deg(1))*n_cells(2)
                call loop2_clamped( i, j, shift, deg, n_cells, ind, spmat )
             end do
          end do
       end do
    end do
    do k = deg(3)+1, n_cells(3)-deg(3) 
       do j=1,n_cells(2)
          do i=1,n_cells(1)+deg(1)
             do k2 = k-deg(3), k+deg(3)
                shift=(k2-1)*(n_cells(1)+deg(1))*n_cells(2)
                call loop2_clamped( i, j, shift, deg, n_cells, ind, spmat )
             end do
          end do
       end do
    end do
    do k = n_cells(3)-deg(3)+1, n_cells(3)
       do j=1,n_cells(2)
          do i=1,n_cells(1)+deg(1)
             do k2 = 1, k+deg(3)-n_cells(3)
                shift=(k2-1)*(n_cells(1)+deg(1))*n_cells(2)
                call loop2_clamped( i, j, shift, deg, n_cells, ind, spmat )
             end do
             do k2 = k-deg(3), n_cells(3)
                shift=(k2-1)*(n_cells(1)+deg(1))*n_cells(2)
                call loop2_clamped( i, j, shift, deg, n_cells, ind, spmat )
             end do
          end do
       end do
    end do
    spmat%arr_a = 0.0_f64

    SLL_ASSERT( maxval(spmat%arr_ja) == spmat%n_cols )
    SLL_ASSERT( ind == spmat%n_nnz+1 )

  end subroutine sll_s_spline_fem_sparsity_mass3d_clamped

  !initialise sparse matrix in second direction
  subroutine loop2_clamped( i, j, shift, deg, n_cells, ind, spmat )
    sll_int32,              intent( in    ) :: i !< row in first direction
    sll_int32,              intent( in    ) :: j !< row in second direction
    sll_int32,              intent( in    ) :: shift !< shift in the third direction
    sll_int32,              intent( in    ) :: deg(3) !< spline degree for every direction
    sll_int32,              intent( in    ) :: n_cells(3) !< number of cells (and grid points)
    sll_int32,              intent( inout ) :: ind   !< index 
    type(sll_t_matrix_csr), intent( inout ) :: spmat !< sparse matrix
    !local variables
    sll_int32 :: j2, shift1

    !in the first deg(2) rows the entries are splitted in the second direction
    if(j<=deg(2))then
       do j2 = 1, j+deg(2)
          shift1=(j2-1)*(n_cells(1)+deg(1))+shift
          call loop1_clamped( i, shift1, deg(1), n_cells(1), ind, spmat )
       end do
       do j2 = n_cells(2)+j-deg(2), n_cells(2)
          shift1=(j2-1)*(n_cells(1)+deg(1))+shift
          call loop1_clamped( i, shift1, deg(1), n_cells(1), ind, spmat )
       end do

       !in rows: deg(2)+1-n_cells(2)-deg(2) we have consecutive entries in the second direction  
    elseif( j > deg(2) .and. j<= n_cells(2)-deg(2) )then
       do j2 = j-deg(2), j+deg(2)
          shift1=(j2-1)*(n_cells(1)+deg(1))+shift
          call loop1_clamped( i, shift1, deg(1), n_cells(1), ind, spmat )
       end do

       !in the last deg(2) rows the entries are splitted in the second direction
    elseif( j > n_cells(2)-deg(2) .and. j<= n_cells(2) )then
       do j2 = 1, j+deg(2)-n_cells(2)
          shift1=(j2-1)*(n_cells(1)+deg(1))+shift
          call loop1_clamped( i, shift1, deg(1), n_cells(1), ind, spmat )
       end do
       do j2 = j-deg(2), n_cells(2)
          shift1=(j2-1)*(n_cells(1)+deg(1))+shift
          call loop1_clamped(i,shift1, deg(1),n_cells(1),ind,spmat)
       end do
    else
       print*, 'error in loop2_clamped'
       stop
    end if

  end subroutine loop2_clamped

  !initialise sparse matrix in first direction
  subroutine loop1_clamped( i, shift, deg, n_cells, ind, spmat )
    sll_int32,              intent( in    ) :: i !< row in first direction
    sll_int32,              intent( in    ) :: shift !< shift in the second direction
    sll_int32,              intent( in    ) :: deg !< spline degree in first direction
    sll_int32,              intent( in    ) :: n_cells !< number of cells in the first direction 
    sll_int32,              intent( inout ) :: ind   !< index    
    type(sll_t_matrix_csr), intent( inout ) :: spmat !< sparse matrix
    !local variables
    sll_int32 :: i2

    !in the first deg(1) rows we have between deg(1)+1 and 2*deg(1) entries 
    if(i <= deg)then
       do i2 = 1, i+deg
          spmat%arr_ja(ind) = i2+shift
          ind = ind+1
       end do
       !in rows: deg(1)+1:n_cells(1) we have 2*deg(1)+1 consecutive entries in the first direction  
    elseif( i > deg .and. i <= n_cells) then
       do i2 = i-deg, i+deg
          spmat%arr_ja(ind) = i2+shift
          ind = ind+1
       end do
       !in the last deg(1) rows we have between 2*deg(1) and deg(1)+1 entries 
    elseif( i > n_cells .and. i <= n_cells + deg ) then
       do i2 = i-deg, n_cells+deg
          spmat%arr_ja(ind) = i2+shift
          ind = ind+1
       end do
    else
       print*, 'error in loop1_clamped'
       stop
    end if

  end subroutine loop1_clamped

  !---------------------------------------------------------------------------!
  !> Helper function to create sparsity pattern of the 3d mass matrix
  subroutine sll_s_spline_fem_sparsity_mixedmass3d( deg1, deg2, n_cells, spmat )
    sll_int32,              intent( in    ) :: deg1(3), deg2(3) !< spline deg for every direction
    sll_int32,              intent( in    ) :: n_cells(3) !< number of cells (and grid points)
    type(sll_t_matrix_csr), intent(   out ) :: spmat     !< sparse matrix

    !local variables
    sll_int32 :: n_nnz
    sll_int32 :: n_nnz_per_row
    sll_int32 :: ind, shift,row,i,j,k,l

    ! Compute number of non-zero elements
    n_nnz_per_row = (deg1(1)+deg2(1)+1)*(deg1(2)+deg2(2)+1)*(deg1(3)+deg2(3)+1)
    n_nnz = n_nnz_per_row*n_cells(1)*n_cells(2)*n_cells(3)
    ! Create the csr matrix
    call spmat%create( n_rows=n_cells(1)*n_cells(2)*n_cells(3), n_cols=n_cells(1)*n_cells(2)*n_cells(3), n_nnz=n_nnz )

    ! Set the row and column vector for the elements (row vector compressed)
    ! There are n_nnz_per_row elements in each of the rows
    spmat%arr_ia(1) = 1
    do row = 2, n_cells(1)*n_cells(2)*n_cells(3)+1
       spmat%arr_ia(row) = spmat%arr_ia(row-1) + n_nnz_per_row
    end do

    SLL_ASSERT( spmat%arr_ia(n_cells(1)*n_cells(2)*n_cells(3)+1) == n_nnz+1 )


    ind = 1
    !initialise mass matrix 
    !loop over number of rows
    do k = 1 ,n_cells(3)
       if(k <= deg2(3)) then  
          do j = 1, n_cells(2)
             do i = 1, n_cells(1)
                do l = 1, k+deg1(3)
                   shift = (l-1)*n_cells(1)*n_cells(2)
                   call mloop2( i, j, shift, deg1, deg2, n_cells, ind, spmat )
                end do
                do l = n_cells(3)+k-deg2(3), n_cells(3)
                   shift = (l-1)*n_cells(1)*n_cells(2)
                   call mloop2( i, j, shift, deg1, deg2, n_cells, ind, spmat )
                end do
                SLL_ASSERT( ind == (i+(j-1)*n_cells(1)+(k-1)*n_cells(1)*n_cells(2))*n_nnz_per_row+1 )
             end do
          end do
       elseif( k > deg2(3) .and. k <= n_cells(3)-deg1(3)) then
          do j = 1, n_cells(2)
             do i = 1, n_cells(1)
                do l = k-deg2(3), k+deg1(3)
                   shift = (l-1)*n_cells(1)*n_cells(2)
                   call mloop2( i, j, shift, deg1, deg2, n_cells, ind, spmat )
                end do
                SLL_ASSERT( ind == (i+(j-1)*n_cells(1)+(k-1)*n_cells(1)*n_cells(2))*n_nnz_per_row+1 )
             end do
          end do
       else
          do j = 1, n_cells(2)
             do i = 1, n_cells(1)
                do l = 1, k+deg1(3)-n_cells(3)
                   shift = (l-1)*n_cells(1)*n_cells(2)
                   call mloop2( i, j, shift, deg1, deg2, n_cells, ind, spmat )
                end do
                do l = k-deg2(3), n_cells(3)
                   shift = (l-1)*n_cells(1)*n_cells(2)
                   call mloop2( i, j, shift, deg1, deg2, n_cells, ind, spmat )
                end do
                SLL_ASSERT( ind == (i+(j-1)*n_cells(1)+(k-1)*n_cells(1)*n_cells(2))*n_nnz_per_row+1 )
             end do
          end do
       end if
    end do
    spmat%arr_a = 0.0_f64

    SLL_ASSERT( maxval(spmat%arr_ja) == spmat%n_cols )
    SLL_ASSERT( ind == spmat%n_nnz+1 )

  end subroutine sll_s_spline_fem_sparsity_mixedmass3d

  !initialise sparse matrix in second direction
  subroutine mloop2( i, j, shift, deg1, deg2, n_cells, ind, spmat )
    sll_int32,              intent( in    ) :: i !< row in first direction
    sll_int32,              intent( in    ) :: j !< row in second direction
    sll_int32,              intent( in    ) :: shift !< shift in the third direction
    sll_int32,              intent( in    ) :: deg1(3), deg2(3) !< spline degrees for every direction
    sll_int32,              intent( in    ) :: n_cells(3) !< number of cells (and grid points)
    sll_int32,              intent( inout ) :: ind   !< index 
    type(sll_t_matrix_csr), intent( inout ) :: spmat !< sparse matrix
    !local variables
    sll_int32 :: l, shift1

    !in the first deg2(2) rows the entries are splitted in the second direction
    if(j <= deg2(2))then
       do l = 1, j+deg1(2)
          shift1 = (l-1)*n_cells(1)+shift
          call mloop1( i, shift1, deg1(1), deg2(1), n_cells(1), ind, spmat )
       end do
       do l = n_cells(2)+j-deg2(2), n_cells(2)
          shift1 = (l-1)*n_cells(1)+shift
          call mloop1( i, shift1, deg1(1), deg2(1), n_cells(1), ind, spmat )
       end do

       !in rows: deg2(2)+1-n_cells(2)-deg1(2) we have consecutive entries in the second direction  
    elseif( j > deg2(2) .and. j <= n_cells(2)-deg1(2))then
       do l = j-deg2(2), j+deg1(2)
          shift1 = (l-1)*n_cells(1)+shift
          call mloop1( i, shift1, deg1(1), deg2(1), n_cells(1), ind, spmat )
       end do

       !in the last deg1(2) rows the entries are splitted in the second direction
    else
       do l = 1, j+deg1(2)-n_cells(2)
          shift1 = (l-1)*n_cells(1)+shift
          call mloop1( i, shift1, deg1(1), deg2(1), n_cells(1), ind, spmat )
       end do
       do l = j-deg2(2), n_cells(2)
          shift1 = (l-1)*n_cells(1)+shift
          call mloop1(i, shift1, deg1(1), deg2(1), n_cells(1), ind,spmat)
       end do
    end if

  end subroutine mloop2

  !initialise sparse matrix in first direction
  subroutine mloop1( i, shift, deg1, deg2, n_cells, ind, spmat )
    sll_int32,              intent( in    ) :: i !< row in first direction
    sll_int32,              intent( in    ) :: shift !< shift in the second direction
    sll_int32,              intent( in    ) :: deg1, deg2 !< spline degrees in first direction
    sll_int32,              intent( in    ) :: n_cells !< number of cells in the first direction 
    sll_int32,              intent( inout ) :: ind   !< index     
    type(sll_t_matrix_csr), intent( inout ) :: spmat !< sparse matrix
    !local variables
    sll_int32 :: l

    !in the first deg2(1) rows the entries are splitted in the first direction
    if(i <= deg2)then
       do l = 1, i+deg1
          spmat%arr_ja(ind) = l+shift
          ind = ind+1
       end do
       do l = n_cells+i-deg2, n_cells 
          spmat%arr_ja(ind) = l+shift
          ind = ind+1
       end do

       !in rows: deg2(1)+1 to n_cells(1)-deg1(1) we have consecutive entries in the first direction  
    elseif( i > deg2 .and. i<= n_cells-deg1) then
       do l = i-deg2, i+deg1  
          spmat%arr_ja(ind) = l+shift
          ind = ind+1
       end do

       !in the last deg1(1) rows the entries are splitted in the first direction
    elseif(i > n_cells-deg1 .and. i <= n_cells) then
       do l = 1, i+deg1-n_cells
          spmat%arr_ja(ind) = l+shift
          ind = ind+1
       end do
       do l = i-deg2, n_cells
          spmat%arr_ja(ind) = l+shift
          ind = ind+1
       end do
    end if
  end subroutine mloop1

  !---------------------------------------------------------------------------!
  !> Helper function to create sparsity pattern of the 3d clamped mixed mass matrix
  subroutine sll_s_spline_fem_sparsity_mixedmass3d_clamped( deg1, deg2, n_cells, spmat )
    sll_int32,              intent( in    ) :: deg1(3), deg2(3) !< spline degrees for every direction
    sll_int32,              intent( in    ) :: n_cells(3) !< number of cells (and grid points)
    type(sll_t_matrix_csr), intent(   out ) :: spmat     !< sparse matrix
    !local variables
    sll_int32 :: n_nnz
    sll_int32 :: n_nnz_per_row
    sll_int32 :: ind, shift,row,i,j,k,l

    ! Compute number of non-zero elements
    n_nnz_per_row = (deg1(1)+deg2(1)+1)*(deg1(2)+deg2(2)+1)*(deg1(3)+deg2(3)+1)
    n_nnz =  (n_nnz_per_row*(n_cells(1)-deg1(1))+deg1(1)*(2*deg2(1)+deg1(1)+1)*(deg1(2)+deg2(2)+1)*(deg1(3)+deg2(3)+1))*n_cells(2)*n_cells(3)
    ! Create the csr matrix
    call spmat%create( n_rows=(n_cells(1)+deg1(1))*n_cells(2)*n_cells(3), n_cols=(n_cells(1)+deg2(1))*n_cells(2)*n_cells(3), n_nnz=n_nnz )

    ! Set the row and column vector for the elements (row vector compressed)
    ! There are n_nnz_per_row elements in each of the rows
    spmat%arr_ia(1) = 1
    row=1
    do k = 1, n_cells(3)
       do j = 1, n_cells(2)
          do i=2, deg1(1)+1
             row=row+1
             spmat%arr_ia(row) = spmat%arr_ia(row-1) + (deg2(1) + i - 1)*(deg1(2)+deg2(2)+1)*(deg1(3)+deg2(3)+1)
          end do
          do i = deg1(1)+2, n_cells(1)+1
             row=row+1
             spmat%arr_ia(row) = spmat%arr_ia(row-1) + n_nnz_per_row
          end do
          do i=n_cells(1)+2, n_cells(1)+deg1(1)+1
             row=row+1
             spmat%arr_ia(row) = spmat%arr_ia(row-1) + (deg1(1)+deg2(1) + n_cells(1)+2 - i)*(deg1(2)+deg2(2)+1)*(deg1(3)+deg2(3)+1)
          end do
       end do
    end do

    SLL_ASSERT( row == (n_cells(1)+deg1(1))*n_cells(2)*n_cells(3)+1)
    SLL_ASSERT( spmat%arr_ia(row) == n_nnz+1 )


    ind = 1
    !initialise mass matrix 
    !loop over number of rows
    do k=1,n_cells(3)
       if(k<=deg2(3)) then  
          do j=1,n_cells(2)
             do i=1,n_cells(1)+deg1(1)
                do l = 1, k+deg1(3)
                   shift=(l-1)*(n_cells(1)+deg2(1))*n_cells(2)
                   call mloop2_clamped( i, j, shift, deg1, deg2, n_cells, ind, spmat )
                end do
                do l = n_cells(3)+k-deg2(3), n_cells(3)
                   shift=(l-1)*(n_cells(1)+deg2(1))*n_cells(2)
                   call mloop2_clamped( i, j, shift, deg1, deg2, n_cells, ind, spmat )
                end do
             end do
          end do
       elseif( k > deg2(3) .and. k<= n_cells(3)-deg1(3)) then
          do j=1,n_cells(2)
             do i=1,n_cells(1)+deg1(1)
                do l = k-deg2(3), k+deg1(3)
                   shift=(l-1)*(n_cells(1)+deg2(1))*n_cells(2)
                   call mloop2_clamped( i, j, shift, deg1, deg2, n_cells, ind, spmat )
                end do
             end do
          end do
       else
          do j=1,n_cells(2)
             do i=1, n_cells(1)+deg1(1)
                do l = 1, k+deg1(3)-n_cells(3)
                   shift=(l-1)*(n_cells(1)+deg2(1))*n_cells(2)
                   call mloop2_clamped( i, j, shift, deg1, deg2, n_cells, ind, spmat )
                end do
                do l = k-deg2(3), n_cells(3)
                   shift=(l-1)*(n_cells(1)+deg2(1))*n_cells(2)
                   call mloop2_clamped( i, j, shift, deg1, deg2, n_cells, ind, spmat )
                end do
             end do
          end do
       end if
    end do
    spmat%arr_a = 0.0_f64

    SLL_ASSERT( maxval(spmat%arr_ja) == spmat%n_cols )
    SLL_ASSERT( ind == spmat%n_nnz+1 )

  end subroutine sll_s_spline_fem_sparsity_mixedmass3d_clamped

  !initialise sparse mass matrix in second direction
  subroutine mloop2_clamped( i, j, shift, deg1, deg2, n_cells, ind, spmat )
    sll_int32,              intent( in    ) :: i !< row in first direction
    sll_int32,              intent( in    ) :: j !< row in second direction
    sll_int32,              intent( in    ) :: shift !< shift in the third direction
    sll_int32,              intent( in    ) :: deg1(3), deg2(3) !< spline degrees for every direction
    sll_int32,              intent( in    ) :: n_cells(3) !< number of cells (and grid points)
    sll_int32,              intent( inout ) :: ind   !< index 
    type(sll_t_matrix_csr), intent( inout ) :: spmat !< sparse matrix
    !local variables
    sll_int32 :: l, shift1

    !in the first deg2(2) rows the entries are splitted in the second direction
    if(j<=deg2(2))then
       do l = 1, j+deg1(2)
          shift1=(l-1)*(n_cells(1)+deg2(1))+shift
          call mloop1_clamped( i, shift1, deg1(1), deg2(1), n_cells(1), ind, spmat )
       end do
       do l = n_cells(2)+j-deg2(2), n_cells(2)
          shift1=(l-1)*(n_cells(1)+deg2(1))+shift
          call mloop1_clamped( i, shift1, deg1(1), deg2(1), n_cells(1), ind, spmat )
       end do

       !in rows: deg2(2)+1-n_cells(2)-deg1(2) we have consecutive entries in the second direction  
    elseif( j > deg2(2) .and. j<= n_cells(2)-deg1(2))then
       do l = j-deg2(2), j+deg1(2)
          shift1=(l-1)*(n_cells(1)+deg2(1))+shift
          call mloop1_clamped( i, shift1, deg1(1), deg2(1), n_cells(1), ind, spmat )
       end do

       !in the last deg1(2) rows the entries are splitted in the second direction
    else
       do l = 1, j+deg1(2)-n_cells(2)
          shift1=(l-1)*(n_cells(1)+deg2(1))+shift
          call mloop1_clamped( i, shift1, deg1(1), deg2(1), n_cells(1), ind, spmat )
       end do
       do l = j-deg2(2), n_cells(2)
          shift1=(l-1)*(n_cells(1)+deg2(1))+shift
          call mloop1_clamped(i,shift1, deg1(1), deg2(1),n_cells(1),ind,spmat)
       end do
    end if

  end subroutine mloop2_clamped

  !initialise sparse mass matrix in first direction
  subroutine mloop1_clamped( i, shift, deg1, deg2, n_cells, ind, spmat )
    sll_int32,              intent( in    ) :: i !< row in first direction
    sll_int32,              intent( in    ) :: shift !< shift in the second direction
    sll_int32,              intent( in    ) :: deg1, deg2 !< spline degrees in first direction
    sll_int32,              intent( in    ) :: n_cells !< number of cells in the first direction 
    sll_int32,              intent( inout ) :: ind   !< index     
    type(sll_t_matrix_csr), intent( inout ) :: spmat !< sparse matrix
    !local variables
    sll_int32 :: l

    !in the first deg1(1) rows we have between deg2+1 and deg2+deg1 entries 
    if(i<= deg1)then
       do l = 1, i+deg2
          spmat%arr_ja(ind) = l+shift
          ind = ind+1
       end do
       !in rows: deg1+1 to n_cells we have deg1+deg2+1 consecutive entries in the first direction  
    elseif( i > deg1 .and. i<= n_cells) then
       do l = i-deg1, i+deg2  
          spmat%arr_ja(ind) = l+shift
          ind = ind+1
       end do
       !in the last deg1 rows we have between deg2+deg1 and deg2+1 entries
    elseif( i > n_cells .and. i <= n_cells + deg1 ) then
       do l = i-deg1, n_cells+deg2
          spmat%arr_ja(ind) = l+shift
          ind = ind+1
       end do
    end if
  end subroutine mloop1_clamped


  !---------------------------------------------------------------------------!
  !>  Assemble the given row of the 3d mass matrix
  subroutine assemble_mass3d ( deg, n_cells, mass_line, matrix, row, ind )
    sll_int32,              intent( in    ) :: deg(3)        !< spline degree for every direction
    sll_int32,              intent( in    ) :: n_cells(3)    !< number of cells (and grid points)
    sll_real64,             intent( in    ) :: mass_line((2*deg(1)+1)*(2*deg(2)+1)*(2*deg(3)+1)) !< massline for one row of the sparse mass matrix
    type(sll_t_matrix_csr), intent( inout ) :: matrix        !< sparse mass matrix
    sll_int32,              intent( in    ) :: row(3)        !< current row of the sparse matrix
    sll_int32,              intent( inout ) :: ind           !< index in the vector containing the entries of the sparse matrix
    !local variables
    sll_int32 ::  column, shift

    !assemble massmatrix in third direction
    ! For the first deg(3) rows we need to put the first part to the back due to periodic boundaries
    if(row(3) <= deg(3)) then
       do column = 2-row(3)+deg(3), deg(3)
          shift = (column-1)*(2*deg(1)+1)*(2*deg(2)+1)
          call assemble2( deg, n_cells, mass_line, matrix, row, ind, shift )
       end do
       do column = deg(3)+1,2*deg(3)+1
          shift = (column-1)*(2*deg(1)+1)*(2*deg(2)+1)
          call assemble2( deg, n_cells, mass_line, matrix, row, ind, shift )
       end do

       do column = 1, deg(3)-row(3)+1
          shift=(column-1)*(2*deg(1)+1)*(2*deg(2)+1)
          call assemble2( deg, n_cells, mass_line, matrix, row, ind, shift )
       end do
    else if(row(3) >= deg(3)+1 .and. row(3)<= n_cells(3)-deg(3)) then
       do column = 1, 2*deg(3)+1
          shift = (column-1)*(2*deg(1)+1)*(2*deg(2)+1)
          call assemble2( deg, n_cells, mass_line, matrix, row, ind, shift )
       end do

       ! For the last deg(3) rows, we need to put the second part to the front due to periodic bounaries
    else if(row(3) >= n_cells(3)-deg(3)+1 .and. row(3) <= n_cells(3)) then

       do column = n_cells(3)-row(3)+deg(3)+2, 2*deg(3)+1
          shift = (column-1)*(2*deg(1)+1)*(2*deg(2)+1)
          call assemble2( deg, n_cells, mass_line, matrix, row, ind, shift )
       end do

       do column = 1, deg(3)+1
          shift = (column-1)*(2*deg(1)+1)*(2*deg(2)+1)
          call assemble2( deg, n_cells, mass_line, matrix, row, ind, shift )
       end do

       do column = 2+deg(3), n_cells(3)-row(3)+deg(3)+1 
          shift = (column-1)*(2*deg(1)+1)*(2*deg(2)+1)
          call assemble2( deg, n_cells, mass_line, matrix, row, ind, shift )
       end do
    else
       SLL_ERROR('assemble','error in row in assemble_mass')
    end if

    SLL_ASSERT( ind == (row(1)+(row(2)-1)*n_cells(1)+(row(3)-1)*n_cells(1)*n_cells(2))*(2*deg(1)+1)*(2*deg(2)+1)*(2*deg(3)+1)+1)

  end subroutine assemble_mass3d

  !assemble mass matrix in second direction
  subroutine assemble2( deg, n_cells, mass_line, matrix, row, ind, shift )
    sll_int32,              intent( in    ) :: deg(3) !< spline degree for every direction
    sll_int32,              intent( in    ) :: n_cells(3) !< number of cells (and grid points) 
    sll_real64,             intent( in    ) :: mass_line((2*deg(1)+1)*(2*deg(2)+1)*(2*deg(3)+1)) !< massline for one row of the sparse mass matrix
    type(sll_t_matrix_csr), intent( inout ) :: matrix !< sparse matrix
    sll_int32,              intent( in    ) :: row(3) !< current row of the sparse matrix
    sll_int32,              intent( inout ) :: ind  !< index 
    sll_int32,              intent( in    ) :: shift !< shift in the third direction
    !local variables
    sll_int32 ::  column, shift1
    ! For the first deg(2) rows we need to put the first part to the back due to periodic boundaries
    if(row(2) <= deg(2)) then
       do column = 2-row(2)+deg(2), deg(2)
          shift1 = (column-1)*(2*deg(1)+1)+shift
          call assemble1( deg, n_cells(1), mass_line, matrix, row(1), ind, shift1 )
       end do
       do column = deg(2)+1, 2*deg(2)+1
          shift1 = (column-1)*(2*deg(1)+1)+shift
          call assemble1( deg, n_cells(1), mass_line, matrix, row(1), ind, shift1 )
       end do

       do column = 1, deg(2)-row(2)+1
          shift1 = (column-1)*(2*deg(1)+1)+shift
          call assemble1( deg, n_cells(1), mass_line, matrix, row(1), ind, shift1 )
       end do
    else if(row(2) >= deg(2)+1 .and. row(2)<= n_cells(2)-deg(2)) then
       do column = 1, 2*deg(2)+1
          shift1 = (column-1)*(2*deg(1)+1)+shift
          call assemble1( deg, n_cells(1), mass_line, matrix, row(1), ind, shift1 )
       end do

       ! For the last deg(2) rows, we need to put the second part to the front due to periodic boundaries
    else if(row(2) >= n_cells(2)-deg(2)+1 .and. row(2) <= n_cells(2)) then

       do column = n_cells(2)-row(2)+deg(2)+2, 2*deg(2)+1
          shift1 = (column-1)*(2*deg(1)+1)+shift
          call assemble1( deg, n_cells(1), mass_line, matrix, row(1), ind, shift1 )
       end do

       do column = 1, deg(2)+1
          shift1 = (column-1)*(2*deg(1)+1)+shift
          call assemble1(deg, n_cells(1), mass_line, matrix, row(1), ind, shift1)
       end do

       do column = 2+deg(2), n_cells(2)-row(2)+deg(2)+1
          shift1 = (column-1)*(2*deg(1)+1)+shift
          call assemble1( deg, n_cells(1), mass_line, matrix, row(1), ind, shift1 )
       end do
    else
       SLL_ERROR('assemble','error in row in assemble_mass')
    end if

  end subroutine assemble2

  !assemble mass matrix in first direction
  subroutine assemble1( deg, n_cells, mass_line, matrix, row, ind, shift )
    sll_int32,              intent( in    ) :: deg(3) !< spline degree for every direction
    sll_int32,              intent( in    ) :: n_cells !< number of cells in the first direction
    sll_real64,             intent( in    ) :: mass_line((2*deg(1)+1)*(2*deg(2)+1)*(2*deg(3)+1)) !< massline for one row of the sparse mass matrix
    type(sll_t_matrix_csr), intent( inout ) :: matrix !< sparse matrix
    sll_int32,              intent( in    ) :: row  !< current row in the first direction 
    sll_int32,              intent( inout ) :: ind !< index
    sll_int32,              intent( in    ) :: shift !< shift in the second direction
    !local variables
    sll_int32 ::  column

    ! For the first deg(1) rows we need to put the first part to the back due to periodic boundaries
    if(row <= deg(1)) then
       do column = 2-row+deg(1), deg(1)
          matrix%arr_a(1,ind) = mass_line(column+shift)

          ind = ind+1
       end do
       matrix%arr_a(1,ind:ind+deg(1)) = mass_line(deg(1)+1+shift:2*deg(1)+1+shift)
       ind = ind+deg(1)+1

       do column = 1, deg(1)-row+1
          matrix%arr_a(1,ind) = mass_line(column+shift) 
          ind = ind+1
       end do
    else if(row > deg(1) .and. row <= n_cells-deg(1)) then
       matrix%arr_a(1,ind:ind+2*deg(1)) = mass_line(1+shift:2*deg(1)+1+shift)

       ind = ind+2*deg(1)+1
       ! For the last deg(1) rows, we need to put the second part to the front due to periodic boundaries
    else if(row >= n_cells-deg(1)+1 .and. row <= n_cells) then

       do column = n_cells-row+deg(1)+2, 2*deg(1)+1
          matrix%arr_a(1,ind) = mass_line(column+shift)
          ind = ind+1
       end do

       matrix%arr_a(1,ind:ind+deg(1)) = mass_line(1+shift:deg(1)+1+shift)
       ind = ind+deg(1)+1

       do column = 2+deg(1), n_cells-row+deg(1)+1
          matrix%arr_a(1,ind) = mass_line(column+shift)
          ind = ind+1
       end do
    else
       SLL_ERROR('assemble','error in row in assemble_mass')
    end if

  end subroutine assemble1

  !---------------------------------------------------------------------------!
  !>  Assemble the given row of the clamped mass matrix
  subroutine assemble_mass3d_clamped ( deg, n_cells, mass_line, matrix, row, ind )
    sll_int32,              intent( in    ) :: deg(3)                                       !< spline degree for every direction
    sll_int32,              intent( in    ) :: n_cells(3)                                    !< number of rows (and grid points)
    sll_real64,             intent( in    ) :: mass_line((2*deg(1)+1)*(2*deg(2)+1)*(2*deg(3)+1)) !< massline for one row of the sparse mass matrix
    type(sll_t_matrix_csr), intent( inout ) :: matrix                                       !< sparse mass matrix
    sll_int32,              intent( in    ) :: row(3)                                       !< current row of the sparse matrix
    sll_int32,              intent( inout ) :: ind                                          !< index in the vector containing the entries of the sparse matrix
    !local variables
    sll_int32 ::  column, shift, begin

    begin= ind

    !assemble massmatrix in third direction
    ! For the first deg(3) rows we need to put the first part to the back due to periodic boundaries
    if(row(3)<=deg(3)) then
       do column = 2-row(3)+deg(3), deg(3)
          shift=(column-1)*(2*deg(1)+1)*(2*deg(2)+1)
          call assemble2_clamped( deg, n_cells, mass_line, matrix, row, ind, shift )
       end do
       do column= deg(3)+1,2*deg(3)+1
          shift=(column-1)*(2*deg(1)+1)*(2*deg(2)+1)
          call assemble2_clamped( deg, n_cells, mass_line, matrix, row, ind, shift )
       end do

       do column = 1, deg(3)-row(3)+1
          shift=(column-1)*(2*deg(1)+1)*(2*deg(2)+1)
          call assemble2_clamped( deg, n_cells, mass_line, matrix, row, ind, shift )
       end do
    else if(row(3) >= deg(3)+1 .and. row(3)<= n_cells(3)-deg(3)) then
       do column=1, 2*deg(3)+1
          shift=(column-1)*(2*deg(1)+1)*(2*deg(2)+1)
          call assemble2_clamped( deg, n_cells, mass_line, matrix, row, ind, shift )
       end do

       ! For the last deg(3) rows, we need to put the second part to the front due to periodic boundaries
    else if(row(3) >= n_cells(3)-deg(3)+1 .and. row(3) <= n_cells(3)) then

       do column = n_cells(3)-row(3)+deg(3)+2, 2*deg(3)+1
          shift=(column-1)*(2*deg(1)+1)*(2*deg(2)+1)
          call assemble2_clamped( deg, n_cells, mass_line, matrix, row, ind, shift )
       end do

       do column=1, deg(3)+1
          shift=(column-1)*(2*deg(1)+1)*(2*deg(2)+1)
          call assemble2_clamped( deg, n_cells, mass_line, matrix, row, ind, shift )
       end do

       do column = 2+deg(3), n_cells(3)-row(3)+deg(3)+1 
          shift=(column-1)*(2*deg(1)+1)*(2*deg(2)+1)
          call assemble2_clamped( deg, n_cells, mass_line, matrix, row, ind, shift )
       end do
    else
       SLL_ERROR('assemble','error in row in assemble3_mass_clamped')
    end if
    SLL_ASSERT( ind == begin+(2*deg(1)+1)*(2*deg(2)+1)*(2*deg(3)+1))

  end subroutine assemble_mass3d_clamped

  !assemble mass matrix in second direction
  subroutine assemble2_clamped( deg, n_cells, mass_line, matrix, row, ind, shift )
    sll_int32,              intent( in    ) :: deg(3) !< spline degree for every direction
    sll_int32,              intent( in    ) :: n_cells(3) !< number of cells (and grid points) 
    sll_real64,             intent( in    ) :: mass_line((2*deg(1)+1)*(2*deg(2)+1)*(2*deg(3)+1)) !< massline for one row of the sparse mass matrix
    type(sll_t_matrix_csr), intent( inout ) :: matrix !< sparse matrix
    sll_int32,              intent( in    ) :: row(3) !< current row of the sparse matrix
    sll_int32,              intent( inout ) :: ind  !< index 
    sll_int32,              intent( in    ) :: shift !< shift in the third direction
    !local variables
    sll_int32 ::  column, shift1
    ! For the first deg(2) rows we need to put the first part to the back due to periodic boundaries
    if(row(2)<=deg(2)) then
       do column = 2-row(2)+deg(2), deg(2)
          shift1=(column-1)*(2*deg(1)+1)+shift
          call assemble1_clamped( deg, n_cells(1), mass_line, matrix, row(1), ind, shift1 )
       end do
       do column= deg(2)+1,2*deg(2)+1
          shift1=(column-1)*(2*deg(1)+1)+shift
          call assemble1_clamped( deg, n_cells(1), mass_line, matrix, row(1), ind, shift1 )
       end do

       do column = 1, deg(2)-row(2)+1
          shift1=(column-1)*(2*deg(1)+1)+shift
          call assemble1_clamped( deg, n_cells(1), mass_line, matrix, row(1), ind, shift1 )
       end do
    else if(row(2) >= deg(2)+1 .and. row(2)<= n_cells(2)-deg(2)) then
       do column=1, 2*deg(2)+1
          shift1=(column-1)*(2*deg(1)+1)+shift
          call assemble1_clamped( deg, n_cells(1), mass_line, matrix, row(1), ind, shift1 )
       end do
       ! For the last deg(2) rows, we need to put the second part to the front due to periodic boundaries
    else if(row(2) >= n_cells(2)-deg(2)+1 .and. row(2) <= n_cells(2)) then

       do column = n_cells(2)-row(2)+deg(2)+2, 2*deg(2)+1
          shift1=(column-1)*(2*deg(1)+1)+shift
          call assemble1_clamped( deg, n_cells(1), mass_line, matrix, row(1), ind, shift1 )
       end do

       do column=1, deg(2)+1
          shift1=(column-1)*(2*deg(1)+1)+shift
          call assemble1_clamped(deg, n_cells(1), mass_line, matrix, row(1), ind, shift1)
       end do

       do column = 2+deg(2), n_cells(2)-row(2)+deg(2)+1
          shift1=(column-1)*(2*deg(1)+1)+shift
          call assemble1_clamped( deg, n_cells(1), mass_line, matrix, row(1), ind, shift1 )
       end do
    else
       SLL_ERROR('assemble','error in row in assemble2_mass_clamped')
    end if

  end subroutine assemble2_clamped

  !assemble mass matrix in first direction
  subroutine assemble1_clamped( deg, n_cells, mass_line, matrix, row, ind, shift )
    sll_int32,              intent( in    ) :: deg(3) !< spline degree for every direction
    sll_int32,              intent( in    ) :: n_cells !< number of cells in the first direction
    sll_real64,             intent( in    ) :: mass_line((2*deg(1)+1)*(2*deg(2)+1)*(2*deg(3)+1)) !< massline for one row of the sparse mass matrix
    type(sll_t_matrix_csr), intent( inout ) :: matrix !< sparse matrix
    sll_int32,              intent( in    ) :: row  !< current row in the first direction 
    sll_int32,              intent( inout ) :: ind !< index
    sll_int32,              intent( in    ) :: shift !< shift in the second direction

    if(row >= 2*deg(1) .and. row<= n_cells-deg(1)+1) then
       matrix%arr_a(1,ind:ind+2*deg(1))=mass_line(1+shift:2*deg(1)+1+shift)

       ind = ind+2*deg(1)+1
    else
       SLL_ERROR('assemble','error in row in assemble1_mass_clamped')
    end if

  end subroutine assemble1_clamped

  !>  Assemble the boundary part of the clamped mass matrix
  subroutine assemble_mass3d_clamped_boundary ( deg, n_cells, mass_line, matrix, row, ind )
    sll_int32,              intent( in    ) :: deg(3)                                       !< spline degree for every direction
    sll_int32,              intent( in    ) :: n_cells(3)                                    !< number of rows (and grid points)
    sll_real64,             intent( in    ) :: mass_line((7*deg(1)**2-deg(1)-2)/2*(2*deg(2)+1)*(2*deg(3)+1)) !< boundary massline for the sparse mass matrix
    type(sll_t_matrix_csr), intent( inout ) :: matrix                                       !< sparse mass matrix
    sll_int32,              intent( in    ) :: row(3)                                      !< current row of the sparse matrix
    sll_int32,              intent( inout ) :: ind                                          !< index in the vector containing the entries of the sparse matrix
    !local variables
    sll_int32 ::  column, shift

    !assemble massmatrix in third direction
    ! For the first deg(3) rows we need to put the first part to the back due to periodic boundaries
    if(row(3)<=deg(3)) then
       do column = 2-row(3)+deg(3), deg(3)
          shift=(column-1)*(7*deg(1)**2-deg(1)-2)/2*(2*deg(2)+1)
          call assemble2_clamped_boundary( deg, n_cells, mass_line, matrix, row, ind, shift )
       end do
       do column= deg(3)+1,2*deg(3)+1
          shift=(column-1)*(7*deg(1)**2-deg(1)-2)/2*(2*deg(2)+1)
          call assemble2_clamped_boundary( deg, n_cells, mass_line, matrix, row, ind, shift )
       end do
       do column = 1, deg(3)-row(3)+1
          shift=(column-1)*(7*deg(1)**2-deg(1)-2)/2*(2*deg(2)+1)
          call assemble2_clamped_boundary( deg, n_cells, mass_line, matrix, row, ind, shift )
       end do
    else if(row(3) >= deg(3)+1 .and. row(3)<= n_cells(3)-deg(3)) then
       do column=1, 2*deg(3)+1
          shift=(column-1)*(7*deg(1)**2-deg(1)-2)/2*(2*deg(2)+1)
          call assemble2_clamped_boundary( deg, n_cells, mass_line, matrix, row, ind, shift )
       end do

       ! For the last deg(3) rows, we need to put the second part to the front due to periodic boundaries
    else if(row(3) >= n_cells(3)-deg(3)+1 .and. row(3) <= n_cells(3)) then
       do column = n_cells(3)-row(3)+deg(3)+2, 2*deg(3)+1
          shift=(column-1)*(7*deg(1)**2-deg(1)-2)/2*(2*deg(2)+1)
          call assemble2_clamped_boundary( deg, n_cells, mass_line, matrix, row, ind, shift )
       end do
       do column=1, deg(3)+1
          shift=(column-1)*(7*deg(1)**2-deg(1)-2)/2*(2*deg(2)+1)
          call assemble2_clamped_boundary( deg, n_cells, mass_line, matrix, row, ind, shift )
       end do
       do column = 2+deg(3), n_cells(3)-row(3)+deg(3)+1 
          shift=(column-1)*(7*deg(1)**2-deg(1)-2)/2*(2*deg(2)+1)
          call assemble2_clamped_boundary( deg, n_cells, mass_line, matrix, row, ind, shift )
       end do
    else
       SLL_ERROR('assemble','error in row in assemble3_mass_clamped_boundary')
    end if

  end subroutine assemble_mass3d_clamped_boundary

  !assemble boundary part in second direction
  subroutine assemble2_clamped_boundary( deg, n_cells, mass_line, matrix, row, ind, shift )
    sll_int32,              intent( in    ) :: deg(3) !< spline degree for every direction
    sll_int32,              intent( in    ) :: n_cells(3) !< number of cells (and grid points) 
    sll_real64,             intent( in    ) :: mass_line((7*deg(1)**2-deg(1)-2)/2*(2*deg(2)+1)*(2*deg(3)+1)) !< boundary massline for the sparse mass matrix
    type(sll_t_matrix_csr), intent( inout ) :: matrix !< sparse matrix
    sll_int32,              intent( in    ) :: row(3) !< current row of the sparse matrix
    sll_int32,              intent( inout ) :: ind  !< index 
    sll_int32,              intent( in    ) :: shift !< shift in the third direction
    !local variables
    sll_int32 ::  column, shift1

    ! For the first deg(2) rows we need to put the first part to the back due to periodic boundaries
    if(row(2)<=deg(2)) then
       do column = 2-row(2)+deg(2), deg(2)
          shift1=(column-1)*(7*deg(1)**2-deg(1)-2)/2+shift
          call assemble1_clamped_boundary( deg, n_cells(1), mass_line, matrix, row(1), ind, shift1 )
       end do
       do column= deg(2)+1,2*deg(2)+1
          shift1=(column-1)*(7*deg(1)**2-deg(1)-2)/2+shift
          call assemble1_clamped_boundary( deg, n_cells(1), mass_line, matrix, row(1), ind, shift1 )
       end do

       do column = 1, deg(2)-row(2)+1
          shift1=(column-1)*(7*deg(1)**2-deg(1)-2)/2+shift
          call assemble1_clamped_boundary( deg, n_cells(1), mass_line, matrix, row(1), ind, shift1 )
       end do
    else if(row(2) >= deg(2)+1 .and. row(2)<= n_cells(2)-deg(2)) then
       do column=1, 2*deg(2)+1
          shift1=(column-1)*(7*deg(1)**2-deg(1)-2)/2+shift
          call assemble1_clamped_boundary( deg, n_cells(1), mass_line, matrix, row(1), ind, shift1 )
       end do

       ! For the last deg(2) rows, we need to put the second part to the front due to periodic boundaries
    else if(row(2) >= n_cells(2)-deg(2)+1 .and. row(2) <= n_cells(2)) then

       do column = n_cells(2)-row(2)+deg(2)+2, 2*deg(2)+1
          shift1=(column-1)*(7*deg(1)**2-deg(1)-2)/2+shift
          call assemble1_clamped_boundary( deg, n_cells(1), mass_line, matrix, row(1), ind, shift1 )
       end do

       do column=1, deg(2)+1
          shift1=(column-1)*(7*deg(1)**2-deg(1)-2)/2+shift
          call assemble1_clamped_boundary( deg, n_cells(1), mass_line, matrix, row(1), ind, shift1)
       end do

       do column = 2+deg(2), n_cells(2)-row(2)+deg(2)+1
          shift1=(column-1)*(7*deg(1)**2-deg(1)-2)/2+shift
          call assemble1_clamped_boundary( deg, n_cells(1), mass_line, matrix, row(1), ind, shift1 )
       end do
    else
       SLL_ERROR('assemble','error in row in assemble2_mass_clamped_boundary')
    end if

  end subroutine assemble2_clamped_boundary

  !assemble boundary part in first direction
  subroutine assemble1_clamped_boundary( deg, n_cells, mass_line, matrix, row, ind, shift )
    sll_int32,              intent( in    ) :: deg(3) !< spline degree for every direction
    sll_int32,              intent( in    ) :: n_cells !< number of cells in the first direction
    sll_real64,             intent( in    ) :: mass_line((7*deg(1)**2-deg(1)-2)/2*(2*deg(2)+1)*(2*deg(3)+1)) !< boundary massline for the sparse mass matrix
    type(sll_t_matrix_csr), intent( inout ) :: matrix !< sparse matrix
    sll_int32,              intent( in    ) :: row  !< current row in the first direction 
    sll_int32,              intent( inout ) :: ind !< index
    sll_int32,              intent( in    ) :: shift !< shift in the second direction
    !local variables
    sll_int32 ::  shift1, column

    !use boundary mass for first deg(1) rows
    if(row <= deg(1)) then
       shift1 = ((row-1)*(2*deg(1)+row))/2 + shift
       do column = 1, deg(1)+row
          matrix%arr_a(1,ind) = mass_line(column+shift1)
          ind = ind+1
       end do
    else if(row > deg(1) .and. row < 2*deg(1)) then
       shift1 = (row-deg(1)-1)*(2*deg(1)+1)+(3*deg(1)**2+deg(1))/2 + shift
       do column = 1, 2*deg(1)+1
          matrix%arr_a(1,ind) = mass_line(column+shift1)
          ind = ind+1
       end do
       !use boundary mass for last deg(1) rows
    else if(row >= n_cells-deg(1)+2 .and. row <= n_cells) then
       shift1 = (n_cells-row)*(2*deg(1)+1)+(3*deg(1)**2+deg(1))/2 + shift
       do column = 2*deg(1)+1,1,-1
          matrix%arr_a(1,ind)=mass_line(column+shift1)
          ind = ind+1
       end do
    else if(row > n_cells .and. row <= n_cells+deg(1)) then
       shift1 = ((n_cells+deg(1)-row)*(2*deg(1)+n_cells+deg(1)+1-row))/2 + shift
       do column = deg(1)+1+n_cells+deg(1)-row,1,-1
          matrix%arr_a(1,ind) = mass_line(column+shift1)
          ind = ind+1
       end do
    else
       SLL_ERROR('assemble','error in row in assemble1_mass_clamped_boundary')
    end if

  end subroutine assemble1_clamped_boundary


  !---------------------------------------------------------------------------!
  !>  Assemble the given row of the mixed mass matrix
  subroutine assemble_mixedmass3d ( deg1, deg2, n_cells, mass_line, matrix, row, ind )
    sll_int32,              intent( in    ) :: deg1(3), deg2(3)                             !< spline deg for every direction
    sll_int32,              intent( in    ) :: n_cells(3)                                    !< number of rows (and grid points)
    sll_real64,             intent( in    ) :: mass_line((deg1(1)+deg2(1)+1)*(deg1(2)+ deg2(2)+1)*(deg1(3)+ deg2(3)+1)) !< massline for one row of the sparse mass matrix
    type(sll_t_matrix_csr), intent( inout ) :: matrix                                       !< sparse mass matrix
    sll_int32,              intent( in    ) :: row(3)                                      !< current row of the sparse matrix
    sll_int32,              intent( inout ) :: ind                                          !< index in the vector containing the entries of the sparse matrix
    !local variables
    sll_int32 ::  column, shift

    !assemble massmatrix in third direction
    ! For the first deg2(3) rows we need to put the first part to the back due to periodic boundaries
    if(row(3) <= deg2(3)) then
       do column = 2-row(3)+deg2(3), deg2(3)
          shift = (column-1)*(deg1(1)+deg2(1)+1)*(deg1(2)+deg2(2)+1)
          call massemble2( deg1, deg2, n_cells, mass_line, matrix, row, ind, shift )
       end do
       do column = deg2(3)+1,deg1(3)+deg2(3)+1
          shift = (column-1)*(deg1(1)+deg2(1)+1)*(deg1(2)+deg2(2)+1)
          call massemble2( deg1, deg2, n_cells, mass_line, matrix, row, ind, shift )
       end do

       do column = 1, deg2(3)-row(3)+1
          shift = (column-1)*(deg1(1)+deg2(1)+1)*(deg1(2)+deg2(2)+1)
          call massemble2( deg1, deg2, n_cells, mass_line, matrix, row, ind, shift )
       end do

       SLL_ASSERT( ind == (row(1)+(row(2)-1)*n_cells(1)+(row(3)-1)*n_cells(1)*n_cells(2))*(deg1(1)+deg2(1)+1)*(deg1(2)+deg2(2)+1)*(deg1(3)+deg2(3)+1)+1 )
    else if(row(3) > deg2(3) .and. row(3) <= n_cells(3)-deg1(3)) then
       do column = 1, deg1(3)+deg2(3)+1
          shift = (column-1)*(deg1(1)+deg2(1)+1)*(deg1(2)+deg2(2)+1)
          call massemble2( deg1, deg2, n_cells, mass_line, matrix, row, ind, shift )
       end do
       SLL_ASSERT( ind == (row(1)+(row(2)-1)*n_cells(1)+(row(3)-1)*n_cells(1)*n_cells(2))*(deg1(1)+deg2(1)+1)*(deg1(2)+deg2(2)+1)*(deg1(3)+deg2(3)+1)+1 )
       ! For the last deg(3) rows, we need to put the second part to the front due to periodic boundaries
    else if(row(3) > n_cells(3)-deg1(3) .and. row(3) <= n_cells(3)) then

       do column = n_cells(3)-row(3)+deg2(3)+2, deg1(3)+deg2(3)+1 
          shift = (column-1)*(deg1(1)+deg2(1)+1)*(deg1(2)+deg2(2)+1)
          call massemble2( deg1, deg2, n_cells, mass_line, matrix, row, ind, shift )
       end do

       do column = 1, deg2(3)+1
          shift = (column-1)*(deg1(1)+deg2(1)+1)*(deg1(2)+deg2(2)+1)
          call massemble2( deg1, deg2, n_cells, mass_line, matrix, row, ind, shift )
       end do

       do column = 2+deg2(3), n_cells(3)-row(3)+deg2(3)+1 
          shift = (column-1)*(deg1(1)+deg2(1)+1)*(deg1(2)+deg2(2)+1)
          call massemble2( deg1, deg2, n_cells, mass_line, matrix, row, ind, shift )
       end do
       SLL_ASSERT( ind == (row(1)+(row(2)-1)*n_cells(1)+(row(3)-1)*n_cells(1)*n_cells(2))*(deg1(1)+deg2(1)+1)*(deg1(2)+deg2(2)+1)*(deg1(3)+deg2(3)+1)+1 )
    else
       SLL_ERROR('assemble','error in row in assemble_mass')
    end if

  end subroutine assemble_mixedmass3d

  !assemble mass matrix in second direction
  subroutine massemble2( deg1, deg2, n_cells, mass_line, matrix, row, ind, shift )
    sll_int32,              intent( in    ) :: deg1(3), deg2(3) !< spline degrees for every direction
    sll_int32,              intent( in    ) :: n_cells(3) !< number of cells (and grid points) 
    sll_real64,             intent( in    ) :: mass_line((deg1(1)+deg2(1)+1)*(deg1(2)+ deg2(2)+1)*(deg1(3)+ deg2(3)+1)) !< massline for one row of the sparse mass matrix
    type(sll_t_matrix_csr), intent( inout ) :: matrix !< sparse matrix
    sll_int32,              intent( in    ) :: row(3) !< current row of the sparse matrix
    sll_int32,              intent( inout ) :: ind !< index 
    sll_int32,              intent( in    ) :: shift  !< shift in the third direction
    !local variables
    sll_int32 ::  column, shift1
    ! For the first deg2(2) rows we need to put the first part to the back due to periodic boundaries
    if(row(2) <= deg2(2)) then
       do column = 2-row(2)+deg2(2), deg2(2)
          shift1 = (column-1)*(deg1(1)+deg2(1)+1)+shift
          call massemble1( deg1, deg2, n_cells(1), mass_line, matrix, row(1), ind, shift1 )
       end do
       do column = deg2(2)+1,deg1(2)+deg2(2)+1
          shift1 = (column-1)*(deg1(1)+deg2(1)+1)+shift
          call massemble1( deg1, deg2, n_cells(1), mass_line, matrix, row(1), ind, shift1 )
       end do

       do column = 1, deg2(2)-row(2)+1
          shift1 = (column-1)*(deg1(1)+deg2(1)+1)+shift
          call massemble1( deg1, deg2, n_cells(1), mass_line, matrix, row(1), ind, shift1 )
       end do
    else if(row(2) > deg2(2) .and. row(2) <= n_cells(2)-deg1(2)) then
       do column = 1, deg1(2)+deg2(2)+1
          shift1 = (column-1)*(deg1(1)+deg2(1)+1)+shift
          call massemble1( deg1, deg2, n_cells(1), mass_line, matrix, row(1), ind, shift1 )
       end do

       ! For the last deg1(2) rows, we need to put the second part to the front due to periodic boundaries
    else if(row(2) > n_cells(2)-deg1(2) .and. row(2) <= n_cells(2)) then

       do column = n_cells(2)-row(2)+deg2(2)+2, deg1(2)+deg2(2)+1
          shift1 = (column-1)*(deg1(1)+deg2(1)+1)+shift
          call massemble1( deg1, deg2, n_cells(1), mass_line, matrix, row(1), ind, shift1 )
       end do

       do column = 1, deg2(2)+1
          shift1 = (column-1)*(deg1(1)+deg2(1)+1)+shift
          call massemble1( deg1, deg2, n_cells(1), mass_line, matrix, row(1), ind, shift1)
       end do

       do column = 2+deg2(2), n_cells(2)-row(2)+deg2(2)+1
          shift1 = (column-1)*(deg1(1)+deg2(1)+1)+shift
          call massemble1( deg1, deg2, n_cells(1), mass_line, matrix, row(1), ind, shift1 )
       end do
    else
       SLL_ERROR('assemble','error in row in assemble_mass')
    end if

  end subroutine massemble2

  !assemble mass matrix in first direction
  subroutine massemble1( deg1, deg2, n_cells, mass_line, matrix, row, ind, shift )
    sll_int32,              intent( in    ) :: deg1(3), deg2(3)  !< spline degrees for every direction
    sll_int32,              intent( in    ) :: n_cells !< number of cells in the first direction
    sll_real64,             intent( in    ) :: mass_line((deg1(1)+deg2(1)+1)*(deg1(2)+ deg2(2)+1)*(deg1(3)+ deg2(3)+1)) !< massline for one row of the sparse mass matrix
    type(sll_t_matrix_csr), intent( inout ) :: matrix !< sparse matrix
    sll_int32,              intent( in    ) :: row !< current row in the first direction
    sll_int32,              intent( inout ) :: ind !< index
    sll_int32,              intent( in    ) :: shift !< shift in the second direction
    !local variables
    sll_int32 ::  column

    ! For the first deg2(1) rows we need to put the first part to the back due to periodic boundaries
    if(row <= deg2(1)) then
       do column = 2-row+deg2(1), deg2(1)
          matrix%arr_a(1,ind) = mass_line(column+shift) 
          ind = ind+1
       end do
       matrix%arr_a(1,ind:ind+deg1(1)) = mass_line(deg2(1)+1+shift:deg1(1)+deg2(1)+1+shift)
       ind = ind+deg1(1)+1

       do column = 1, deg2(1)-row+1
          matrix%arr_a(1,ind) = mass_line(column+shift) 
          ind = ind+1
       end do
    else if(row > deg2(1) .and. row <= n_cells-deg1(1)) then
       matrix%arr_a(1,ind:ind+deg1(1)+deg2(1))=mass_line(1+shift:deg1(1)+deg2(1)+1+shift)

       ind = ind+deg1(1)+deg2(1)+1
       ! For the last deg1(1) rows, we need to put the second part to the front due to periodic boundaries
    else if(row > n_cells-deg1(1) .and. row <= n_cells) then

       do column = n_cells-row+deg2(1)+2, deg1(1)+deg2(1)+1
          matrix%arr_a(1,ind) = mass_line(column+shift)
          ind = ind+1
       end do

       matrix%arr_a(1,ind:ind+deg2(1)) = mass_line(1+shift:deg2(1)+1+shift)
       ind = ind+deg2(1)+1

       do column = 2+deg2(1), n_cells-row+deg2(1)+1
          matrix%arr_a(1,ind) = mass_line(column+shift)
          ind = ind+1
       end do
    else
       SLL_ERROR('assemble','error in row in assemble_mass')
    end if

  end subroutine massemble1

  !---------------------------------------------------------------------------!
  !>  Assemble the given row of the clamped mixed mass matrix
  subroutine assemble_mixedmass3d_clamped ( deg1, deg2, n_cells, mass_line, matrix, row, ind )
    sll_int32,              intent( in    ) :: deg1(3), deg2(3)                             !< spline degrees for every direction
    sll_int32,              intent( in    ) :: n_cells(3)                                    !< number of rows (and grid points)
    sll_real64,             intent( in    ) :: mass_line((deg1(1)+deg2(1)+1)*(deg1(2)+ deg2(2)+1)*(deg1(3)+ deg2(3)+1)) !< massline for one row of the sparse mass matrix
    type(sll_t_matrix_csr), intent( inout ) :: matrix                                       !< sparse mass matrix
    sll_int32,              intent( in    ) :: row(3)                                      !< current row of the sparse matrix
    sll_int32,              intent( inout ) :: ind                                          !< index in the vector containing the entries of the sparse matrix
    !local variables
    sll_int32 ::  column, shift, begin

    begin= ind

    !assemble massmatrix in third direction
    ! For the first deg2(3) rows we need to put the first part to the back due to periodic boundaries
    if(row(3)<=deg2(3)) then
       do column = 2-row(3)+deg2(3), deg2(3)
          shift=(column-1)*(deg1(1)+deg2(1)+1)*(deg1(2)+deg2(2)+1)
          call massemble2_clamped( deg1, deg2, n_cells, mass_line, matrix, row, ind, shift )
       end do
       do column= deg2(3)+1,deg1(3)+deg2(3)+1
          shift=(column-1)*(deg1(1)+deg2(1)+1)*(deg1(2)+deg2(2)+1)
          call massemble2_clamped( deg1, deg2, n_cells, mass_line, matrix, row, ind, shift )
       end do

       do column = 1, deg2(3)-row(3)+1
          shift=(column-1)*(deg1(1)+deg2(1)+1)*(deg1(2)+deg2(2)+1)
          call massemble2_clamped( deg1, deg2, n_cells, mass_line, matrix, row, ind, shift )
       end do
    else if(row(3) > deg2(3) .and. row(3)<= n_cells(3)-deg1(3)) then
       do column=1, deg1(3)+deg2(3)+1
          shift=(column-1)*(deg1(1)+deg2(1)+1)*(deg1(2)+deg2(2)+1)
          call massemble2_clamped( deg1, deg2, n_cells, mass_line, matrix, row, ind, shift )
       end do
       ! For the last deg(3) rows, we need to put the second part to the front due to periodic boundaries
    else if(row(3) > n_cells(3)-deg1(3) .and. row(3) <= n_cells(3)) then

       do column = n_cells(3)-row(3)+deg2(3)+2, deg1(3)+deg2(3)+1 
          shift=(column-1)*(deg1(1)+deg2(1)+1)*(deg1(2)+deg2(2)+1)
          call massemble2_clamped( deg1, deg2, n_cells, mass_line, matrix, row, ind, shift )
       end do

       do column=1, deg2(3)+1
          shift=(column-1)*(deg1(1)+deg2(1)+1)*(deg1(2)+deg2(2)+1)
          call massemble2_clamped( deg1, deg2, n_cells, mass_line, matrix, row, ind, shift )
       end do

       do column = 2+deg2(3), n_cells(3)-row(3)+deg2(3)+1 
          shift=(column-1)*(deg1(1)+deg2(1)+1)*(deg1(2)+deg2(2)+1)
          call massemble2_clamped( deg1, deg2, n_cells, mass_line, matrix, row, ind, shift )
       end do
    else
       SLL_ERROR('assemble','error in row in assemble_mass')
    end if
    SLL_ASSERT( ind == begin+(deg1(1)+deg2(1)+1)*(deg1(2)+deg2(2)+1)*(deg1(3)+deg2(3)+1) ) 

  end subroutine assemble_mixedmass3d_clamped

  !assemble mass matrix in second direction
  subroutine massemble2_clamped( deg1, deg2, n_cells, mass_line, matrix, row, ind, shift )
    sll_int32,              intent( in    ) :: deg1(3), deg2(3) !< spline degrees for every direction
    sll_int32,              intent( in    ) :: n_cells(3) !< number of cells (and grid points) 
    sll_real64,             intent( in    ) :: mass_line((deg1(1)+deg2(1)+1)*(deg1(2)+ deg2(2)+1)*(deg1(3)+ deg2(3)+1)) !< massline for one row of the sparse mass matrix
    type(sll_t_matrix_csr), intent( inout ) :: matrix !< sparse matrix
    sll_int32,              intent( in    ) :: row(3) !< current row of the sparse matrix
    sll_int32,              intent( inout ) :: ind !< index 
    sll_int32,              intent( in    ) :: shift  !< shift in the third direction
    !local variables
    sll_int32 ::  column, shift1

    ! For the first deg2(2) rows we need to put the first part to the back due to periodic boundaries
    if(row(2)<=deg2(2)) then
       do column = 2-row(2)+deg2(2), deg2(2)
          shift1=(column-1)*(deg1(1)+deg2(1)+1)+shift
          call massemble1_clamped( deg1, deg2, n_cells(1), mass_line, matrix, row(1), ind, shift1 )
       end do
       do column= deg2(2)+1,deg1(2)+deg2(2)+1
          shift1=(column-1)*(deg1(1)+deg2(1)+1)+shift
          call massemble1_clamped( deg1, deg2, n_cells(1), mass_line, matrix, row(1), ind, shift1 )
       end do

       do column = 1, deg2(2)-row(2)+1
          shift1=(column-1)*(deg1(1)+deg2(1)+1)+shift
          call massemble1_clamped( deg1, deg2, n_cells(1), mass_line, matrix, row(1), ind, shift1 )
       end do
    else if(row(2) > deg2(2) .and. row(2)<= n_cells(2)-deg1(2)) then
       do column=1, deg1(2)+deg2(2)+1
          shift1=(column-1)*(deg1(1)+deg2(1)+1)+shift
          call massemble1_clamped( deg1, deg2, n_cells(1), mass_line, matrix, row(1), ind, shift1 )
       end do
       ! For the last deg1(2) rows, we need to put the second part to the front due to periodic boundaries
    else if(row(2) > n_cells(2)-deg1(2) .and. row(2) <= n_cells(2)) then

       do column = n_cells(2)-row(2)+deg2(2)+2, deg1(2)+deg2(2)+1
          shift1=(column-1)*(deg1(1)+deg2(1)+1)+shift
          call massemble1_clamped( deg1, deg2, n_cells(1), mass_line, matrix, row(1), ind, shift1 )
       end do

       do column=1, deg2(2)+1
          shift1=(column-1)*(deg1(1)+deg2(1)+1)+shift
          call massemble1_clamped( deg1, deg2, n_cells(1), mass_line, matrix, row(1), ind, shift1 )
       end do

       do column = 2+deg2(2), n_cells(2)-row(2)+deg2(2)+1
          shift1=(column-1)*(deg1(1)+deg2(1)+1)+shift
          call massemble1_clamped( deg1, deg2, n_cells(1), mass_line, matrix, row(1), ind, shift1 )
       end do
    else
       SLL_ERROR('assemble','error in row in assemble_mass')
    end if

  end subroutine massemble2_clamped

  !assemble mass matrix in first direction
  subroutine massemble1_clamped( deg1, deg2, n_cells, mass_line, matrix, row, ind, shift )
    sll_int32,              intent( in    ) :: deg1(3), deg2(3)  !< spline degrees for every direction
    sll_int32,              intent( in    ) :: n_cells !< number of cells in the first direction
    sll_real64,             intent( in    ) :: mass_line((deg1(1)+deg2(1)+1)*(deg1(2)+ deg2(2)+1)*(deg1(3)+ deg2(3)+1)) !< massline for one row of the sparse mass matrix
    type(sll_t_matrix_csr), intent( inout ) :: matrix !< sparse matrix
    sll_int32,              intent( in    ) :: row !< current row in the first direction
    sll_int32,              intent( inout ) :: ind !< index
    sll_int32,              intent( in    ) :: shift !< shift in the second direction

    if(row > deg1(1)+deg2(1) .and. row <= n_cells-deg2(1)) then
       matrix%arr_a(1,ind:ind+deg1(1)+deg2(1))=mass_line(1+shift:deg1(1)+deg2(1)+1+shift)
       ind = ind+deg1(1)+deg2(1)+1
    else
       SLL_ERROR('assemble','error in row in assemble1_mixedmass_clamped')
    end if

  end subroutine massemble1_clamped

  !---------------------------------------------------------------------------!
  !>  Assemble the boundary part of the clamped mixed mass matrix
  subroutine assemble_mixedmass3d_clamped_boundary ( deg1, deg2, n_cells, mass_line, matrix, row, ind )
    sll_int32,              intent( in    ) :: deg1(3), deg2(3)                             !< spline degrees for every direction
    sll_int32,              intent( in    ) :: n_cells(3)                                    !< number of rows (and grid points)
    sll_real64,             intent( in    ) :: mass_line(((deg1(1)+deg2(1))**2+(2*deg2(1)+deg1(1)-deg1(1)**2)/2)*(deg1(2)+deg2(2)+1)*(deg1(3)+deg2(3)+1)) !< boundary massline for the sparse mass matrix
    type(sll_t_matrix_csr), intent( inout ) :: matrix                                       !< sparse mass matrix
    sll_int32,              intent( in    ) :: row(3)                                      !< current row of the sparse matrix
    sll_int32,              intent( inout ) :: ind                                          !< index in the vector containing the entries of the sparse matrix
    !local variables
    sll_int32 ::  column, shift

    !assemble massmatrix in third direction
    ! For the first deg2(3) rows we need to put the first part to the back due to periodic boundaries
    if(row(3)<=deg2(3)) then
       do column = 2-row(3)+deg2(3), deg2(3)
          shift=(column-1)*((deg1(1)+deg2(1))**2+(2*deg2(1)+deg1(1)-deg1(1)**2)/2)*(deg1(2)+deg2(2)+1)
          call massemble2_clamped_boundary( deg1, deg2, n_cells, mass_line, matrix, row, ind, shift )
       end do
       do column= deg2(3)+1,deg1(3)+deg2(3)+1
          shift=(column-1)*((deg1(1)+deg2(1))**2+(2*deg2(1)+deg1(1)-deg1(1)**2)/2)*(deg1(2)+deg2(2)+1)
          call massemble2_clamped_boundary( deg1, deg2, n_cells, mass_line, matrix, row, ind, shift )
       end do

       do column = 1, deg2(3)-row(3)+1
          shift=(column-1)*((deg1(1)+deg2(1))**2+(2*deg2(1)+deg1(1)-deg1(1)**2)/2)*(deg1(2)+deg2(2)+1)
          call massemble2_clamped_boundary( deg1, deg2, n_cells, mass_line, matrix, row, ind, shift )
       end do
    else if(row(3) > deg2(3) .and. row(3)<= n_cells(3)-deg1(3)) then
       do column=1, deg1(3)+deg2(3)+1
          shift=(column-1)*((deg1(1)+deg2(1))**2+(2*deg2(1)+deg1(1)-deg1(1)**2)/2)*(deg1(2)+deg2(2)+1)
          call massemble2_clamped_boundary( deg1, deg2, n_cells, mass_line, matrix, row, ind, shift )
       end do
       ! For the last deg(3) rows, we need to put the second part to the front due to periodic boundaries
    else if(row(3) > n_cells(3)-deg1(3) .and. row(3) <= n_cells(3)) then

       do column = n_cells(3)-row(3)+deg2(3)+2, deg1(3)+deg2(3)+1 
          shift=(column-1)*((deg1(1)+deg2(1))**2+(2*deg2(1)+deg1(1)-deg1(1)**2)/2)*(deg1(2)+deg2(2)+1)
          call massemble2_clamped_boundary( deg1, deg2, n_cells, mass_line, matrix, row, ind, shift )
       end do

       do column=1, deg2(3)+1
          shift=(column-1)*((deg1(1)+deg2(1))**2+(2*deg2(1)+deg1(1)-deg1(1)**2)/2)*(deg1(2)+deg2(2)+1)
          call massemble2_clamped_boundary( deg1, deg2, n_cells, mass_line, matrix, row, ind, shift )
       end do

       do column = 2+deg2(3), n_cells(3)-row(3)+deg2(3)+1 
          shift=(column-1)*((deg1(1)+deg2(1))**2+(2*deg2(1)+deg1(1)-deg1(1)**2)/2)*(deg1(2)+deg2(2)+1)
          call massemble2_clamped_boundary( deg1, deg2, n_cells, mass_line, matrix, row, ind, shift )
       end do
    else
       SLL_ERROR('assemble','error in row in assemble_mass')
    end if

  end subroutine assemble_mixedmass3d_clamped_boundary

  !assemble boundary part in second direction
  subroutine massemble2_clamped_boundary( deg1, deg2, n_cells, mass_line, matrix, row, ind, shift )
    sll_int32,              intent( in    ) :: deg1(3), deg2(3) !< spline degrees for every direction
    sll_int32,              intent( in    ) :: n_cells(3) !< number of cells (and grid points) 
    sll_real64,             intent( in    ) :: mass_line(((deg1(1)+deg2(1))**2+(2*deg2(1)+deg1(1)-deg1(1)**2)/2)*(deg1(2)+ deg2(2)+1)*(deg1(3)+ deg2(3)+1)) !< boundary massline for the sparse mass matrix
    type(sll_t_matrix_csr), intent( inout ) :: matrix !< sparse matrix
    sll_int32,              intent( in    ) :: row(3) !< current row of the sparse matrix
    sll_int32,              intent( inout ) :: ind  !< index 
    sll_int32,              intent( in    ) :: shift !< shift in the third direction
    !local variables
    sll_int32 ::  column, shift1

    ! For the first deg2(2) rows we need to put the first part to the back due to periodic boundaries
    if(row(2)<=deg2(2)) then
       do column = 2-row(2)+deg2(2), deg2(2)
          shift1=(column-1)*((deg1(1)+deg2(1))**2+(2*deg2(1)+deg1(1)-deg1(1)**2)/2)+shift
          call massemble1_clamped_boundary( deg1, deg2, n_cells(1), mass_line, matrix, row(1), ind, shift1 )
       end do
       do column= deg2(2)+1,deg1(2)+deg2(2)+1
          shift1=(column-1)*((deg1(1)+deg2(1))**2+(2*deg2(1)+deg1(1)-deg1(1)**2)/2)+shift
          call massemble1_clamped_boundary( deg1, deg2, n_cells(1), mass_line, matrix, row(1), ind, shift1 )
       end do

       do column = 1, deg2(2)-row(2)+1
          shift1=(column-1)*((deg1(1)+deg2(1))**2+(2*deg2(1)+deg1(1)-deg1(1)**2)/2)+shift
          call massemble1_clamped_boundary( deg1, deg2, n_cells(1), mass_line, matrix, row(1), ind, shift1 )
       end do
    else if(row(2) > deg2(2) .and. row(2)<= n_cells(2)-deg1(2)) then
       do column=1, deg1(2)+deg2(2)+1
          shift1=(column-1)*((deg1(1)+deg2(1))**2+(2*deg2(1)+deg1(1)-deg1(1)**2)/2)+shift
          call massemble1_clamped_boundary( deg1, deg2, n_cells(1), mass_line, matrix, row(1), ind, shift1 )
       end do
       ! For the last deg1(2) rows, we need to put the second part to the front due to periodic boundaries
    else if(row(2) > n_cells(2)-deg1(2) .and. row(2) <= n_cells(2)) then

       do column = n_cells(2)-row(2)+deg2(2)+2, deg1(2)+deg2(2)+1
          shift1=(column-1)*((deg1(1)+deg2(1))**2+(2*deg2(1)+deg1(1)-deg1(1)**2)/2)+shift
          call massemble1_clamped_boundary( deg1, deg2, n_cells(1), mass_line, matrix, row(1), ind, shift1 )
       end do

       do column=1, deg2(2)+1
          shift1=(column-1)*((deg1(1)+deg2(1))**2+(2*deg2(1)+deg1(1)-deg1(1)**2)/2)+shift
          call massemble1_clamped_boundary( deg1, deg2, n_cells(1), mass_line, matrix, row(1), ind, shift1)
       end do

       do column = 2+deg2(2), n_cells(2)-row(2)+deg2(2)+1
          shift1=(column-1)*((deg1(1)+deg2(1))**2+(2*deg2(1)+deg1(1)-deg1(1)**2)/2)+shift
          call massemble1_clamped_boundary( deg1, deg2, n_cells(1), mass_line, matrix, row(1), ind, shift1 )
       end do
    else
       SLL_ERROR('assemble','error in row in assemble_mass')
    end if

  end subroutine massemble2_clamped_boundary

  !assemble boundary part in first direction
  subroutine massemble1_clamped_boundary( deg1, deg2, n_cells, mass_line, matrix, row, ind, shift )
    sll_int32,              intent( in    ) :: deg1(3), deg2(3) !< spline degrees for every direction
    sll_int32,              intent( in    ) :: n_cells !< number of cells in the first direction
    sll_real64,             intent( in    ) :: mass_line(((deg1(1)+deg2(1))**2+(2*deg2(1)+deg1(1)-deg1(1)**2)/2)*(deg1(2)+ deg2(2)+1)*(deg1(3)+ deg2(3)+1)) !< boundary massline for the sparse mass matrix
    type(sll_t_matrix_csr), intent( inout ) :: matrix !< sparse matrix
    sll_int32,              intent( in    ) :: row  !< current row in the first direction
    sll_int32,              intent( inout ) :: ind !< index
    sll_int32,              intent( in    ) :: shift !< shift in the second direction
    !local variables
    sll_int32 ::  column, shift1

    !use boundary mass for first deg1(1) rows
    if(row<=deg1(1)) then
       shift1 = ((row-1)*(2*deg2(1)+row))/2 +shift
       do column = 1, deg2(1)+row
          matrix%arr_a(1,ind) = mass_line(column+shift1)
          ind = ind+1
       end do
    else if(row > deg1(1) .and. row <= deg1(1)+deg2(1)) then
       shift1 = (row-deg1(1)-1)*(deg1(1)+deg2(1)+1)+(deg1(1)*(2*deg2(1)+deg1(1)+1))/2 + shift
       do column = 1, deg1(1)+deg2(1)+1
          matrix%arr_a(1,ind) = mass_line(column+shift1)
          ind = ind+1
       end do

       !use boundary mass for last deg1(1) rows
    else if(row >= n_cells-deg2(1)+1 .and. row <= n_cells) then
       shift1 = (n_cells-row)*(deg1(1)+deg2(1)+1)+(deg1(1)*(2*deg2(1)+deg1(1)+1))/2 + shift
       do column = deg1(1)+deg2(1)+1,1,-1
          matrix%arr_a(1,ind)=mass_line(column+shift1)
          ind = ind+1
       end do
    else if(row > n_cells .and. row <= n_cells+deg1(1)) then
       shift1 = ((n_cells+deg1(1)-row)*(2*deg2(1)+n_cells+deg1(1)+1-row))/2 + shift
       do column = deg2(1)+1+n_cells+deg1(1)-row,1,-1
          matrix%arr_a(1,ind) = mass_line(column+shift1)
          ind = ind+1
       end do
    else
       SLL_ERROR('assemble','error in row in assemble_mass')
    end if

  end subroutine massemble1_clamped_boundary


end module sll_m_spline_fem_utilities_3d_helper
