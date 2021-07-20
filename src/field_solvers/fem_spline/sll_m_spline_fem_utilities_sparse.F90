!> @ingroup maxwell_solvers
!> @brief
!> Utilites for Maxwell solver's with spline finite elements using sparse matrix
!> linear solvers
!> @details
!> 
!> @author
!> Katharina Kormann

module sll_m_spline_fem_utilities_sparse
  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#include "sll_assert.h"
#include "sll_errors.h"
#include "sll_memory.h"
#include "sll_working_precision.h"

  use sll_m_matrix_csr, only: &
       sll_t_matrix_csr

  use sll_m_mapping_3d, only: &
       sll_t_mapping_3d

  use sll_m_splines_pp, only: &
       sll_t_spline_pp_1d

  use sll_m_spline_fem_utilities, only: &
       sll_i_profile_function_1d, &
       sll_s_spline_fem_mass_line_boundary_full, &
       sll_s_spline_fem_mixedmass_line_boundary, &
       sll_s_spline_fem_mass_line_full, &
       sll_s_spline_fem_mixedmass_line_full

  implicit none

  public :: &
       sll_s_spline_fem_mass1d, &
       sll_s_spline_fem_mass1d_clamped, &
       sll_s_spline_fem_mass1d_clamped_full, &
       sll_s_spline_fem_mass1d_full, &
       sll_s_spline_fem_sparsity_mass, &
       sll_s_spline_fem_mixedmass1d, &
       sll_s_spline_fem_mixedmass1d_full, &
       sll_s_spline_fem_mixedmass1d_clamped_full      

  private
  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

contains

  !---------------------------------------------------------------------------!
  !> Helper function to create sparsity pattern of the mass matrix
  subroutine sll_s_spline_fem_sparsity_mass( degree, n_cells, spmat )
    sll_int32,              intent( in    ) :: degree   !< spline degree
    sll_int32,              intent( in    ) :: n_cells  !< no. of grid cells
    type(sll_t_matrix_csr), intent(   out ) :: spmat    !< sparse matrix

    !local variables
    sll_int32 :: n_nnz
    sll_int32 :: n_nnz_per_row
    sll_int32 :: row, column, ind

    ! Compute number of non-zero elements
    n_nnz_per_row = (2*degree+1)
    n_nnz = n_nnz_per_row*n_cells

    if (n_cells<2*degree) then
       SLL_ERROR("sll_s_spline_fem_sparsity_mass", "The number of grid cells is smaller than twice the degree. This function is not written for this case")
    end if

    ! Create the csr matrix
    call spmat%create( n_rows=n_cells, n_cols=n_cells, n_nnz=n_nnz )

    ! Set the row and column vector for the elements (row vector compressed)
    ! There are n_nnz_per_row elements in each of the rows
    spmat%arr_ia(1) = 1
    do row = 2, n_cells+1
       spmat%arr_ia(row) = spmat%arr_ia(row-1) + n_nnz_per_row
    end do
    ! For row i the columns are i-deg, i-deg+1, ..., i, i+1, ... i+deg with periodic boundaries

    ind = 1
!!$    do row = 1, n_cells
!!$       do column = row-degree, row+degree
!!$          spmat%arr_ja(ind) = modulo(column-1, n_cells)+1
!!$          ind = ind+1
!!$       end do
!!$    end do

    do row = 1, degree
       do column = 1, row+degree
          spmat%arr_ja(ind) = column
          ind = ind+1
       end do
       do column = n_cells+row-degree, n_cells
          spmat%arr_ja(ind) = column
          ind = ind+1
       end do
    end do

    do row = degree+1, n_cells-degree
       do column = row-degree, row+degree
          spmat%arr_ja(ind) = column
          ind = ind+1
       end do
    end do

    do row = n_cells-degree+1, n_cells
       do column = 1, row+degree-n_cells
          spmat%arr_ja(ind) = column
          ind = ind+1
       end do
       do column = row-degree, n_cells
          spmat%arr_ja(ind) = column
          ind = ind+1
       end do
    end do

    spmat%arr_a = 0.0_f64

    SLL_ASSERT( ind == spmat%n_nnz+1 )


  end subroutine sll_s_spline_fem_sparsity_mass


  !> Assemble sparse matrix for mass matrix
  subroutine assemble_mass ( n_cells, degree, mass_line, matrix )
    sll_int32,  intent( in    ) :: n_cells  !< no. of grid cells
    sll_int32,  intent( in    ) :: degree   !< spline degree
    sll_real64, intent( in    ) :: mass_line(degree+1) !< mass line
    type(sll_t_matrix_csr), intent( inout ) :: matrix !< sparse matrix (mass matrix on output)

    !local variables
    sll_int32 :: row, column, ind

    ind = 1
    ! For the first degree rows we need to put the first part to the back due to periodic boundaries
    do row = 1, degree

       do column = -row+1, -1
          matrix%arr_a(1,ind) = mass_line(1-column) 
          ind = ind+1
       end do
       matrix%arr_a(1,ind:ind+degree) = mass_line
       ind = ind+degree+1

       do column = -degree, -row
          matrix%arr_a(1,ind) = mass_line(1-column) 
          ind = ind+1
       end do

    end do

    do row = degree+1, n_cells-degree

       do column = -degree, -1
          matrix%arr_a(1,ind) = mass_line(1-column) 
          ind = ind+1
       end do
       matrix%arr_a(1,ind:ind+degree) = mass_line
       ind = ind+degree+1
    end do

    ! For the last degree rows, we need to put the second part to the front due to periodic boundaries
    do row = n_cells-degree+1, n_cells

       do column = n_cells-row+1, degree
          matrix%arr_a(1,ind) = mass_line(column+1)
          ind = ind+1
       end do

       do column = -degree, -1
          matrix%arr_a(1,ind) = mass_line(1-column)
          ind = ind+1
       end do

       do column = 0, n_cells-row
          matrix%arr_a(1,ind) = mass_line(column+1)
          ind = ind+1
       end do
    end do

    SLL_ASSERT( ind == matrix%n_nnz+1)


  end subroutine assemble_mass


  !> Set up 1d mass matrix (sparsity pattern and values)
  subroutine sll_s_spline_fem_mass1d( n_cells, s_deg, mass_line, matrix )
    sll_int32, intent( in ) :: n_cells        !< no. of grid cells 
    sll_int32, intent( in ) :: s_deg          !< spline degree
    sll_real64, intent( in ) :: mass_line(:)  !< mass line
    type(sll_t_matrix_csr) :: matrix          !< sparse matrix (containing mass matrix on output)

    call sll_s_spline_fem_sparsity_mass( s_deg, n_cells, matrix )
    call assemble_mass( n_cells, s_deg, mass_line, matrix )
    ! write(24,*) matrix%arr_a(1,:)
  end subroutine sll_s_spline_fem_mass1d


  !>  Assemble the given row of a sparse matrix for a full mass line
  subroutine assemble_mass_full ( n_cells, degree, mass_line, matrix, row, ind )
    sll_int32,  intent( in    ) :: n_cells     !< no. of grid cells
    sll_int32,  intent( in    ) :: degree      !< spline degree
    sll_real64, intent( in    ) :: mass_line(2*degree+1)  !< mass line (without symmetry)
    type(sll_t_matrix_csr), intent( inout ) :: matrix !< sparse matrix 
    sll_int32, intent ( in    ) :: row               !< current row of the sparse matrix
    sll_int32, intent(inout) :: ind                  !< index in the vector containing the entries of the sparse matrix
    !local variables
    sll_int32 ::  column
    ! For the first degree rows we need to put the first part to the back due to periodic boundaries
    if(row<=degree) then
       do column = 2-row+degree, degree
          matrix%arr_a(1,ind) = mass_line(column) 
          ind = ind+1
       end do
       matrix%arr_a(1,ind:ind+degree) = mass_line(degree+1:2*degree+1)
       ind = ind+degree+1

       do column = 1, degree-row+1
          matrix%arr_a(1,ind) = mass_line(column) 
          ind = ind+1
       end do
    else if(row >= degree+1 .and. row<= n_cells-degree) then
       matrix%arr_a(1,ind:ind+2*degree)=mass_line

       ind = ind+2*degree+1
       ! For the last degree rows, we need to put the second part to the front due to periodic boundaries
    else if(row >= n_cells-degree+1 .and. row <= n_cells) then

       do column = n_cells-row+degree+2, 2*degree+1
          matrix%arr_a(1,ind) = mass_line(column)
          ind = ind+1
       end do

       matrix%arr_a(1,ind:ind+degree) = mass_line(1:degree+1)
       ind = ind+degree+1


       do column = 2+degree, n_cells-row+degree+1
          matrix%arr_a(1,ind) = mass_line(column)
          ind = ind+1
       end do
    else
       print*, 'error in row in assemble_mass'
    end if

    SLL_ASSERT( ind == row*(2*degree+1)+1)


  end subroutine assemble_mass_full


  !> Set up 1d mass matrix for specific spline degree and profile function
  Subroutine sll_s_spline_fem_mass1d_full( n_cells, s_deg, matrix, profile)
    sll_int32, intent( in ) :: n_cells !< no. of grid cells 
    sll_int32, intent( in ) :: s_deg !< spline degree
    type(sll_t_matrix_csr), intent(out) :: matrix !< sparse matrix (containing mass matrix on output)
    procedure(sll_i_profile_function_1d) :: profile !< profile function
    !local variables
    sll_int32 :: row, ind
    sll_real64 :: mass_line(2*s_deg+1)
    call sll_s_spline_fem_sparsity_mass( s_deg, n_cells, matrix )
    ind=1
    do row = 1, n_cells
       call sll_s_spline_fem_mass_line_full( s_deg, profile, mass_line, row, n_cells )
       !scale with delta_x=1/n_cells
       mass_line=mass_line/real(n_cells,f64)
       call assemble_mass_full( n_cells, s_deg, mass_line, matrix, row, ind )
    end do
  end subroutine sll_s_spline_fem_mass1d_full

  !---------------------------------------------------------------------------!
  !> Helper function to create sparsity pattern of the mass matrix with mixed degree
  !> \a deg1 and \a deg2
  subroutine sll_s_maxwell_spline_fem_sparsity_mixedmass( deg1, deg2, n_cells, spmat )
    sll_int32,              intent( in    ) :: deg1, deg2  !< degree test and basis functions
    sll_int32,              intent( in    ) :: n_cells     !< no of grid cells
    type(sll_t_matrix_csr), intent(   out ) :: spmat       !< sparse matrix

    !local variables
    sll_int32 :: n_nnz
    sll_int32 :: n_nnz_per_row
    sll_int32 :: row, column, ind

    ! Compute number of non-zero elements
    n_nnz_per_row = deg1+deg2+1
    n_nnz = n_nnz_per_row*n_cells

    ! Create the csr matrix
    call spmat%create( n_rows=n_cells, n_cols=n_cells, n_nnz=n_nnz )

    ! Set the row and column vector for the elements (row vector compressed)
    ! There are n_nnz_per_row elements in each of the rows
    spmat%arr_ia(1) = 1
    do row = 2, n_cells+1
       spmat%arr_ia(row) = spmat%arr_ia(row-1) + n_nnz_per_row
    end do
    ! For row i the columns are i-deg, i-deg+1, ..., i, i+1, ... i+deg with periodic boundaries

    ind = 1

    do row = 1, deg2
       do column = 1, row+deg1
          spmat%arr_ja(ind) = column
          ind = ind+1
       end do
       do column = n_cells+row-deg2, n_cells
          spmat%arr_ja(ind) = column
          ind = ind+1
       end do
    end do

    do row = deg2+1, n_cells-deg1
       do column = row-deg2, row+deg1
          spmat%arr_ja(ind) = column
          ind = ind+1
       end do
    end do

    do row = n_cells-deg1+1, n_cells
       do column = 1, row+deg1-n_cells
          spmat%arr_ja(ind) = column
          ind = ind+1
       end do
       do column = row-deg2, n_cells
          spmat%arr_ja(ind) = column
          ind = ind+1
       end do
    end do

    spmat%arr_a = 0.0_f64

    SLL_ASSERT( ind == spmat%n_nnz+1 )


  end subroutine sll_s_maxwell_spline_fem_sparsity_mixedmass


  !> Assemble the mixed mass matrix
  subroutine assemble_mixedmass ( n_cells, degree, mass, matrix )
    sll_int32,  intent( in    ) :: n_cells
    sll_int32,  intent( in    ) :: degree
    sll_real64, intent( in    ) :: mass(degree*2)
    type(sll_t_matrix_csr), intent( inout ) :: matrix

    !local variables
    sll_int32 :: row, ind

    ind = 1
    ! For the first degree rows we need to put the first part to the back due to periodic boundaries
    do row = 1, degree-1

       matrix%arr_a(1,ind:ind+degree+row-1) = mass(degree-row+1:2*degree)
       ind = ind + degree+row
       matrix%arr_a(1,ind:ind+degree-row-1 ) = mass(1:degree-row)
       ind = ind + degree-row
    end do

    do row = degree, n_cells-degree
       matrix%arr_a(1,ind:ind+degree*2-1) = mass
       ind = ind + degree*2
    end do

    ! For the last degree rows, we need to put the second part to the front due to periodic boundaries
    do row = n_cells-degree+1, n_cells
       matrix%arr_a(1,ind:ind+degree+row-n_cells-1 ) = mass(degree+n_cells-row+1:2*degree)
       ind = ind + degree-n_cells+row

       matrix%arr_a(1,ind:ind+degree+n_cells-row-1 ) = mass(1:degree+n_cells-row)
       ind = ind + degree+n_cells-row
    end do

    SLL_ASSERT( ind == matrix%n_nnz+1)


  end subroutine assemble_mixedmass


  !> Assemble the given row of a mixed mass matrix for a full mass_line
  subroutine assemble_mixedmass_full ( n_cells, deg1, deg2, mass_line, matrix, row, ind )
    sll_int32,  intent( in    ) :: n_cells !< no. of grid cells
    sll_int32,  intent( in    ) :: deg1, deg2 !< spline degrees
    sll_real64, intent( in    ) :: mass_line(deg1+deg2+1)  !< mass line 
    type(sll_t_matrix_csr), intent( inout ) :: matrix !< sparse matrix 
    sll_int32, intent ( in    ) :: row !< row of the sparse matrix
    sll_int32, intent(inout) :: ind !< index in the vector containing the entries of the sparse matrix
    !local variables
    sll_int32 ::  column

    ! For the first degree rows we need to put the first part to the back due to periodic boundaries

    if(row<=deg2) then
       do column = 2-row+deg2, deg2
          matrix%arr_a(1,ind) = mass_line(column) 
          ind = ind+1
       end do
       matrix%arr_a(1,ind:ind+deg1) = mass_line(deg2+1:deg1+deg2+1)
       ind = ind+deg1+1

       do column = 1, deg2-row+1
          matrix%arr_a(1,ind) = mass_line(column) 
          ind = ind+1
       end do

    else if( row > deg2 .and. row <= n_cells-deg1) then
       matrix%arr_a(1,ind:ind+deg1+deg2) = mass_line
       ind = ind + deg1+deg2+1

       ! For the last degree rows, we need to put the second part to the front due to periodic boundaries
    else if ( row > n_cells-deg1 .and. row <=n_cells) then
       do column = n_cells-row+deg2+2, deg1+deg2+1
          matrix%arr_a(1,ind) = mass_line(column)
          ind = ind+1
       end do

       matrix%arr_a(1,ind:ind+deg2) = mass_line(1:deg2+1)
       ind = ind+deg2+1

       do column = 2+deg2, n_cells-row+deg2+1
          matrix%arr_a(1,ind) = mass_line(column)
          ind = ind+1
       end do
    else
       print*, 'error in row in assemble_mixedmass'
    end if

    SLL_ASSERT( ind == row*(deg1+deg2+1)+1)


  end subroutine assemble_mixedmass_full


  subroutine sll_s_spline_fem_mixedmass1d( n_cells, s_deg, mass_line, matrix )
    sll_int32, intent( in ) :: n_cells !< no. of grid cells 
    sll_int32, intent( in ) :: s_deg !< spline degree
    sll_real64, intent( in ) :: mass_line(:) !< mass line
    type(sll_t_matrix_csr) :: matrix !< sparse matrix (containing mass matrix on output)


    call sll_s_maxwell_spline_fem_sparsity_mixedmass( s_deg, s_deg-1, n_cells, matrix )
    call assemble_mixedmass( n_cells, s_deg, mass_line, matrix )

  end subroutine sll_s_spline_fem_mixedmass1d


  subroutine sll_s_spline_fem_mixedmass1d_full( n_cells, deg1, deg2, matrix, profile )
    sll_int32, intent( in ) :: n_cells !< no. of grid cells 
    sll_int32, intent( in ) :: deg1, deg2 !< spline degrees
    type(sll_t_matrix_csr), intent(out):: matrix !< sparse matrix (containing mass matrix on output)
    procedure(sll_i_profile_function_1d) :: profile !< profile function
    !local variables
    sll_int32 :: row, ind, s_deg
    sll_real64 :: mass_line(deg1+deg2+1)
    sll_real64 :: delta_x
    delta_x=1._f64/real(n_cells,f64)
    s_deg = max(deg1,deg2)

    call  sll_s_maxwell_spline_fem_sparsity_mixedmass( deg1, deg2, n_cells, matrix )
    ind = 1
    do row = 1, n_cells
       call sll_s_spline_fem_mixedmass_line_full( deg1, deg2, profile, mass_line, row, n_cells)
       mass_line=mass_line*delta_x
       call assemble_mixedmass_full( n_cells, deg1, deg2, mass_line, matrix, row, ind )
    end do
  end subroutine sll_s_spline_fem_mixedmass1d_full


  !> Helper function to create sparsity pattern of the clamped mass matrix
  subroutine sll_s_spline_fem_sparsity_mass_clamped( degree, n_cells, spmat )
    sll_int32,              intent( in    ) :: degree !< spline degree
    sll_int32,              intent( in    ) :: n_cells !< no. of grid cells
    type(sll_t_matrix_csr), intent(   out ) :: spmat !< sparse matrix
    !local variables
    sll_int32 :: n_nnz
    sll_int32 :: n_nnz_per_row
    sll_int32 :: row, column, ind

    ! Compute number of non-zero elements
    n_nnz_per_row = (2*degree+1)
    n_nnz = n_nnz_per_row*(n_cells-degree)+(3*degree**2+degree)

    ! Create the csr matrix
    call spmat%create( n_rows=n_cells+degree, n_cols=n_cells+degree, n_nnz=n_nnz )

    ! Set the row and column vector for the elements (row vector compressed)
    ! There are n_nnz_per_row elements in each of the rows of the mainblock
    spmat%arr_ia(1) = 1
    do row=2, degree+1
       spmat%arr_ia(row) = spmat%arr_ia(row-1) + degree + row - 1
    end do
    do row = degree+2, n_cells+1
       spmat%arr_ia(row) = spmat%arr_ia(row-1) + n_nnz_per_row
    end do
    do row=n_cells+2, n_cells+degree+1
       spmat%arr_ia(row) = spmat%arr_ia(row-1) + 2*degree + n_cells+2-row
    end do
    SLL_ASSERT( spmat%arr_ia(n_cells+degree+1) == n_nnz+1 )


    ! For row i the columns are i-deg, i-deg+1, ..., i, i+1, ... i+deg with clamped boundaries
    ind = 1
    do row = 1, degree
       do column = 1, row+degree
          spmat%arr_ja(ind) = column
          ind = ind+1
       end do
    end do

    do row = degree+1, n_cells
       do column = row-degree, row+degree
          spmat%arr_ja(ind) = column
          ind = ind+1
       end do
    end do

    do row = n_cells+1, n_cells+degree
       do column = row-degree, n_cells+degree
          spmat%arr_ja(ind) = column
          ind = ind+1
       end do
    end do

    spmat%arr_a = 0.0_f64
    SLL_ASSERT( maxval(spmat%arr_ja) == spmat%n_cols )
    SLL_ASSERT( ind == spmat%n_nnz+1 )


  end subroutine sll_s_spline_fem_sparsity_mass_clamped


  !> Assemble the clamped mass matrix
  subroutine assemble_mass_clamped ( n_cells, degree, mass_line, mass_line_b, matrix )
    sll_int32,  intent( in    ) :: n_cells !< no. of grid cells
    sll_int32,  intent( in    ) :: degree !< spline degree
    sll_real64, intent( in    ) :: mass_line(degree+1) !< mass line
    sll_real64, intent( in    ) :: mass_line_b((degree+1)*degree) !< boundary mass line
    type(sll_t_matrix_csr), intent( inout ) :: matrix !< sparse matrix (mass matrix on output)

    !local variables
    sll_int32 :: row, column, ind, q

    q=degree+1

    ind = 1
    ! For the first degree rows we need to use the boundary mass_line
    do row = 1, degree
       do column = q-row, degree-1
          matrix%arr_a(1,ind) = mass_line_b(row+(row+column-q)*degree) 
          ind=ind+1
       end do
       matrix%arr_a(1,ind:ind+degree) = mass_line_b((row-1)*q+1:row*q)
       ind = ind+q
    end do

    !first part from boundary mass, rest from normal mass
    do row = degree+1, 2*degree-1
       do column = 1, degree+q-row
          matrix%arr_a(1,ind) = mass_line_b(row+(column-1+row-q)*degree) 
          ind=ind+1
       end do
       do column = 2*q-row, degree
          matrix%arr_a(1,ind) = mass_line(degree+2-column) 
          ind=ind+1
       end do
       matrix%arr_a(1,ind:ind+degree) = mass_line
       ind = ind+q
    end do

    do row = 2*degree, n_cells-degree+1
       do column = -degree, -1
          matrix%arr_a(1,ind) = mass_line(1-column) 
          ind = ind+1
       end do
       matrix%arr_a(1,ind:ind+degree) = mass_line
       ind = ind+q
    end do

    !first part from normal mass, rest from boundary mass
    do row = 1, degree-1
       do column = 1, q
          matrix%arr_a(1,ind) = mass_line(q+1-column)
          ind = ind+1
       end do
       do column = row, degree-2
          matrix%arr_a(1,ind) = mass_line(2+(column-row))
          ind=ind+1
       end do
       do column = 0, row
          matrix%arr_a(1,ind) = mass_line_b(q*degree-row-degree*column)
          ind=ind+1
       end do
    end do


    ! For the last degree rows, we need to use the boundary mass_line
    do row = 1, degree
       do column = 1, degree
          matrix%arr_a(1,ind) = mass_line_b((q-row)*q+1-column) 
          ind=ind+1
       end do
       do column = 1, q-row  
          matrix%arr_a(1,ind) = mass_line_b((degree-row)*q+1-(column-1)*degree) !:ind+degree-row
          ind = ind+1
       end do
    end do

    SLL_ASSERT( ind == matrix%n_nnz+1)

  end subroutine assemble_mass_clamped

  !> Set up 1d mass matrix (sparsity pattern and values)
  subroutine sll_s_spline_fem_mass1d_clamped( n_cells, s_deg, mass_line, mass_line_b, matrix )
    sll_int32, intent( in ) :: n_cells !< no. of grid cells 
    sll_int32, intent( in ) :: s_deg !< spline degree
    sll_real64, intent( in ) :: mass_line(:) !< mass line
    sll_real64, intent( in ) :: mass_line_b(:) !< boundary mass line
    type(sll_t_matrix_csr) :: matrix !< sparse matrix (containing mass matrix on output)

    SLL_ASSERT( n_cells >= 2*s_deg)

    call sll_s_spline_fem_sparsity_mass_clamped( s_deg, n_cells, matrix )
    call assemble_mass_clamped( n_cells, s_deg, mass_line, mass_line_b, matrix )

  end subroutine sll_s_spline_fem_mass1d_clamped

  !>  Assemble the given row of a sparse matrix for a full mass line
  subroutine assemble_mass_clamped_full ( n_cells, deg, mass, mass_b, matrix, row, ind )
    sll_int32,  intent( in    ) :: n_cells !< no. of grid cells
    sll_int32,  intent( in    ) :: deg !< spline degree
    sll_real64, intent( in    ) :: mass(2*deg+1) !< mass line (without symmetry)
    sll_real64, intent( in    ) :: mass_b((7*deg**2-deg-2)/2) !< boundary mass line 
    type(sll_t_matrix_csr), intent( inout ) :: matrix !< sparse matrix
    sll_int32, intent ( in    ) :: row  !< row of the sparse matrix
    sll_int32, intent(inout) :: ind !< index in the vector containing the entries of the sparse matrix

    !local variables
    sll_int32 :: column, shift1

    !use boundary mass for first deg rows
    if(row <= deg) then
       shift1 = ((row-1)*(2*deg+row))/2 
       do column = 1, deg+row
          matrix%arr_a(1,ind) = mass_b(column+shift1)
          ind = ind+1
       end do
    else if(row > deg .and. row < 2*deg) then
       shift1 = (row-deg-1)*(2*deg+1)+(3*deg**2+deg)/2 
       do column = 1, 2*deg+1
          matrix%arr_a(1,ind) = mass_b(column+shift1)
          ind = ind+1
       end do
    else if( row >= 2*deg .and. row<= n_cells-deg+1) then
       matrix%arr_a(1,ind:ind+2*deg) = mass
       ind = ind+2*deg+1
       !use boundary mass for last deg rows
    else if(row >= n_cells-deg+2 .and. row <= n_cells) then
       shift1 = (n_cells-row)*(2*deg+1)+(3*deg**2+deg)/2 
       do column = 2*deg+1,1,-1
          matrix%arr_a(1,ind)=mass_b(column+shift1)
          ind = ind+1
       end do
    else if(row > n_cells .and. row <= n_cells+deg) then
       shift1 = ((n_cells+deg-row)*(2*deg+n_cells+deg+1-row))/2
       do column = deg+1+n_cells+deg-row,1,-1
          matrix%arr_a(1,ind) = mass_b(column+shift1)
          ind = ind+1
       end do
    else
       print*, 'ERROR in assemble mass_full_clamped'
    end if

  end subroutine assemble_mass_clamped_full


  !> Set up 1d mass matrix for specific spline degree and profile function
  subroutine sll_s_spline_fem_mass1d_clamped_full( n_cells, deg, spline, matrix, profile)
    sll_int32, intent( in ) :: n_cells !< no. of grid cells 
    sll_int32, intent( in ) :: deg !< spline degree
    type(sll_t_spline_pp_1d), intent(in) :: spline !< pp-spline
    type(sll_t_matrix_csr), intent(out) :: matrix !< sparse matrix (containing mass matrix on output)
    procedure(sll_i_profile_function_1d) ::  profile !< profile function
    !local variables
    sll_int32 :: row, ind, begin
    sll_real64 :: mass_line(2*deg+1)
    sll_real64 :: mass_line_b((7*deg**2-deg-2)/2)

    call sll_s_spline_fem_sparsity_mass_clamped( deg, n_cells, matrix )

    call sll_s_spline_fem_mass_line_boundary_full ( deg, profile, spline, 1, n_cells, mass_line_b )
    mass_line_b = mass_line_b/real(n_cells,f64)

    ind=1
    begin=ind
    do row = 1, 2*deg-1
       call assemble_mass_clamped_full( n_cells, deg, mass_line, mass_line_b, matrix, row, ind )
    end do
    SLL_ASSERT( ind == begin+(7*deg**2-deg-2)/2 )
    do row = 2*deg, n_cells-deg+1
       call sll_s_spline_fem_mass_line_full( deg, profile, mass_line, row-deg, n_cells )
       !scale with delta_x=1/n_cells
       mass_line = mass_line/real(n_cells,f64)
       begin=ind
       call assemble_mass_clamped_full( n_cells, deg, mass_line, mass_line_b, matrix, row, ind )
       SLL_ASSERT( ind == begin+(2*deg+1) )
    end do

    call sll_s_spline_fem_mass_line_boundary_full ( deg, profile, spline, n_cells+1, n_cells, mass_line_b )
    mass_line_b = mass_line_b/real(n_cells,f64)
    begin=ind
    do row = n_cells-deg+2, n_cells+deg
       call assemble_mass_clamped_full( n_cells, deg, mass_line, mass_line_b, matrix, row, ind )
    end do
    SLL_ASSERT( ind == begin+(7*deg**2-deg-2)/2 )
    SLL_ASSERT( ind == matrix%n_nnz+1)
  end subroutine sll_s_spline_fem_mass1d_clamped_full


  !---------------------------------------------------------------------------!
  !> Helper function to create sparsity pattern of the mass matrix with mixed degree
  !> \a deg1 and \a deg2
  subroutine sll_s_maxwell_spline_fem_sparsity_mixedmass_clamped( deg1, deg2, n_cells, spmat )
    sll_int32,              intent( in    ) :: deg1, deg2  !< degree test and basis functions
    sll_int32,              intent( in    ) :: n_cells     !< no of grid cells
    type(sll_t_matrix_csr), intent(   out ) :: spmat       !< sparse matrix
    !local variables
    sll_int32 :: n_nnz
    sll_int32 :: n_nnz_per_row
    sll_int32 :: row, column, ind

    ! Compute number of non-zero elements
    n_nnz_per_row = deg1+deg2+1
    n_nnz = n_nnz_per_row*(n_cells-deg1)+deg1*(2*deg2+deg1+1)

    ! Create the csr matrix
    call spmat%create( n_rows=n_cells+deg1, n_cols=n_cells+deg2, n_nnz=n_nnz )

    ! Set the row and column vector for the elements (row vector compressed)
    spmat%arr_ia(1) = 1
    do row = 2, deg1+1
       spmat%arr_ia(row) = spmat%arr_ia(row-1) + deg2 + row - 1
    end do
    do row = deg1+2, n_cells+1
       spmat%arr_ia(row) = spmat%arr_ia(row-1) + n_nnz_per_row
    end do
    do row = n_cells+2, n_cells+deg1+1
       spmat%arr_ia(row) = spmat%arr_ia(row-1) + deg1+deg2 + n_cells+2-row
    end do
    SLL_ASSERT( spmat%arr_ia(n_cells+deg1+1) == n_nnz+1 )

    ind = 1
    do row = 1, deg1
       do column = 1, row+deg2
          spmat%arr_ja(ind) = column
          ind = ind+1
       end do
    end do
    do row = deg1+1, n_cells
       do column = row-deg1, row+deg2  
          spmat%arr_ja(ind) = column
          ind = ind+1
       end do
    end do
    do row = n_cells+1, n_cells +deg1
       do column = row-deg1, n_cells+deg2
          spmat%arr_ja(ind) = column
          ind = ind+1
       end do
    end do

    spmat%arr_a = 0.0_f64

    SLL_ASSERT( ind == spmat%n_nnz+1 )


  end subroutine sll_s_maxwell_spline_fem_sparsity_mixedmass_clamped


  !> Assemble the given row of a mixed mass matrix for a full mass_line
  subroutine assemble_mixedmass_clamped_full( n_cells, deg1, deg2, mass_line, matrix, row, ind )
    sll_int32,  intent( in    ) :: n_cells !< no. of grid cells
    sll_int32,  intent( in    ) :: deg1, deg2 !< spline degrees
    sll_real64, intent( in    ) :: mass_line(deg1+deg2+1)  !< mass line 
    type(sll_t_matrix_csr), intent( inout ) :: matrix !< sparse matrix 
    sll_int32, intent ( in    ) :: row !< row of the sparse matrix
    sll_int32, intent(inout) :: ind !< index in the vector containing the entries of the sparse matrix

    if(row > deg1+deg2 .and. row <= n_cells-deg2) then
       matrix%arr_a(1,ind:ind+deg1+deg2)=mass_line
       ind = ind+deg1+deg2+1
    else
       SLL_ERROR('assemble','error in row in assemble1_mixedmass_clamped')
    end if

  end subroutine assemble_mixedmass_clamped_full


  !> Assemble the boundary part of a mixed mass matrix
  subroutine assemble_mixedmass_clamped_boundary( n_cells, deg1, deg2, mass_line, matrix, row, ind )
    sll_int32,  intent( in    ) :: n_cells !< no. of grid cells
    sll_int32,  intent( in    ) :: deg1, deg2 !< spline degrees
    sll_real64, intent( in    ) :: mass_line( (deg1+deg2)**2+(2*deg2+deg1-deg1**2)/2 )  !< mass line 
    type(sll_t_matrix_csr), intent( inout ) :: matrix !< sparse matrix 
    sll_int32, intent ( in    ) :: row !< row of the sparse matrix
    sll_int32, intent(inout) :: ind !< index in the vector containing the entries of the sparse matrix
    !local variables
    sll_int32 ::  column

    !use boundary mass for first deg1 rows
    if(row<=deg1) then
       do column = 1, deg2+row
          matrix%arr_a(1,ind) = mass_line(column)
          ind = ind+1
       end do
    else if(row > deg1 .and. row <= deg1+deg2) then
       do column = 1, deg1+deg2+1
          matrix%arr_a(1,ind) = mass_line(column)
          ind = ind+1
       end do

       !use boundary mass for last deg1 rows
    else if(row >= n_cells-deg2+1 .and. row <= n_cells) then
       do column = deg1+deg2+1,1,-1
          matrix%arr_a(1,ind)=mass_line(column)
          ind = ind+1
       end do
    else if(row > n_cells .and. row <= n_cells+deg1) then
       do column = deg2+1+n_cells+deg1-row,1,-1
          matrix%arr_a(1,ind) = mass_line(column)
          ind = ind+1
       end do
    end if

  end subroutine assemble_mixedmass_clamped_boundary


  subroutine sll_s_spline_fem_mixedmass1d_clamped_full( n_cells, deg1, deg2, matrix, profile, spline1, spline2 )
    sll_int32, intent( in ) :: n_cells !< no. of grid cells 
    sll_int32, intent( in ) :: deg1, deg2 !< spline degrees
    type(sll_t_matrix_csr), intent(out):: matrix !< sparse matrix (containing mass matrix on output)
    procedure(sll_i_profile_function_1d) :: profile !< profile function
    type(sll_t_spline_pp_1d), intent( in    ) :: spline1 !< 1D pp-spline
    type(sll_t_spline_pp_1d), intent( in    ) :: spline2 !< 1D pp-spline
    !local variables
    sll_int32 :: row, ind
    sll_real64 :: mass_line(deg1+deg2+1)
    sll_real64 :: mass_line_b( (deg1+deg2)**2+(2*deg2+deg1-deg1**2)/2 )
    sll_real64 :: delta_x

    delta_x = 1._f64/real(n_cells,f64)

    call sll_s_maxwell_spline_fem_sparsity_mixedmass_clamped( deg1, deg2, n_cells, matrix )

    call sll_s_spline_fem_mixedmass_line_boundary( deg1, deg2, profile, spline1, spline2, 1, n_cells, mass_line_b )
    !scale with delta_x=1/n_cells
    mass_line_b = mass_line_b*delta_x

    ind = 1
    do row = 1, deg1+deg2
       call assemble_mixedmass_clamped_boundary( n_cells, deg1, deg2, mass_line_b, matrix, row, ind )
    end do
    do row = deg1+deg2+1, n_cells-deg2
       call sll_s_spline_fem_mixedmass_line_full( deg1, deg2, profile, mass_line, row-deg1, n_cells)
       !scale with delta_x=1/n_cells
       mass_line=mass_line*delta_x
       call assemble_mixedmass_clamped_full( n_cells, deg1, deg2, mass_line, matrix, row, ind )
    end do

    call sll_s_spline_fem_mixedmass_line_boundary ( deg1, deg2, profile, spline1, spline2, n_cells+1, n_cells, mass_line_b )
    !scale with delta_x=1/n_cells
    mass_line_b = mass_line_b*delta_x
    do row = n_cells-deg2+1, n_cells+deg1
       call assemble_mixedmass_clamped_boundary( n_cells, deg1, deg2, mass_line_b, matrix, row, ind )
    end do
    SLL_ASSERT( ind == matrix%n_nnz+1)
  end subroutine sll_s_spline_fem_mixedmass1d_clamped_full


end module sll_m_spline_fem_utilities_sparse
