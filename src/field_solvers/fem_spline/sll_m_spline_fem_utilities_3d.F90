!> @ingroup maxwell_solvers
!> @brief
!> Utilites for 3D Maxwell solvers with spline finite elements
!> @details
!> 
!> @author
!> Benedikt Perse

module sll_m_spline_fem_utilities_3d
  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#include "sll_assert.h"
#include "sll_errors.h"
#include "sll_memory.h"
#include "sll_working_precision.h"

  use sll_m_matrix_csr, only: &
       sll_t_matrix_csr

  use sll_m_low_level_bsplines, only: &
       sll_s_uniform_bsplines_eval_basis

  use sll_m_gauss_legendre_integration, only: &
       sll_f_gauss_legendre_points_and_weights

  use sll_m_spline_fem_utilities_3d_helper

  use sll_m_mapping_3d, only: &
       sll_t_mapping_3d

  implicit none

  public :: &
       sll_s_spline_fem_mass3d, &
       sll_s_spline_fem_mass_line, &
       sll_s_spline_fem_mixedmass3d, &
       sll_s_spline_fem_mixedmass_line, &
       sll_i_profile_function


  private
  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  abstract interface
     function sll_i_profile_function( x, component ) result(res)
       use sll_m_working_precision
       sll_real64, intent(in) :: x(3)
       sll_int32, optional, intent(in)  :: component(:)
       sll_real64             :: res
     end function sll_i_profile_function
  end interface

contains

  !> Set up 3d mass matrix for specific spline degrees and profile function
  subroutine sll_s_spline_fem_mass3d( n_cells, deg, component, matrix, profile, n_cells_min, n_cells_max)
    sll_int32,              intent( in    ) :: n_cells(3) !< number of cells (and grid points)
    sll_int32,              intent( in    ) :: deg(3)     !< highest spline degree
    sll_int32,              intent( in    ) :: component !< Specify the component 
    type(sll_t_matrix_csr), intent(   out ) :: matrix    !< sparse mass matrix
    procedure(sll_i_profile_function) :: profile !< profile function
    sll_int32, optional,   intent( in    ) :: n_cells_min(3) !< minimal cell number for mpi process
    sll_int32, optional,   intent( in    ) :: n_cells_max(3) !< maximal cell number for mpi process
    !local variables
    sll_int32  :: i1, i2, i3, ind, k, begin(3), limit(3)
    sll_int32  :: q(3)
    sll_real64 :: mass_line( 1:(2*deg(1)+1)*(2*deg(2)+1)*(2*deg(3)+1))
    sll_real64, allocatable :: xw_gauss_d1(:,:), xw_gauss_d2(:,:), xw_gauss_d3(:,:)
    sll_real64, allocatable :: bspl_d1(:,:), bspl_d2(:,:), bspl_d3(:,:)

    q = 2*deg+1!max(2*deg-1,1) !number of quadrature points

    allocate( xw_gauss_d1(1:2, 1:q(1)) )
    allocate( xw_gauss_d2(1:2, 1:q(2)) )
    allocate( xw_gauss_d3(1:2, 1:q(3)) )
    allocate( bspl_d1(1:deg(1)+1, 1:q(1)) )
    allocate( bspl_d2(1:deg(2)+1, 1:q(2)) )
    allocate( bspl_d3(1:deg(3)+1, 1:q(3)) )

    if( present(n_cells_min) .and. present(n_cells_max) )then
       ind = 1+((n_cells_min(2)-1)*n_cells(1)+(n_cells_min(3)-1)*n_cells(2)*n_cells(1))*&
            (2*deg(1)+1)*(2*deg(2)+1)*(2*deg(3)+1)
       begin = n_cells_min
       limit = n_cells_max
    else
       ind = 1
       begin = 1
       limit = n_cells
    end if

    call sll_s_spline_fem_sparsity_mass3d( deg(:), n_cells, matrix )
    bspl_d1 = 0._f64
    ! take enough Gauss points so that projection is exact for splines of deg deg
    ! rescale on [0,1] for compatibility with B-splines
    xw_gauss_d1 = sll_f_gauss_legendre_points_and_weights( q(1), 0._f64, 1._f64 )
    ! Compute bsplines at gauss_points
    do k = 1, q(1)
       call sll_s_uniform_bsplines_eval_basis( deg(1), xw_gauss_d1(1,k), bspl_d1(:, k) )
    end do
    bspl_d2 = 0._f64
    xw_gauss_d2 = sll_f_gauss_legendre_points_and_weights( q(2), 0._f64, 1._f64 )
    ! Compute bsplines at gauss_points
    do k=1, q(2)
       call sll_s_uniform_bsplines_eval_basis( deg(2), xw_gauss_d2(1,k), bspl_d2(:,k) )
    end do
    bspl_d3 = 0._f64
    xw_gauss_d3 = sll_f_gauss_legendre_points_and_weights( q(3), 0._f64, 1._f64 )
    ! Compute bsplines at gauss_points
    do k = 1, q(3)
       call sll_s_uniform_bsplines_eval_basis( deg(3), xw_gauss_d3(1,k), bspl_d3(:,k) )
    end do
    !loop over rows of the mass matrix to compute and assemble the massline
    do i3 = begin(3), limit(3)
       do i2 = begin(2), limit(2)
          do i1 = 1, n_cells(1)
             call sll_s_spline_fem_mass_line( q, deg, profile, [i1,i2,i3], n_cells, component, mass_line, xw_gauss_d1, xw_gauss_d2, xw_gauss_d3, bspl_d1, bspl_d2, bspl_d3  )
             !scale with 1/n_cells
             mass_line = mass_line/real(n_cells(1)*n_cells(2)*n_cells(3),f64)
             call assemble_mass3d( deg(:), n_cells, mass_line(:), matrix, [i1,i2,i3], ind )
          end do
       end do
    end do

    !SLL_ASSERT( ind == matrix%n_nnz+1 )

  end subroutine sll_s_spline_fem_mass3d

  !---------------------------------------------------------------------------!
  !> Computes the mass line for a mass matrix with \a degree splines
  subroutine sll_s_spline_fem_mass_line ( q, deg, profile, row, n_cells, component, mass_line,  xw_gauss_d1,  xw_gauss_d2,  xw_gauss_d3, bspl_d1, bspl_d2, bspl_d3  )
    sll_int32,              intent( in    ) :: q(3)   !< number of quadrature points 
    sll_int32,              intent( in    ) :: deg(3) !< spline degree in every direction
    procedure(sll_i_profile_function) ::  profile !< profile function
    sll_int32,              intent( in    ) :: row(3)                                           !< current row in the matrix
    sll_int32,              intent( in    ) :: n_cells(3)                                         !< gridpoints in every direction
    sll_int32,              intent( in    ) :: component                                         !< specify component of the form for which the massline is computed
    sll_real64,             intent(   out ) :: mass_line((2*deg(1)+1)*(2*deg(2)+1)*(2*deg(3)+1)) !< massline for one row of the sparse mass matrix
    sll_real64, intent ( in ) :: xw_gauss_d1(2, q(1)), xw_gauss_d2(2, q(2)), xw_gauss_d3(2, q(3)) !< quadrature gauss points
    sll_real64, intent(  in ) :: bspl_d1(deg(1)+1, q(1)), bspl_d2(deg(2)+1, q(2)), bspl_d3(deg(3)+1, q(3)) !< spline values
    ! local variables
    sll_int32  :: int1, int2, int3, j1, j2, j3, k1, k2, k3, l1, l2, l3
    sll_real64 :: c(3)
    sll_int32  :: limit(3), ind1(3), ind2(3), ind3(3), ind4(3)

    mass_line = 0._f64    
    !loop over massline entries in third dimension
    do j3 = 1, 2*deg(3)+1 
       if(j3 <= deg(3))then
          limit(3) = deg(3)
          ind1(3)  = deg(3)+1-j3 !position in massline
       else
          limit(3) = 2*deg(3)+1
          ind1(3)  = j3 !position in massline
       end if
       !loop over spline cells
       do int3 = j3, limit(3)
          if(j3 <= deg(3))then
             ind2(3) = deg(3)+1-int3+j3 !spline1
             ind3(3) = deg(3)+1-int3 !spline2
             ind4(3) = int3-j3 !left cellborder
          else
             ind2(3) = int3-j3+1 !spline1
             ind3(3) = int3-deg(3) !spline2
             ind4(3) = deg(3)+j3-int3 !left cellborder
          end if
          !loop over massline entries in second dimension
          do j2 = 1, 2*deg(2)+1
             if(j2 <= deg(2))then
                limit(2) = deg(2)
                ind1(2)  = deg(2)+1-j2
             else
                limit(2) = 2*deg(2)+1
                ind1(2)  = j2
             end if
             !loop over spline cells
             do int2 = j2, limit(2)
                if(j2 <= deg(2))then
                   ind2(2) = deg(2)+1-int2+j2 !spline1
                   ind3(2) = deg(2)+1-int2 !spline2
                   ind4(2) = int2-j2 !left cellborder
                else
                   ind2(2) = int2-j2+1 !spline1
                   ind3(2) = int2-deg(2) !spline2
                   ind4(2) = deg(2)+j2-int2 !left cellborder
                end if
                !loop over massline entries in first dimension
                do j1 = 1, 2*deg(1)+1
                   if(j1 <= deg(1))then
                      limit(1) = deg(1)
                      ind1(1)  = deg(1)+1-j1
                   else
                      limit(1) = 2*deg(1)+1
                      ind1(1)  = j1
                   end if
                   !loop over spline cells
                   do int1 = j1, limit(1)
                      if(j1 <= deg(1))then
                         ind2(1) = deg(1)+1-int1+j1 !spline1
                         ind3(1) = deg(1)+1-int1 !spline2
                         ind4(1) = int1-j1 !left cellborder
                      else
                         ind2(1) = int1-j1+1 !spline1
                         ind3(1) = int1-deg(1) !spline2
                         ind4(1) = deg(1)+j1-int1 !left cellborder
                      end if
                      ! loop over Gauss points
                      do l3 = 1, q(3)
                         k3 = 1 - (-1)**l3 * l3/2 + q(3)*(l3+1-2*((l3+1)/2) )!modulo(l3+1,2)
                         c(3) = (xw_gauss_d3(1,k3)+ real( row(3)-1+ind4(3)-((row(3)-1+ind4(3))/n_cells(3))*n_cells(3),f64 ))/real(n_cells(3),f64)
                         do l2 = 1, q(2)
                            k2 = 1 - (-1)**l2 * l2/2 + q(2)*(l2+1-2*((l2+1)/2) )!modulo(l2+1,2)
                            c(2) = (xw_gauss_d2(1,k2)+ real( row(2)-1+ind4(2)-((row(2)-1+ind4(2))/n_cells(2))*n_cells(2),f64 ))/real(n_cells(2),f64)
                            do l1 = 1, q(1)
                               k1 = 1 - (-1)**l1 * l1/2 + q(1)*(l1+1-2*((l1+1)/2) )!modulo(l1+1,2)
                               c(1) = (xw_gauss_d1(1,k1)+ real( row(1)-1+ind4(1)-((row(1)-1+ind4(1))/n_cells(1))*n_cells(1),f64 ))/real(n_cells(1),f64)
                               mass_line(ind1(1)+(ind1(2)-1)*(2*deg(1)+1)+(ind1(3)-1)*(2*deg(1)+1)*(2*deg(2)+1)) = &
                                    mass_line(ind1(1)+(ind1(2)-1)*(2*deg(1)+1)+(ind1(3)-1)*(2*deg(1)+1)*(2*deg(2)+1)) + &
                                    xw_gauss_d1(2,k1)* xw_gauss_d2(2,k2)* xw_gauss_d3(2,k3) *&
                                    bspl_d1(ind2(1), k1)*&
                                    bspl_d2(ind2(2), k2)*&
                                    bspl_d3(ind2(3), k3)*&
                                    bspl_d1(ind3(1), k1)*&
                                    bspl_d2(ind3(2), k2)*&
                                    bspl_d3(ind3(3), k3)*&
                                    profile( c, [component] )
                            enddo
                         enddo
                      end do
                   end do
                end do
             end do
          enddo
       end do
    end do

  end subroutine sll_s_spline_fem_mass_line


  !> Set up 3d mixed mass matrix for specific spline degree and profile function
  Subroutine sll_s_spline_fem_mixedmass3d( n_cells, deg1, deg2, component,  matrix, profile, n_cells_min, n_cells_max )
    sll_int32,              intent( in    ) :: n_cells(3)  !< number of cells (and grid points)
    sll_int32,              intent( in    ) :: deg1(3), deg2(3)      !< spline degrees
    sll_int32,              intent( in    ) :: component(2)  !< Specify the component 
    procedure(sll_i_profile_function) :: profile !< profile function
    type(sll_t_matrix_csr), intent(   out ) :: matrix     !< sparse mass matrix
    sll_int32, optional,   intent( in    ) :: n_cells_min(3) !< minimal cell number for mpi process
    sll_int32, optional,   intent( in    ) :: n_cells_max(3) !< maximal cell number for mpi process
    !local variables
    sll_int32  :: i1, i2, i3, ind, k, begin(3), limit(3), q(3)
    sll_real64 :: mass_line(1:(deg1(1)+deg2(1)+1)*(deg1(2)+deg2(2)+1)*(deg1(3)+deg2(3)+1))
    sll_real64, allocatable :: xw_gauss_d1(:,:), xw_gauss_d2(:,:), xw_gauss_d3(:,:)
    sll_real64, allocatable :: bspl_d1a(:,:), bspl_d2a(:,:), bspl_d3a(:,:)
    sll_real64, allocatable :: bspl_d1b(:,:), bspl_d2b(:,:), bspl_d3b(:,:)

    !quadrature points
    q = 2*max(deg1,deg2)+1!max(2*max(deg1,deg2), 1)
    
    if( present(n_cells_min) .and. present(n_cells_max) )then
       ind = 1+((n_cells_min(2)-1)*n_cells(1)+(n_cells_min(3)-1)*n_cells(2)*n_cells(1))*&
            (deg1(1)+deg2(1)+1)*(deg1(2)+deg2(2)+1)*(deg1(3)+deg2(3)+1)
       begin = n_cells_min
       limit = n_cells_max
    else
       ind = 1
       begin = 1
       limit = n_cells
    end if

    call sll_s_spline_fem_sparsity_mixedmass3d( deg1(:), deg2(:), n_cells, matrix ) 
    allocate( xw_gauss_d1(1:2, 1:q(1)) )
    allocate( xw_gauss_d2(1:2, 1:q(2)) )
    allocate( xw_gauss_d3(1:2, 1:q(3)) )
    allocate( bspl_d1a(1:deg1(1)+1, 1:q(1)) )
    allocate( bspl_d2a(1:deg1(2)+1, 1:q(2)) )
    allocate( bspl_d3a(1:deg1(3)+1, 1:q(3)) )
    allocate( bspl_d1b(1:deg2(1)+1, 1:q(1)) )
    allocate( bspl_d2b(1:deg2(2)+1, 1:q(2)) )
    allocate( bspl_d3b(1:deg2(3)+1, 1:q(3)) )

    bspl_d1a = 0._f64
    bspl_d1b = 0._f64
    ! take enough Gauss points so that projection is exact for splines of deg deg
    ! rescale on [0,1] for compatibility with B-splines
    xw_gauss_d1 = sll_f_gauss_legendre_points_and_weights( q(1), 0._f64, 1._f64 )
    ! Compute bsplines at gauss_points
    do k = 1, q(1)
       call sll_s_uniform_bsplines_eval_basis( deg1(1),xw_gauss_d1(1,k), bspl_d1a(1:deg1(1)+1, k) )
       call sll_s_uniform_bsplines_eval_basis( deg2(1),xw_gauss_d1(1,k), bspl_d1b(1:deg2(1)+1, k) )
    end do
    bspl_d2a=0._f64
    bspl_d2b=0._f64
    xw_gauss_d2 = sll_f_gauss_legendre_points_and_weights( q(2), 0._f64, 1._f64 )
    ! Compute bsplines at gauss_points
    do k=1, q(2)
       call sll_s_uniform_bsplines_eval_basis( deg1(2),xw_gauss_d2(1,k), bspl_d2a(1:deg1(2)+1, k) )
       call sll_s_uniform_bsplines_eval_basis( deg2(2),xw_gauss_d2(1,k), bspl_d2b(1:deg2(2)+1, k) )
    end do
    bspl_d3a = 0._f64
    bspl_d3b = 0._f64
    xw_gauss_d3 = sll_f_gauss_legendre_points_and_weights( q(3), 0._f64, 1._f64 )
    ! Compute bsplines at gauss_points
    do k = 1, q(3)
       call sll_s_uniform_bsplines_eval_basis( deg1(3),xw_gauss_d3(1,k), bspl_d3a(1:deg1(3)+1, k) )
       call sll_s_uniform_bsplines_eval_basis( deg2(3),xw_gauss_d3(1,k), bspl_d3b(1:deg2(3)+1, k) )
    end do
    !loop over rows of the mass matrix to compute and assemble the massline
    do i3 = begin(3), limit(3)
       do i2 = begin(2), limit(2)
          do i1 = 1, n_cells(1)
             call sll_s_spline_fem_mixedmass_line( q(:), deg1(:), deg2(:), profile, [i1,i2,i3], n_cells, component, mass_line(:),  xw_gauss_d1, xw_gauss_d2, xw_gauss_d3, bspl_d1a, bspl_d2a, bspl_d3a, bspl_d1b, bspl_d2b, bspl_d3b )
             !scale with delta_x=1/n_cells
             mass_line = mass_line/real(n_cells(1)*n_cells(2)*n_cells(3),f64)
             call assemble_mixedmass3d( deg1(:), deg2(:), n_cells, mass_line(:), matrix, [i1,i2,i3], ind )
          end do
       end do
    end do

    !    SLL_ASSERT( ind == matrix%n_nnz+1 )

    deallocate( xw_gauss_d1 )
    deallocate( xw_gauss_d2 )
    deallocate( xw_gauss_d3 )
    deallocate( bspl_d1a )
    deallocate( bspl_d2a )
    deallocate( bspl_d3a )
    deallocate( bspl_d1b )
    deallocate( bspl_d2b )
    deallocate( bspl_d3b )

  end subroutine sll_s_spline_fem_mixedmass3d


  ! Computes the mixed mass line for a mass matrix with \a degree splines
  subroutine sll_s_spline_fem_mixedmass_line( q, deg1, deg2, profile, row, n_cells, component, mass_line,  xw_gauss_d1,  xw_gauss_d2,  xw_gauss_d3, bspl_d1a, bspl_d2a, bspl_d3a, bspl_d1b, bspl_d2b, bspl_d3b )
    sll_int32,              intent( in    ) :: q(3)   !< number of quadrature points 
    sll_int32,              intent( in    ) :: deg1(3), deg2(3)                           !< spline degrees in every direction
    procedure(sll_i_profile_function) :: profile !< profile function
    sll_int32,              intent( in    ) :: row(3)                                           !< current row in the mass matrix
    sll_int32,              intent( in    ) :: n_cells(3)                                         !< gridpoints in every direction
    sll_int32,              intent( in    ) :: component(2)                                         !< specify component of the form for which the massline is computed
    sll_real64,             intent(   out ) :: mass_line((deg1(1)+deg2(1)+1)*(deg1(2)+deg2(2)+1)*(deg1(3)+deg2(3)+1)) !< massline for one row of the sparse mass matrix

    sll_real64, intent( in ) :: xw_gauss_d1(2, q(1)), xw_gauss_d2(2, q(2)), xw_gauss_d3(2, q(3))  !< quadrature gauss points
    sll_real64, intent( in ) :: bspl_d1a(deg1(1)+1, q(1)), bspl_d2a(deg1(2)+1, q(2)), bspl_d3a(deg1(3)+1, q(3)) !< spline values
    sll_real64, intent( in ) :: bspl_d1b(deg2(1)+1, q(1)), bspl_d2b(deg2(2)+1, q(2)), bspl_d3b(deg2(3)+1, q(3)) !< spline values
    ! local variables
    sll_int32  :: int1, int2, int3, j1, j2, j3, k1, k2, k3, l1, l2, l3
    sll_real64 :: c(3)
    sll_int32  :: limit(3), ind1(3), ind2(3), ind3(3), ind4(3)

    mass_line = 0._f64

    !loop over massline entries in third dimension
    do j3 = 1, deg1(3)+deg2(3)+1 
       if(j3 <= deg2(3))then
          limit(3) = deg2(3)
          ind1(3)  = deg2(3)+1-j3 !position in massline
       else
          limit(3) = deg1(3)+deg2(3)+1 
          ind1(3) = j3!position in massline
       end if
       ! loop over spline cells
       do int3 = j3, min(limit(3), j3+deg2(3)) 
          if(j3 <= deg2(3))then
             ind2(3) = deg1(3)+1-int3+j3 !spline1
             ind3(3) = deg2(3)+1-int3 !spline2
             ind4(3) = int3-j3 !left cellborder
          else
             ind2(3) = limit(3)-int3+1 !spline1
             ind3(3) = deg2(3)+1-int3+j3 !spline2
             ind4(3) = int3-deg2(3)-1 !left cellborder
          end if
          !loop over massline entries in second dimension
          do j2 = 1, deg1(2)+deg2(2)+1    
             if(j2 <= deg2(2))then
                limit(2) = deg2(2)
                ind1(2)  = deg2(2)+1-j2
             else
                limit(2) = deg1(2)+deg2(2)+1
                ind1(2)  = j2
             end if
             ! loop over spline cells
             do int2 = j2,  min(limit(2), j2+deg2(2))
                if(j2 <= deg2(2))then
                   ind2(2) = deg1(2)+1-int2+j2 !spline1
                   ind3(2) = deg2(2)+1-int2 !spline2
                   ind4(2) = int2-j2 !left cellborder
                else
                   ind2(2) = limit(2)-int2+1 !spline1
                   ind3(2) = deg2(2)+1-int2+j2 !spline2
                   ind4(2) = int2-deg2(2)-1 !left cellborder
                end if
                !loop over massline entries in first dimension
                do j1 = 1, deg1(1)+deg2(1)+1 
                   if(j1 <= deg2(1))then
                      limit(1) = deg2(1)
                      ind1(1)  = deg2(1)+1-j1
                   else
                      limit(1) = deg1(1)+deg2(1)+1
                      ind1(1)  = j1
                   end if
                   ! loop over spline cells
                   do int1 = j1, min(limit(1), j1+deg2(1))
                      if(j1 <= deg2(1))then
                         ind2(1) = deg1(1)+1-int1+j1 !spline1
                         ind3(1) = deg2(1)+1-int1 !spline2
                         ind4(1) = int1-j1 !left cellborder
                      else
                         ind2(1) = limit(1)-int1+1 !spline1
                         ind3(1) = deg2(1)+1-int1+j1 !spline2
                         ind4(1) = int1-deg2(1)-1 !left cellborder
                      end if
                      ! loop over Gauss points
                      do l3 = 1, q(3)
                         k3 = 1 - (-1)**l3 * l3/2 + q(3)*(l3+1-2*((l3+1)/2) )!modulo(l3+1,2)
                         c(3) = (xw_gauss_d3(1,k3)+ real( row(3)-1+ind4(3)-((row(3)-1+ind4(3))/n_cells(3))*n_cells(3),f64 ))/real(n_cells(3),f64)
                         do l2 = 1, q(2)
                            k2 = 1 - (-1)**l2 * l2/2 + q(2)*(l2+1-2*((l2+1)/2) )!modulo(l2+1,2)
                            c(2) = (xw_gauss_d2(1,k2)+ real( row(2)-1+ind4(2)-((row(2)-1+ind4(2))/n_cells(2))*n_cells(2),f64 ))/real(n_cells(2),f64)
                            do l1 = 1, q(1)
                               k1 = 1 - (-1)**l1 * l1/2 + q(1)*(l1+1-2*((l1+1)/2) )!modulo(l1+1,2)
                               c(1) = (xw_gauss_d1(1,k1)+ real( row(1)-1+ind4(1)-((row(1)-1+ind4(1))/n_cells(1))*n_cells(1),f64 ))/real(n_cells(1),f64)

                               mass_line(ind1(1)+(ind1(2)-1)*(deg1(1)+deg2(1)+1)+(ind1(3)-1)*(deg1(1)+deg2(1)+1)*(deg1(2)+deg2(2)+1)) = &
                                    mass_line(ind1(1)+(ind1(2)-1)*(deg1(1)+deg2(1)+1)+(ind1(3)-1)*(deg1(1)+deg2(1)+1)*(deg1(2)+deg2(2)+1)) + &
                                    xw_gauss_d1(2,k1)* xw_gauss_d2(2,k2)* xw_gauss_d3(2,k3) *&
                                    bspl_d1a(ind2(1), k1)*&
                                    bspl_d2a(ind2(2), k2)*&
                                    bspl_d3a(ind2(3), k3)*&
                                    bspl_d1b(ind3(1), k1)*&
                                    bspl_d2b(ind3(2), k2)*&
                                    bspl_d3b(ind3(3), k3)*&
                                    profile( c, component )
                            enddo
                         enddo
                      end do
                   end do
                end do
             end do
          enddo
       end do
    end do


  end subroutine sll_s_spline_fem_mixedmass_line


end module sll_m_spline_fem_utilities_3d
