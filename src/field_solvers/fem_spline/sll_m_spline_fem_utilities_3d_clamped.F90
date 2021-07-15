!> @ingroup maxwell_solvers
!> @brief
!> Utilites for Maxwell solver's with spline finite elements using sparse matrices
!> @details
!> 
!> @author
!> Benedikt Perse

module sll_m_spline_fem_utilities_3d_clamped
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

  use sll_m_mapping_3d, only: &
       sll_t_mapping_3d

  use sll_m_spline_fem_utilities_3d

  use sll_m_spline_fem_utilities_3d_helper

  use sll_m_splines_pp

  implicit none

  public :: &
       sll_s_spline_fem_mass3d_clamped, &
       sll_s_spline_fem_mixedmass3d_clamped, &
       sll_s_multiply_g_clamped, &
       sll_s_multiply_gt_clamped, &
       sll_s_multiply_c_clamped, &
       sll_s_multiply_ct_clamped, &
       sll_s_multiply_g_clamped_1d, &
       sll_s_multiply_gt_clamped_1d


  private
  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

contains
  
  !> Set up 3d clamped mass matrix for specific spline degrees and profile function
  subroutine sll_s_spline_fem_mass3d_clamped( n_cells, deg, component, matrix, profile, spline, n_cells_min, n_cells_max)
    sll_int32,                intent( in    ) :: n_cells(3) !< number of cells (and grid points)
    sll_int32,                intent( in    ) :: deg(3)     !< highest spline degree
    sll_int32,                intent( in    ) :: component !< Specify the component 
    type(sll_t_matrix_csr),   intent(   out ) :: matrix    !< sparse mass matrix
    procedure(sll_i_profile_function) :: profile !< profile function
    type(sll_t_spline_pp_1d), intent( in    ) :: spline !< 1D pp-spline
    sll_int32, optional,   intent( in    ) :: n_cells_min(3) !< minimal cell number for mpi process 
    sll_int32, optional,   intent( in    ) :: n_cells_max(3) !< maximal cell number for mpi process
    !local variables
    sll_int32  :: i1, i2, i3, ind, k, j, begin(3), limit(3)
    sll_int32  :: q(3), start
    sll_real64 :: mass_line_b( 1:(7*deg(1)**2-deg(1)-2)/2*(2*deg(2)+1)*(2*deg(3)+1))
    sll_real64 :: mass_line( 1:(2*deg(1)+1)*(2*deg(2)+1)*(2*deg(3)+1))
    sll_real64, allocatable :: xw_gauss_d1(:,:), xw_gauss_d2(:,:), xw_gauss_d3(:,:)
    sll_real64, allocatable :: bspl_d1(:,:), bspl_d2(:,:), bspl_d3(:,:)

    q = deg+1 !number of quadrature points
    allocate( xw_gauss_d1(1:2, 1:q(1)) )
    allocate( xw_gauss_d2(1:2, 1:q(2)) )
    allocate( xw_gauss_d3(1:2, 1:q(3)) )
    allocate( bspl_d1(1:deg(1)+1, 1:q(1)) )
    allocate( bspl_d2(1:deg(2)+1, 1:q(2)) )
    allocate( bspl_d3(1:deg(3)+1, 1:q(3)) )

    if( present(n_cells_min) .and. present(n_cells_max) )then
       ind=1+((n_cells_min(2)-1)+(n_cells_min(3)-1)*n_cells(2))*&
            ( ((n_cells(1)-deg(1))*(2*deg(1)+1)+(3*deg(1)**2+deg(1)))*&
            (2*deg(2)+1)*(2*deg(3)+1) )
       begin=n_cells_min
       limit=n_cells_max
    else
       ind=1
       begin=1
       limit=n_cells
    end if

    call sll_s_spline_fem_sparsity_mass3d_clamped( deg(:), n_cells, matrix )

    bspl_d1=0._f64
    ! take enough Gauss points so that projection is exact for splines of deg deg
    ! rescale on [0,1] for compatibility with B-splines
    xw_gauss_d1 = sll_f_gauss_legendre_points_and_weights( q(1), 0._f64, 1._f64 )
    ! Compute bsplines at gauss_points
    do k=1, q(1)
       call sll_s_uniform_bsplines_eval_basis( deg(1), xw_gauss_d1(1,k), bspl_d1(:, k) )
    end do
    bspl_d2=0._f64
    xw_gauss_d2 = sll_f_gauss_legendre_points_and_weights( q(2), 0._f64, 1._f64 )
    ! Compute bsplines at gauss_points
    do k=1, q(2)
       call sll_s_uniform_bsplines_eval_basis( deg(2), xw_gauss_d2(1,k), bspl_d2(:,k) )
    end do
    bspl_d3=0._f64
    xw_gauss_d3 = sll_f_gauss_legendre_points_and_weights( q(3), 0._f64, 1._f64 )
    ! Compute bsplines at gauss_points
    do k=1, q(3)
       call sll_s_uniform_bsplines_eval_basis( deg(3), xw_gauss_d3(1,k), bspl_d3(:,k) )
    end do

    !loop over rows of the mass matrix to compute and assemble the massline
    do i3=begin(3), limit(3)
       do i2=begin(2), limit(2)
          call sll_s_spline_fem_mass_line_boundary( q, deg, profile, [1,i2,i3], n_cells, component, mass_line_b, xw_gauss_d1, xw_gauss_d2, xw_gauss_d3, bspl_d1, bspl_d2, bspl_d3, spline )

          !scale with delta_x=1/n_cells
          mass_line_b=mass_line_b/real(n_cells(1)*n_cells(2)*n_cells(3),f64)

          start=ind
          do i1=1, 2*deg(1)-1
             call assemble_mass3d_clamped_boundary( deg, n_cells, mass_line_b, matrix, [i1,i2,i3], ind )
          end do
          SLL_ASSERT( ind == start+(7*deg(1)**2-deg(1)-2)/2*(2*deg(2)+1)*(2*deg(3)+1))
          do i1 = 2*deg(1), n_cells(1)-deg(1)+1
             call sll_s_spline_fem_mass_line( q, deg, profile, [i1-deg(1),i2,i3], n_cells, component, mass_line, xw_gauss_d1, xw_gauss_d2, xw_gauss_d3, bspl_d1, bspl_d2, bspl_d3 )
             !scale with delta_x=1/n_cells
             mass_line=mass_line/real(n_cells(1)*n_cells(2)*n_cells(3),f64)
             call assemble_mass3d_clamped( deg, n_cells, mass_line, matrix, [i1,i2,i3], ind )
          end do
          call sll_s_spline_fem_mass_line_boundary( q, deg, profile, [n_cells(1)+1,i2,i3], n_cells, component, mass_line_b, xw_gauss_d1, xw_gauss_d2, xw_gauss_d3, bspl_d1, bspl_d2, bspl_d3, spline )

          !scale with delta_x=1/n_cells
          mass_line_b=mass_line_b/real(n_cells(1)*n_cells(2)*n_cells(3),f64)

          start=ind
          do i1= n_cells(1)-deg(1)+2, n_cells(1)+deg(1)
             call assemble_mass3d_clamped_boundary( deg, n_cells, mass_line_b, matrix, [i1,i2,i3], ind )
          end do
          SLL_ASSERT( ind == start+(7*deg(1)**2-deg(1)-2)/2*(2*deg(2)+1)*(2*deg(3)+1))
       end do
    end do

  end subroutine sll_s_spline_fem_mass3d_clamped

  ! Computes the mass line for a mass matrix with \a degree splines
  subroutine sll_s_spline_fem_mass_line_boundary ( q, deg, profile, row, n_cells, component, mass_line,  xw_gauss_d1,  xw_gauss_d2,  xw_gauss_d3, bspl_d1, bspl_d2, bspl_d3, spline  )
    sll_int32,              intent( in    ) :: q(3)   !< number of quadrature points  
    sll_int32,              intent( in    ) :: deg(3) !< spline degree in every direction
    procedure(sll_i_profile_function) ::  profile !< profile function 
    sll_int32,              intent( in    ) :: row(3)        !< current row in the mass matrix
    sll_int32,              intent( in    ) :: n_cells(3)    !< gridpoints in every direction
    sll_int32,              intent( in    ) :: component     !< specify component of the form for which the massline is computed
    sll_real64,             intent(   out ) :: mass_line((7*deg(1)**2-deg(1)-2)/2*(2*deg(2)+1)*(2*deg(3)+1)) !< massline for one row of the sparse mass matrix

    sll_real64, intent ( in ) :: xw_gauss_d1(2, q(1)), xw_gauss_d2(2, q(2)), xw_gauss_d3(2, q(3)) !< quadrature gauss points
    sll_real64, intent(  in ) :: bspl_d1(deg(1)+1, q(1)), bspl_d2(deg(2)+1, q(2)), bspl_d3(deg(3)+1, q(3)) !< spline values
    type(sll_t_spline_pp_1d), intent( in    ) :: spline !< 1D pp-spline
    ! local variables
    sll_int32  :: int1, int2, int3, j1, j2, j3, k1, k2, k3, interval
    sll_real64 :: c(3)
    sll_int32  :: limit(3), ind1(3), ind2(3), ind3(3), ind4(3)
    sll_real64 :: bspl_b(1:deg(1)+1, 1:q(1),1:2*deg(1)-1)

    mass_line=0._f64

    do interval = 1, deg(1)-1
       do k1 = 1, q(1)
          do j1 = 1, deg(1)+1
             bspl_b(j1,k1,interval) =   sll_f_spline_pp_horner_1d( deg(1), spline%poly_coeffs_boundary_left(:,:,interval), xw_gauss_d1(1,k1), j1)
          end do
       end do
    end do
    do interval = deg(1), 2*deg(1)-1
       bspl_b(:,:,interval)=bspl_d1
    end do
    if( row(1) < 2*deg(1) )then
       ind4(1)= 1
    else if ( row(1) > n_cells(1) - deg(1)+1 .and. row(1) <= n_cells(1)+deg(1) )then
       ind4(1)= -1
    end if

    !loop over massline entries in third dimension
    do j3 = 1, 2*deg(3)+1 
       if(j3<=deg(3))then
          limit(3)=deg(3)
          ind1(3)=deg(3)+1-j3 !position in massline
       else
          limit(3)=2*deg(3)+1
          ind1(3)=j3 !position in massline
       end if
       !loop over spline cells
       do int3=j3, limit(3)
          if(j3<=deg(3))then
             ind2(3)=deg(3)+1-int3+j3 !spline1
             ind3(3)=deg(3)+1-int3 !spline2
             ind4(3)=int3-j3 !left cellborder
          else
             ind2(3)=int3-j3+1 !spline1
             ind3(3)=int3-deg(3) !spline2
             ind4(3)=deg(3)+j3-int3 !left cellborder
          end if
          !loop over massline entries in second dimension
          do j2 = 1, 2*deg(2)+1
             if(j2<=deg(2))then
                limit(2)=deg(2)
                ind1(2)=deg(2)+1-j2
             else
                limit(2)=2*deg(2)+1
                ind1(2)=j2
             end if
             !loop over spline cells
             do int2= j2, limit(2)
                if(j2<=deg(2))then
                   ind2(2)=deg(2)+1-int2+j2 !spline1
                   ind3(2)=deg(2)+1-int2 !spline2
                   ind4(2)=int2-j2 !left cellborder
                else
                   ind2(2)=int2-j2+1 !spline1
                   ind3(2)=int2-deg(2) !spline2
                   ind4(2)=deg(2)+j2-int2 !left cellborder
                end if
                !loop over massline entries in first dimension
                !interval
                do interval = 1, 2*deg(1)-1
                   !row
                   do j1 = interval, min(interval+deg(1),2*deg(1)-1)
                      !column
                      do int1 = interval, deg(1)+interval
                         if(j1 <= deg(1)+1)then
                            ind1(1) = int1+((j1-1)*(2*deg(1)+j1))/2 !position in first dimension in massline
                         else if(j1 > deg(1)+1 .and. j1 < 2*deg(1))then
                            ind1(1) = int1 + (j1-deg(1)-1)*2*deg(1)+(3*deg(1)**2+deg(1))/2 !position in first dimension in massline
                         end if
                         ind2(1) = j1-interval+1 !spline1
                         ind3(1) = int1-interval+1 !spline2

                         ! loop over Gauss points
                         do k3=1, q(3)
                            c(3)=(xw_gauss_d3(1,k3)+ real( row(3)-1+ind4(3)-((row(3)-1+ind4(3))/n_cells(3))*n_cells(3),f64 ))/real(n_cells(3),f64)
                            do k2=1, q(2)
                               c(2)=(xw_gauss_d2(1,k2)+ real( row(2)-1+ind4(2)-((row(2)-1+ind4(2))/n_cells(2))*n_cells(2),f64 ))/real(n_cells(2),f64)
                               do k1=1, q(1)
                                  c(1)=(real(ind4(1),f64)*xw_gauss_d1(1,k1)+real(row(1)-1+ind4(1)*(interval-1),f64))/real(n_cells(1),f64) 

                                  mass_line(ind1(1)+(ind1(2)-1)*(7*deg(1)**2-deg(1)-2)/2+(ind1(3)-1)*(7*deg(1)**2-deg(1)-2)/2*(2*deg(2)+1)) = &
                                       mass_line(ind1(1)+(ind1(2)-1)*(7*deg(1)**2-deg(1)-2)/2+(ind1(3)-1)*(7*deg(1)**2-deg(1)-2)/2*(2*deg(2)+1)) + &
                                       xw_gauss_d1(2,k1)* xw_gauss_d2(2,k2)* xw_gauss_d3(2,k3) *&
                                       bspl_b(ind2(1), k1,interval)*&
                                       bspl_d2(ind2(2), k2)*&
                                       bspl_d3(ind2(3), k3)*&
                                       bspl_b(ind3(1), k1,interval)*&
                                       bspl_d2(ind3(2), k2)*&
                                       bspl_d3(ind3(3), k3)*&
                                       profile( c, [component] )
                               enddo
                            enddo
                         enddo
                      end do
                   end do
                end do
             end do
          enddo
       end do
    end do

  end subroutine sll_s_spline_fem_mass_line_boundary

  !> Set up 3d clamped mixed mass matrix for specific spline degree and profile function
  Subroutine sll_s_spline_fem_mixedmass3d_clamped( n_cells, deg1, deg2, component, matrix, profile, spline1, spline2, n_cells_min, n_cells_max )
    sll_int32,                intent( in    ) :: n_cells(3)  !< number of cells (and grid points)
    sll_int32,                intent( in    ) :: deg1(3), deg2(3)      !< spline degrees
    sll_int32,                intent( in    ) :: component(2)  !< Specify the component 
    type(sll_t_matrix_csr),   intent(   out ) :: matrix     !< sparse mass matrix
    procedure(sll_i_profile_function) :: profile !< profile function
    type(sll_t_spline_pp_1d), intent( in    ) :: spline1 !< 1D pp-spline
    type(sll_t_spline_pp_1d), intent( in    ) :: spline2 !< 1D pp-spline
    sll_int32, optional,   intent( in    ) :: n_cells_min(3) !< minimal cell number for mpi process
    sll_int32, optional,   intent( in    ) :: n_cells_max(3) !< maximal cell number for mpi process
    !local variables
    sll_int32 :: q(3),  begin(3), limit(3)
    sll_int32  :: i1, i2, i3, ind, k, j, start
    sll_real64 :: mass_line_b( 1:((deg1(1)+deg2(1))**2+(2*deg2(1)+deg1(1)-deg1(1)**2)/2)*(deg1(2)+deg2(2)+1)*(deg1(3)+deg2(3)+1) ) 
    sll_real64 :: mass_line( 1:(deg1(1)+deg2(1)+1)*(deg1(2)+deg2(2)+1)*(deg1(3)+deg2(3)+1) )
    sll_real64, allocatable :: xw_gauss_d1(:,:), xw_gauss_d2(:,:), xw_gauss_d3(:,:)
    sll_real64, allocatable :: bspl_d1a(:,:), bspl_d2a(:,:), bspl_d3a(:,:)
    sll_real64, allocatable :: bspl_d1b(:,:), bspl_d2b(:,:), bspl_d3b(:,:)

    !quadrature points
    q = max(deg1,deg2)+1

    allocate( xw_gauss_d1(1:2, 1:q(1)) )
    allocate( xw_gauss_d2(1:2, 1:q(2)) )
    allocate( xw_gauss_d3(1:2, 1:q(3)) )
    allocate( bspl_d1a(1:deg1(1)+1, 1:q(1)) )
    allocate( bspl_d2a(1:deg1(2)+1, 1:q(2)) )
    allocate( bspl_d3a(1:deg1(3)+1, 1:q(3)) )
    allocate( bspl_d1b(1:deg2(1)+1, 1:q(1)) )
    allocate( bspl_d2b(1:deg2(2)+1, 1:q(2)) )
    allocate( bspl_d3b(1:deg2(3)+1, 1:q(3)) )

    if( present(n_cells_min) .and. present(n_cells_max) )then
       ind=1+( (n_cells_min(2)-1)+(n_cells_min(3)-1)*n_cells(2) )*&
            ( ((n_cells(1)-deg1(1))*(deg1(1)+deg2(1)+1)+deg1(1)*(2*deg2(1)+deg1(1)+1))*&
            (deg1(2)+deg2(2)+1)*(deg1(3)+deg2(3)+1) )
       begin=n_cells_min
       limit=n_cells_max
    else
       ind=1
       begin=1
       limit=n_cells
    end if

    call sll_s_spline_fem_sparsity_mixedmass3d_clamped( deg1(:), deg2(:), n_cells, matrix )
    bspl_d1a=0._f64
    bspl_d1b=0._f64
    ! take enough Gauss points so that projection is exact for splines of deg deg
    ! rescale on [0,1] for compatibility with B-splines
    xw_gauss_d1 = sll_f_gauss_legendre_points_and_weights( q(1), 0._f64, 1._f64 )
    ! Compute bsplines at gauss_points
    do k=1, q(1)
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
    bspl_d3a=0._f64
    bspl_d3b=0._f64
    xw_gauss_d3 = sll_f_gauss_legendre_points_and_weights( q(3), 0._f64, 1._f64 )
    ! Compute bsplines at gauss_points
    do k=1, q(3)
       call sll_s_uniform_bsplines_eval_basis( deg1(3),xw_gauss_d3(1,k), bspl_d3a(1:deg1(3)+1, k) )
       call sll_s_uniform_bsplines_eval_basis( deg2(3),xw_gauss_d3(1,k), bspl_d3b(1:deg2(3)+1, k) )
    end do


    !loop over rows of the mass matrix to compute assemble the massline
    do i3=begin(3), limit(3)
       do i2=begin(2), limit(2)
          call sll_s_spline_fem_mixedmass_line_boundary( q, deg1, deg2, profile, [1,i2,i3], n_cells, component, mass_line_b,  xw_gauss_d1, xw_gauss_d2, xw_gauss_d3, bspl_d1a, bspl_d2a, bspl_d3a, bspl_d1b, bspl_d2b, bspl_d3b, spline1, spline2 )

          !scale with delta_x=1/n_cells
          mass_line_b=mass_line_b/real(n_cells(1)*n_cells(2)*n_cells(3),f64)
          start=ind
          do i1=1, deg1(1)+deg2(1)
             call assemble_mixedmass3d_clamped_boundary( deg1(:), deg2(:), n_cells, mass_line_b, matrix, [i1,i2,i3], ind )
          end do
          SLL_ASSERT( ind == start+((deg1(1)+deg2(1))**2+(2*deg2(1)+deg1(1)-deg1(1)**2)/2)*(deg1(2)+deg2(2)+1)*(deg1(3)+deg2(3)+1) )
          do i1=deg1(1)+deg2(1)+1, n_cells(1)-deg2(1)
             call sll_s_spline_fem_mixedmass_line( q, deg1, deg2, profile, [i1-deg1(1),i2,i3], n_cells, component, mass_line,  xw_gauss_d1, xw_gauss_d2, xw_gauss_d3, bspl_d1a, bspl_d2a, bspl_d3a, bspl_d1b, bspl_d2b, bspl_d3b )
             !scale with delta_x=1/n_cells
             mass_line=mass_line/real(n_cells(1)*n_cells(2)*n_cells(3),f64)
             call assemble_mixedmass3d_clamped( deg1(:), deg2(:), n_cells, mass_line(:), matrix, [i1,i2,i3], ind )
          end do
          call sll_s_spline_fem_mixedmass_line_boundary( q, deg1, deg2, profile, [n_cells(1)+1,i2,i3], n_cells, component, mass_line_b,  xw_gauss_d1, xw_gauss_d2, xw_gauss_d3, bspl_d1a, bspl_d2a, bspl_d3a, bspl_d1b, bspl_d2b, bspl_d3b, spline1, spline2 )

          !scale with delta_x=1/n_cells
          mass_line_b=mass_line_b/real(n_cells(1)*n_cells(2)*n_cells(3),f64)
          start=ind
          do i1 = n_cells(1)-deg2(1)+1, n_cells(1)+deg1(1)
             call assemble_mixedmass3d_clamped_boundary( deg1, deg2, n_cells, mass_line_b, matrix, [i1,i2,i3], ind )
          end do
          SLL_ASSERT( ind == start+((deg1(1)+deg2(1))**2+(2*deg2(1)+deg1(1)-deg1(1)**2)/2)*(deg1(2)+deg2(2)+1)*(deg1(3)+deg2(3)+1) )
       end do
    end do

    deallocate( xw_gauss_d1 )
    deallocate( xw_gauss_d2 )
    deallocate( xw_gauss_d3 )
    deallocate( bspl_d1a )
    deallocate( bspl_d2a )
    deallocate( bspl_d3a )
    deallocate( bspl_d1b )
    deallocate( bspl_d2b )
    deallocate( bspl_d3b )

  end subroutine sll_s_spline_fem_mixedmass3d_clamped


   !> Computes the mixed mass line for a mass matrix with \a degree splines
  subroutine sll_s_spline_fem_mixedmass_line_boundary( q, deg1, deg2, profile, row, n_cells, component, mass_line,  xw_gauss_d1,  xw_gauss_d2,  xw_gauss_d3, bspl_d1a, bspl_d2a, bspl_d3a, bspl_d1b, bspl_d2b, bspl_d3b, spline1, spline2 )
    sll_int32,              intent( in    ) :: q(3)   !< number of quadrature points 
    sll_int32,              intent( in    ) :: deg1(3), deg2(3)      !< spline degrees in every direction
    procedure(sll_i_profile_function) ::  profile !< profile function
    sll_int32,              intent( in    ) :: row(3)                !< current row in the mass matrix
    sll_int32,              intent( in    ) :: n_cells(3)            !< gridpoints in every direction
    sll_int32,              intent( in    ) :: component(2)          !< specify component of the form for which the massline is computed
    sll_real64,             intent(   out ) :: mass_line( 1:((deg1(1)+deg2(1))**2+(2*deg2(1)+deg1(1)-deg1(1)**2)/2)*(deg1(2)+deg2(2)+1)*(deg1(3)+deg2(3)+1) ) !< massline for one row of the sparse mass matrix
    sll_real64, intent( in ) :: xw_gauss_d1(2, q(1)), xw_gauss_d2(2, q(2)), xw_gauss_d3(2, q(3)) !< quadrature gauss points
    sll_real64, intent( in ) :: bspl_d1a(deg1(1)+1, q(1)), bspl_d2a(deg1(2)+1, q(2)), bspl_d3a(deg1(3)+1, q(3)) !< spline values
    sll_real64, intent( in ) :: bspl_d1b(deg2(1)+1, q(1)), bspl_d2b(deg2(2)+1, q(2)), bspl_d3b(deg2(3)+1, q(3)) !< spline values
    type(sll_t_spline_pp_1d), intent( in    ) :: spline1 !< 1D pp-spline
    type(sll_t_spline_pp_1d), intent( in    ) :: spline2 !< 1D pp-spline
    ! local variables
    sll_int32  :: int1, int2, int3, j1, j2, j3, k1, k2, k3, interval
    sll_real64 :: c(3)
    sll_int32  :: limit(3), ind1(3), ind2(3), ind3(3), ind4(3)
    sll_real64 :: bspl_b1(1:deg1(1)+1, 1:q(1),1:deg1(1)+deg2(1))
    sll_real64 :: bspl_b2(1:deg2(1)+1, 1:q(1),1:deg1(1)+deg2(1))

    mass_line=0._f64

    do interval = 1, deg1(1)-1
       do k1 = 1, q(1)
          do j1 = 1, deg1(1)+1
             bspl_b1(j1,k1,interval) =   sll_f_spline_pp_horner_1d( deg1(1), spline1%poly_coeffs_boundary_left(:,:,interval), xw_gauss_d1(1,k1), j1)
          end do
       end do
    end do
    if( deg1(1) == 0 ) then
       do interval = 1, deg2(1)
          bspl_b1(:,:,interval)=bspl_d1a
       end do
    else
       do interval = deg1(1), deg1(1)+deg2(1)
          bspl_b1(:,:,interval)=bspl_d1a
       end do
    end if

    do interval = 1, deg2(1)-1
       do k1 = 1, q(1)
          do j1 = 1, deg2(1)+1
             bspl_b2(j1,k1,interval) =   sll_f_spline_pp_horner_1d( deg2(1), spline2%poly_coeffs_boundary_left(:,:,interval), xw_gauss_d1(1,k1), j1)
          end do
       end do
    end do
    if( deg2(1) == 0 ) then
       do interval = 1, deg1(1)
          bspl_b2(:,:,interval)=bspl_d1b
       end do
    else
       do interval = deg2(1), deg1(1)+deg2(1)
          bspl_b2(:,:,interval)=bspl_d1b
       end do
    end if

    if( row(1) <= deg1(1)+deg2(1) )then
       ind4(1)= 1
    else if ( row(1) > n_cells(1) - deg2(1) .and. row(1) <= n_cells(1)+deg1(1) )then
       ind4(1)= -1
    end if


    !loop over massline entries in third dimension
    do j3 = 1, deg1(3)+deg2(3)+1 
       if(j3<=deg2(3))then
          limit(3)=deg2(3)
          ind1(3)=deg2(3)+1-j3 !position in massline
       else
          limit(3)=deg1(3)+deg2(3)+1 
          ind1(3)=j3!position in massline
       end if
       ! loop over spline cells
       do int3=j3,  min(limit(3),j3+deg2(3)) 
          if(j3<=deg2(3))then
             ind2(3)=deg1(3)+1-int3+j3 !spline1
             ind3(3)=deg2(3)+1-int3 !spline2
             ind4(3)=int3-j3 !left cellborder
          else
             ind2(3)=limit(3)-int3+1 !spline1
             ind3(3)=deg2(3)+1-int3+j3 !spline2
             ind4(3)=int3-deg2(3)-1 !left cellborder
          end if
          !loop over massline entries in second dimension
          do j2 = 1, deg1(2)+deg2(2)+1    
             if(j2<=deg2(2))then
                limit(2)=deg2(2)
                ind1(2)=deg2(2)+1-j2
             else
                limit(2)=deg1(2)+deg2(2)+1
                ind1(2)=j2
             end if
             ! loop over spline cells
             do int2=j2,  min(limit(2),j2+deg2(2))
                if(j2<=deg2(2))then
                   ind2(2)=deg1(2)+1-int2+j2 !spline1
                   ind3(2)=deg2(2)+1-int2 !spline2
                   ind4(2)=int2-j2 !left cellborder
                else
                   ind2(2)=limit(2)-int2+1 !spline1
                   ind3(2)=deg2(2)+1-int2+j2 !spline2
                   ind4(2)=int2-deg2(2)-1 !left cellborder
                end if
                !loop over massline entries in first dimension
                !interval
                do interval = 1, deg1(1)+deg2(1)
                   !row
                   do j1 = interval, min(interval+deg1(1),deg1(1)+deg2(1))
                      !column
                      do int1 = interval, deg2(1)+interval
                         if(j1 <= deg1(1)+1)then
                            ind1(1) = int1 + ((j1-1)*(2*deg2(1)+j1))/2
                         else if(j1 > deg1(1)+1 .and. j1 <= deg1(1)+deg2(1))then
                            ind1(1) = int1 + (j1-deg1(1)-1)*(deg1(1)+deg2(1)) + (deg1(1)*(2*deg2(1)+deg1(1)+1))/2!position in first dimension in massline
                         end if
                         ind2(1) = j1-interval+1 !spline1
                         ind3(1) = int1-interval+1 !spline2
                         ! loop over Gauss points
                         do k3=1, q(3)
                            c(3)=(xw_gauss_d3(1,k3)+ real( row(3)-1+ind4(3)-((row(3)-1+ind4(3))/n_cells(3))*n_cells(3),f64 ))/real(n_cells(3),f64)
                            do k2=1, q(2)
                               c(2)=(xw_gauss_d2(1,k2)+ real( row(2)-1+ind4(2)-((row(2)-1+ind4(2))/n_cells(2))*n_cells(2),f64 ))/real(n_cells(2),f64)
                               do k1=1, q(1)
                                  c(1)=(real(ind4(1),f64)*xw_gauss_d1(1,k1)+real(row(1)-1+ind4(1)*(interval-1),f64))/real(n_cells(1),f64)

                                  mass_line(ind1(1)+(ind1(2)-1)*((deg1(1)+deg2(1))**2+(2*deg2(1)+deg1(1)-deg1(1)**2)/2)+(ind1(3)-1)*((deg1(1)+deg2(1))**2+(2*deg2(1)+deg1(1)-deg1(1)**2)/2)*(deg1(2)+deg2(2)+1)) = &
                                       mass_line(ind1(1)+(ind1(2)-1)*((deg1(1)+deg2(1))**2+(2*deg2(1)+deg1(1)-deg1(1)**2)/2)+(ind1(3)-1)*((deg1(1)+deg2(1))**2+(2*deg2(1)+deg1(1)-deg1(1)**2)/2)*(deg1(2)+deg2(2)+1)) + &
                                       xw_gauss_d1(2,k1)* xw_gauss_d2(2,k2)* xw_gauss_d3(2,k3) *&
                                       bspl_b1(ind2(1), k1,interval)*&
                                       bspl_d2a(ind2(2), k2)*&
                                       bspl_d3a(ind2(3), k3)*&
                                       bspl_b2(ind3(1), k1,interval)*&
                                       bspl_d2b(ind3(2), k2)*&
                                       bspl_d3b(ind3(3), k3)*&
                                       profile( c, component )
                               end do
                            end do
                         enddo
                      end do
                   end do
                enddo
             end do
          enddo
       end do
    enddo

  end subroutine sll_s_spline_fem_mixedmass_line_boundary


  !> Multiply by dicrete gradient matrix
  subroutine sll_s_multiply_g_clamped( n_dofs, delta_x, s_deg_0, field_in, field_out )
    sll_int32,  intent( in    )   :: n_dofs(3) !< number of cells (and grid points)
    sll_real64, intent( in    )   :: delta_x(3) !< grid spacing
    sll_int32,  intent( in    )   :: s_deg_0(3) !< highest spline degree
    sll_real64, intent( in    )   :: field_in(:) !< field_in 
    sll_real64, intent(   out )   :: field_out(:) !< G*field_in
    !local variables
    sll_real64 :: coef
    sll_int32  :: stride, stride_end
    sll_int32  :: ind3d, ind1d
    sll_int32  :: i,j,k

    coef = 1.0_f64/ delta_x(1)

    ind1d = 0
    ind3d = 0
    do k = 1, n_dofs(3)
       do j = 1, n_dofs(2)
          ind3d = ind3d + 1
          ind1d = ind1d + 1
          field_out(ind3d) = real(s_deg_0(1),f64)*coef*( field_in(ind1d+1) )
          do i = 2, s_deg_0(1)-1  !1, s_deg_0 -1 perfect condutctor BC
             ind3d = ind3d + 1
             ind1d = ind1d + 1
             field_out(ind3d) = real(s_deg_0(1),f64)/real(i,f64)*coef*( field_in(ind1d+1) - field_in(ind1d) )
          end do
          do i = s_deg_0(1), n_dofs(1)-s_deg_0(1)
             ind3d = ind3d + 1
             ind1d = ind1d + 1
             field_out(ind3d) = coef*( field_in(ind1d+1) - field_in(ind1d) )
          end do
          do i = n_dofs(1)-s_deg_0(1)+1, n_dofs(1)-2 !n_dofs-1 perfect conductor BC
             ind3d = ind3d + 1
             ind1d = ind1d + 1
             field_out(ind3d) = real(s_deg_0(1),f64)/real(n_dofs(1)-i,f64)*coef*( field_in(ind1d+1) - field_in(ind1d) )
          end do
          ind3d = ind3d + 1
          ind1d = ind1d + 1
          field_out(ind3d) = real(s_deg_0(1),f64)*coef*( - field_in(ind1d) )
          ind1d = ind1d + 1
       end do
    end do

    coef = 1.0_f64/ delta_x(2)
    stride = n_dofs(1)
    stride_end = (1-n_dofs(2))*n_dofs(1)

    ind1d = 0
    do k = 1, n_dofs(3)
       ind3d = ind3d + 1
       ind1d = ind1d + 1
       field_out(ind3d) =  0._f64
       do i = 2, n_dofs(1)-1 !1: n_dofs(1) perfect conductor BC
          ind3d = ind3d + 1
          ind1d = ind1d + 1
          field_out(ind3d) =  coef * ( field_in(ind1d) - field_in(ind1d-stride_end) )
       end do
       ind3d = ind3d + 1
       ind1d = ind1d + 1
       field_out(ind3d) =  0._f64

       do j = 2, n_dofs(2)
          ind3d = ind3d + 1
          ind1d = ind1d + 1
          field_out(ind3d) =  0._f64
          do i = 2, n_dofs(1)-1 !1: n_dofs(1) perfect conductor BC
             ind3d = ind3d + 1
             ind1d = ind1d + 1
             field_out(ind3d) =  coef * ( field_in(ind1d) - field_in(ind1d-stride) )
          end do
          ind3d = ind3d + 1
          ind1d = ind1d + 1
          field_out(ind3d) =  0._f64
       end do
    end do

    coef = 1.0_f64/ delta_x(3)
    stride = n_dofs(1)*n_dofs(2)
    stride_end = (1-n_dofs(3))*n_dofs(1)*n_dofs(2)

    ind1d = 0
    do j = 1, n_dofs(2)
       ind3d = ind3d + 1
       ind1d = ind1d + 1
       field_out(ind3d) =  0._f64
       do i = 2, n_dofs(1)-1 !1: n_dofs(1) perfect conductor BC 
          ind3d = ind3d + 1
          ind1d = ind1d + 1
          field_out(ind3d) = coef * ( field_in(ind1d) - field_in(ind1d-stride_end) )
       end do
       ind3d = ind3d + 1
       ind1d = ind1d + 1
       field_out(ind3d) =  0._f64
    end do
    do k = 2, n_dofs(3)
       do j = 1, n_dofs(2)
          ind3d = ind3d + 1
          ind1d = ind1d + 1
          field_out(ind3d) =  0._f64
          do i = 2, n_dofs(1)-1 !1: n_dofs(1) perfect conductor BC
             ind3d = ind3d + 1
             ind1d = ind1d + 1
             field_out(ind3d) = coef * ( field_in(ind1d) - field_in(ind1d-stride) )
          end do
          ind3d = ind3d + 1
          ind1d = ind1d + 1
          field_out(ind3d) =  0._f64
       end do
    end do

  end subroutine sll_s_multiply_g_clamped


  !> Multiply by transpose of dicrete gradient matrix
  subroutine sll_s_multiply_gt_clamped( n_dofs, delta_x, s_deg_0, field_in, field_out )
    sll_int32,  intent( in    )  :: n_dofs(3) !< number of cells (and grid points)
    sll_real64, intent( in    )  :: delta_x(3) !< grid spacing
    sll_int32,  intent( in    )  :: s_deg_0(3) !< highest spline degree
    sll_real64, intent( in    )  :: field_in(:)!< field_in
    sll_real64, intent(   out )  :: field_out(:) !< Gt*field_in
    !local variables
    sll_real64 :: coef
    sll_int32  :: stride, stride_end
    sll_int32  :: ind3d, ind1d
    sll_int32  :: i,j,k


    coef = 1.0_f64/ delta_x(1)

    ind1d = 0
    ind3d = 0
    do k = 1, n_dofs(3)
       do j = 1, n_dofs(2)
          ind3d = ind3d + 1
          ind1d = ind1d + 1
          field_out(ind1d) =  0._f64!- real(s_deg_0(1),f64)*coef*field_in(ind3d) perfect condutctor BC
          do i = 2, s_deg_0(1)
             ind3d = ind3d + 1
             ind1d = ind1d + 1
             field_out(ind1d) = real(s_deg_0(1),f64)*coef*( field_in(ind3d-1)/real(i-1,f64) - field_in(ind3d)/real(i,f64) )
          end do
          do i = s_deg_0(1)+1, n_dofs(1)-s_deg_0(1)
             ind3d = ind3d + 1
             ind1d = ind1d + 1
             field_out(ind1d) = coef*( field_in(ind3d-1) - field_in(ind3d) )
          end do
          do i = n_dofs(1)-s_deg_0(1)+1, n_dofs(1)-1
             ind3d = ind3d + 1
             ind1d = ind1d + 1
             field_out(ind1d) = real(s_deg_0(1),f64)*coef*( field_in(ind3d-1)/real(n_dofs(1)+1-i,f64) - field_in(ind3d)/real(n_dofs(1)-i,f64) )
          end do
          !ind3d = ind3d + 1
          ind1d = ind1d + 1
          field_out(ind1d) =  0._f64!real(s_deg_0(1),f64)*coef*field_in(ind3d)!-1) perfect condutctor BC
          !ind3d = ind3d - 1
       end do
    end do


    coef = 1.0_f64/ delta_x(2)
    stride = n_dofs(1)
    stride_end = (1-n_dofs(2))*n_dofs(1)

    ind1d = 0
    do k = 1, n_dofs(3)
       do j = 1, n_dofs(2)-1
          ind3d = ind3d + 1
          ind1d = ind1d + 1
          field_out(ind1d) = 0._f64
          do i = 2, n_dofs(1)-1 !1, n_dofs(1) perfect conductor BC
             ind3d = ind3d + 1
             ind1d = ind1d + 1

             field_out(ind1d) =  field_out(ind1d) + &
                  coef * ( field_in(ind3d) - field_in(ind3d+stride) )

          end do
          ind3d = ind3d + 1
          ind1d = ind1d + 1
          field_out(ind1d) = 0._f64
       end do
       ind3d = ind3d + 1
       ind1d = ind1d + 1
       field_out(ind1d) = 0._f64
       do i= 2, n_dofs(1)-1 !1, n_dofs(1) perfect conductor BC
          ind3d = ind3d + 1
          ind1d = ind1d + 1

          field_out(ind1d) =  field_out(ind1d) + &
               coef * ( field_in(ind3d) - field_in(ind3d+stride_end) )

       end do
       ind3d = ind3d + 1
       ind1d = ind1d + 1
       field_out(ind1d) = 0._f64
    end do
    coef = 1.0_f64/ delta_x(3)
    stride = n_dofs(1)*n_dofs(2)
    stride_end = (1-n_dofs(3))*n_dofs(1)*n_dofs(2)

    ind1d = 0
    do k = 1, n_dofs(3)-1
       do j = 1, n_dofs(2)
          ind3d = ind3d + 1
          ind1d = ind1d + 1
          field_out(ind1d) = 0._f64
          do i = 2, n_dofs(1)-1 !1, n_dofs(1) perfect conductor BC
             ind3d = ind3d + 1
             ind1d = ind1d + 1

             field_out(ind1d) =  field_out(ind1d) + &
                  coef * ( field_in(ind3d) - field_in(ind3d+stride) )
          end do
          ind3d = ind3d + 1
          ind1d = ind1d + 1
          field_out(ind1d) = 0._f64
       end do
    end do
    do j = 1, n_dofs(2)
       ind3d = ind3d + 1
       ind1d = ind1d + 1
       field_out(ind1d) = 0._f64
       do i = 2, n_dofs(1)-1  !1, n_dofs(1) perfect conductor BC
          ind3d = ind3d + 1
          ind1d = ind1d + 1

          field_out(ind1d) =   field_out(ind1d) + &
               coef * ( field_in(ind3d) - field_in(ind3d+stride_end) )
       end do
       ind3d = ind3d + 1
       ind1d = ind1d + 1
       field_out(ind1d) = 0._f64
    end do

  end subroutine sll_s_multiply_gt_clamped


  !> Multiply by discrete curl matrix
  subroutine sll_s_multiply_c_clamped( n_dofs, delta_x, s_deg_0, field_in, field_out )
    sll_int32,  intent( in    )  :: n_dofs(3)    !< number of cells (and grid points) 
    sll_real64, intent( in    )  :: delta_x(3)   !< grid spacing
    sll_int32,  intent( in    )  :: s_deg_0(3)   !< highest spline degree
    sll_real64, intent( in    )  :: field_in(:)  !< field_in
    sll_real64, intent(   out )  :: field_out(:) !< C*field_in
    ! local variables
    sll_real64 :: coef(2)
    sll_int32 :: stride(2), jump(2), indp(2), n_total0, n_total1
    sll_int32 :: i,j,k, ind3d, ind3d_1, ind3d_2, ind

    n_total0 = product(n_dofs)
    n_total1 = (n_dofs(1)-1)*n_dofs(2)*n_dofs(3)

    ! TODO: Avoid the IF for periodic boundaries
    ! First component
    coef(1) = 1.0_f64/ delta_x(2)
    coef(2) = -1.0_f64/ delta_x(3)

    stride(1) = n_dofs(1)
    stride(2) = n_dofs(1)*n_dofs(2)

    jump(1) = n_total1+n_total0
    jump(2) = n_total1


    ind3d = 0
    do k = 1, n_dofs(3)
       if (k == 1) then
          indp(2) = stride(2)*(n_dofs(3)-1)
       else
          indp(2) = - stride(2)
       end if
       do j = 1, n_dofs(2)
          if ( j==1 ) then
             indp(1) = stride(1)*(n_dofs(2)-1)
          else
             indp(1) = - stride(1)
          end if
          ind3d = ind3d + 1
          field_out(ind3d) = 0._f64 
          do i = 2, n_dofs(1)-1  !1, n_dofs(1) perfect conductor BC
             ind3d = ind3d + 1

             ind3d_1 = ind3d +indp(1)+jump(1)
             ind3d_2 = ind3d +indp(2)+jump(2)

             field_out(ind3d) =  &
                  coef(1) * ( field_in( ind3d+jump(1) ) -&
                  field_in( ind3d_1 )) + &
                  coef(2) * ( field_in(ind3d+jump(2) ) - &
                  field_in( ind3d_2 ))

          end do
          ind3d = ind3d + 1
          field_out(ind3d) =  0._f64
       end do
    end do


    ! Second component
    coef(1) = 1.0_f64/ delta_x(3)
    coef(2) = -1.0_f64/ delta_x(1)

    stride(1) = (n_dofs(1)-1)*n_dofs(2)
    stride(2) = 1

    jump(1) = -n_total0 
    jump(2) = n_total1

    ind = ind3d
    do k = 1, n_dofs(3)
       if (k == 1) then
          indp(1) = stride(1)*(n_dofs(3)-1)
       else
          indp(1) = - stride(1)
       end if
       do j = 1, n_dofs(2)
          ind3d = ind3d + 1
          ind = ind +1
          ind3d_1 = ind3d +indp(1)+jump(1)
          ind3d_2 = ind -stride(2)+jump(2)

          field_out(ind3d) = &
               coef(1) * ( field_in( ind3d+jump(1) ) -&
               field_in( ind3d_1 ) )+ &  
               real(s_deg_0(1),f64)*coef(2)* &
               ( field_in(ind+jump(2)+1 ) )
          do i = 2, s_deg_0(1)-1 !1, s_deg_0(1) -1 perfect condutctor BC
             ind3d = ind3d + 1
             ind = ind +1
             ind3d_1 = ind3d +indp(1)+jump(1)
             ind3d_2 = ind -stride(2)+jump(2)

             field_out(ind3d) = &
                  coef(1) * ( field_in( ind3d+jump(1) ) -&
                  field_in( ind3d_1 ) )+ &  
                  real(s_deg_0(1),f64)/real(i,f64)*coef(2)* &
                  ( field_in(ind+jump(2)+1 ) - &
                  field_in( ind3d_2+1 ) )
          end do
          do i = s_deg_0(1), n_dofs(1)-s_deg_0(1)
             ind3d = ind3d + 1
             ind = ind +1
             ind3d_1 = ind3d +indp(1)+jump(1)
             ind3d_2 = ind -stride(2)+jump(2)

             field_out(ind3d) = &
                  coef(1) * ( field_in( ind3d+jump(1) ) -&
                  field_in( ind3d_1 ) )+ &
                  coef(2) * ( field_in(ind+jump(2)+1 ) - &
                  field_in( ind3d_2+1 ) )
          end do
          do i = n_dofs(1)-s_deg_0(1)+1, n_dofs(1)-2 !n_dofs(1)-1 perfect conductor BC
             ind3d = ind3d + 1
             ind = ind +1
             ind3d_1 = ind3d +indp(1)+jump(1)
             ind3d_2 = ind -stride(2)+jump(2)

             field_out(ind3d) = &
                  coef(1) * ( field_in( ind3d+jump(1) ) -&
                  field_in( ind3d_1 ) )+ &
                  real(s_deg_0(1),f64)/real(n_dofs(1)-i,f64) * &
                  coef(2) * ( field_in(ind+jump(2)+1 ) - &
                  field_in( ind3d_2+1 ) )
          end do
          ind3d = ind3d + 1
          ind = ind +1
          ind3d_1 = ind3d +indp(1)+jump(1)
          ind3d_2 = ind -stride(2)+jump(2)

          field_out(ind3d) = &
               coef(1) * ( field_in( ind3d+jump(1) ) -&
               field_in( ind3d_1 ) )+ &
               real(s_deg_0(1),f64) * &
               coef(2) * ( - field_in( ind3d_2+1 ) )
          ind = ind + 1
       end do
    end do

    ! Third component
    coef(1) = 1.0_f64/ delta_x(1)
    coef(2) = -1.0_f64/ delta_x(2)

    stride(1) = 1
    stride(2) = n_dofs(1)-1

    jump(1) = -n_total0
    jump(2) = -(n_total0+n_total1)


    ind = ind3d
    do k = 1, n_dofs(3)
       do j = 1, n_dofs(2)
          if (j == 1) then
             indp(2) = stride(2)*(n_dofs(2)-1)
          else
             indp(2) = - stride(2)
          end if
          ind3d = ind3d + 1
          ind = ind + 1
          ind3d_1 = ind -stride(1)+jump(1)
          ind3d_2 = ind3d +indp(2)+jump(2)

          field_out(ind3d) = &
               real(s_deg_0(1),f64) * &
               coef(1) * ( field_in( ind+jump(1)+1 ) )+ &
               coef(2)*( field_in(ind3d+jump(2) ) - &
               field_in( ind3d_2 ) )
          do i = 2, s_deg_0(1)-1 !1, s_deg_0(1) -1 perfect condutctor BC
             ind3d = ind3d + 1
             ind = ind + 1
             ind3d_1 = ind -stride(1)+jump(1)
             ind3d_2 = ind3d +indp(2)+jump(2)

             field_out(ind3d) = &
                  real(s_deg_0(1),f64)/real(i,f64) * &
                  coef(1) * ( field_in( ind+jump(1)+1 ) -&
                  field_in( ind3d_1+1 ) )+ &
                  coef(2)*( field_in(ind3d+jump(2) ) - &
                  field_in( ind3d_2 ) )
          end do
          do i = s_deg_0(1), n_dofs(1)-s_deg_0(1)
             ind3d = ind3d + 1
             ind = ind + 1
             ind3d_1 = ind -stride(1)+jump(1)
             ind3d_2 = ind3d +indp(2)+jump(2)

             field_out(ind3d) = &
                  coef(1) * ( field_in( ind+jump(1)+1 ) -&
                  field_in( ind3d_1+1 ) )+ &
                  coef(2) * ( field_in(ind3d+jump(2) ) - &
                  field_in( ind3d_2 ) )
          end do
          do i = n_dofs(1)-s_deg_0(1)+1, n_dofs(1)-2  !n_dofs(1)-1 perfect conductor BC
             ind3d = ind3d + 1
             ind = ind + 1
             ind3d_1 = ind -stride(1)+jump(1)
             ind3d_2 = ind3d +indp(2)+jump(2)

             field_out(ind3d) = &
                  real(s_deg_0(1),f64)/real(n_dofs(1)-i,f64) * &
                  coef(1) * ( field_in( ind+jump(1)+1 ) -&
                  field_in( ind3d_1+1 ) )+ &
                  coef(2) * ( field_in(ind3d+jump(2) ) - &
                  field_in( ind3d_2 ) )
          end do
          ind3d = ind3d + 1
          ind = ind + 1
          ind3d_1 = ind -stride(1)+jump(1)
          ind3d_2 = ind3d +indp(2)+jump(2)

          field_out(ind3d) = &
               real(s_deg_0(1),f64) * &
               coef(1) * ( - field_in( ind3d_1+1 ) )+ &
               coef(2) * ( field_in(ind3d+jump(2) ) - &
               field_in( ind3d_2 ) )

          ind = ind + 1
       end do
    end do

  end subroutine sll_s_multiply_c_clamped


  !> Multiply by transpose of discrete curl matrix
  subroutine sll_s_multiply_ct_clamped( n_dofs, delta_x, s_deg_0, field_in, field_out )
    sll_int32,  intent( in    )  :: n_dofs(3)    !< number of cells (and grid points) 
    sll_real64, intent( in    )  :: delta_x(3)   !< grid spacing
    sll_int32,  intent( in    )  :: s_deg_0(3)   !< highest spline degree
    sll_real64, intent( in    )  :: field_in(:)  !< field_in
    sll_real64, intent(   out )  :: field_out(:) !< C^T*field_in
    ! local variables
    sll_real64 :: coef(2)
    sll_int32 :: stride(2), jump(2), indp(2), n_total0, n_total1
    sll_int32 :: i,j,k, ind3d, ind3d_1, ind3d_2, ind 

    n_total0 = product(n_dofs)
    n_total1 = (n_dofs(1)-1)*n_dofs(2)*n_dofs(3)

    ! TODO: Avoid the IF for periodic boundaries
    ! First component
    coef(1) = -1.0_f64/ delta_x(2)
    coef(2) = 1.0_f64/ delta_x(3)

    stride(1) = n_dofs(1)-1
    stride(2) = (n_dofs(1)-1)*n_dofs(2)

    jump(1) = n_total0
    jump(2) = n_total0+n_total1


    ind3d = 0
    do k = 1, n_dofs(3)
       if (k == n_dofs(3)) then
          indp(2) = -stride(2)*(n_dofs(3)-1)
       else
          indp(2) = stride(2)
       end if
       do j = 1, n_dofs(2)
          if ( j== n_dofs(2)) then
             indp(1) = -stride(1)*(n_dofs(2)-1)
          else
             indp(1) = stride(1)
          end if
          do i = 1, n_dofs(1)-1
             ind3d = ind3d + 1

             ind3d_1 = ind3d +indp(1)+jump(2)
             ind3d_2 = ind3d +indp(2)+jump(1)

             field_out(ind3d) =  &
                  coef(1) * ( field_in( ind3d+jump(2) ) -&
                  field_in( ind3d_1 )) + &
                  coef(2) * ( field_in(ind3d+jump(1) ) - &
                  field_in( ind3d_2 ))

          end do
       end do
    end do


    ! Second component
    coef(1) = -1.0_f64/ delta_x(3)
    coef(2) = 1.0_f64/ delta_x(1)

    stride(1) = n_dofs(1)*n_dofs(2)
    stride(2) = 1


    jump(1) = n_total0
    jump(2) = -n_total1

    ind = ind3d
    do k = 1, n_dofs(3)
       if (k == n_dofs(3)) then
          indp(1) = -stride(1)*(n_dofs(3)-1)
       else
          indp(1) = stride(1)
       end if
       do j = 1, n_dofs(2)
          ind3d = ind3d + 1
          ind = ind + 1
          ind3d_1 = ind3d +indp(1)+jump(2)

          field_out(ind3d) = 0._f64 
!!$               coef(1) * ( field_in( ind3d+jump(2) ) -&
!!$               field_in( ind3d_1 ) )+ &
!!$               real(s_deg_0(1),f64) * coef(2) * &
!!$               (- field_in( ind+jump(1) ) )
          do i = 2, s_deg_0(1)
             ind3d = ind3d + 1 
             ind = ind + 1
             ind3d_1 = ind3d +indp(1)+jump(2)
             ind3d_2 = ind -stride(2)+jump(1)
             field_out(ind3d) = &
                  coef(1) * ( field_in( ind3d+jump(2) ) -&
                  field_in( ind3d_1 ) )+ &
                  real(s_deg_0(1),f64)*coef(2) * &
                  ( field_in( ind3d_2 )/real(i-1,f64) - &
                  field_in( ind+jump(1) )/real(i,f64) )
          end do
          do i= s_deg_0(1)+1, n_dofs(1)-s_deg_0(1)
             ind3d = ind3d + 1
             ind = ind + 1
             ind3d_1 = ind3d +indp(1)+jump(2)
             ind3d_2 = ind -stride(2)+jump(1)

             field_out(ind3d) = &
                  coef(1) * ( field_in( ind3d+jump(2) ) -&
                  field_in( ind3d_1 ) )+ &
                  coef(2) * ( field_in(ind3d_2 ) - &
                  field_in( ind+jump(1) ))
          end do
          do i = n_dofs(1)-s_deg_0(1)+1, n_dofs(1)-1
             ind3d = ind3d + 1
             ind = ind + 1
             ind3d_1 = ind3d +indp(1)+jump(2)
             ind3d_2 = ind -stride(2)+jump(1)

             field_out(ind3d) = &
                  coef(1) * ( field_in( ind3d+jump(2) ) -&
                  field_in( ind3d_1 ) )+ &
                  real(s_deg_0(1),f64)*coef(2) * &
                  ( field_in(ind3d_2 )/real(n_dofs(1)+1-i,f64) - &
                  field_in( ind+jump(1) )/real(n_dofs(1)-i,f64) )
          end do
          ind3d = ind3d + 1
          ind = ind + 1
          ind3d_1 = ind3d +indp(1)+jump(2)
          ind3d_2 = ind -stride(2)+jump(1)
          field_out(ind3d) = 0._f64
!!$               coef(1) * ( field_in( ind3d+jump(2) ) -&
!!$               field_in( ind3d_1 ) )+ &
!!$               real(s_deg_0(1),f64)*coef(2) * &
!!$               field_in(ind3d_2 )
          ind = ind - 1
       end do
    end do

    ! Third component
    coef(1) = -1.0_f64/ delta_x(1)
    coef(2) = 1.0_f64/ delta_x(2)

    stride(1) = 1
    stride(2) = n_dofs(1)

    jump(1) = -(n_total1+n_total0)
    jump(2) = -n_total1

    ind = ind3d
    do k = 1, n_dofs(3)
       do j = 1, n_dofs(2)
          if (j == n_dofs(2)) then
             indp(2) = -stride(2)*(n_dofs(2)-1)
          else
             indp(2) = stride(2)
          end if
          ind3d = ind3d + 1
          ind = ind + 1
          ind3d_2 = ind3d +indp(2)+jump(1)

          field_out(ind3d) =  0._f64
!!$               real(s_deg_0(1),f64) * coef(1) * &
!!$               ( -field_in( ind+jump(2) ) )+ &
!!$               coef(2) * ( field_in(ind3d+jump(1) ) - &
!!$               field_in( ind3d_2 ))
          do i = 2, s_deg_0(1)
             ind3d = ind3d + 1
             ind = ind + 1
             ind3d_1 = ind -stride(1)+jump(2)
             ind3d_2 = ind3d +indp(2)+jump(1)

             field_out(ind3d) = &
                  real(s_deg_0(1),f64) * coef(1) * &
                  ( field_in( ind3d_1 )/real(i-1,f64) -&
                  field_in( ind+jump(2) )/real(i,f64) )+ &
                  coef(2) * ( field_in(ind3d+jump(1) ) - &
                  field_in( ind3d_2 ))
          end do
          do i = s_deg_0(1)+1, n_dofs(1)-s_deg_0(1)
             ind3d = ind3d + 1
             ind = ind + 1
             ind3d_1 = ind -stride(1)+jump(2)
             ind3d_2 = ind3d +indp(2)+jump(1)

             field_out(ind3d) = &
                  coef(1) * ( field_in( ind3d_1 ) -&
                  field_in( ind+jump(2)) )+ &
                  coef(2) * ( field_in(ind3d+jump(1) ) - &
                  field_in( ind3d_2 ))

          end do
          do i = n_dofs(1)-s_deg_0(1)+1, n_dofs(1)-1
             ind3d = ind3d + 1
             ind = ind + 1
             ind3d_1 = ind -stride(1)+jump(2)
             ind3d_2 = ind3d +indp(2)+jump(1)

             field_out(ind3d) = &
                  real(s_deg_0(1),f64) * coef(1) * &
                  ( field_in( ind3d_1 )/real(n_dofs(1)+1-i,f64) -&
                  field_in( ind+jump(2) )/real(n_dofs(1)-i,f64) )+ &
                  coef(2) * ( field_in(ind3d+jump(1) ) - &
                  field_in( ind3d_2 ))
          end do
          ind3d = ind3d + 1
          ind = ind + 1
          ind3d_1 = ind -stride(1)+jump(2)
          ind3d_2 = ind3d +indp(2)+jump(1)

          field_out(ind3d) =  0._f64
!!$               real(s_deg_0(1),f64) * coef(1) * &
!!$               field_in( ind3d_1 ) + &
!!$               coef(2) * ( field_in(ind3d+jump(1) ) - &
!!$               field_in( ind3d_2 ))
          ind = ind - 1 
       end do
    end do

  end subroutine sll_s_multiply_ct_clamped


  !> Multiplication of the input vector \a in by the clamped derivative matrix G
  subroutine sll_s_multiply_g_clamped_1d( n_dofs, delta_x, s_deg_0, in, out )
    sll_int32,  intent( in    )  :: n_dofs  !< number of cells (and grid points)
    sll_real64, intent( in    )  :: delta_x !< grid spacing
    sll_int32,  intent( in    )  :: s_deg_0 !< highest spline degree 
    sll_real64, intent( in    )  :: in(:)   !< field_in
    sll_real64, intent(   out )  :: out(:)  !< G*field_in
    !local variable
    sll_int32 :: i

    out(1) = real(s_deg_0,f64)*( in(2) )/delta_x
    do i = 2, s_deg_0-1  !1, s_deg_0 -1 perfect condutctor BC
       out(i) = real(s_deg_0,f64)*( in(i+1) - in(i) )/(real(i,f64)*delta_x)
    end do
    do i = s_deg_0, n_dofs-s_deg_0
       out(i) = ( in(i+1) - in(i) )/delta_x
    end do
    do i = n_dofs-s_deg_0+1, n_dofs-2 !n_dofs-1 perfect conductor BC
       out(i) = real(s_deg_0,f64)*( in(i+1) - in(i) )/(real(n_dofs-i,f64)*delta_x)
    end do
    out( n_dofs-1) = real(s_deg_0,f64)*( - in( n_dofs-1) )/delta_x

  end subroutine sll_s_multiply_g_clamped_1d

  !> Multiplication of the input vector \a in by the transposed clamped derivative matrix G^T 
  subroutine sll_s_multiply_gt_clamped_1d( n_dofs, delta_x, s_deg_0, in, out )
    sll_int32,  intent( in    )  :: n_dofs  !< number of cells (and grid points)
    sll_real64, intent( in    )  :: delta_x !< grid spacing
    sll_int32,  intent( in    )  :: s_deg_0 !< highest spline degree 
    sll_real64, intent( in    )  :: in(:)   !< field_in
    sll_real64, intent(   out )  :: out(:)  !< G*field_in
    !local variable
    sll_int32 :: i

    out(1) = 0._f64!- real(s_deg_0,f64)*in(1)/delta_x perfect condutctor BC
    do i = 2, s_deg_0
       out(i) = real(s_deg_0,f64)*( in(i-1)/real(i-1,f64) - in(i)/real(i,f64) )/delta_x
    end do
    do i = s_deg_0+1, n_dofs-s_deg_0
       out(i) = ( in(i-1) - in(i) )/delta_x
    end do
    do i = n_dofs-s_deg_0+1, n_dofs-1
       out(i) = real(s_deg_0,f64)*( in(i-1)/real(n_dofs+1-i,f64) - in(i)/real(n_dofs-i,f64) )/delta_x
    end do
    out(n_dofs) = 0._f64!real(s_deg_0,f64)*in(n_dofs-1)/delta_x perfect condutctor BC

  end subroutine sll_s_multiply_gt_clamped_1d


end module sll_m_spline_fem_utilities_3d_clamped

