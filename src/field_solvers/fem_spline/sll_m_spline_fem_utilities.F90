!> @ingroup maxwell_solvers
!> @brief
!> Utilites for Maxwell solver's with spline finite elements
!> @details
!> 
!> @author
!> Katharina Kormann

module sll_m_spline_fem_utilities
  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#include "sll_assert.h"
#include "sll_memory.h"
#include "sll_working_precision.h"

  use sll_m_low_level_bsplines, only: &
       sll_s_uniform_bsplines_eval_basis

  use sll_m_splines_pp

  use sll_m_constants, only : &
       sll_p_twopi

  use sll_m_gauss_legendre_integration, only: &
       sll_f_gauss_legendre_points_and_weights
  
  implicit none

  public :: &
       sll_s_spline_fem_mass_line, &
       sll_s_spline_fem_mass_line_boundary, &
       sll_s_spline_fem_mass_line_boundary_full, &
       sll_s_spline_fem_mass_line_full, &
       sll_s_spline_fem_mixedmass_line, &
       sll_s_spline_fem_mixedmass_line_full, &
       sll_s_spline_fem_compute_mass_eig, &
       sll_s_spline_fem_multiply_mass, &
       sll_s_spline_fem_multiply_massmixed, &
       sll_s_spline_fem_interpolation_eigenvalues, &
       sll_s_multiply_g_1d, &
       sll_s_multiply_gt_1d, &
       sll_s_multiply_g, &
       sll_s_multiply_gt, &
       sll_s_multiply_c, &
       sll_s_multiply_ct, &
       sll_i_profile_function_1d


  private
  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  abstract interface
     function sll_i_profile_function_1d( x ) result(res)
       use sll_m_working_precision
       sll_real64, intent(in) :: x
       sll_real64             :: res
     end function sll_i_profile_function_1d
  end interface
  
contains

  !---------------------------------------------------------------------------!
  !> Computes the mass line for a mass matrix with \a degree splines
  subroutine sll_s_spline_fem_mass_line ( degree, mass_line )
    sll_int32,  intent(in   ) :: degree !< spline degree
    sll_real64, intent(  out) :: mass_line(degree+1) !< values of the righter half (including diagonal) of the mass line 

    !local variables
    sll_int32 :: q, j, int, quad
    sll_real64, allocatable :: quad_xw(:,:), spline_val(:,:)

    q = degree+1

    allocate(quad_xw(2, q))
    allocate(spline_val( degree+1, q ))

    quad_xw = sll_f_gauss_legendre_points_and_weights( q , 0.0_f64, 1.0_f64 )


    do j=1,q
       call sll_s_uniform_bsplines_eval_basis( degree, quad_xw(1,j), spline_val(:,j) )
    end do

    mass_line = 0.0_f64
    do j=1,degree+1
       do int= j, degree+1
          do quad = 1, q
             mass_line(j) = mass_line(j) + &
                  spline_val(int,quad) * spline_val(int-j+1,quad)*quad_xw(2,quad)
          end do
       end do
    end do


  end subroutine sll_s_spline_fem_mass_line

  !> Function computing the boundary part for a mass matrix with clamped splines of \a degree
  subroutine sll_s_spline_fem_mass_line_boundary ( degree, spline, mass_line )
    sll_int32,  intent(in   ) :: degree  !< Spline degree
    type(sll_t_spline_pp_1d), intent(in) :: spline !< pp spline class
    sll_real64, intent(  out) :: mass_line(:) !< Values of the boundary part of the mass line on output

    !local variables
    sll_int32 :: q, j, int, quad
    sll_real64, allocatable :: quad_xw(:,:), spline_val(:,:),spline_val_b(:,:)

    q = degree+1

    allocate(quad_xw(2, q))
    allocate(spline_val( degree+1, q ))
    allocate(spline_val_b( degree+1, q ))
    spline_val_b=0._f64
    spline_val=0._f64

    quad_xw = sll_f_gauss_legendre_points_and_weights( q , 0.0_f64, 1.0_f64 )


    select case(degree)
    case(0)
       print*, 'fem_utilities:degree 0 not implemented'
    case(1)
       do quad = 1, q
          do int = 1, degree+1
             spline_val(int,quad) = sll_f_spline_pp_horner_1d( degree, spline%poly_coeffs, quad_xw(1,quad), int)
          end do
       end do

       !first interval
       mass_line = 0.0_f64
       do j = 1, degree
          do int = j, degree+1
             do quad = 1, q
                mass_line(int) = mass_line(int) + &
                     spline_val(j,quad) * spline_val(int,quad) * quad_xw(2,quad)
             end do
          end do
       end do
    case(2)
       do quad = 1, q
          do int = 1, degree+1
             spline_val_b(int,quad) = sll_f_spline_pp_horner_1d( degree, spline%poly_coeffs_boundary_left(:,:,1), quad_xw(1,quad), int)
          end do
       end do

       mass_line = 0.0_f64
       !first interval
       do j = 1, degree
          do int = j, degree+1
             do quad = 1, q
                mass_line(int-j+1+q*(j-1)) = mass_line(int-j+1+q*(j-1)) + &
                     spline_val_b(j,quad) * spline_val_b(int,quad) * quad_xw(2,quad)
             end do
          end do
       end do


       do quad = 1, q
          do int = 1, degree+1
             spline_val(int,quad) = sll_f_spline_pp_horner_1d( degree, spline%poly_coeffs, quad_xw(1,quad), int)
          end do
       end do

       !last interval
       do j = degree, degree
          do int = 1, degree+1
             do quad = 1, q
                mass_line(int+q*(j-1)) = mass_line(int+q*(j-1)) + &
                     spline_val(1,quad) * spline_val(int,quad) * quad_xw(2,quad)
             end do
          end do
       end do

    case(3)
       do quad = 1, q
          do int = 1, degree+1
             spline_val_b(int,quad) = sll_f_spline_pp_horner_1d( degree, spline%poly_coeffs_boundary_left(:,:,1), quad_xw(1,quad), int)
          end do
       end do

       mass_line = 0.0_f64
       !first interval
       do j = 1, degree
          do int = j, degree+1
             do quad = 1, q
                mass_line(int-j+1+q*(j-1)) = mass_line(int-j+1+q*(j-1)) + &
                     spline_val_b(j,quad) * spline_val_b(int,quad) * quad_xw(2,quad)
             end do
          end do
       end do

       spline_val_b = 0._f64
       do quad = 1, q
          do int = 1, degree+1
             spline_val_b(int,quad) = sll_f_spline_pp_horner_1d( degree, spline%poly_coeffs_boundary_left(:,:,2), quad_xw(1,quad), int)
          end do
       end do

       !second interval
       do j = 1, degree-1
          do int = j, degree+1
             do quad = 1, q
                mass_line(int-j+1+q*j)=mass_line(int-j+1+q*j) + &
                     spline_val_b(j,quad) * spline_val_b(int,quad) * quad_xw(2,quad)
             end do
          end do
       end do

       do quad = 1, q
          do int = 1, degree+1
             spline_val(int,quad) = sll_f_spline_pp_horner_1d( degree, spline%poly_coeffs, quad_xw(1,quad), int)
          end do
       end do

       !last interval
       do j = degree, degree
          do int = 1, degree+1
             do quad = 1, q
                mass_line(int+q*(j-1)) = mass_line(int+q*(j-1)) + &
                     spline_val(1,quad) * spline_val(int,quad) * quad_xw(2,quad)
             end do
          end do
       end do

    case default
       print*, "for chosen spline degree, no boundary mass_line is implemented"
    end select


  end subroutine sll_s_spline_fem_mass_line_boundary

  !---------------------------------------------------------------------------!
  !for use of affine transformation
  !> Function computing the boundary part for a mass matrix with clamped splines of \a degree (full version without symmetry part)
  subroutine sll_s_spline_fem_mass_line_boundary_full ( deg, profile, spline, row, n_cells, mass_line )
    sll_int32, intent( in   ) :: deg  !< spline degree
    procedure(sll_i_profile_function_1d) ::  profile !< profile function
    type(sll_t_spline_pp_1d), intent(in) :: spline   !< pp spline class
    sll_int32, intent( in   ) :: row  !< row of the mass line in the matrix
    sll_int32, intent( in   ) :: n_cells !< no of grid cells
    sll_real64, intent(  out) :: mass_line((7*deg**2-deg-2)/2) !< Values of the boundary part of the mass line on output
    !local variables
    sll_int32 :: q, j, int, ind1, quad, interval, ind4
    sll_real64, allocatable :: quad_xw(:,:),spline_val_b(:,:,:), spline_val(:,:)
    sll_real64 :: c

    q = min(3*deg+1,10)

    if(deg==0)then
       print*, 'fem_utilities:degree 0 not implemented'
       stop
    end if

    allocate(quad_xw(2, q))
    allocate(spline_val( deg+1, q ))
    allocate(spline_val_b( deg+1, q, 2*deg-1 ))

    quad_xw = sll_f_gauss_legendre_points_and_weights( q , 0.0_f64, 1.0_f64 )
    spline_val=0._f64
    do quad = 1, q
       call sll_s_uniform_bsplines_eval_basis( deg, quad_xw(1,quad), spline_val(:, quad) )
    end do

    spline_val_b=0._f64
    do interval = 1, deg-1
       do quad = 1, q
          do int = 1, deg+1
             spline_val_b(int,quad,interval) =   sll_f_spline_pp_horner_1d( deg, spline%poly_coeffs_boundary_left(:,:,interval), quad_xw(1,quad), int)
          end do
       end do
    end do
    do interval = deg, 2*deg-1
       spline_val_b(:,:,interval) = spline_val
    end do

    if(row<=2*deg)then
       ind4 = 1
    else if(row>n_cells-deg+1 .and. row<= n_cells+deg)then
       ind4 = -1
    end if
    mass_line = 0.0_f64

    !interval
    do interval = 1, 2*deg-1
       !row
       do j = interval, min(interval+deg,2*deg-1)
          !column
          do int = interval, deg+interval
             if(j <= deg+1)then
                ind1 = int + ((j-1)*(2*deg+j))/2 !position in first dimension in massline
             else if( j > deg+1 .and. j < 2*deg)then
                ind1 = int + (j-deg-1)*2*deg + (3*deg**2+deg)/2 !position in first dimension in massline
             end if
             do quad = 1, q
                c = (real(ind4,f64)*quad_xw(1,quad)+real(row-1+ind4*(interval-1),f64))/real(n_cells,f64)
                mass_line(ind1) = mass_line(ind1) + &
                     spline_val_b(j-interval+1,quad,interval) * spline_val_b(int+1-interval,quad,interval) * quad_xw(2,quad)*profile(c)
             end do
          end do
       end do
    end do

  end subroutine sll_s_spline_fem_mass_line_boundary_full



  !---------------------------------------------------------------------------!
  !for use of affine transformation
  !> Helper function to find line of mass matrix (full version without symmetry part)
  subroutine sll_s_spline_fem_mass_line_full ( deg, profile, mass_line, row, n_cells )
    sll_int32,  intent(in   ) :: deg !< spline degree
    procedure(sll_i_profile_function_1d) ::  profile !< profile function
    sll_real64, intent(  out) :: mass_line(2*deg+1)   !< full mass line
    sll_int32, intent(in) :: row        !< row of the mass line in the matrix
    sll_int32, intent(in) :: n_cells      !< no of grid cells
    !local variables
    sll_int32 :: j, int, quad, q, ind2, ind3, ind4
    sll_real64, allocatable :: quad_xw(:,:), spline_val(:,:)
    sll_real64 :: c

    q = min(3*deg+1,10)

    allocate(quad_xw(2, q))
    allocate(spline_val( deg+1, q ))

    quad_xw = sll_f_gauss_legendre_points_and_weights( q , 0.0_f64, 1.0_f64 )
    spline_val=0._f64
    do quad = 1, q
       call sll_s_uniform_bsplines_eval_basis( deg, quad_xw(1,quad), spline_val(:,quad) )
    end do
    
    mass_line = 0.0_f64
    do j = 1, deg
       do int = j, deg
          ind2 = deg+1-int+j !spline1
          ind3 = deg+1-int !spline2
          ind4 = int-j !left cellborder
          do quad = 1, q
             c = (quad_xw(1,quad)+ real( row-1+ind4-((row-1+ind4)/n_cells)*n_cells,f64 ))/real(n_cells,f64)    

             mass_line(deg+1-j)=mass_line(deg+1-j)+ &
                  spline_val(ind2,quad)*spline_val(ind3,quad)*quad_xw(2,quad)*profile(c)
          end do
       end do
    end do
    do j=1,deg+1
       do int= j, deg+1
          ind2 = int-j+1 !spline1
          ind3 = int !spline2
          ind4 = deg+j-int !left cellborder
          do quad = 1, q
             c = (quad_xw(1,quad)+ real( row-1+ind4-((row-1+ind4)/n_cells)*n_cells,f64 ))/real(n_cells,f64)    

             mass_line(j+deg) = mass_line(j+deg) + &
                  spline_val(ind2,quad)*spline_val(ind3,quad)*quad_xw(2,quad)*profile(c)
          end do
       end do
    end do

  end subroutine sll_s_spline_fem_mass_line_full



  !---------------------------------------------------------------------------!
!!$  !> Helper function to find line of mass matrix ( N_i^p N_j^{p-1}))
!!$  subroutine sll_s_spline_fem_mixedmass_line_full ( deg, mass_line, cell, n_cells) 
!!$    sll_int32,  intent(in   ) :: deg
!!$    sll_real64, intent(  out) :: mass_line(deg*2)
!!$    sll_int32, intent(in) :: cell
!!$    sll_int32, intent(in) :: n_cells
!!$    !local variables
!!$    sll_int32 :: q, i, int, quad
!!$    sll_real64, allocatable :: quad_xw(:,:), spline_val_0(:,:), spline_val_1(:,:)
!!$    ! sll_real64 :: jacobian
!!$
!!$    q = min(2*deg+1,10)
!!$
!!$    allocate(quad_xw(2, q))
!!$    allocate(spline_val_0( deg+1, q ))
!!$    allocate(spline_val_1( deg, q ))
!!$
!!$    quad_xw = sll_f_gauss_legendre_points_and_weights( q , 0.0_f64, 1.0_f64 )
!!$
!!$
!!$    do quad = 1, q
!!$       call sll_s_uniform_bsplines_eval_basis( deg, quad_xw(1,quad), spline_val_0(:,quad) )
!!$       call sll_s_uniform_bsplines_eval_basis( deg-1, quad_xw(1,quad), spline_val_1(:,quad) )
!!$    end do
!!$
!!$    mass_line = 0.0_f64
!!$    do i= 2, deg+1
!!$       do int= i, deg+1
!!$          do quad = 1, q
!!$             !jacobian= map%jacobian((quad_xw(1,quad)+cell+int-1)*delta_x, 0.0_f64)
!!$             mass_line(i+deg-1) = mass_line(i+deg-1) + &
!!$                  spline_val_0(int,quad) * spline_val_1(int-i+1,quad)*quad_xw(2,quad)!*jacobian
!!$          end do
!!$       end do
!!$    end do
!!$
!!$    do  i= -deg+1, 0
!!$       do int = 1, deg+i
!!$          do quad = 1, q
!!$             ! jacobian= map%jacobian((quad_xw(1,quad)+cell+int-1)*delta_x, 0.0_f64)
!!$             mass_line(i+deg) = mass_line(i+deg) + &
!!$                  spline_val_0(int,quad) * spline_val_1(int-i,quad)*quad_xw(2,quad)!*jacobian
!!$          end do
!!$       end do
!!$    end do
!!$
!!$
!!$  end subroutine sll_s_spline_fem_mixedmass_line_full



    !> Helper function to find line of mass matrix ( N_i^p N_j^{p-1}))
  subroutine sll_s_spline_fem_mixedmass_line_full ( deg1, deg2, profile, mass_line, cell, n_cells) 
    sll_int32,  intent(in   ) :: deg1, deg2
    procedure(sll_i_profile_function_1d) ::  profile !< profile function
    sll_real64, intent(  out) :: mass_line(deg1+deg2+1)
    sll_int32, intent(in) :: cell
    sll_int32, intent(in) :: n_cells
    !local variables
    sll_int32 :: q, j, int, k, l, limit, ind1, ind2, ind3, ind4
    sll_real64, allocatable :: quad_xw(:,:), spline_val_1(:,:), spline_val_2(:,:)
    sll_real64 :: c

    q = min(deg1+deg2+1,10)

    allocate(quad_xw(2, q))
    allocate(spline_val_1( deg1+1, q ))
    allocate(spline_val_2( deg2+1, q ))

    quad_xw = sll_f_gauss_legendre_points_and_weights( q , 0.0_f64, 1.0_f64 )


    do k = 1, q
       call sll_s_uniform_bsplines_eval_basis( deg1, quad_xw(1,k), spline_val_1(:,k) )
       call sll_s_uniform_bsplines_eval_basis( deg2, quad_xw(1,k), spline_val_2(:,k) )
    end do

    mass_line = 0.0_f64
    do j = 1, deg1+deg2+1 
       if(j <= deg2)then
          limit = deg2
          ind1  = deg2+1-j
       else
          limit = deg1+deg2+1
          ind1  = j
       end if
       ! loop over spline cells
       do int = j, min(limit, j+deg2)
          if(j <= deg2)then
             ind2 = deg1+1-int+j !spline1
             ind3 = deg2+1-int !spline2
             ind4 = int-j !left cellborder
          else
             ind2 = limit-int+1 !spline1
             ind3 = deg2+1-int+j !spline2
             ind4 = int-deg2-1 !left cellborder
          end if
          do l = 1, q
             k = 1 - (-1)**l * l/2 + modulo(l+1,2)*q
             c = (quad_xw(1,k)+ real( cell-1+ind4-((cell-1+ind4)/n_cells)*n_cells,f64 ))/real(n_cells,f64)
             !jacobian= map%jacobian( [c, 0._f64, 0._f64] )

             mass_line(ind1)=mass_line(ind1)+ &
                  spline_val_1(ind2,k)*spline_val_2(ind3,k)*quad_xw(2,k)*profile(c)!*abs(jacobian)
          end do
       end do
    end do
    
  end subroutine sll_s_spline_fem_mixedmass_line_full

  !---------------------------------------------------------------------------!
  !> Helper function to find line of mass matrix ( N_i^p N_j^{p-1}))
  subroutine sll_s_spline_fem_mixedmass_line ( deg, mass_line )
    sll_int32,  intent(in   ) :: deg
    sll_real64, intent(  out) :: mass_line(deg*2)

    !local variables
    sll_int32 :: q, j, int, quad
    sll_real64, allocatable :: quad_xw(:,:), spline_val_0(:,:), spline_val_1(:,:)

    q = min(3*deg+1,10)

    allocate(quad_xw(2, q))
    allocate(spline_val_0( deg+1, q ))
    allocate(spline_val_1( deg, q ))

    quad_xw = sll_f_gauss_legendre_points_and_weights( q , 0.0_f64, 1.0_f64 )


    do j=1,q
       call sll_s_uniform_bsplines_eval_basis( deg, quad_xw(1,j), spline_val_0(:,j) )
       call sll_s_uniform_bsplines_eval_basis( deg-1, quad_xw(1,j), spline_val_1(:,j) )
    end do

    mass_line = 0.0_f64
    do j=2,deg+1
       do int= j, deg+1
          do quad = 1, q
             mass_line(j+deg-1) = mass_line(j+deg-1) + &
                  spline_val_0(int,quad) * spline_val_1(int-j+1,quad)*quad_xw(2,quad)
          end do
       end do
    end do

    do j=-deg+1,0
       do int= 1, deg+j
          do quad = 1, q
             mass_line(j+deg) = mass_line(j+deg) + &
                  spline_val_0(int,quad) * spline_val_1(int-j,quad)*quad_xw(2,quad)
          end do
       end do
    end do


  end subroutine sll_s_spline_fem_mixedmass_line


  !> Compute eigenvalues of mass matrix with given line \a mass_line
  subroutine sll_s_spline_fem_compute_mass_eig( n_cells, degree, mass_line, eig_values_mass )
    sll_int32, intent( in ) :: n_cells   !< Number of cells
    sll_int32, intent( in ) :: degree    !< Spline degree 
    sll_real64, intent( in ) :: mass_line(0:degree)  !< mass line 
    sll_real64, intent( out ) :: eig_values_mass(n_cells)  !< eigenvalues of the mass matrix

    sll_int32 :: j,k
    sll_real64 :: factor


    factor = sll_p_twopi/real(n_cells,f64)
    do k=0, n_cells-1
       eig_values_mass(k+1) = mass_line(0)
       do j=1,degree
          eig_values_mass(k+1) = eig_values_mass(k+1)+ &
               mass_line(j) * 2.0_f64 * cos( factor*real(k*j,f64))
       end do
    end do

  end subroutine sll_s_spline_fem_compute_mass_eig

!> Multiply the vector \a invec with the spline FEM mass matrix with mass line \a mass 
  subroutine sll_s_spline_fem_multiply_mass ( n_cells, degree, mass, invec, outvec )
    sll_int32,  intent( in    ) :: n_cells  !< no. of cells in the grid
    sll_int32,  intent( in    ) :: degree   !< spline degree
    sll_real64, intent( in    ) :: mass(degree+1)   !< mass line
    sll_real64, intent( in    ) :: invec(:)  !< input vector
    sll_real64, intent( out   ) :: outvec(:) !< output vector

    !local variables
    sll_int32 :: row, column, ind

    ind = 1
    ! For the first degree rows we need to put the first part to the back due to periodic boundaries
    do row = 1, degree
       outvec(row) = mass(1)*invec(row)
       do column = 1,row-1
          outvec(row) = outvec(row) + mass(column+1)*(invec(row+column) + &
               invec(row-column))
       end do
       do column = row,degree
          outvec(row)  = outvec(row) + mass(column+1)*(invec(row+column) + &
               invec(row-column+n_cells))
       end do
    end do

    do row = degree+1, n_cells-degree
       outvec(row) = mass(1)*invec(row)
       do column = 1, degree
          outvec(row) = outvec(row) + mass(column+1)*(invec(row+column) + &
               invec(row-column))
       end do
    end do

    ! For the last degree rows, we need to put the second part to the front due to periodic boundaries
    do row = n_cells-degree+1, n_cells
       outvec(row) = mass(1)*invec(row)
       do column = 1,n_cells-row
          outvec(row) = outvec(row) + mass(column+1)*(invec(row+column) + &
               invec(row-column))
       end do
       do column = n_cells-row+1,degree
          outvec(row)  = outvec(row) + mass(column+1)*(invec(row+column-n_cells) + &
               invec(row-column))
       end do
    end do


  end subroutine sll_s_spline_fem_multiply_mass


!> Multiplication of the mix mass matrix given by a mass line \a mass
  subroutine sll_s_spline_fem_multiply_massmixed ( n_cells, degree, mass, invec, outvec )
    sll_int32,  intent( in    ) :: n_cells  !< no. of grid cells
    sll_int32,  intent( in    ) :: degree   !< spline degree
    sll_real64, intent( in    ) :: mass(degree*2) !< values of the mass line 
    sll_real64, intent( in    ) :: invec(:) !< input vector 
    sll_real64, intent( out   ) :: outvec(:)!< output vector

    !local variables
    sll_int32 :: row, column, ind

    ind = 1
    outvec = 0.0_f64
    ! For the first degree rows we need to put the first part to the back due to periodic boundaries
    do row = 1, degree-1
       do column=-row+1,degree
          outvec(row) = outvec(row) + mass(degree+column)*invec(row+column)
       end do
       do column = -degree+1,-row
          outvec(row) = outvec(row) + mass(degree+column)*invec(row+column+n_cells)
       end do
    end do

    do row = degree, n_cells-degree
       do column = -degree+1, degree
          outvec(row) = outvec(row) + mass(degree+column)*invec(row+column)
       end do
    end do

    ! For the last degree rows, we need to put the second part to the front due to periodic boundaries
    do row = n_cells-degree+1, n_cells
       do column = -degree+1, n_cells-row
          outvec(row) = outvec(row) + mass(degree+column)*invec(row+column)
       end do
       do column = n_cells-row+1, degree
          outvec(row) = outvec(row) + mass(degree+column)*invec(row+column-n_cells)
       end do
    end do


  end subroutine sll_s_spline_fem_multiply_massmixed


  !> Compute the eigenvalues of the interpolation matrix for splines of degree \a degree (with first grid point as starting point of the first spline)
  !> The interpolation is for the Greville points, i.e. the grid points for odd degree and the midpoints for
  !> even degree. The grid points are numbered as they appear in the grid and the first mid point is the one
  !> in the first cell.
  subroutine sll_s_spline_fem_interpolation_eigenvalues(degree, ndofs, eig)
    sll_int32, intent(in) :: degree   !< spline degree
    sll_int32, intent(in) :: ndofs    !< no. of degrees of freedom 
    sll_real64, intent(out) :: eig(ndofs) !< eigenvalues

    ! local variables
    sll_int32 :: k
    sll_int32 :: p
    sll_real64 :: spline_coeff(degree+1)
     sll_real64 :: ni !1/n_dofs in double precision
     ni=1.0_f64/real(ndofs,f64)

     if ( modulo( degree, 2) == 0 ) then
        call sll_s_uniform_bsplines_eval_basis( degree, 0.5_f64, spline_coeff )
     else
        call sll_s_uniform_bsplines_eval_basis( degree, 0.0_f64, spline_coeff )
     end if
     print*, 'spline_coeffs', spline_coeff

    eig(1) = 1.0_f64
    
    do k=1,(ndofs+1)/2-1
       ! real part
       eig(k+1) = spline_coeff(degree+1)
       do p=1,degree
          eig(k+1) = eig(k+1) + cos(sll_p_twopi*real(k*p,f64)*ni)*spline_coeff(degree+1-p)
       end do
       ! imaginary part
       eig(ndofs-k+1) = 0.0_f64
       do p=1,degree
          eig(ndofs-k+1) = eig(ndofs-k+1) - sin(sll_p_twopi*real(k*p,f64)*ni)*spline_coeff(degree+1-p)
       end do
    end do

    if ( modulo( ndofs, 2) == 0 ) then
       eig(ndofs/2+1) = spline_coeff(degree+1)
       do p=1,degree
          eig(ndofs/2+1) = eig(ndofs/2+1) + cos(sll_p_twopi*real(p,f64))*spline_coeff(degree+1-p)
       end do
    end if
     
  end subroutine sll_s_spline_fem_interpolation_eigenvalues

  !> Multiplication of the input vector \a in by the derivative matrix G
  subroutine sll_s_multiply_g_1d( n_dofs, delta_x,  in, out )
    sll_int32,  intent( in    ) :: n_dofs   !< no. of degrees of freedom
    sll_real64, intent( in    ) :: delta_x  !< grid spacing
    sll_real64, intent( in    ) :: in(:) !< Coefficient for each DoF (input vector)
    sll_real64, intent(   out ) :: out(:) !< Coefficient for each DoF (output vector)
    !local variables
    sll_int32 :: i
    
    ! treat Periodic point
    out(1) = ( in(1) - in(n_dofs)  )/delta_x
    do i = 2, n_dofs
       out(i) =  ( in(i) - in(i-1) )/delta_x
    end do
    
  end subroutine sll_s_multiply_g_1d



  !> Multiplication of the input vector \a in by the transposed derivative matrix G^T  
  subroutine sll_s_multiply_gt_1d( n_dofs, delta_x, in, out)
    sll_int32,  intent( in    ) :: n_dofs  !< no. of degrees of freedom
    sll_real64, intent( in    ) :: delta_x !< grid spacing
    sll_real64, intent( in    ) :: in(:) !< Coefficient for each DoF (input vector)
    sll_real64, intent(   out ) :: out(:) !< Coefficient for each DoF (output vector)
    !local variables
    sll_int32 :: i
    
    do i= 1, n_dofs-1
       out(i) =  ( in(i) - in(i+1)  )/delta_x
    end do
    ! treat Periodic point
    out(n_dofs) = ( in(n_dofs) - in(1) )/delta_x

  end subroutine sll_s_multiply_gt_1d


  
   !> Multiply by dicrete gradient matrix (3d version)
  subroutine sll_s_multiply_g( n_dofs, delta_x, field_in, field_out )
    sll_int32,  intent( in    )   :: n_dofs(3)   !< no. of degrees of freedom in each direction
    sll_real64, intent( in    )   :: delta_x(3)  !< grid spacing in each direction
    sll_real64, intent( in    )   :: field_in(:)  !< input vector  (size n1*n2*n3)
    sll_real64, intent(   out )   :: field_out(:)  !< output vector (size 3*n1*n2*n3)

    sll_real64 :: coef
    sll_int32  :: jump, jump_end
    sll_int32  :: ind3d, ind1d
    sll_int32  :: i,j,k

    coef = 1.0_f64/ delta_x(1)
    jump = 1
    jump_end = 1-n_dofs(1)

    ind3d = 0
    do k=1,n_dofs(3)
       do j=1,n_dofs(2)
          ind3d = ind3d + 1
          field_out(ind3d) = coef * ( field_in(ind3d) - field_in(ind3d-jump_end) )
          do i=2,n_dofs(1)
             ind3d = ind3d + 1
             field_out(ind3d) =  coef * ( field_in(ind3d) - field_in(ind3d-jump) )
          end do
       end do
    end do

    coef = 1.0_f64/ delta_x(2)
    jump = n_dofs(1)
    jump_end = (1-n_dofs(2))*n_dofs(1)

    ind1d = 0
    do k=1,n_dofs(3)
       do i=1,n_dofs(1)
          ind3d = ind3d + 1
          ind1d = ind1d + 1
          field_out(ind3d) =  coef * ( field_in(ind1d) - field_in(ind1d-jump_end) )
       end do
       do j=2,n_dofs(2)
          do i=1,n_dofs(1)
             ind3d = ind3d + 1
             ind1d = ind1d + 1
             field_out(ind3d) =  coef * ( field_in(ind1d) - field_in(ind1d-jump) )
          end do
       end do
    end do

    coef = 1.0_f64/ delta_x(3)
    jump = n_dofs(1)*n_dofs(2)
    jump_end = (1-n_dofs(3))*n_dofs(1)*n_dofs(2)

    ind1d = 0
    do j=1,n_dofs(2)
       do i=1,n_dofs(1)
          ind3d = ind3d + 1
          ind1d = ind1d + 1
          field_out(ind3d) = coef * ( field_in(ind1d) - field_in(ind1d-jump_end) )
       end do
    end do
    do k=2,n_dofs(3)
       do j=1,n_dofs(2)
          do i=1,n_dofs(1)
             ind3d = ind3d + 1
             ind1d = ind1d + 1
             field_out(ind3d) = coef * ( field_in(ind1d) - field_in(ind1d-jump) )
          end do
       end do
    end do


  end subroutine sll_s_multiply_g


  !> Multiply by transpose of dicrete gradient matrix (3d version)
  subroutine sll_s_multiply_gt( n_dofs, delta_x, field_in, field_out )
    sll_int32,  intent( in    )   :: n_dofs(3)  !< no. of degrees of freedom per direction
    sll_real64, intent( in    )   :: delta_x(3) !< grid spacing in each direction
    sll_real64, intent( in    )   :: field_in(:)  !< input vector
    sll_real64, intent(   out )   :: field_out(:)  !< output vector

    sll_real64 :: coef
    sll_int32  :: jump, jump_end
    sll_int32  :: ind3d, ind1d
    sll_int32  :: i,j,k

    coef = 1.0_f64/ delta_x(1)
    jump = 1
    jump_end = 1-n_dofs(1)

    ind3d = 0
    do k=1,n_dofs(3)
       do j=1,n_dofs(2)
          do i=1,n_dofs(1)-1
             ind3d = ind3d + 1
             field_out(ind3d) =  &
                  coef * ( field_in(ind3d) - field_in(ind3d+jump) )
          end do
          ind3d = ind3d + 1
          field_out(ind3d) = coef * ( field_in(ind3d) - field_in(ind3d+jump_end) )
       end do
    end do

    coef = 1.0_f64/ delta_x(2)
    jump = n_dofs(1)
    jump_end = (1-n_dofs(2))*n_dofs(1)

    ind1d = 0
    do k=1,n_dofs(3)
       do j=1,n_dofs(2)-1
          do i=1,n_dofs(1)
             ind3d = ind3d + 1
             ind1d = ind1d + 1

             field_out(ind1d) =  field_out(ind1d) + &
                  coef * ( field_in(ind3d) - field_in(ind3d+jump) )
          end do
       end do
       do i=1,n_dofs(1)
          ind3d = ind3d + 1
          ind1d = ind1d + 1

          field_out(ind1d) =  field_out(ind1d) + &
               coef * ( field_in(ind3d) - field_in(ind3d+jump_end) )
       end do
    end do

    coef = 1.0_f64/ delta_x(3)
    jump = n_dofs(1)*n_dofs(2)
    jump_end = (1-n_dofs(3))*n_dofs(1)*n_dofs(2)

    ind1d = 0
    do k=1,n_dofs(3)-1
       do j=1,n_dofs(2)
          do i=1,n_dofs(1)
             ind3d = ind3d + 1
             ind1d = ind1d + 1

             field_out(ind1d) =  field_out(ind1d) + &
                  coef * ( field_in(ind3d) - field_in(ind3d+jump) )
          end do
       end do
    end do
    do j=1,n_dofs(2)
       do i=1,n_dofs(1)
          ind3d = ind3d + 1
          ind1d = ind1d + 1

          field_out(ind1d) =   field_out(ind1d) + &
               coef * ( field_in(ind3d) - field_in(ind3d+jump_end) )
       end do
    end do

  end subroutine sll_s_multiply_gt

  
  !> Multiply by discrete curl matrix
  subroutine sll_s_multiply_c( n_dofs, delta_x, field_in, field_out )
    sll_int32,  intent( in    )   :: n_dofs(3)   !< no. of degrees of freedom per direction
    sll_real64, intent( in    )   :: delta_x(3)  !< grid spacing per direction
    sll_real64, intent( in    )   :: field_in(:) !< input vector
    sll_real64, intent(   out )   :: field_out(:)!< output vector  

    ! local variables
    sll_real64 :: coef(2)
    sll_int32  :: stride(2), jump(2), indp(2), n_total
    sll_int32  :: i,j,k, ind3d, ind3d_1, ind3d_2

    n_total = product(n_dofs)

    ! TODO: Avoid the IF for periodic boundaries
    ! First component
    coef(1) = 1.0_f64/ delta_x(2)
    coef(2) = -1.0_f64/ delta_x(3)

    stride(1) = n_dofs(1)
    stride(2) = n_dofs(1)*n_dofs(2)

    jump(1) = n_total
    jump(2) = 2*n_total

    ind3d = 0
    do k=1,n_dofs(3)
       if (k == 1) then
          indp(2) = stride(2)*(n_dofs(3)-1)
       else
          indp(2) = - stride(2)
       end if
       do j=1,n_dofs(2)
          if ( j==1 ) then
             indp(1) = stride(1)*(n_dofs(2)-1)
          else
             indp(1) = - stride(1)
          end if
          do i=1,n_dofs(1)
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
    coef(1) = 1.0_f64/ delta_x(3)
    coef(2) = -1.0_f64/ delta_x(1)

    stride(2) = 1
    stride(1) = n_dofs(1)*n_dofs(2)

    jump(1) = n_total
    jump(2) = -n_total

    do k=1,n_dofs(3)
       if (k == 1) then
          indp(1) = stride(1)*(n_dofs(3)-1)
       else
          indp(1) = - stride(1)
       end if
       do j=1,n_dofs(2)
          do i=1,n_dofs(1)
             if ( i==1 ) then
                indp(2) = stride(2)*(n_dofs(1)-1)
             else
                indp(2) = - stride(2)
             end if
             ind3d = ind3d + 1

             ind3d_1 = ind3d +indp(1)+jump(2)
             ind3d_2 = ind3d +indp(2)+jump(1)

             field_out(ind3d) = &
                  coef(1) * ( field_in( ind3d+jump(2) ) -&
                  field_in( ind3d_1 ) )+ &
                  coef(2) * ( field_in(ind3d+jump(1) ) - &
                  field_in( ind3d_2 ))
          end do
       end do
    end do

    ! Third component
    coef(1) = 1.0_f64/ delta_x(1)
    coef(2) = -1.0_f64/ delta_x(2)

    stride(1) = 1
    stride(2) = n_dofs(1)

    jump(1) = -2*n_total
    jump(2) = -n_total

    do k=1,n_dofs(3)
       do j=1,n_dofs(2)
          if (j == 1) then
             indp(2) = stride(2)*(n_dofs(2)-1)
          else
             indp(2) = - stride(2)
          end if
          do i=1,n_dofs(1)
             if ( i==1 ) then
                indp(1) = stride(1)*(n_dofs(1)-1)
             else
                indp(1) = - stride(1)
             end if
             ind3d = ind3d + 1

             ind3d_1 = ind3d +indp(1)+jump(2)
             ind3d_2 = ind3d +indp(2)+jump(1)

             field_out(ind3d) =  &
                  coef(1) * ( field_in( ind3d+jump(2) ) -&
                  field_in( ind3d_1 ) )+ &
                  coef(2) * ( field_in(ind3d+jump(1) ) - &
                  field_in( ind3d_2 ) )
          end do
       end do
    end do

  end subroutine sll_s_multiply_c


  !> Multiply by transpose of discrete curl matrix
  subroutine sll_s_multiply_ct( n_dofs, delta_x, field_in, field_out )
    sll_int32,  intent( in    )   :: n_dofs(3)    !< no. of degrees of freedom per direction
    sll_real64, intent( in    )   :: delta_x(3)   !< grid spacing per direction
    sll_real64, intent( in    )   :: field_in(:)  !< Matrix to be multiplied
    sll_real64, intent(   out )   :: field_out(:) !< C*field_in

    ! local variables
    sll_real64 :: coef(2)
    sll_int32  :: stride(2), jump(2), indp(2), n_total
    sll_int32  :: i,j,k, ind3d, ind3d_1, ind3d_2

    n_total = product(n_dofs)
    ! TODO: Avoid the IF for periodic boundaries
    ! First component
    coef(1) = -1.0_f64/ delta_x(2)
    coef(2) = 1.0_f64/ delta_x(3)

    stride(1) = n_dofs(1)
    stride(2) = n_dofs(1)*n_dofs(2)

    jump(1) = n_total
    jump(2) = 2*n_total

    ind3d = 0
    do k=1,n_dofs(3)
       if (k == n_dofs(3)) then
          indp(2) = -stride(2)*(n_dofs(3)-1)
       else
          indp(2) = stride(2)
       end if
       do j=1,n_dofs(2)
          if ( j== n_dofs(2)) then
             indp(1) = -stride(1)*(n_dofs(2)-1)
          else
             indp(1) = stride(1)
          end if
          do i=1,n_dofs(1)
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

    stride(2) = 1
    stride(1) = n_dofs(1)*n_dofs(2)

    jump(1) = n_total
    jump(2) = -n_total

    do k=1,n_dofs(3)
       if (k == n_dofs(3)) then
          indp(1) = -stride(1)*(n_dofs(3)-1)
       else
          indp(1) = stride(1)
       end if
       do j=1,n_dofs(2)
          do i=1,n_dofs(1)
             if ( i==n_dofs(1) ) then
                indp(2) = -stride(2)*(n_dofs(1)-1)
             else
                indp(2) = stride(2)
             end if
             ind3d = ind3d + 1

             ind3d_1 = ind3d +indp(1)+jump(2)
             ind3d_2 = ind3d +indp(2)+jump(1)

             field_out(ind3d) = &
                  coef(1) * ( field_in( ind3d+jump(2) ) -&
                  field_in( ind3d_1 ) )+ &
                  coef(2) * ( field_in(ind3d+jump(1) ) - &
                  field_in( ind3d_2 ))
          end do
       end do
    end do

    ! Third component
    coef(1) = -1.0_f64/ delta_x(1)
    coef(2) = 1.0_f64/ delta_x(2)

    stride(1) = 1
    stride(2) = n_dofs(1)

    jump(1) = -2*n_total
    jump(2) = -n_total

    do k=1,n_dofs(3)
       do j=1,n_dofs(2)
          if (j == n_dofs(2)) then
             indp(2) = -stride(2)*(n_dofs(2)-1)
          else
             indp(2) = stride(2)
          end if
          do i=1,n_dofs(1)
             if ( i==n_dofs(1) ) then
                indp(1) = -stride(1)*(n_dofs(1)-1)
             else
                indp(1) = stride(1)
             end if
             ind3d = ind3d + 1

             ind3d_1 = ind3d +indp(1)+jump(2)
             ind3d_2 = ind3d +indp(2)+jump(1)

             field_out(ind3d) =  &
                  coef(1) * ( field_in( ind3d+jump(2) ) -&
                  field_in( ind3d_1 ) )+ &
                  coef(2) * ( field_in(ind3d+jump(1) ) - &
                  field_in( ind3d_2 ) )
          end do
       end do
    end do

  end subroutine sll_s_multiply_ct


end module sll_m_spline_fem_utilities
