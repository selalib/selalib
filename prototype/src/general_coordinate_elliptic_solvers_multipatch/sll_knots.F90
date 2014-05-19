module sll_knots
  implicit none

  integer, parameter :: KNOTS_PERIODIC = 0, KNOTS_DIRICHLET = 1

contains

  ! SUBROUTINE ASSEMBLING KNOTS ARRAY
  subroutine initialize_knots(&
       spline_degree,&
       num_cells, &
       eta_min, &
       eta_max, &
       bc_left, &
       bc_right, &
       knots)
    
    integer, intent(in)   :: spline_degree 
    integer, intent(in)   :: num_cells
    real(8), intent(in)   :: eta_min
    real(8), intent(in)   :: eta_max
    integer, intent(in)   :: bc_left
    integer, intent(in)   :: bc_right
    real(8), dimension(:) :: knots
    integer               :: i, k 
    real(8)               :: eta
    real(8)               :: delta

    ! This routine needs some error-checking...
    delta = (eta_max - eta_min)/num_cells

    if ( (bc_left  == KNOTS_PERIODIC) .and. &
         (bc_right == KNOTS_PERIODIC) ) then 
       ! it is intentional that eta_min is not used, one intends to consider
       ! only the [0,1] interval...
       do k = -spline_degree, spline_degree+1
          knots(k+spline_degree+1) = delta*k 
       end do
    else if ( (bc_left  == KNOTS_DIRICHLET) .and. &
              (bc_right == KNOTS_DIRICHLET) ) then 
       do i = 1, spline_degree + 1
          knots(i) = eta_min
       enddo
       eta = eta_min
       do i = spline_degree + 2, num_cells + 1 + spline_degree
          eta = eta + delta
          knots(i) = eta
       enddo
       do i = num_cells + spline_degree + 1, num_cells + 1 + 2*spline_degree
          knots(i) = eta_max
       enddo
    end if

  end subroutine initialize_knots


  subroutine initialize_knots_per(&
       spline_degree,&
       num_cells, &
       eta_min, &
       eta_max, &
       knots)

    integer, intent(in)   :: spline_degree 
    integer, intent(in)   :: num_cells
    real(8), intent(in)   :: eta_min
    real(8), intent(in)   :: eta_max
    real(8), dimension(:) :: knots
    integer               :: i
    real(8)               :: delta

    delta = (eta_max - eta_min)/num_cells
    
    do i = - spline_degree, num_cells + spline_degree
       
       knots( i+ spline_degree + 1 ) = eta_min + i* delta
    end do
  end subroutine initialize_knots_per

  subroutine initialize_knots_dir(&
       spline_degree,&
       num_cells, &
       eta_min, &
       eta_max, &
       knots)
    
    integer, intent(in)   :: spline_degree 
    integer, intent(in)   :: num_cells
    real(8), intent(in)   :: eta_min
    real(8), intent(in)   :: eta_max
    real(8), dimension(:) :: knots
    integer               :: i
    real(8)               :: delta
    real(8)               :: eta
    
    delta = (eta_max - eta_min)/num_cells
    
    do i = 1, spline_degree + 1
       knots(i) = eta_min
    enddo
    eta = eta_min
    do i = spline_degree + 2, num_cells + 1 + spline_degree
       eta = eta + delta
       knots(i) = eta
    enddo
    do i = num_cells + spline_degree + 1, num_cells + 1 + 2*spline_degree
       knots(i) = eta_max
    enddo
  end subroutine initialize_knots_dir


  subroutine initialize_knots_all(&
       spline_degree,&
       num_cells, &
       eta_min, &
       eta_max, &
       bc_left, &
       bc_right, &
       knots)
    
    integer, intent(in)   :: spline_degree 
    integer, intent(in)   :: num_cells
    real(8), intent(in)   :: eta_min
    real(8), intent(in)   :: eta_max
    integer, intent(in)   :: bc_left
    integer, intent(in)   :: bc_right
    real(8), dimension(:) :: knots
    integer               :: i, k 
    real(8)               :: eta
    real(8)               :: delta
    
    ! This routine needs some error-checking...
    delta = (eta_max - eta_min)/num_cells
    
    if ( (bc_left  == KNOTS_PERIODIC) .and. &
         (bc_right == KNOTS_PERIODIC) ) then 
      
       call initialize_knots_per(&
            spline_degree,&
            num_cells, &
            eta_min, &
            eta_max, &
            knots)
    else if ( (bc_left  == KNOTS_DIRICHLET) .and. &
         (bc_right == KNOTS_DIRICHLET) ) then 
       call initialize_knots_dir(&
            spline_degree,&
            num_cells, &
            eta_min, &
            eta_max, &
            knots)
       
    end if
    
  end subroutine initialize_knots_all
  
end module sll_knots
