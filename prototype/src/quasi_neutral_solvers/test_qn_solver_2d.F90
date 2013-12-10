program test_qn_solver_2d

#include "sll_working_precision.h"
#include "sll_memory.h"
#include "sll_assert.h"
#include "sll_constants.h"

  use sll_qn_solver_2d

  implicit none
  
  sll_int32                    :: BC ! Boundary_conditions
  sll_int32                             :: NP_r, NP_theta
  ! NP_r and NP_theta are the numbers of points in directions r and 
  ! theta respectively
  sll_real64                            :: rmin, rmax, Zi
  sll_real64, dimension(:), allocatable :: Te
  sll_int32                             :: ierr, i

  NP_r = 256
  NP_theta = 256
  rmin = 1.d0
  rmax = 10.d0
  Zi = 1.d0

  do i=1,2 ! 2 test functions

     if (i==1) then
        BC = SLL_NEUMANN
     else
        BC = SLL_DIRICHLET
     endif
     print*, ' '
     print*, 'Testing sll_qns2d_angular_spect_method_seq with ', BC
     print*, ' '
     SLL_ALLOCATE(Te(NP_r), ierr)
     Te = 1.d0
     call test_process(BC, NP_r, NP_theta, rmin, rmax, Te, Zi)
     SLL_DEALLOCATE_ARRAY(Te, ierr)

  enddo

  print*, 'test_sll_qns2d_angular_spect_method_seq: PASSED'
  print*, ' '

contains


  subroutine test_process(BC, NP_r, NP_theta, rmin, rmax, Te, Zi)

    sll_int32                               :: BC ! Boundary_conditions
    sll_int32                               :: NP_r, NP_theta
    ! NP_r and NP_theta are the numbers of points in directions r and 
    ! theta respectively
    sll_real64                              :: rmin, rmax, Zi
    sll_real64, dimension(:)                :: Te
    sll_real64, dimension(:),   allocatable :: c, f, g
    sll_int32                               :: ierr
    sll_real64                              :: dr, dtheta
    sll_real64                              :: r, theta
    sll_real64, dimension(NP_r,NP_theta)    :: rho, phi
    sll_real64, dimension(NP_r,NP_theta)    :: phi_exact
    sll_int32                               :: i, j, i_test
    type (qn_solver_2d), pointer            :: plan
    sll_real64                              :: average_err
    sll_real64                              :: average_err_bound
    sll_real64                              :: Mr, Mtheta

    if (BC==SLL_NEUMANN) then
       dr = (rmax-rmin)/(NP_r-1)
    else 
       dr = (rmax-rmin)/(NP_r+1)
    endif

    dtheta = 2*sll_pi/NP_theta

    SLL_ALLOCATE(c(NP_r), ierr)
    SLL_ALLOCATE(f(NP_theta), ierr)
    SLL_ALLOCATE(g(NP_theta), ierr)

    plan => new(BC, rmin, rmax, NP_r, NP_theta)

    do i_test=1,2

       f = 0.d0
       average_err_bound  = 0.d0

       do j=1,NP_theta

          theta = (j-1)*dtheta
          Mr = 4*abs(cos(theta))

          if (BC==SLL_NEUMANN) then
             if (i_test==1) then
                f(j) = sin(rmax-rmin)*cos(theta)
             else
                f(j)= sin(rmax-rmin) * exp(-.5*(theta-sll_pi)**2)/sqrt(2*sll_pi)
             endif
          endif

          do i=1,NP_r
             if (BC==SLL_NEUMANN) then
                r = rmin + (i-1)*dr
                c(i) = 2/r
             else ! 'dirichlet'
                r = rmin + i*dr
                c(i) = (rmax+rmin-2*r) / ( (rmax-r)*(r-rmin) )
             endif
          
             if (i_test==1) then

                phi_exact(i,j)  = sin(r-rmin)*sin(rmax-r)*cos(theta)

                rho(i,j) = cos(theta) * ( 2*cos(rmin+rmax-2*r) - &
                                      c(i)* sin(rmin+rmax-2*r) + &
                    (1/r**2+1/(Zi*Te(i)))*sin(rmax-r)*sin(r-rmin))
             else
               phi_exact(i,j)  = sin(r-rmin)*sin(rmax-r) * &
                   exp(-.5*(theta-sll_pi)**2)/sqrt(2*sll_pi)

               rho(i,j) = ( 2*cos(rmax+rmin-2*r) - c(i)*sin(rmax+rmin-2*r) ) * &
                          exp(-.5*(theta-sll_pi)**2)/sqrt(2*sll_pi) + phi_exact(i,j) &
                                       * ( 1/(Zi*Te(i)) - ((theta-sll_pi)**2-1)/r**2 )
             endif

             Mtheta = abs(sin(r-rmin)*sin(rmax-r))
             average_err_bound = average_err_bound + &
             Mr*dr**2/12 + abs(c(i))*Mr*dr**2/6 + Mtheta*dtheta**2/(r**2*12)
          enddo

       enddo

       g = -f

       call solve(plan, rho, c, Te, f, g, Zi, phi)

       average_err = sum(abs(phi_exact-phi))/(NP_r*NP_theta)
       average_err_bound = average_err_bound/(NP_r*NP_theta)

       print*, 'Error =', average_err
       print*, 'Boundary error =', average_err_bound
       print*, ' '

       if ( average_err > average_err_bound ) then
          print*, 'test_sll_qns2d_angular_spect_method_seq: FAILED'
          print*, ' '
          stop
       endif
    
    enddo

    call delete(plan)

    SLL_DEALLOCATE_ARRAY(c, ierr)
    SLL_DEALLOCATE_ARRAY(f, ierr)
    SLL_DEALLOCATE_ARRAY(g, ierr)
 
  end subroutine test_process


end program test_qn_solver_2d
