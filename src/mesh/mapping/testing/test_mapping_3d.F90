!>$L_x$ domain dimensions and M is an integer.
!>$
!>$
!>$
!>B_z(x,y,t) =   \cos(\frac{2 M \pi}{L_x} x)  \cos(\frac{2 M \pi}{L_x} t)
!>$
!>$
!>E_y(x,y,t) = \sin(\frac{2 M \pi}{L_x} x)  \sin(\frac{2 M \pi}{L_x} t)
!>$
!
!  Contact : Benedikt Perse
!
program test_mapping_3d
  !------------------------------------------------------------------------
  !  test 3D Maxwell spline finite element solver on a periodic grid
  !------------------------------------------------------------------------
  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#include "sll_memory.h"
#include "sll_working_precision.h"
#include "sll_maxwell_solvers_macros.h"

  use sll_m_constants, only: &
       sll_p_pi, sll_p_twopi, sll_p_fourpi

  use sll_m_mapping_3d, only: &
       sll_t_mapping_3d

  use sll_m_gauss_legendre_integration, only: &
       sll_f_gauss_legendre_points_and_weights

  use sll_m_3d_coordinate_transformations

  use sll_m_timer, only: &
       sll_s_set_time_mark, &
       sll_f_time_elapsed_between, &
       sll_t_time_mark

  implicit none
  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  !  type(sll_t_arbitrary_degree_spline_1d), pointer :: aspl
  !  sll_real64, dimension(:), allocatable :: knots

  sll_real64                      :: eta1_max, eta1_min
  sll_real64                      :: delta_eta(3)
  sll_int32                       :: nc_eta(3), nc_total
  type(sll_t_mapping_3d)          :: map
  sll_real64        :: x, y, z
  sll_real64         :: error1, error2, error3, error4
  sll_int32                       :: i, j, k, c1, c2, c3
  sll_int32                       :: deg(3)
  sll_real64                      :: matrix(3,3)
  sll_real64                      :: val(3)
  sll_real64                      :: params(6), x1, x2, maxerror, maxerror1, maxerror2, maxerror3, maxerror4
  sll_real64, allocatable :: g1(:,:), g2(:,:), g3(:,:)
  type(sll_t_time_mark) :: start, end, start1, end1

  params = 0._f64
  params(1) = 0._f64
  params(2) = 1._f64!sll_p_twopi!
  params(3) = 1._f64!sll_p_twopi!

!!$  params(1)= sll_p_twopi
!!$  params(2)= sll_p_twopi
!!$  params(3)= sll_p_twopi
!!$  params(4)= 0.1_f64
!!$  params(5)= 0.1_f64


  ! Define computational domain in physical coordinates
  eta1_min = 0._f64
  eta1_max = 1._f64
  nc_eta = [32, 512, 2 ]
  nc_total = product(nc_eta+1)
  delta_eta = 1._f64/real(nc_eta,f64)
  ! Set spline degree of 0-forms
  deg = [3, 3, 1]

  call sll_s_set_time_mark( start )
  call map%init(params,&
       sll_f_cylindrical_x1,&
       sll_f_cylindrical_x2,&
       sll_f_cylindrical_x3,&
       sll_f_cylindrical_jac11,&
       sll_f_cylindrical_jac12,&
       sll_f_cylindrical_jac13,&
       sll_f_cylindrical_jac21,&
       sll_f_cylindrical_jac22,&
       sll_f_cylindrical_jac23,&
       sll_f_cylindrical_jac31,&
       sll_f_cylindrical_jac32,&
       sll_f_cylindrical_jac33,&
       sll_f_cylindrical_jacobian, &
       flag2d = .true., n_cells=nc_eta, deg=deg )
  print*,'Initialisation finished'
  call sll_s_set_time_mark( end )
  write(*, "(A, F10.3)") "Mapping init run time [s] = ", sll_f_time_elapsed_between( start, end)
  call sll_s_set_time_mark( start )


  error1 = 0._f64
  error2 = 0._f64
  error3 = 0._f64
  error4 = 0._f64

  allocate(g3(2,deg(3)+1))
  allocate(g2(2,deg(2)+1))
  allocate(g1(2,deg(1)+1))



  g1 = sll_f_gauss_legendre_points_and_weights(deg(1)+1, 0.0_f64, 1.0_f64)
  g2 = sll_f_gauss_legendre_points_and_weights(deg(2)+1, 0.0_f64, 1.0_f64)
  g3 = sll_f_gauss_legendre_points_and_weights(deg(3)+1, 0.0_f64, 1.0_f64)

  call sll_s_set_time_mark( start1 )
  do k = 1, nc_eta(3)
     do j = 1, nc_eta(2)  
        do i = 1, nc_eta(1)
           do c3 = 1, deg(3)+1
              do c2 = 1, deg(2)+1
                 do c1 = 1, deg(1)+1
                    z = delta_eta(3) * (g3(1,c3) + real(k-1,f64))
                    y = delta_eta(2) * (g2(1,c2) + real(j-1,f64))
                    x = delta_eta(1) * (g1(1,c1) + real(i-1,f64))


                    x1 = map%get_x1( [x, y, z] )
                    x2 = map%x1_func( [x, y, z], map%params )
                    error1 = max(abs(x1 - x2), error1)
                    !print*, x1, x2, 'error1', x1-x2

                    x1 = map%get_x2( [x, y, z] )
                    x2 = map%x2_func( [x, y, z], map%params )
                    error2 = max(abs(x1 - x2), error2)
                    !print*,  x1, x2,'error2', abs(x1-x2)

                    x1 = map%get_x3( [x, y, z] )
                    x2 = map%x3_func( [x, y, z], map%params )
                    error3 = max(abs(x1 - x2), error3)
                    !print*, x1, x2, 'error3', x1-x2

                    x1 = map%jacob( [x, y, z], map%params )
                    x2 = map%jacobian( [x, y, z] )
                    error4 = max(abs(x1 - x2), error4)
                    !print*, x1, x2,  'error4', abs(x1-x2) 
                 end do
              end do
           end do
        end do
     end do
  end do
  call sll_s_set_time_mark( end1 )
  write(*, "(A, F10.3)") "Error run time [s] = ", sll_f_time_elapsed_between( start1, end1)
  maxerror1 = (error1)
  maxerror2 = (error2)
  maxerror3 = (error3)
  maxerror4 = (error4)
  maxerror= maxval( [maxerror1, maxerror2, maxerror3, maxerror4] )
  print*, 'interpolation error', maxerror1, maxerror2,  maxerror3,  maxerror4
  call sll_s_set_time_mark( start1 )
  do k = 1, 2
     do j = 1, 10  
        do i = 1, 10 
           val = map%get_x( [x, y, z] )
        end do
     end do
  end do
  call sll_s_set_time_mark( end1 )
  write(*, "(A, F10.3)") "get_x run time [s] = ", sll_f_time_elapsed_between( start1, end1)

  call sll_s_set_time_mark( start1 )
  do k = 1, nc_eta(3)
     do j = 1, nc_eta(2)  
        do i = 1, nc_eta(1)
           do c3 = 1, deg(3)+1
              do c2 = 1, deg(2)+1
                 do c1 = 1, deg(1)+1
                    z = delta_eta(3) * (g3(1,c3) + real(k-1,f64))
                    y = delta_eta(2) * (g2(1,c2) + real(j-1,f64))
                    x = delta_eta(1) * (g1(1,c1) + real(i-1,f64))
                    matrix =  map%jacobian_matrix( [x, y, z] )
!!$           x1 = map%j_matrix(1,1)%f( [x, y, z], map%params )
!!$           x2 = map%j_matrix(2,2)%f( [x, y, z], map%params )
!!$           print*, 'matrix11', matrix(1,1), x1, 'error' , abs(matrix(1,1)- x1)
!!$           print*, 'matrix22', matrix(2,2), x2, 'error' , abs(matrix(2,2)- x2)
                 end do
              end do
           end do
        end do
     end do
  end do
  call sll_s_set_time_mark( end1 )
  write(*, "(A, F10.3)") "jacobian matrix run time [s] = ", sll_f_time_elapsed_between( start1, end1)

  call sll_s_set_time_mark( start1 )
  do k = 1, nc_eta(3)
     do j = 1, nc_eta(2)  
        do i = 1, nc_eta(1)
           do c3 = 1, deg(3)+1
              do c2 = 1, deg(2)+1
                 do c1 = 1, deg(1)+1
                    z = delta_eta(3) * (g3(1,c3) + real(k-1,f64))
                    y = delta_eta(2) * (g2(1,c2) + real(j-1,f64))
                    x = delta_eta(1) * (g1(1,c1) + real(i-1,f64))
                    val(1) =  map%jacobian( [x, y, z] )
                 end do
              end do
           end do
        end do
     end do
  end do
  call sll_s_set_time_mark( end1 )
  write(*, "(A, F10.3)") "jacobian run time [s] = ", sll_f_time_elapsed_between( start1, end1)


  call sll_s_set_time_mark( start1 )
  do k = 1, nc_eta(3)
     do j = 1, nc_eta(2)  
        do i = 1, nc_eta(1)
           do c3 = 1, deg(3)+1
              do c2 = 1, deg(2)+1
                 do c1 = 1, deg(1)+1
                    z = delta_eta(3) * (g3(1,c3) + real(k-1,f64))
                    y = delta_eta(2) * (g2(1,c2) + real(j-1,f64))
                    x = delta_eta(1) * (g1(1,c1) + real(i-1,f64))
                    matrix =  map%jacobian_matrix_transposed( [x, y, z] )
                 end do
              end do
           end do
        end do
     end do
  end do
  call sll_s_set_time_mark( end1 )
  write(*, "(A, F10.3)") "jacobian matrix transposed run time [s] = ", sll_f_time_elapsed_between( start1, end1)

  call sll_s_set_time_mark( start1 )
  do k = 1, nc_eta(3)
     do j = 1, nc_eta(2)  
        do i = 1, nc_eta(1)
           do c3 = 1, deg(3)+1
              do c2 = 1, deg(2)+1
                 do c1 = 1, deg(1)+1
                    z = delta_eta(3) * (g3(1,c3) + real(k-1,f64))
                    y = delta_eta(2) * (g2(1,c2) + real(j-1,f64))
                    x = delta_eta(1) * (g1(1,c1) + real(i-1,f64))
                    matrix = map%jacobian_matrix_inverse( [x, y, z] )
                 end do
              end do
           end do
        end do
     end do
  end do
  call sll_s_set_time_mark( end1 )
  write(*, "(A, F10.3)") "jacobian matrix inverse run time [s] = ", sll_f_time_elapsed_between( start1, end1)

  call sll_s_set_time_mark( start1 )
  do k = 1, nc_eta(3)
     do j = 1, nc_eta(2)  
        do i = 1, nc_eta(1)
           do c3 = 1, deg(3)+1
              do c2 = 1, deg(2)+1
                 do c1 = 1, deg(1)+1
                    z = delta_eta(3) * (g3(1,c3) + real(k-1,f64))
                    y = delta_eta(2) * (g2(1,c2) + real(j-1,f64))
                    x = delta_eta(1) * (g1(1,c1) + real(i-1,f64))
                    matrix = map%jacobian_matrix_inverse_transposed( [x, y, z] )
                 end do
              end do
           end do
        end do
     end do
  end do
  call sll_s_set_time_mark( end1 )
  write(*, "(A, F10.3)") "jacobian matrix inverse transposed run time [s] = ", sll_f_time_elapsed_between( start1, end1)

  call sll_s_set_time_mark( start1 )
  do k = 1, nc_eta(3)
     do j = 1, nc_eta(2)  
        do i = 1, nc_eta(1)
           do c3 = 1, deg(3)+1
              do c2 = 1, deg(2)+1
                 do c1 = 1, deg(1)+1
                    z = delta_eta(3) * (g3(1,c3) + real(k-1,f64))
                    y = delta_eta(2) * (g2(1,c2) + real(j-1,f64))
                    x = delta_eta(1) * (g1(1,c1) + real(i-1,f64))
                    matrix = map%metric( [x, y, z] )/map%jacobian( [x, y, z] )
                 end do
              end do
           end do
        end do
     end do
  end do
  call sll_s_set_time_mark( end1 )
  write(*, "(A, F10.3)") "metric run time [s] = ", sll_f_time_elapsed_between( start1, end1)

  call sll_s_set_time_mark( start1 )
  do k = 1, nc_eta(3)
     do j = 1, nc_eta(2)  
        do i = 1, nc_eta(1)
           do c3 = 1, deg(3)+1
              do c2 = 1, deg(2)+1
                 do c1 = 1, deg(1)+1
                    z = delta_eta(3) * (g3(1,c3) + real(k-1,f64))
                    y = delta_eta(2) * (g2(1,c2) + real(j-1,f64))
                    x = delta_eta(1) * (g1(1,c1) + real(i-1,f64))
                    val(1) = map%metric_single_jacobian( [x, y, z], 2, 2 )
                 end do
              end do
           end do
        end do
     end do
  end do
  call sll_s_set_time_mark( end1 )
  write(*, "(A, F10.3)") "metric_single run time [s] = ", sll_f_time_elapsed_between( start1, end1)

  call sll_s_set_time_mark( start1 )
  do k = 1, nc_eta(3)
     do j = 1, nc_eta(2)  
        do i = 1, nc_eta(1)
           do c3 = 1, deg(3)+1
              do c2 = 1, deg(2)+1
                 do c1 = 1, deg(1)+1
                    z = delta_eta(3) * (g3(1,c3) + real(k-1,f64))
                    y = delta_eta(2) * (g2(1,c2) + real(j-1,f64))
                    x = delta_eta(1) * (g1(1,c1) + real(i-1,f64))
                    matrix = map%metric_inverse( [x, y, z] )*map%jacobian( [x, y, z] )
                 end do
              end do
           end do
        end do
     end do
  end do
  call sll_s_set_time_mark( end1 )
  write(*, "(A, F10.3)") "metric_inverse run time [s] = ", sll_f_time_elapsed_between( start1, end1)

  call sll_s_set_time_mark( start1 )
  do k = 1, nc_eta(3)
     do j = 1, nc_eta(2)  
        do i = 1, nc_eta(1)
           do c3 = 1, deg(3)+1
              do c2 = 1, deg(2)+1
                 do c1 = 1, deg(1)+1
                    z = delta_eta(3) * (g3(1,c3) + real(k-1,f64))
                    y = delta_eta(2) * (g2(1,c2) + real(j-1,f64))
                    x = delta_eta(1) * (g1(1,c1) + real(i-1,f64))
                    val(1) = map%metric_inverse_single_jacobian( [x, y, z], 1, 1 )
                 end do
              end do
           end do
        end do
     end do
  end do
  call sll_s_set_time_mark( end1 )
  write(*, "(A, F10.3)") "metric_inverse_single run time [s] = ", sll_f_time_elapsed_between( start1, end1)

  ! Clean up
  call map%free()

  if (  maxerror<6.d-9) then
     print*, 'PASSED.'
  else
     print*, 'FAILED.'
  end if
  call sll_s_set_time_mark( end )
  write(*, "(A, F10.3)") "Main part run time [s] = ", sll_f_time_elapsed_between( start, end)

end program test_mapping_3d
