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
!  Contact : Eric Sonnendrucker, Katharina Kormann
!
program test_curl_curl_part
  !------------------------------------------------------------------------
  !  test 3D Maxwell spline finite element solver on a periodic grid
  !------------------------------------------------------------------------
  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#include "sll_memory.h"
#include "sll_working_precision.h"
#include "sll_maxwell_solvers_macros.h"

  use sll_m_low_level_bsplines, only: &
       sll_s_eval_uniform_periodic_spline_curve_with_zero

  use sll_m_constants, only: &
       sll_p_pi

  use sll_m_maxwell_3d_fem, only: &
       sll_t_maxwell_3d_fem

  use sll_m_maxwell_1d_base, only: &
       sll_s_plot_two_fields_1d

  use sll_m_constants, only: sll_p_pi, sll_p_twopi

  use sll_m_timer, only: &
       sll_s_set_time_mark, &
       sll_f_time_elapsed_between, &
       sll_t_time_mark

  implicit none
  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  !  type(sll_t_arbitrary_degree_spline_1d), pointer :: aspl
  !  sll_real64, dimension(:), allocatable :: knots

  sll_real64 :: eta1_max, eta1_min
  sll_real64 :: delta_eta(3)

  sll_int32  :: nc_eta(3), nc_total

  type(sll_t_maxwell_3d_fem)  :: maxwell_3d
  sll_real64, allocatable :: x(:,:,:), y(:,:,:), z(:,:,:)
  sll_real64, allocatable :: efield(:), bfield(:)
  sll_real64, allocatable :: efield_ref(:), bfield_ref(:)
  sll_real64, allocatable :: efield_val(:), bfield_val(:)
  sll_real64, allocatable :: rho(:), rho_ref(:), current(:)

  sll_int32                               :: i, j, k, istep, nsteps, ind
  sll_real64                                :: w1, w2
  sll_real64                                :: time
  sll_real64                                :: delta_t

  sll_real64                              :: Lx(3)
  sll_real64, dimension(3,2)              :: domain
  sll_int32                               :: deg(3)

  sll_real64                              :: error(6)
  sll_real64                              :: energy(2)
  type(sll_t_time_mark) :: start, end

  call sll_s_set_time_mark( start )

  ! Define computational domain
  eta1_min = .0_f64; eta1_max = 2.0_f64*sll_p_pi
  nc_eta = [32, 16, 32]
  nc_total = product(nc_eta)
  Lx(1) = eta1_max-eta1_min
  Lx(2) = Lx(1); Lx(3) = Lx(2)
  delta_eta = Lx/real(nc_eta,f64)
  domain(1,:) = [eta1_min, eta1_max]
  domain(2,:) = [eta1_min, eta1_max]
  domain(3,:) = [eta1_min, eta1_max]
  ! Set spline degree of 0-forms
  deg =3! [2,2,3]
  ! Time loop
  delta_t = 0.01_f64
  nsteps = 3
  w1 = sqrt(3.0_f64)
  w2 = + sqrt(3.0_f64)
  time = 0.0_f64
  ! Initialise maxwell FEM object
  call maxwell_3d%init(domain, nc_eta, deg)


  allocate(x(1:nc_eta(1)+1,1:nc_eta(2)+1,1:nc_eta(3)+1) )
  allocate(y(1:nc_eta(1)+1,1:nc_eta(2)+1,1:nc_eta(3)+1) )
  allocate(z(1:nc_eta(1)+1,1:nc_eta(2)+1,1:nc_eta(3)+1) )
  allocate(efield(1:nc_total*3))
  allocate(bfield(1:nc_total*3))
  allocate(efield_ref(1:nc_total*3))
  allocate(bfield_ref(1:nc_total*3))
  allocate(efield_val(1:nc_total*3))
  allocate(bfield_val(1:nc_total*3))
  allocate( rho(1:nc_total) )
  allocate( rho_ref( 1:nc_total) )
  allocate(current(1:nc_total*3))

  do k = 1, nc_eta(3)+1
     do j = 1, nc_eta(2)+1
        do i = 1, nc_eta(1)+1
           x(i,j,k) = eta1_min + real(i-1,f64)*delta_eta(1)
           y(i,j,k) = eta1_min + real(j-1,f64)*delta_eta(2)
           z(i,j,k) = eta1_min + real(k-1,f64)*delta_eta(3)
        end do
     end do
  end do


  current = 0._f64
!!$  call maxwell_3d%compute_rhs_from_function( 1, 1, current(1:nc_total), cos_k )
!!$  call maxwell_3d%compute_rhs_from_function( 1, 3, current(nc_total*2+1:nc_total*3), cos_k )
!!$  current(1:nc_total) = 3._f64* current(1:nc_total)
!!$  current(nc_total*2+1:nc_total*3) = -3._f64*current(nc_total*2+1:nc_total*3)
  call maxwell_3d%compute_rhs_from_function( 1, 1, current(1:nc_total), jx )
  call maxwell_3d%compute_rhs_from_function( 1, 2, current(nc_total+1:nc_total*2), jy )
  call maxwell_3d%compute_rhs_from_function( 1, 3, current(nc_total*2+1:nc_total*3), jz )
!!$  call random_number( efield )
!!$  call maxwell_3d%multiply_ct( efield, current )
  !call random_number( current )
  rho = 0._f64
  call maxwell_3d%multiply_gt( current, rho )
  print*, 'rhs divergence', maxval(abs(rho))

 ! maxwell_3d%curl_matrix%epsilon = delta_t
  !call maxwell_3d%curl_solver%solve( current, efield )
  call maxwell_3d%uzawa_iterator%solve( current, efield )

  efield_ref = 0._f64
!!$  call maxwell_3d%L2projection( 1, 1, efield_ref(1:nc_total), cos_k )
!!$  call maxwell_3d%L2projection( 1, 3, efield_ref(nc_total*2+1:nc_total*3), cos_k )
!!$  efield_ref(nc_total*2+1:nc_total*3) = -efield_ref(nc_total*2+1:nc_total*3)

  call maxwell_3d%L2projection( 1, 1, efield_ref(1:nc_total), ax )
  call maxwell_3d%L2projection( 1, 2, efield_ref(nc_total+1:nc_total*2), ay )
  call maxwell_3d%L2projection( 1, 3, efield_ref(nc_total*2+1:nc_total*3), az )
  
  rho = 0._f64
  call maxwell_3d%compute_rho_from_E( efield_ref, rho )
  print*, 'lhs divergence', maxval(abs(rho))

  call maxwell_3d%L2projection( 0, 1, rho_ref, cos_k )
   
  call maxwell_3d%curl_operator%dot(efield_ref, efield_val )
  call maxwell_3d%MG_operator%dot(rho_ref, bfield_val )

  efield_val = efield_val + bfield_val
  
  print*, 'error operator', maxval(abs(efield_val(1:nc_total) - current(1:nc_total) ) ), maxval(abs(efield_val(nc_total+1:nc_total*2) - current(nc_total+1:nc_total*2) ) ), maxval(abs(efield_val(nc_total*2+1:nc_total*3) - current(nc_total*2+1:nc_total*3) ) )

!!$  call sll_s_plot_two_fields_1d('currentx',nc_total,efield_val(1:nc_total),current(1:nc_total),istep,time)
!!$  call sll_s_plot_two_fields_1d('currenty',nc_total,efield_val(nc_total+1:nc_total*2),current(nc_total+1:nc_total*2),istep,time)
!!$  call sll_s_plot_two_fields_1d('currentz',nc_total,efield_val(nc_total*2+1:nc_total*3),current(nc_total*2+1:nc_total*3),istep,time)

  print*, 'error solver', maxval(abs(efield(1:nc_total) -efield_ref(1:nc_total) ) ), maxval(abs(efield(nc_total+1:nc_total*2) -efield_ref(nc_total+1:nc_total*2) ) ), maxval(abs(efield(nc_total*2+1:nc_total*3) -efield_ref(nc_total*2+1:nc_total*3) ) ), maxval(abs(rho_ref- maxwell_3d%uzawa_iterator%x_0))





!!$  call maxwell_3d%L2projection( 1, 1, efield(1:nc_total), cos_k )
!!$  call maxwell_3d%L2projection( 1, 2, efield(nc_total+1:nc_total*2), cos_k )
!!$  call maxwell_3d%L2projection( 1, 3, efield(nc_total*2+1:nc_total*3), cos_k )
!!$  efield(nc_total+1:nc_total*2) = - 2.0_f64 *  efield(nc_total+1:nc_total*2)
!!$
!!$  call  evaluate_spline_3d ( nc_eta, [deg(1)-1,deg(2),deg(3)], efield(1:nc_total), efield_val(1:nc_total) )
!!$  call  evaluate_spline_3d ( nc_eta, [deg(1),deg(2)-1,deg(3)], efield(1+nc_total:2*nc_total), efield_val(1+nc_total:2*nc_total) )
!!$  call  evaluate_spline_3d ( nc_eta, [deg(1),deg(2),deg(3)-1], efield(1+nc_total*2:3*nc_total), efield_val(1+nc_total*2:3*nc_total) )
!!$
!!$
!!$  ! Reference solutions
!!$  time = real(nsteps,f64)*delta_t
!!$  ind = 1
!!$  do k = 1, nc_eta(3)
!!$     do j = 1, nc_eta(2)
!!$        do i = 1, nc_eta(1)
!!$           efield_ref(ind) = cos_k([x(i,j,k), y(i,j,k), z(i,j,k)])
!!$           efield_ref(ind+nc_total) = -2.0_f64*efield_ref(ind)
!!$           efield_ref(ind+nc_total*2) = efield_ref(ind)
!!$           ind = ind+1
!!$        end do
!!$     end do
!!$  end do


  ! Clean up
  call maxwell_3d%free()
  deallocate(x)
  deallocate(y)
  deallocate(z)
  deallocate(efield)
  deallocate(bfield)
  deallocate(efield_ref)
  deallocate(bfield_ref)
  deallocate(efield_val)
  deallocate(bfield_val)
  deallocate( rho )
  deallocate( rho_ref )


  call sll_s_set_time_mark( end )
  write(*, "(A, F10.3)") "Main part run time [s] = ", sll_f_time_elapsed_between( start, end)



contains
  function cos_k(x)
    sll_real64             :: cos_k
    sll_real64, intent(in) :: x(3)

    cos_k = cos((x(1)+x(2)+x(3))-w1*time) 
  end function cos_k


  function sin_k(x)
    sll_real64             :: sin_k
    sll_real64, intent(in) :: x(3)

    sin_k = sin((x(1)+x(2)+x(3))-w1*time) 
  end function sin_k

  function jx(x)
    sll_real64             :: jx
    sll_real64, intent(in) :: x(3)

    jx = 3._f64*( sin(x(1)+x(2)+x(3)) + cos(x(1)+x(2)+x(3)) ) -sin(x(1)+x(2)+x(3))
  end function jx

  function jy(x)
    sll_real64             :: jy
    sll_real64, intent(in) :: x(3)

    jy = -3._f64*( cos(x(1)+x(2)+x(3)) - 4._f64* cos(2._f64*(x(1)+x(2)+x(3)) ) )  -sin(x(1)+x(2)+x(3))
  end function jy

  function jz(x)
    sll_real64             :: jz
    sll_real64, intent(in) :: x(3)

    jz = -3._f64*( sin(x(1)+x(2)+x(3)) + 4._f64* cos(2._f64*(x(1)+x(2)+x(3)) ) )  -sin(x(1)+x(2)+x(3))
  end function jz

  function ax(x)
    sll_real64             :: ax
    sll_real64, intent(in) :: x(3)

    ax =  sin(x(1)+x(2)+x(3)) + cos(x(1)+x(2)+x(3)) 
  end function ax

  function ay(x)
    sll_real64             :: ay
    sll_real64, intent(in) :: x(3)

    ay = -sin(x(1)+x(2)+x(3))**2 + cos(x(1)+x(2)+x(3))**2 -cos(x(1)+x(2)+x(3)) 
  end function ay

  function az(x)
    sll_real64             :: az
    sll_real64, intent(in) :: x(3)

    az = ( sin(x(1)+x(2)+x(3)) -1._f64) * sin(x(1)+x(2)+x(3)) - cos(x(1)+x(2)+x(3))**2 
  end function az

  subroutine evaluate_spline_3d ( ndofs, deg, dofs, vals )
    sll_int32, intent( in ) :: ndofs(3)
    sll_int32, intent( in ) :: deg(3)
    sll_real64, intent( in ) :: dofs(:)
    sll_real64, intent( out ) :: vals(:)

    sll_int32 :: i,j,k,istart,iend
    sll_real64 :: a_in(ndofs(2)), a_out(ndofs(2)),b_in(ndofs(3)), b_out(ndofs(3))


    istart = 1
    iend = ndofs(1)
    do k=1,ndofs(3)
       do j=1,ndofs(2)
          call sll_s_eval_uniform_periodic_spline_curve_with_zero(deg(1), dofs(istart:iend), vals(istart:iend))
          istart = iend+1
          iend = iend + ndofs(1)
       end do
    end do

    do k=1,ndofs(3)
       do i=1,ndofs(1)
          istart = (k-1)*ndofs(2)*ndofs(1)+i
          do j=1,ndofs(2)
             a_in(j) = vals(istart+(j-1)*ndofs(1))
          end do
          call sll_s_eval_uniform_periodic_spline_curve_with_zero(deg(2), a_in, a_out)
          do j=1,ndofs(2)
             vals(istart+(j-1)*ndofs(1)) = a_out(j)
          end do
       end do
    end do

    do j=1,ndofs(2)
       do i=1,ndofs(1)
          istart = (j-1)*ndofs(1)+i
          do k=1,ndofs(3)
             b_in(k) = vals(istart+(k-1)*ndofs(1)*ndofs(2))
          end do
          call sll_s_eval_uniform_periodic_spline_curve_with_zero(deg(3), b_in, b_out)
          do k=1,ndofs(3)
             vals(istart+(k-1)*ndofs(1)*ndofs(2)) = b_out(k)
          end do
       end do
    end do


  end subroutine evaluate_spline_3d

end program test_curl_curl_part
