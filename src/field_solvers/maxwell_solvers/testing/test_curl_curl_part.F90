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
       sll_p_pi, sll_p_twopi

  use sll_m_maxwell_3d_fem, only: &
       sll_t_maxwell_3d_fem

  use sll_m_maxwell_clamped_3d_fem

  use sll_m_maxwell_1d_base, only: &
       sll_s_plot_two_fields_1d

  use sll_m_constants, only: sll_p_pi, sll_p_twopi

  use sll_m_splines_pp

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

  sll_int32  :: nc_eta(3), nc_total, nc_total0, nc_total1

  type(sll_t_maxwell_3d_fem)  :: maxwell_3d
  !type(sll_t_maxwell_clamped_3d_fem)  :: maxwell_3d
  sll_real64, allocatable :: x(:,:,:), y(:,:,:), z(:,:,:)
  sll_real64, allocatable :: afield(:)
  sll_real64, allocatable :: afield_ref(:)
  sll_real64, allocatable :: afield_val(:), scratch(:)
  sll_real64, allocatable :: rho(:), rho_ref(:), current(:), noise(:)

  sll_int32                               :: i, j, k, ind
  sll_real64                                :: delta_t
  sll_real64                              :: eps
  sll_real64                              :: Lx(3)
  sll_real64, dimension(3,2)              :: domain
  sll_int32                               :: deg(3), boundary(3)

  type(sll_t_time_mark) :: start, end

  call sll_s_set_time_mark( start )

  ! Define computational domain
  eta1_min = .0_f64; eta1_max = 2.0_f64*sll_p_pi
  nc_eta = 16![8, 16, 32]
  nc_total = product(nc_eta)
  Lx(1) = eta1_max-eta1_min
  Lx(2) = Lx(1); Lx(3) = Lx(2)
  delta_eta = Lx/real(nc_eta,f64)
  domain(1,:) = [eta1_min, eta1_max]
  domain(2,:) = [eta1_min, eta1_max]
  domain(3,:) = [eta1_min, eta1_max]
  ! Set spline degree of 0-forms
  deg = 3![2,2,3]
  ! Time loop
  delta_t = 0.01_f64


  if( deg(1) == 2 ) then
     boundary = [ sll_p_boundary_clamped_square, sll_p_boundary_periodic, sll_p_boundary_periodic]
  else if( deg(1) == 3 ) then
     boundary = [ sll_p_boundary_clamped_cubic, sll_p_boundary_periodic, sll_p_boundary_periodic]
  else
     boundary = [ sll_p_boundary_clamped, sll_p_boundary_periodic, sll_p_boundary_periodic]
  end if

  ! Initialise maxwell FEM object
  call maxwell_3d%init(domain, nc_eta, deg)
  !call maxwell_3d%init(domain, nc_eta, deg, boundary)

  nc_total0 = maxwell_3d%n_total0
  nc_total1 = maxwell_3d%n_total1

  allocate(x(1:nc_eta(1)+1,1:nc_eta(2)+1,1:nc_eta(3)+1) )
  allocate(y(1:nc_eta(1)+1,1:nc_eta(2)+1,1:nc_eta(3)+1) )
  allocate(z(1:nc_eta(1)+1,1:nc_eta(2)+1,1:nc_eta(3)+1) )
  allocate( afield(1:nc_total1+nc_total0*2) )
  allocate( afield_ref(1:nc_total1+nc_total0*2) )
  allocate( afield_val(1:nc_total1+nc_total0*2) )
  allocate( scratch(1:nc_total1+nc_total0*2) )

  allocate( rho( nc_total0 ) )
  allocate( rho_ref( nc_total0 ) )
  allocate( current(1:nc_total1+nc_total0*2) )
  allocate( noise(1:nc_total1+nc_total0*2) )

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
!!$  call maxwell_3d%compute_rhs_from_function( 1, 1, current(1:nc_total1), cos_k )
!!$  call maxwell_3d%compute_rhs_from_function( 1, 3, current(nc_total1+nc_total0+1:nc_total1+nc_total0*2), cos_k )
!!$  current(1:nc_total1) = 3._f64* current(1:nc_total1)
!!$  current(nc_total1+nc_total0+1:nc_total1+nc_total0*2) = -3._f64*current(nc_total1+nc_total0+1:nc_total1+nc_total0*2)
  call maxwell_3d%compute_rhs_from_function( 1, 1, current(1:nc_total1), jx )
  call maxwell_3d%compute_rhs_from_function( 1, 2, current(nc_total1+1:nc_total1+nc_total0), jy )
  call maxwell_3d%compute_rhs_from_function( 1, 3, current(nc_total1+nc_total0+1:nc_total1+nc_total0*2), jz )
!!$  call random_number( afield )
!!$  call maxwell_3d%multiply_ct( afield, noise )
  call random_number( noise )
  eps = 1d-6
  current = (1._f64-eps) * current +  eps* noise
!!$  rho = 0._f64
  call maxwell_3d%multiply_gt( current, rho )
  print*, 'rhs divergence', maxval(abs(rho))

 ! maxwell_3d%curl_matrix%epsilon = 6.0_f64!1d-0
  !call maxwell_3d%curl_solver%solve( current, afield )
  !call maxwell_3d%preconditioner_curl_fft%solve( current, afield )
  call maxwell_3d%uzawa_iterator%solve( current, afield )

  !print*, 'afield', afield

  !print*, 'p',maxwell_3d%uzawa_iterator%x_0

  afield_ref = 0._f64
!!$  call maxwell_3d%L2projection( 1, 1, afield_ref(1:nc_total1), cos_k )
!!$  call maxwell_3d%L2projection( 1, 3, afield_ref(nc_total1+nc_total0+1:nc_total1+nc_total0*2), cos_k )
!!$  afield_ref(nc_total1+nc_total0+1:nc_total1+nc_total0*2) = -afield_ref(nc_total1+nc_total0+1:nc_total1+nc_total0*2)

  call maxwell_3d%L2projection( 1, 1, afield_ref(1:nc_total1), ax )
  call maxwell_3d%L2projection( 1, 2, afield_ref(nc_total1+1:nc_total1+nc_total0), ay )
  call maxwell_3d%L2projection( 1, 3, afield_ref(nc_total1+nc_total0+1:nc_total1+nc_total0*2), az )
  
  rho = 0._f64
  call maxwell_3d%compute_rho_from_E( afield_ref, rho )
!!$  print*, 'lhs divergence', maxval(abs(rho))
!!$
  rho_ref = 0._f64
  call maxwell_3d%L2projection( 0, 1, rho_ref, cos_k )
  maxwell_3d%curl_matrix%epsilon = 0._f64
  call maxwell_3d%curl_matrix%dot(afield, afield_val )
  call maxwell_3d%MG_operator%dot(maxwell_3d%uzawa_iterator%x_0, scratch )

  afield_val = afield_val + scratch
  
  print*, 'error operator', maxval(abs(afield_val(1:nc_total1) - current(1:nc_total1) ) ), maxval(abs(afield_val(nc_total1+1:nc_total1+nc_total0) - current(nc_total1+1:nc_total1+nc_total0) ) ), maxval(abs(afield_val(nc_total1+nc_total0+1:nc_total0+nc_total*2) - current(nc_total1+nc_total0+1:nc_total0+nc_total*2) ) )

!!$  call sll_s_plot_two_fields_1d('currentx',nc_total,afield_val(1:nc_total),current(1:nc_total),0._f64,0._f64)
!!$  call sll_s_plot_two_fields_1d('currenty',nc_total,afield_val(nc_total+1:nc_total*2),current(nc_total+1:nc_total*2),0._f64,0._f64)
!!$  call sll_s_plot_two_fields_1d('currentz',nc_total,afield_val(nc_total*2+1:nc_total*3),current(nc_total*2+1:nc_total*3),0._f64,0._f64)

  print*, 'error solver', maxval(abs(afield(1:nc_total1) -afield_ref(1:nc_total1) ) ), maxval(abs(afield(nc_total1+1:nc_total1+nc_total0) -afield_ref(nc_total1+1:nc_total1+nc_total0) ) ), maxval(abs(afield(nc_total1+nc_total0+1:nc_total0+nc_total*2) -afield_ref(nc_total1+nc_total0+1:nc_total0+nc_total*2) ) ), maxval(abs(rho_ref- maxwell_3d%uzawa_iterator%x_0))





!!$  call maxwell_3d%L2projection( 1, 1, afield(1:nc_total), cos_k )
!!$  call maxwell_3d%L2projection( 1, 2, afield(nc_total+1:nc_total*2), cos_k )
!!$  call maxwell_3d%L2projection( 1, 3, afield(nc_total*2+1:nc_total*3), cos_k )
!!$  afield(nc_total+1:nc_total*2) = - 2.0_f64 *  afield(nc_total+1:nc_total*2)
!!$
!!$  call  evaluate_spline_3d ( nc_eta, [deg(1)-1,deg(2),deg(3)], afield(1:nc_total), afield_val(1:nc_total) )
!!$  call  evaluate_spline_3d ( nc_eta, [deg(1),deg(2)-1,deg(3)], afield(1+nc_total:2*nc_total), afield_val(1+nc_total:2*nc_total) )
!!$  call  evaluate_spline_3d ( nc_eta, [deg(1),deg(2),deg(3)-1], afield(1+nc_total*2:3*nc_total), afield_val(1+nc_total*2:3*nc_total) )
!!$
!!$
!!$  ! Reference solutions
!!$  ind = 1
!!$  do k = 1, nc_eta(3)
!!$     do j = 1, nc_eta(2)
!!$        do i = 1, nc_eta(1)
!!$           afield_ref(ind) = cos_k([x(i,j,k), y(i,j,k), z(i,j,k)])
!!$           afield_ref(ind+nc_total) = -2.0_f64*afield_ref(ind)
!!$           afield_ref(ind+nc_total*2) = afield_ref(ind)
!!$           ind = ind+1
!!$        end do
!!$     end do
!!$  end do


  ! Clean up
  call maxwell_3d%free()
  deallocate(x)
  deallocate(y)
  deallocate(z)
  deallocate(afield)
  deallocate(afield_ref)
  deallocate(afield_val)
  deallocate(scratch)
  deallocate( rho )
  deallocate( rho_ref )


  call sll_s_set_time_mark( end )
  write(*, "(A, F10.3)") "Main part run time [s] = ", sll_f_time_elapsed_between( start, end)



contains
  function cos_k(x)
    sll_real64             :: cos_k
    sll_real64, intent(in) :: x(3)

    cos_k = cos((x(1)+x(2)+x(3))) 
  end function cos_k


  function sin_k(x)
    sll_real64             :: sin_k
    sll_real64, intent(in) :: x(3)

    sin_k = sin((x(1)+x(2)+x(3))) 
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
