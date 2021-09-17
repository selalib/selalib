! Contact: Katharina Kormann


!------------------------------------------------------------------------
!  test 3D Poisson spline finite element solver on a periodic grid
!------------------------------------------------------------------------
program test_poisson_3d_fem_fft
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#include "sll_memory.h"
#include "sll_working_precision.h"
#include "sll_maxwell_solvers_macros.h"

  use sll_m_low_level_bsplines, only: &
    sll_s_eval_uniform_periodic_spline_curve

  use sll_m_constants, only: &
    sll_p_pi

  use sll_m_poisson_3d_fem_fft, only: &
       sll_t_poisson_3d_fem_fft
  
  use sll_m_constants, only: sll_p_pi, sll_p_twopi

  implicit none
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


  sll_real64 :: eta1_max, eta1_min
  sll_real64 :: delta_eta(3)

  sll_int32  :: nc_eta(3), nc_total

  type(sll_t_poisson_3d_fem_fft)  :: poisson_3d
  sll_real64, allocatable :: x(:,:,:), y(:,:,:), z(:,:,:)
  sll_real64, allocatable :: efield(:), bfield(:)
  sll_real64, allocatable :: efield_ref(:), bfield_ref(:)
  sll_real64, allocatable :: efield_val(:), bfield_val(:)
  sll_real64, allocatable :: rho(:), rho_ref(:)
  
  sll_int32                               :: i, j, k, nsteps, ind
  sll_real64                                :: time
  sll_real64                                :: delta_t
  
  sll_real64                              :: Lx(3)
  sll_real64, dimension(3,2)              :: domain
  sll_int32                               :: deg

  sll_real64                              :: error

  ! Define computational domain
  eta1_min = .0_f64; eta1_max = 2.0_f64*sll_p_pi
  nc_eta = [16, 8, 32]
  nc_total = product(nc_eta)
  Lx(1) = eta1_max-eta1_min
  Lx(2) = Lx(1); Lx(3) = Lx(2)
  delta_eta = Lx/nc_eta
  domain(1,:) = [eta1_min, eta1_max]
  domain(2,:) = [eta1_min, eta1_max]
  domain(3,:) = [eta1_min, eta1_max]
  ! Set spline degree of 0-forms
  deg = 3
  ! Time loop
  delta_t = 0.01_f64
  nsteps = 3
  
  ! Initialise maxwell FEM object
  call poisson_3d%init(nc_eta, [deg,deg,deg], delta_eta)


  allocate(x(1:nc_eta(1)+1,1:nc_eta(2)+1,1:nc_eta(3)+1) )
  allocate(y(1:nc_eta(1)+1,1:nc_eta(2)+1,1:nc_eta(3)+1) )
  allocate(z(1:nc_eta(1)+1,1:nc_eta(2)+1,1:nc_eta(3)+1) )
  allocate(efield(1:nc_total*3))
  allocate(bfield(1:nc_total*3))
  allocate(efield_ref(1:nc_total*3))
  allocate(bfield_ref(1:nc_total*3))
  allocate(efield_val(1:nc_total*3))
  allocate(bfield_val(1:nc_total*3))
  
  do k = 1, nc_eta(3)+1
     do j = 1, nc_eta(2)+1
        do i = 1, nc_eta(1)+1
           x(i,j,k) = eta1_min + (i-1)*delta_eta(1)
           y(i,j,k) = eta1_min + (j-1)*delta_eta(2)
           z(i,j,k) = eta1_min + (k-1)*delta_eta(3)
        end do
     end do
  end do



  ! Poisson problem
  allocate( rho( nc_total) )
  allocate( rho_ref( nc_total) )
  time = 0.0_f64
  
  call poisson_3d%compute_rhs_from_function( cos_k, rho )
  call poisson_3d%compute_e_from_rho( rho, efield )
  call  evaluate_spline_3d ( nc_eta, [deg-1,deg,deg], efield(1:nc_total),  &
       efield_val(1:nc_total))
  call  evaluate_spline_3d ( nc_eta, [deg,deg-1,deg], efield(1+nc_total:nc_total*2),  &
       efield_val(1+nc_total:nc_total*2))
  call  evaluate_spline_3d ( nc_eta, [deg,deg,deg-1], efield(1+nc_total*2:nc_total*3),  &
       efield_val(1+nc_total*2:nc_total*3))
  ! Reference solution
  ind = 1
  do k = 1, nc_eta(3)
     do j = 1, nc_eta(2)
        do i = 1, nc_eta(1)
           efield_ref(ind) = sin_k([x(i,j,k), y(i,j,k), z(i,j,k)])/3.0_f64
           efield_ref(ind+nc_total) = efield_ref(ind)
           efield_ref(ind+nc_total*2) = efield_ref(ind)
           ind = ind+1
        end do
     end do
  end do
  error = maxval(abs(efield_val-efield_ref))
  print*, 'Error Poisson:', error

  write(11,*) efield_val
  write(12,*) efield_ref
  
  call poisson_3d%free()


  if ( error < 5.5E-4 ) then
     print*, 'PASSED.'
  else
     print*, 'FAILED.'
  end if


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
          call sll_s_eval_uniform_periodic_spline_curve(deg(1), dofs(istart:iend), vals(istart:iend))
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
          call sll_s_eval_uniform_periodic_spline_curve(deg(2), a_in, a_out)
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
          call sll_s_eval_uniform_periodic_spline_curve(deg(3), b_in, b_out)
          do k=1,ndofs(3)
             vals(istart+(k-1)*ndofs(1)*ndofs(2)) = b_out(k)
          end do
       end do
    end do
    
    
  end subroutine evaluate_spline_3d
  

end program test_poisson_3d_fem_fft
