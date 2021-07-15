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
!  Contact : Katharina Kormann
!
program test_maxwell_2d_fem_fft
  !------------------------------------------------------------------------
  !  test 3D Maxwell spline finite element solver on a periodic grid
  !------------------------------------------------------------------------
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#include "sll_memory.h"
#include "sll_working_precision.h"
#include "sll_maxwell_solvers_macros.h"

  use sll_m_low_level_bsplines, only: &
    sll_s_eval_uniform_periodic_spline_curve

  use sll_m_constants, only: &
    sll_p_pi

  use sll_m_maxwell_2d_fem_fft, only: &
       sll_t_maxwell_2d_fem_fft
  
  use sll_m_constants, only: sll_p_pi, sll_p_twopi

  implicit none
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  !  type(sll_t_arbitrary_degree_spline_1d), pointer :: aspl
  !  sll_real64, dimension(:), allocatable :: knots

  sll_real64 :: eta1_max, eta1_min
  sll_real64 :: delta_eta(2)

  sll_int32  :: nc_eta(2), nc_total

  type(sll_t_maxwell_2d_fem_fft)  :: maxwell
  sll_real64, allocatable :: x(:,:), y(:,:)
  sll_real64, allocatable :: efield(:), bfield(:)
  sll_real64, allocatable :: efield_ref(:), bfield_ref(:)
  sll_real64, allocatable :: efield_val(:), bfield_val(:)
  sll_real64, allocatable :: rho(:), rho_ref(:)
  
  sll_int32                               :: i, j, istep, nsteps, ind
  sll_real64                                :: w1, w2
  sll_real64                                :: time
  sll_real64                                :: delta_t
  
  sll_real64                              :: Lx(2)
  sll_real64, dimension(2,2)              :: domain
  sll_int32                               :: deg

  sll_real64                              :: error(6)

  ! Define computational domain
  eta1_min = .0_f64; eta1_max = 2.0_f64*sll_p_pi
  nc_eta = [16,32]
  nc_total = product(nc_eta)
  Lx(1) = eta1_max-eta1_min
  Lx(2) = Lx(1); 
  delta_eta = Lx/nc_eta
  domain(1,:) = [eta1_min, eta1_max]
  domain(2,:) = [eta1_min, eta1_max]
  ! Set spline degree of 0-forms
  deg = 3
  ! Time loop
  delta_t = 0.01_f64
  nsteps = 30
  
  ! Initialise maxwell FEM object
  call maxwell%init(domain, nc_eta, deg)


  allocate(x(1:nc_eta(1)+1,1:nc_eta(2)+1) )
  allocate(y(1:nc_eta(1)+1,1:nc_eta(2)+1) )
  allocate(efield(1:nc_total*3))
  allocate(bfield(1:nc_total*3))
  allocate(efield_ref(1:nc_total*3))
  allocate(bfield_ref(1:nc_total*3))
  allocate(efield_val(1:nc_total*3))
  allocate(bfield_val(1:nc_total*3))
  
  do j = 1, nc_eta(2)+1
     do i = 1, nc_eta(1)+1
        x(i,j) = eta1_min + (i-1)*delta_eta(1)
        y(i,j) = eta1_min + (j-1)*delta_eta(2)
     end do
  end do

  w1 = sqrt(3.0_f64)
  w2 = + sqrt(3.0_f64)
  
  ! Poisson problem
  allocate( rho( nc_total) )
  allocate( rho_ref( nc_total) )
  time = 0.0_f64
  
  call maxwell%compute_rhs_from_function( cos_k, 1, 0, rho )

  
  call maxwell%compute_e_from_rho( rho, efield )
  call  evaluate_spline_2d ( nc_eta, [deg-1,deg], efield(1:nc_total),  &
       efield_val(1:nc_total))
  call  evaluate_spline_2d ( nc_eta, [deg,deg-1], efield(1+nc_total:nc_total*2),  &
       efield_val(1+nc_total:nc_total*2))
  call  evaluate_spline_2d ( nc_eta, [deg,deg], efield(1+nc_total*2:nc_total*3),  &
       efield_val(1+nc_total*2:nc_total*3))
  ! Reference solution
  ind = 1
  do j = 1, nc_eta(2)
     do i = 1, nc_eta(1)
        efield_ref(ind) = sin_k([x(i,j), y(i,j)])/2.0_f64
        efield_ref(ind+nc_total) = efield_ref(ind)
        efield_ref(ind+nc_total*2) = 0.0_f64
        ind = ind+1
     end do
  end do
  error(1) = maxval(abs(efield_val-efield_ref))
  print*, 'Error Poisson:', error(1)

  ! Now assemble initial efield and bfield
  time = 0.0_f64

  time = -0.5_f64*delta_t

  call maxwell%L2projection( e1, 1, 2, bfield(1:nc_total) )
  call maxwell%L2projection( e2, 2, 2, bfield(nc_total+1:nc_total*2) )
  call maxwell%L2projection( b3, 3, 2, bfield(nc_total*2+1:nc_total*3) )

  time = 0.0_f64
  
  call maxwell%L2projection( e1, 1, 1, efield(1:nc_total) )
  call maxwell%L2projection( e2, 2, 1, efield(nc_total+1:nc_total*2) )
  call maxwell%L2projection( b3, 3, 1, efield(nc_total*2+1:nc_total*3) )
  efield(nc_total*2+1:nc_total*3) = -efield(nc_total*2+1:nc_total*3)
  
  ! Time stepping
  do istep = 1, nsteps
     call maxwell%compute_b_from_e( delta_t, efield, bfield )
     call maxwell%compute_e_from_b( delta_t, bfield, efield )
  end do

  ! Evaluate E and B at the grid points
  call  evaluate_spline_2d ( nc_eta, [deg,deg-1], bfield(1:nc_total), bfield_val(1:nc_total) )
  call  evaluate_spline_2d ( nc_eta, [deg-1,deg], bfield(1+nc_total:2*nc_total), bfield_val(1+nc_total:2*nc_total)  )
  call  evaluate_spline_2d ( nc_eta, [deg-1,deg-1], bfield(1+nc_total*2:3*nc_total), bfield_val(1+nc_total*2:3*nc_total) )
  call  evaluate_spline_2d ( nc_eta, [deg-1,deg], efield(1:nc_total), efield_val(1:nc_total) )
  call  evaluate_spline_2d ( nc_eta, [deg,deg-1], efield(1+nc_total:2*nc_total), efield_val(1+nc_total:2*nc_total) )
  call  evaluate_spline_2d ( nc_eta, [deg,deg], efield(1+nc_total*2:3*nc_total), efield_val(1+nc_total*2:3*nc_total) )

  ! Reference solutions
  time = (nsteps)*delta_t
  ind = 1
  do j = 1, nc_eta(2)
     do i = 1, nc_eta(1)
        efield_ref(ind) = e1([x(i,j), y(i,j)])
        efield_ref(ind+nc_total) = e2([x(i,j), y(i,j)])
        efield_ref(ind+nc_total*2) = -b3([x(i,j), y(i,j)])
        ind = ind+1
     end do
  end do
  time = (nsteps-0.5_f64)*delta_t
  ind = 1
  do j = 1, nc_eta(2)
     do i = 1, nc_eta(1)
        bfield_ref(ind) = e1([x(i,j), y(i,j)])
        bfield_ref(ind+nc_total) = e2([x(i,j), y(i,j)])
        bfield_ref(ind+nc_total*2) = b3([x(i,j), y(i,j)])
        ind = ind+1
     end do
  end do
  error(3) = maxval(abs(efield_val-efield_ref))
  error(4) = maxval(abs(bfield_val-bfield_ref))
  print*, 'Error efield:', error(3)
  print*, 'Error bfield:', error(4)


  ! Test compute_e_from_j
  call maxwell%L2projection( cos_k, 1, 1, efield(1:nc_total) )
  
  error(2) = maxwell%inner_product( efield(1:nc_total), efield(1:nc_total), 1, 1 ) - 2.0_f64*sll_p_pi**2
  print*, 'Error in L2 norm squared:', error(2)
  
  call maxwell%compute_rhs_from_function( sin_k, 1, 1, rho )

  call maxwell%compute_E_from_j( rho, 1, efield(1:nc_total) )
  call  evaluate_spline_2d ( nc_eta, [deg-1,deg,deg], efield(1:nc_total),  &
       efield_val(1:nc_total))
  ! Reference solution
  ind = 1
  do j = 1, nc_eta(2)
     do i = 1, nc_eta(1)
        efield_ref(ind) = cos_k([x(i,j), y(i,j)])-&
             sin_k([x(i,j), y(i,j)])
        ind = ind+1
     end do
  end do
  error(5) = maxval(abs(efield_val(1:nc_total)-efield_ref(1:nc_total)))
  print*, 'Error compute_e_from_j:', error(5)

  ! Test compute_rho_from_e
  call maxwell%compute_rhs_from_function( cos_k, 1, 0, rho_ref )
  rho_ref = 2.0_f64*rho_ref
  call maxwell%L2projection( sin_k, 1, 1, efield(1:nc_total) )
  call maxwell%L2projection( sin_k, 2, 1, efield(nc_total+1:nc_total*2) )
  call maxwell%L2projection( sin_k, 3, 1, efield(nc_total*2+1:nc_total*3) )

  call maxwell%compute_rho_from_e( efield, rho )
  
  error(6) =  maxval( abs( rho - rho_ref ) )
  print*, 'Error compute_rho_from_e:', error(6)
  

  ! Clean up
  call maxwell%free()
  deallocate(x)
  deallocate(y)
  deallocate(efield)
  deallocate(bfield)
  deallocate(efield_ref)
  deallocate(bfield_ref)
  deallocate(efield_val)
  deallocate(bfield_val)


  if ( error(1) < 5.0E-5 .AND. error(2) < 2.0E-6 .AND. error(3) < 3.2E-5 .AND. error(4) < 9.9E-5 .AND. error(5)<1.4E-4 .AND. error(6)<1.4E-9) then
     print*, 'PASSED.'
  else
     print*, 'FAILED.'
  end if

contains
  function cos_k(x)
    sll_real64             :: cos_k
    sll_real64, intent(in) :: x(2)

    cos_k = cos((x(1)+x(2))-w1*time) 
  end function cos_k

  function b3(x)
    sll_real64             :: b3
    sll_real64, intent(in) :: x(2)

    b3 = - cos(x(1))*cos(x(2))*cos(sqrt(2.0_f64)*time)
    
  end function b3

  
  function e1(x)
    sll_real64             :: e1
    sll_real64, intent(in) :: x(2)

    e1 = cos(x(1))*sin(x(2))*sin(sqrt(2.0_f64)*time)/sqrt(2.0_f64)
    
  end function e1
  
  function e2(x)
    sll_real64             :: e2
    sll_real64, intent(in) :: x(2)

    e2 = - sin(x(1))*cos(x(2))*sin(sqrt(2.0_f64)*time)/sqrt(2.0_f64)
    
  end function e2


  function sin_k(x)
    sll_real64             :: sin_k
    sll_real64, intent(in) :: x(2)

    sin_k = sin((x(1)+x(2))-w1*time) 
  end function sin_k

  subroutine evaluate_spline_2d ( ndofs, deg, dofs, vals )
    sll_int32, intent( in ) :: ndofs(2)
    sll_int32, intent( in ) :: deg(2)
    sll_real64, intent( in ) :: dofs(:)
    sll_real64, intent( out ) :: vals(:)

    sll_int32 :: i,j,istart,iend
    sll_real64 :: a_in(ndofs(2)), a_out(ndofs(2))


    istart = 1
    iend = ndofs(1)
    do j=1,ndofs(2)
       call sll_s_eval_uniform_periodic_spline_curve(deg(1), dofs(istart:iend), vals(istart:iend))
       istart = iend+1
       iend = iend + ndofs(1)
    end do

    do i=1,ndofs(1)
       istart = i
       do j=1,ndofs(2)
          a_in(j) = vals(istart+(j-1)*ndofs(1))
       end do
       call sll_s_eval_uniform_periodic_spline_curve(deg(2), a_in, a_out)
       do j=1,ndofs(2)
          vals(istart+(j-1)*ndofs(1)) = a_out(j)
       end do
    end do
    
    
  end subroutine evaluate_spline_2d
  

end program test_maxwell_2d_fem_fft
