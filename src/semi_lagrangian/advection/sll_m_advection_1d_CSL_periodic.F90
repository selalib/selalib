!**************************************************************
!  Copyright INRIA
!  Authors : 
!     CALVI project team
!  
!  This code SeLaLib (for Semi-Lagrangian-Library) 
!  is a parallel library for simulating the plasma turbulence 
!  in a tokamak.
!  
!  This software is governed by the CeCILL-B license 
!  under French law and abiding by the rules of distribution 
!  of free software.  You can  use, modify and redistribute 
!  the software under the terms of the CeCILL-B license as 
!  circulated by CEA, CNRS and INRIA at the following URL
!  "http://www.cecill.info". 
!**************************************************************

! in development; should be at least cubic splines
! attached with computation of characteristics


module sll_m_advection_1d_CSL_periodic
#include "sll_working_precision.h"
#include "sll_memory.h"
#include "sll_assert.h"
#include "sll_errors.h"
use sll_m_boundary_condition_descriptors
use sll_m_advection_1d_base
use sll_m_characteristics_1d_base
use sll_m_interpolators_1d_base
use sll_m_hermite_interpolation_1d
use sll_m_gauss_legendre_integration
use sll_m_lagrange_interpolation
implicit none

  type,extends(sll_advection_1d_base) :: CSL_periodic_1d_advector
  
    class(sll_interpolator_1d_base), pointer  :: interp
    class(sll_characteristics_1d_base), pointer  :: charac
    sll_real64, dimension(:), pointer :: eta_coords
    sll_real64, dimension(:), pointer :: charac_feet
    sll_real64, dimension(:), pointer :: charac_feet_inside
    !sll_real64, dimension(:), pointer :: csl_mat_init
    !sll_real64, dimension(:), pointer :: csl_mat
    !sll_real64, dimension(:), pointer :: fft_buf
    sll_real64, dimension(:,:), pointer :: deriv
    sll_int32 :: Npts
    sll_int32 :: csl_degree
  contains
    procedure, pass(adv) :: initialize => &
       initialize_CSL_periodic_1d_advector
    procedure, pass(adv) :: advect_1d => &
      CSL_periodic_advect_1d
    procedure, pass(adv) :: advect_1d_constant => &
      CSL_periodic_advect_1d_constant
    procedure, pass(adv) :: delete => delete_CSL_periodic_1d_adv  
  end type CSL_periodic_1d_advector
   




contains
  function new_CSL_periodic_1d_advector( &
    interp, &
    charac, &
    Npts, &
    eta_min, &
    eta_max, &
    eta_coords, &
    csl_degree) &  
    result(adv)      
    type(CSL_periodic_1d_advector), pointer :: adv
    class(sll_interpolator_1d_base), pointer :: interp
    class(sll_characteristics_1d_base), pointer  :: charac
    sll_int32, intent(in) :: Npts
    sll_real64, intent(in), optional :: eta_min
    sll_real64, intent(in), optional :: eta_max
    sll_real64, dimension(:), pointer, optional :: eta_coords
    sll_int32, intent(in), optional :: csl_degree
    sll_int32 :: ierr
    
    SLL_ALLOCATE(adv,ierr)
        
    call initialize_CSL_periodic_1d_advector(&
      adv, &
      interp, &
      charac, &
      Npts, &
      eta_min, &
      eta_max, &
      eta_coords, &
      csl_degree)    
    
  end function  new_CSL_periodic_1d_advector


  subroutine initialize_CSL_periodic_1d_advector(&
    adv, &
    interp, &
    charac, &
    Npts, &
    eta_min, &
    eta_max, &
    eta_coords, &
    csl_degree)    
    class(CSL_periodic_1d_advector), intent(inout) :: adv
    class(sll_interpolator_1d_base), pointer :: interp
    class(sll_characteristics_1d_base), pointer  :: charac
    sll_int32, intent(in) :: Npts
    sll_real64, intent(in), optional :: eta_min
    sll_real64, intent(in), optional :: eta_max
    sll_real64, dimension(:), pointer, optional :: eta_coords
    sll_int32, optional :: csl_degree
    sll_int32 :: ierr
    sll_int32 :: i
    sll_real64 :: delta_eta
    !sll_real64 :: hermite_inversibility     
    
    !hermite_inversibility = 1._f64
    
    if(present(csl_degree))then
      adv%csl_degree = csl_degree
    else
      adv%csl_degree = 3  
    endif
    
    adv%Npts = Npts
    adv%interp => interp
    adv%charac => charac
    SLL_ALLOCATE(adv%eta_coords(Npts),ierr)

    SLL_ALLOCATE(adv%charac_feet(Npts),ierr)
    SLL_ALLOCATE(adv%charac_feet_inside(Npts),ierr)
    !SLL_ALLOCATE(adv%csl_mat(Npts-1),ierr)
    !SLL_ALLOCATE(adv%csl_mat_init(Npts-1),ierr)
    !SLL_ALLOCATE(adv%fft_buf(2*(Npts-1)+15),ierr)
    SLL_ALLOCATE(adv%deriv(2,Npts-1),ierr)
    
    !call dffti(Npts-1,adv%fft_buf)
    
    if(present(eta_min).and.present(eta_max))then
      if(present(eta_coords))then
        print *,'#provide either eta_coords or eta_min and eta_max'
        print *,'#and not both in subroutine initialize_BSL_1d_advector'
        stop
      else
        delta_eta = (eta_max-eta_min)/real(Npts-1,f64)
        do i=1,Npts
          adv%eta_coords(i) = eta_min+real(i-1,f64)*delta_eta
        enddo
      endif
    else if(present(eta_coords))then
      if(size(eta_coords)<Npts)then
        print *,'#bad size for eta_coords in initialize_BSL_1d_advector'
        stop
      else
        adv%eta_coords(1:Npts) = eta_coords(1:Npts)
      endif     
    else
      print *,'#Warning, we assume eta_min = 0._f64 eta_max = 1._f64'
      delta_eta = 1._f64/real(Npts-1,f64)
      do i=1,Npts
          adv%eta_coords(i) = real(i-1,f64)*delta_eta
      enddo                      
    endif
    
!    call compute_csl_mat( &
!      interp, &
!      adv%eta_coords(1), &
!      adv%eta_coords(Npts), &
!      Npts, &
!      adv%csl_mat)
    
    
    !if(present(hermite_degree))then
    !  call compute_csl_hermite_mat( &
    !    hermite_degree, &
    !    Npts-1,adv%csl_mat)
      !print *,'#hermite_degree=',hermite_degree
      !print *,'#csl_mat='
      !do i=1,Npts-1
      !  print *,i,adv%csl_mat(i)
      !enddo
      !stop
      !do i=1,Npts-1
      !  adv%csl_mat_init(i) = adv%csl_mat(i)
      !enddo
      
      
      
      !call dfftf(Npts-1,adv%csl_mat,adv%fft_buf)
      !hermite_inversibility = abs(adv%csl_mat(1))
      !print *,0,adv%csl_mat(1)
      !do i=1,(Npts-1)/2
        !print *,i-1,adv%csl_mat(2*(i-1)),adv%csl_mat(2*(i-1)+1), &
        !hermite_inversibility = min(hermite_inversibility, &
        !  sqrt(adv%csl_mat(2*i)**2+adv%csl_mat(2*i+1)**2))
      !enddo
      !if(mod(Npts-1,2)==0)then
      !  print *,(Npts-1)/2,adv%csl_mat(Npts-1)
      !endif
      
      !print *,'#hermite_inversibility=',hermite_inversibility
      !if(hermite_inversibility<1.-12)then
      !  print *,'#csl_mat not invertible',hermite_inversibility
      !  stop
      !endif
      !call compute_inverse_csl_hermite_mat(Npts-1,adv%csl_mat)



      !call check_solve_csl_mat( &
      !  adv%csl_mat, &
      !  adv%csl_mat_init, &
      !  Npts-1, &
      !  adv%fft_buf)
      
      
      
      !print *,'#inverse mat'
      !do i=1,Npts-1
      !  print *,i,adv%csl_mat(i)*sqrt(real(Npts-1,f64))
      !enddo
      !stop
           
  end subroutine initialize_CSL_periodic_1d_advector

  subroutine CSL_periodic_advect_1d(&
    adv, &
    A, &
    dt, &
    input, &
    output)
    class(CSL_periodic_1d_advector) :: adv
    sll_real64, dimension(:), intent(in) :: A
    sll_real64, intent(in) :: dt 
    sll_real64, dimension(:), intent(in) :: input
    sll_real64, dimension(:), intent(out) :: output      
    sll_int32 :: i
    
    call adv%charac%compute_characteristics( &
      A, &
      dt, &
      adv%eta_coords, &
      adv%charac_feet)
    
    call check_charac_feet( &
      adv%eta_coords, &
      adv%charac_feet, &
      adv%Npts)
    
    do i=1,adv%Npts
      adv%charac_feet_inside(i) = &
        process_outside_point_periodic( &
          adv%charac_feet(i), &
          adv%eta_coords(1), &
          adv%eta_coords(adv%Npts))
    enddo
    call adv%interp%compute_interpolants(input)

    !compute the derivatives at time tn
    !for future use in Hermite form
    call compute_node_derivative_order3( &
      adv%interp, &
      adv%deriv, &
      adv%Npts-1, &
      adv%eta_coords(1), &
      adv%eta_coords(adv%Npts))
    
    call update_solution_csl_periodic( &
      adv%interp, &
      input, &
      adv%deriv, &
      adv%charac_feet, &
      adv%Npts, &
      adv%eta_coords(1), &
      adv%eta_coords(adv%Npts), &
      output, &
      adv%csl_degree)

!    call compute_csl_integral( &
!      adv%Npts, &
!      adv%interp, &
!      adv%eta_coords, &
!      adv%charac_feet, &
!      output(1:adv%Npts-1))
!    if(maxval(input)>1.e-1)then
!      print *,'before'
!      do i=1,adv%Npts-1
!        print *,i,input(i),output(i),0.5_f64*(input(i)+input(i+1))
!      enddo
!!      stop
!    endif             
    
!    call mult_circulant_mat( &
!      adv%Npts-1, &
!      adv%csl_mat, &
!      output(1:adv%Npts-1), &
!      adv%fft_buf)
!    
!    output(adv%Npts) = output(1)
!    if(maxval(input)>1.e-1)then
!      print *,'after'
!      do i=1,adv%Npts-1
!        print *,i,input(i),output(i),0.5_f64*(output(i)+output(i+1))
!        !adv%csl_mat_init(i)
!      enddo
!      stop
!    endif             
!    call adv%interp%compute_interpolants( &
!      input, &
!      adv%eta1_coords, &
!      adv%Npts1, &
!      adv%eta2_coords, &
!      adv%Npts2 )

!    output = adv%interp%interpolate_array( &
!      adv%Npts, &
!      input, &
!      adv%charac_feet_inside)      
          
  end subroutine CSL_periodic_advect_1d


  subroutine CSL_periodic_advect_1d_constant(&
    adv, &
    A, &
    dt, &
    input, &
    output)
    class(CSL_periodic_1d_advector) :: adv
    sll_real64, intent(in) :: A
    sll_real64, intent(in) :: dt 
    sll_real64, dimension(:), intent(in) :: input
    sll_real64, dimension(:), intent(out) :: output      
    sll_real64, dimension(:), allocatable :: A1
    sll_int32 :: ierr
    
    !this version is not optimized
    
    SLL_ALLOCATE(A1(adv%Npts),ierr)
    
    A1 = A
    
    call adv%charac%compute_characteristics( &
      A1, &
      dt, &
      adv%eta_coords, &
      adv%charac_feet)

!    call adv%interp%compute_interpolants( &
!      input, &
!      adv%eta1_coords, &
!      adv%Npts1, &
!      adv%eta2_coords, &
!      adv%Npts2 )

    output = adv%interp%interpolate_array( &
      adv%Npts, &
      input, &
      adv%charac_feet)      

    SLL_DEALLOCATE_ARRAY(A1,ierr)

          
  end subroutine CSL_periodic_advect_1d_constant


  subroutine delete_CSL_periodic_1d_adv(adv)
    class(CSL_periodic_1d_advector), intent(inout) :: adv
    sll_int32 :: ierr
    SLL_DEALLOCATE(adv%eta_coords,ierr)
    SLL_DEALLOCATE(adv%charac_feet,ierr)
  end subroutine delete_CSL_periodic_1d_adv


  subroutine check_charac_feet( &
    origin, &
    feet, &
    Npts, &
    epsilon)
    sll_real64, dimension(:), intent(in) :: origin
    sll_real64, dimension(:), intent(inout) :: feet
    sll_int32, intent(in) :: Npts
    sll_real64, intent(in), optional :: epsilon
    sll_real64 :: length
    sll_real64 :: gap_min
    sll_real64 :: gap_max
    !sll_int32 :: i
    sll_real64 :: eps
    
    if(size(origin)<Npts)then
      SLL_ERROR("check_charac_feet","bad size for origin")
    endif
    if(size(feet)<Npts)then
      SLL_ERROR("check_charac_feet","bad size for feet")
    endif
    length = origin(Npts)-origin(1)
    feet(Npts) = feet(1) +length !added because errors
    if(abs(length-(feet(Npts)-feet(1)))>1.e-12_f64)then
      print *,'origin',origin
      print *,'feet',feet
      print *,feet(Npts)-feet(1),origin(Npts)-origin(1)      
      SLL_ERROR("check_charac_feet","bad length")
    endif
    
    if(present(epsilon))then
      eps = epsilon
    else
      eps = 1.e-12_f64  
    endif
    
    gap_min = compute_gap_min(origin,Npts)
    gap_max = compute_gap_max(origin,Npts)
    
    !print *,'#gap1',gap_min,gap_max
    if(gap_min<=eps)then
      SLL_ERROR("check_charac_feet","bad gap_min")
    endif    

    gap_min = compute_gap_min(feet,Npts)
    gap_max = compute_gap_max(feet,Npts)
    !print *,'#gap2',gap_min,gap_max
    if(gap_min<=eps)then
      SLL_ERROR("check_charac_feet","pb for charac_feet")
    endif    
    
    
    
  end subroutine check_charac_feet

  function compute_gap_min(input,Npts) result(res)
    sll_real64, dimension(:), intent(in) :: input
    sll_int32, intent(in) :: Npts
    sll_real64 :: res
    sll_int32 :: i
    sll_real64 :: val
    
    if(size(input)<=1)then
      SLL_ERROR("compute_gap_min","bad size for input")      
    endif
    res = input(2)-input(1)
    do i=1,Npts-1
      val = input(i+1)-input(i)
      if(val<res)then
        res = val
      endif
    enddo
    
  end function compute_gap_min

  function compute_gap_max(input,Npts) result(res)
    sll_real64, dimension(:), intent(in) :: input
    sll_int32, intent(in) :: Npts
    sll_real64 :: res
    sll_int32 :: i
    sll_real64 :: val
    
    if(size(input)<=1)then
      SLL_ERROR("compute_gap_max","bad size for input")      
    endif
    res = input(2)-input(1)
    do i=1,Npts-1
      val = input(i+1)-input(i)
      if(val>res)then
        res = val
      endif
    enddo
    
  end function compute_gap_max

  subroutine compute_csl_hermite_mat(d,N,output)
    sll_int32, intent(in) :: d
    sll_int32, intent(in) :: N
    sll_real64, dimension(:), intent(out) :: output
    sll_real64 :: w_left(-d/2:(d+1)/2)
    sll_real64 :: w_right((-d+1)/2:d/2+1)
    !sll_real64 :: tmp
    sll_int32 :: r_left
    sll_int32 :: r_right
    sll_int32 :: s_left
    sll_int32 :: s_right 
    sll_int32 :: i
    sll_int32 :: ind
    
    if(N<1)then
      SLL_ERROR("csl_hermite_mat", "bad size of N")
    endif

    if(d<0)then
      SLL_ERROR("csl_hermite_mat", "bad size of d")
    endif

    
    if(size(output)<N)then
      SLL_ERROR("csl_hermite_mat", "bad size of output")
    endif
    output(1:N) = 0._f64
    output(1) = 0.5_f64
    !output(2) = 0.5_f64
    output(N) = 0.5_f64 !because we have to take a_{-i}
    !instead of a_i, if we want to perform further fft of a

    r_left=-d/2
    s_left=(d+1)/2
    r_right=(-d+1)/2
    s_right=d/2+1
    
    call compute_w_hermite_1d(w_left,r_left,s_left)
    if((2*(d/2)-d)==0)then
      w_right(r_right:s_right) = w_left(r_left:s_left)
    else
      w_right(r_right:s_right) = -w_left(s_left:r_left:-1)
    endif    
    
    do i=r_left,s_left
      ind = 1+modulo(i,N)!1+modulo(-i+N,N)
      ind = 1+modulo(N-ind,N)
      output(ind) = output(ind)+w_left(i)/12._f64
    enddo
    do i=r_right,s_right
      ind = 1+modulo(i,N)!1+modulo(-i+N,N)
      ind = 1+modulo(N-ind,N)
      output(ind) = output(ind)-w_right(i)/12._f64
    enddo
    

    
  end subroutine compute_csl_hermite_mat
  
  subroutine compute_inverse_csl_hermite_mat(N,mat)
    sll_int32, intent(in) :: N
    sll_real64, dimension(:), intent(inout) :: mat
    sll_int32 :: i
    sll_real64 :: rea
    sll_real64 :: ima
    sll_real64 :: val
    
    
    if(N<1)then
      print *,'#N=',N
      SLL_ERROR("csl_hermite_mat", "bad size of N")
    endif
    if(size(mat)<N)then
      print *,'#N=',N
      print *,'#size(mat)=',size(mat)
      SLL_ERROR("csl_hermite_mat", "bad size of mat")
    endif
    
    mat(1) = 1._f64/(mat(1)*real(N,f64))
    do i=1,(N-1)/2
      rea=mat(2*i)
      ima=mat(2*i+1)
      val = 1._f64/((rea**2+ima**2)*real(N,f64))
      mat(2*i) = rea*val
      mat(2*i+1) = -ima*val      
    enddo
    if(mod(N,2)==0)then
      !the matrix is not invertible
      if(abs(mat(N))<1.e-12)then
        mat(N) = 0._f64
      !the matrix is invertible
      else
        if(abs(mat(N))<1.e-6)then
          print *,'#mat(N)=',mat(N)
          SLL_WARNING("compute_inverse_csl_hermite_mat", "mat(N) small")
        endif
        mat(N) = 1._f64/(mat(N)*real(N,f64))        
      endif  
    endif
  
  end subroutine compute_inverse_csl_hermite_mat

  subroutine mult_circulant_mat(N,mat,f,fft_buf)
    sll_int32, intent(in) :: N
    sll_real64, dimension(:), intent(in) :: mat
    sll_real64, dimension(:), intent(inout) :: f
    sll_real64, dimension(:), intent(in) :: fft_buf
    sll_int32 :: i
    sll_real64 :: rea 
    sll_real64 :: ima 
    sll_real64 :: reb 
    sll_real64 :: imb 
    
    call dfftf(N,f,fft_buf)
    f(1)=f(1)*mat(1)
    do i=1,(N-1)/2
      rea=f(2*i)
      ima=f(2*i+1)
      reb=mat(2*i)
      imb=mat(2*i+1)
      f(2*i)=rea*reb-ima*imb
      f(2*i+1)=rea*imb+reb*ima
    enddo
    if(mod(N,2)==0)f(N)=f(N)*mat(N)
    call dfftb(N,f,fft_buf)
        
  end subroutine mult_circulant_mat
  
  subroutine compute_csl_integral( &
    Npts, &
    interp, &
    origin, &
    feet, &
    output)
    sll_int32, intent(in) :: Npts
    class(sll_interpolator_1d_base), pointer :: interp 
    sll_real64, dimension(:), intent(in) :: origin
    sll_real64, dimension(:), intent(in) :: feet
    sll_real64, dimension(:), intent(out) :: output
    sll_real64 :: eta_min
    sll_real64 :: eta_max
    sll_real64 :: delta
    sll_real64 :: val
    sll_real64 :: a
    sll_real64 :: b
    sll_int32 :: i     
    sll_int32 :: ii
    sll_int32 :: i_min     
    sll_int32 :: i_max     
         
    
    if(size(output)<Npts-1)then
      print *,'Npts-1=',Npts-1
      print *,'size(output)=',size(output)
      SLL_ERROR("compute_csl_integral","bad size of output")
    endif
    if(size(origin)<Npts)then
      print *,'Npts=',Npts
      print *,'size(origin)=',size(origin)
      SLL_ERROR("compute_csl_integral","bad size of origin")
    endif
    if(size(feet)<Npts)then
      print *,'Npts=',Npts
      print *,'size(feet)=',size(feet)
      SLL_ERROR("compute_csl_integral","bad size of feet")
    endif
    
    !we assume uniform mesh for origin
    eta_min = origin(1)
    eta_max = origin(Npts)
    delta = (eta_max-eta_min)/real(Npts-1,f64)
    
    do i=1,Npts-1
      i_min = floor((feet(i)-eta_min)/delta)
      i_max = floor((feet(i+1)-eta_min)/delta)
      if(i_min==i_max)then
        a = feet(i)
        b = feet(i+1)
        val = compute_simpson_contribution_csl_periodic(interp,a,b,eta_min,eta_max)
      else
        if(i_min>i_max)then
          print *,'i_min,i_max=',i_min,i_max
          SLL_ERROR("compute_csl_integral","bad value of i_min/i_max")
        endif
        a = feet(i)
        b = eta_min+real(i_min+1,f64)*delta
        val = compute_simpson_contribution_csl_periodic(interp,a,b,eta_min,eta_max)
        do ii=i_min+1,i_max-1
          a = eta_min+real(ii,f64)*delta
          b = eta_min+real(ii+1,f64)*delta  
          val = val+compute_simpson_contribution_csl_periodic(interp,a,b,eta_min,eta_max)
        enddo
        a = eta_min+real(i_max,f64)*delta
        b = feet(i+1)
        val = val+compute_simpson_contribution_csl_periodic(interp,a,b,eta_min,eta_max)           
      endif
      output(i) = val/delta
    enddo
    
    !print *,minval(output),maxval(output)
    
  end subroutine compute_csl_integral

  function compute_simpson_contribution_csl_periodic( &
    interp, &
    a, &
    b, &
    eta_min, &
    eta_max) &
    result(res)
    class(sll_interpolator_1d_base), pointer :: interp
    sll_real64, intent(in) :: a
    sll_real64, intent(in) :: b
    sll_real64, intent(in) :: eta_min
    sll_real64, intent(in) :: eta_max
    sll_real64 :: eta
    sll_real64 :: res
    sll_real64 :: nodes(3,2)
    sll_int32 :: j
    
    nodes(1,1) = a
    nodes(3,1) = b
    nodes(2,1) = 0.5_f64*(a+b)
    do j=1,3
      eta = process_outside_point_periodic( &
        nodes(j,1), &
        eta_min, &
        eta_max)
      nodes(j,2) = interp%interpolate_value(eta)
    enddo
    res = nodes(1,2)+4._f64*nodes(2,2)+nodes(3,2) 
    res = res*(b-a)/6._f64
    
  end function compute_simpson_contribution_csl_periodic

  subroutine compute_csl_mat( &
    interp, &
    eta_min, &
    eta_max, &
    Npts, &
    output)
    class(sll_interpolator_1d_base), pointer :: interp
    sll_real64, intent(in) :: eta_min
    sll_real64, intent(in) :: eta_max
    sll_int32, intent(in) :: Npts
    sll_real64, dimension(:), intent(out) :: output
    sll_int32 :: i
    !sll_int32 :: j
    sll_real64, dimension(:), allocatable :: f
    sll_int32 :: ierr
    sll_real64 :: a    
    sll_real64 :: b
    sll_real64 :: delta
    sll_real64 :: df0    
    sll_real64 :: df1
    sll_real64 :: w_deriv_left(4)    
    sll_real64 :: w_deriv_right(4)    
    sll_int32 :: N
    sll_int32 :: ind
    
    w_deriv_left = (/-5.5_f64,9._f64,-4.5_f64,1._f64 /)
    w_deriv_right = (/-1._f64,4.5_f64,-9._f64,5.5_f64 /)
    
    delta = (eta_max-eta_min)/real(Npts-1,f64) 
    
    N = Npts-1
    
    !we suppose periodic boundary conditions
    SLL_ALLOCATE(f(Npts),ierr)
    
    output(1:N) = 0._f64
    output(1) = 0.5_f64
    output(N) = 0.5_f64 !because we have to take a_{-i}
    !instead of a_i, if we want to perform further fft of a

    f = 0._f64
    call interp%compute_interpolants(f)    
    do i=1,Npts-1
      f = 0._f64
      f(i) = 1._f64
      if(i==1)then
        f(Npts) = 1._f64
      endif
      !we first compute the interpolants
      !for this basis function f
      call interp%compute_interpolants(f)
      !a = eta_min+real(i-1,f64)*delta 
      !b = eta_min+real(i,f64)*delta 
      a = eta_min 
      b = eta_min+delta 
      df0 = compute_quadrature(interp,a,b,w_deriv_left)
      df1 = compute_quadrature(interp,a,b,w_deriv_right)
      ind = i !1+modulo(N-i,N)
      output(ind) = output(ind)+(df0-df1)/12._f64      
      
    enddo
    do i=1,N
      print *,i,output(i)
    enddo
    
  end subroutine compute_csl_mat   


  function compute_quadrature(interp,a,b,w) result(res)
    class(sll_interpolator_1d_base), pointer :: interp
    sll_real64, intent(in) :: a
    sll_real64, intent(in) :: b
    sll_real64, dimension(4), intent(in) :: w
    sll_real64 :: res
    sll_real64 :: eta
    sll_int32 :: i
    
    res = 0._f64
    do i=1,4
      eta = a+(real(i-1,f64)/3._f64)*(b-a)
      res = res+ w(i)*interp%interpolate_value(eta)
    enddo
      
  end function compute_quadrature
  
  subroutine check_solve_csl_mat( &
    csl_mat, &
    csl_mat_init, &
    N, &
    fft_buf)
    sll_real64, dimension(:), intent(in) :: csl_mat
    sll_real64, dimension(:), intent(in) :: csl_mat_init
    sll_int32, intent(in) :: N
    sll_real64, dimension(:), intent(in) :: fft_buf
    sll_real64, dimension(:), allocatable :: f
    sll_real64, dimension(:), allocatable :: g
    sll_real64, dimension(:), allocatable :: h
    sll_int32 :: ierr
    sll_real64 :: err
    sll_int32 :: i
    sll_int32 :: j
    
    SLL_ALLOCATE(f(N),ierr)
    SLL_ALLOCATE(g(N),ierr)
    SLL_ALLOCATE(h(N),ierr)
    err = 0._f64
    
    do i=1,N
      f = 0._f64
      f(i) = 1._f64
      h = f
      call mult_circulant_mat( &
        N, &
        csl_mat, &
        f, &
        fft_buf)
      call circ_mat_mul_direct( &
        csl_mat_init, &
        f, &
        g, &
        N)
      err = max(err,maxval(abs(h-g)))
      print *,'#i='
      do j=1,N
        print *,j,f(j),g(j),h(j)
      enddo    
    enddo
    
    print *,'#err=',err
    stop
  
  end subroutine check_solve_csl_mat  
  
  subroutine circ_mat_mul_direct(a,input,output,N)
    sll_real64, dimension(:), intent(in) :: a
    sll_real64, dimension(:), intent(in) :: input
    sll_real64, dimension(:), intent(out) :: output
    sll_int32, intent(in) :: N
    sll_int32 :: i
    sll_int32 :: j
    sll_int32 :: j2
    sll_int32 :: ind
    
    do i=1,N
      output(i) = 0._f64
      do j=1,N
        ind = j!1+modulo(i+j-2,N)
        j2 = 1+modulo(N-i+1,N) !1+modulo(N-j,N)
        output(i) = output(i)+a(j2)*input(ind)
      enddo
    enddo
    
  end subroutine circ_mat_mul_direct

  !we suppose that compute_interpolant is already done
  subroutine compute_node_derivative_order3( &
    interp, &
    deriv, &
    N, &
    eta_min, &
    eta_max)
    class(sll_interpolator_1d_base), pointer :: interp
    sll_real64, dimension(:,:), intent(out) :: deriv
    sll_int32, intent(in) :: N
    sll_real64, intent(in) :: eta_min    
    sll_real64, intent(in) :: eta_max
    sll_int32 :: i
    sll_real64 :: w_deriv_left(4)    
    sll_real64 :: w_deriv_right(4)
    sll_real64 :: a    
    sll_real64 :: b
    sll_real64 :: delta    
    
    w_deriv_left = (/-5.5_f64,9._f64,-4.5_f64,1._f64 /)
    w_deriv_right = (/-1._f64,4.5_f64,-9._f64,5.5_f64 /)

    
    if(N<1)then
      print *,'N=',N
      SLL_ERROR('compute_node_derivative','bad size of N')
    endif
    
    if((size(deriv,1)<2).or.(size(deriv,2)<N))then
      print *,'#size=',size(deriv)
      print *,'#expected size=',2,N
      SLL_ERROR('compute_node_derivative','bad size of deriv')
    endif
    
    delta= (eta_max-eta_min)/real(N,f64)
        
    do i=1,N
      a = eta_min+real(i-1,f64)*delta 
      b = eta_min+real(i,f64)*delta 
      deriv(1,i) = compute_quadrature(interp,a,b,w_deriv_left)
      deriv(2,i) = compute_quadrature(interp,a,b,w_deriv_right)      
    enddo
      
  end subroutine compute_node_derivative_order3

  subroutine update_solution_csl_periodic( &
    interp, &
    input, &
    deriv, &
    charac, &
    Npts, &
    eta_min, &
    eta_max, &
    output, &
    csl_degree)
    class(sll_interpolator_1d_base), pointer :: interp
    sll_real64, dimension(:), intent(in) :: input
    sll_real64, dimension(:,:), intent(inout) :: deriv
    sll_real64, dimension(:), intent(in) :: charac
    sll_int32, intent(in) :: Npts
    sll_real64, intent(in) :: eta_min
    sll_real64, intent(in) :: eta_max
    sll_real64, dimension(:), intent(out) :: output
    sll_int32, intent(in), optional :: csl_degree
    sll_real64, dimension(:), allocatable :: output_bsl
    sll_real64, dimension(:), allocatable :: flux
    sll_int32, dimension(:), allocatable :: jstar
    sll_real64, dimension(:), allocatable :: alpha
    sll_int32 :: ierr
    sll_real64 :: eta
    sll_int32 :: i
    sll_int32 :: ind
    sll_int32 :: N
    sll_real64 :: dof(4)
    sll_real64 :: res
    sll_real64 :: delta
    sll_real64 :: err
    sll_real64, dimension(:), allocatable :: xi
    sll_real64, dimension(:), allocatable  :: fxi
    sll_real64, dimension(:), allocatable :: xval
    sll_real64, dimension(:), allocatable :: fval
    sll_int32 :: ii
    sll_real64 :: a
    sll_real64 :: b
    sll_real64 :: xstarj
    sll_real64 :: xj
    sll_real64 :: xstarj1
    sll_int32 :: i1
    sll_int32 :: ind1
    sll_real64, dimension(:), allocatable :: ww_tmp
    sll_real64, dimension(:), allocatable :: ww
    sll_int32 :: r
    sll_int32 :: s
    sll_int32 :: d
    sll_int32 :: r1
    sll_int32 :: s1
    sll_int32 :: num_gauss_points
    sll_real64, dimension(:,:), allocatable :: xw
    !sll_real64 :: val
    
    
    
    if(present(csl_degree))then
      d = csl_degree
    else
      d = 3  
    endif
    
    r = -d/2
    s = (d+1)/2
    
    
    SLL_ALLOCATE(ww_tmp(r:s-1),ierr)
    
    call compute_csl_ww(ww_tmp,r,s)
    
    !print *,'d=',d
    !print *,'r=,s=',r,s
    !print *,ww_tmp(r:s-1)

    r1 = -d/2+1
    s1 = (d+1)/2
    SLL_ALLOCATE(ww(r1:s1),ierr)
    
    ww(r1:s1) = ww_tmp(r:s-1)
    
    num_gauss_points = (d+1)/2
    SLL_ALLOCATE(xw(2,num_gauss_points),ierr)
    xw = gauss_legendre_points_and_weights( &
      num_gauss_points, &
      0._f64, &
      1._f64 )
    !print *,'points=',xw(1,:)
    !print *,'weights=',xw(2,:)
    
    !stop
    !print *,'r1,s1=',r1,s1
    !print *,ww(r1:s1)
    
    !stop
    SLL_ALLOCATE(xi(r1:s1),ierr)
    SLL_ALLOCATE(fxi(r1:s1),ierr)
    SLL_ALLOCATE(xval(r1:s1),ierr)
    SLL_ALLOCATE(fval(r1:s1),ierr)
    do ii=r1,s1
      xval(ii) = real(ii,f64)
    enddo
    
    !here is the main work for CSL
    !(Jf)(tn+dt,xj) = Jf(tn,xj)+G'(xj)
    !with G(x) = F1(tn,H(tn;x,tn+dt))-F1(tn,x)
    !we have
    !G(xj) = F1(tn,xstarj)-F1(tn,xj)= int(Jf(tn,x),x=xj..xstarj)
    !we write G(x) = int(G1(s),s=x-dx/2..x+dx/2)/dx
    !so, we get
    !G'(xj) = (G1(xj+dx/2)-G1(xj-dx/2))/dx
    !and thus the scheme is automatically conservative
    !we thus have to compute G1(xj+dx/2)
    !we know that 
    !int(G1(s),s=xk-dx/2..xk+dx/2)/dx = int(Jf(tn,x),x=xk..xstark)
    !we first consider the case where the displacement
    !is less than a cell
    !that is xj <= xstarj < xj+dx
    !or xj-dx <=xstarj < xj
    !if xj <= xstarj < xj+dx
    ! val(ell) = int(G1(s),s=x(j+ell)-dx/2..x(j+ell)+dx/2)/dx = 
    !  = int(Jf(tn,x),x=x(j+ell)..xstar(j+ell)), ell=-1,0,1,2
    !  = (1/6)*(xstar(j+ell)-x(j+ell))
    !  *(Jf(tn,x(j+ell))
    !    +4Jf(tn,(x(j+ell)+xstar(j+ell))/2)
    !     +Jf(tn,xstar(j+ell)))
    ! to compute the integral, we use polynomial of degree <= 3
    ! on [xjstar, xjstar+dx]
    !where xjstar<=xstarj<xjstar+dx
    ! we then have G1(xj+dx/2) 
    !  = (7/12)*(val(0)+val(1))-(1/12)*(val(-1)+val(2))
    !
    ! -1 val(-1) 0 val(0) 1 val(1) 2 val(2) 3
    !  |_________|________|________|_________|
    ! -3/2      -1/2     1/2      3/2       5/2 
    !further step is to generalize for
    ! a displacement bigger than one cell
    !G(xj) = int(Jf(tn,x),x=xj..xjstar)+int(Jf(tn,x),x=xjstar..xstarj)
    !G1(xj+dx/2) = A+G2(x0+dx/2)
    !for G2 everything is shifted to xjstar instead of xj
    !and A = sum(Jf(tn,xk),k=j+1..jstar), if jstar>j
    !and A = -sum(Jf(tn,xk),k=jstar..j-1), if jstar<j
    ! so we know
    ! int(Jf(tn,x),x=xstarj..xjstar)
    ! val(ell) = int(G2(s),s=xell-dx/2..xell+dx/2)/dx = 
    !  = int(Jf(tn,x),x=x(jstar+ell)..xstar(j+ell)), ell=-1,0,1,2
    ! to compute the integral, we use polynomial of degree <= 3
    ! on [xjstar, xjstar+dx]
    ! we then have G2(x0+dx/2) 
    !  = (7/12)*(val(0)+val(1))-(1/12)*(val(-1)+val(2))
    !
    ! -1 val(-1) 0 val(0) 1 val(1) 2 val(2) 3
    !  |_________|________|________|_________|
    ! -3/2      -1/2     1/2      3/2       5/2 
    !
    !we compute the flux at 1/2 (we take x0=0,dx=1)
    !to begin we do a first version that just reproduces BSL
    SLL_ALLOCATE(output_bsl(Npts),ierr)
    SLL_ALLOCATE(flux(0:Npts),ierr)
    SLL_ALLOCATE(jstar(1+r1-10:Npts+s1+10),ierr)
    SLL_ALLOCATE(alpha(1+r1-10:Npts+s1+10),ierr)
    
    N = Npts-1
    
    delta = (eta_max-eta_min)/real(N,f64)
    
    
    output_bsl = input
    
    !added; otherwise small diff with bsl
    output_bsl(Npts) = output_bsl(1)
              
    call interp%compute_interpolants(output_bsl)
    call compute_node_derivative_order3( &
      interp, &
      deriv, &
      Npts-1, &
      eta_min, &
      eta_max)
    
    err = 0._f64
    
    
    !print *,'r1=',r1,s1
    !stop
    
    do i=1,Npts
      eta = charac(i)
      eta = process_outside_point_periodic(eta,eta_min,eta_max)
      eta = real(N,f64)*(eta-eta_min)/(eta_max-eta_min)
      ind = floor(eta)
      eta = eta-ind
      if((eta<0).or.(eta>=1))then
        print *,'eta=',eta
        SLL_ERROR('update_solution_csl_periodic','bad value of eta')
      endif            
      if((ind<0).or.(ind>=N))then
        print *,'ind=',eta
        SLL_ERROR('update_solution_csl_periodic','bad value of ind')
      endif
      ind = ind+1            
      dof(1) = output_bsl(ind)
      dof(2) = output_bsl(ind+1)
      dof(3) = deriv(1,ind)
      dof(4) = deriv(2,ind)
      output(i) = evaluate_hermite_1d(eta,dof)
      eta = charac(i)
      eta = process_outside_point_periodic(eta,eta_min,eta_max)
      res = output(i)-interp%interpolate_value(eta)
      if(abs(res)>1.e-10)then
        print *,'#problem detected'
        print *,'#dof=',dof
        print *,'ind=',ind
        eta = charac(i)
        eta = process_outside_point_periodic(eta,eta_min,eta_max)
        eta = real(N,f64)*(eta-eta_min)/(eta_max-eta_min)
        ind = floor(eta)
        eta = eta-ind
        print *,'#eta=',eta
        stop
      endif
      err = max(err,abs(res))
    enddo
    
    !check that the hermite  version of bsl
    !leads to the same result 
    if(err>1.e-14)then
      print *,'#err bsl=',err
    endif
    output_bsl = output
    output = output_bsl
    output_bsl = input
    !now we construct the new scheme
    
    !do i=1,Npts
    !  output_bsl(i) = 1._f64
    !enddo

    call interp%compute_interpolants(output_bsl)
    call compute_node_derivative_order3( &
      interp, &
      deriv, &
      Npts-1, &
      eta_min, &
      eta_max)


    
    !we compute jstar
    do i=1,Npts
      eta = charac(i)
      eta = real(N,f64)*(eta-eta_min)/(eta_max-eta_min)
      jstar(i) = floor(eta)
      alpha(i) = eta-jstar(i)
      jstar(i) = jstar(i)+1
      if((alpha(i)<0._f64).or.(alpha(i)>=1.))then
        print *,'alpha(i)=',i,alpha(i)
        SLL_ERROR('update_solution_csl_periodic','bad value of alpha(i)')
      endif
      !print *,i,jstar(i)-i,jstar(i)
    enddo
    !stop
    err = abs(real(jstar(Npts)-jstar(1)-N,f64))
    err = max(err,abs(alpha(N+1)-alpha(1)))
    if(err>1.e-13)then
      print *,'#err periodicity=',err
    endif
    
    !enforce periodicity
    jstar(N+1) = jstar(1)+N
    alpha(N+1) = alpha(1)
    
    !boundary conditions
    do ii=r1,0
      jstar(ii) = jstar(N+ii)-N
      alpha(ii) = alpha(N+ii)
    enddo
    do ii=1,1+s1
      jstar(N+ii) = jstar(ii)+N
      alpha(N+ii) = alpha(ii)
    enddo
!    jstar(0) = jstar(N)-N
!    jstar(-1) = jstar(N-1)-N
!    jstar(-2) = jstar(N-2)-N
!    jstar(-3) = jstar(N-3)-N
!    jstar(N+2) = jstar(2)+N
!    jstar(N+3) = jstar(3)+N
!    jstar(N+4) = jstar(4)+N
!    alpha(0) = alpha(N)
!    alpha(-1) = alpha(N-1)
!    alpha(N+2) = alpha(2)
!    alpha(N+3) = alpha(3)
    
    !we then check that the jstar and alpha 
    !are well computed
    
    err = 0._f64
    do i=1,Npts
      eta = real(jstar(i)-1,f64)+alpha(i)
      eta = eta_min + eta*delta
      err = max(err,abs(eta-charac(i)))
    enddo    
    
    if(err>1.e-14)then
      print *,'#err charac=',err
    endif
    

    
    !now we try order 3 (?)

    do i=1,N
      !xj = eta_min+real(i-1,f64)*delta
      !xstarj = charac(i)      
      !xstarj = eta_min+(jstar(i)+alpha(i)-1)*delta

      !we now use lagrange interpolation of degree d-1
      !for d=4
      !-1,0,1,2
      !for d=3
      !0,1,2
      !->r1,s1
      !for the moment we use Hermite

      ind = 1+modulo(N+jstar(i)-1,N)
      ind1 = 1+modulo(N+jstar(i),N)
      dof(1) = output_bsl(ind)
      dof(2) = output_bsl(ind1)
      dof(3) = deriv(1,ind)
      dof(4) = deriv(2,ind)
      
      do ii=r1,s1
        fval(ii) = output_bsl(1+modulo(N+jstar(i)+ii-1,N))
      enddo
       

      do ii=r1,s1
        xi(ii) = jstar(i+ii)+alpha(i+ii)-jstar(i)
      enddo

      do ii=r1,s1
        a = real(ii,f64)
        b = xi(ii)
        fxi(ii) = contribution_gauss_lagrange( &
          a, &
          b, &
          xval, &
          fval, &
          r1, &
          s1, &
          xw, &
          num_gauss_points)
!        val = contribution_simpson_hermite(a,b,dof)
!        if(abs(val-fxi(ii))>1.e-12)then
!          print *,fxi(ii),val
!          print *,'a,b=',a,b
!          print *,'dof=',dof
!          print *,'xval=',xval
!          print *,'fval=',fval
!          stop
!        endif
      enddo
      
      flux(i) = 0._f64
      do ii=r1,s1
        flux(i) = flux(i)+ww(ii)*fxi(ii)
      enddo    
      !flux(i) = (7._f64/12._f64)*(fxi(0)+fxi(1))
      !flux(i) = flux(i)-(1._f64/12._f64)*(fxi(-1)+fxi(2))
      
      !flux(i) = fxi(0)

      !flux(i) = (37._f64/60._f64)*(fxi(0)+fxi(1))
      !flux(i) = flux(i)-(8._f64/60._f64)*(fxi(-1)+fxi(2))
      !flux(i) = flux(i)+(1._f64/60._f64)*(fxi(-2)+fxi(3))
      
      
      
      !flux(i) = (6._f64/12._f64)*(fxi(0)+fxi(1))
      !flux(i) = flux(i)-(1._f64/12._f64)*(fxi(-1)+fxi(2))

      do ii=i+1,jstar(i)
        ind1 = 1+modulo(ii+N-1,N)
        flux(i) = flux(i)+output_bsl(ind1)
      enddo  
      do ii=jstar(i)+1,i
        ind1 = 1+modulo(ii+N-1,N)
        flux(i) = flux(i)-output_bsl(ind1)
      enddo
        
    enddo

    flux(N+1) = flux(1)
    
    
    do i=1,N
      ind1 = 1+modulo(i-1+N-1,N)
      output(i) = output_bsl(i)+(flux(i)-flux(ind1))
    enddo
    output(N+1)=output(1)

    return



    !we do a first version that is subject to CFL condition
    ! this is the simplest scheme that should always work
    
    do i=1,Npts
      xj = eta_min+real(i-1,f64)*delta
      xstarj = charac(i)      
      if(xstarj.ge.xj)then
        if(xstarj.ge.(xj+delta))then
          print *,'xstarj>=xj+delta'
          print *,'xstarj=',xstarj
          print *,'xj+delta=',xj+delta
          print *,'i=',i
          SLL_ERROR('update_solution_csl_periodic','dt may be too big')
        endif
        flux(i) = output_bsl(1+modulo(i,N))*(xstarj-xj)/delta !to begin low order approx
      else
        if(xstarj.le.(xj-delta))then
          print *,'xstarj<=xj-delta'
          print *,'xstarj=',xstarj
          print *,'xj-delta=',xj-delta
          print *,'i=',i
          SLL_ERROR('update_solution_csl_periodic','dt may be too big')
        endif
        flux(i) = output_bsl(i)*(xstarj-xj)/delta        
      endif  
    enddo
    
    do i=1,N
      output(i) = output_bsl(i)+(flux(i)-flux(1+modulo(i-1+N-1,N)))
    enddo
    output(N+1)=output(1)
    
    
    
    !now we try something of order 1 for the function
    !previous version was order 0 for function


    do i=1,N
      xj = eta_min+real(i-1,f64)*delta
      xstarj = charac(i)      
      if(xstarj.ge.xj)then
        if(xstarj.ge.(xj+delta))then
          print *,'xstarj>=xj+delta'
          print *,'xstarj=',xstarj
          print *,'xj+delta=',xj+delta
          print *,'i=',i
          SLL_ERROR('update_solution_csl_periodic','dt may be too big')
        endif
        i1 = 1+modulo(i,N)
        a = 0.5_f64*(xj+xstarj)
        fxi(0) =  output_bsl(i)*(xj+delta-a)/delta
        fxi(0) = fxi(0)+output_bsl(i1)*(1._f64-(xj+delta-a)/delta)
        fxi(0) = fxi(0)*(xstarj-xj)/delta        
        xstarj1 = real(jstar(i+1)-1,f64)+alpha(i+1)
        xstarj1 = eta_min + xstarj1*delta
        !do not use charac(i1) because charac is not periodic
        a = 0.5_f64*(xj+delta+xstarj1)
        fxi(1) =  output_bsl(i)*(xj+delta-a)/delta
        fxi(1) = fxi(1)+output_bsl(i1)*(1._f64-(xj+delta-a)/delta)
        fxi(1) = fxi(1)*(xstarj1-xj-delta)/delta
        flux(i) =  0.5_f64*(fxi(0)+fxi(1))
      else
        !SLL_ERROR('update_solution_csl_periodic','temporary we do not want this case') 
        if(xstarj.le.(xj-delta))then
          print *,'xstarj<=xj-delta'
          print *,'xstarj=',xstarj
          print *,'xj-delta=',xj-delta
          print *,'i=',i
          SLL_ERROR('update_solution_csl_periodic','dt may be too big')
        endif
        flux(i) = output_bsl(i)*(xstarj-xj)/delta        
      endif  
    enddo
    
    flux(N+1) = flux(1)
    
    
    do i=1,N
      output(i) = output_bsl(i)+(flux(i)-flux(1+modulo(i-1+N-1,N)))
    enddo
    output(N+1)=output(1)



!
!
!    
!      !if(xstarj.ge.xj)then
!      if(jstar(i).ge.i)then
!        if(jstar(i)>i)then
!        !if(xstarj.ge.(xj+delta))then
!          print *,'xstarj>=xj+delta'
!          print *,'xstarj=',xstarj
!          print *,'xj+delta=',xj+delta
!          print *,'i=',i
!          SLL_ERROR('update_solution_csl_periodic','dt may be too big')
!        endif
!        ind = 1+modulo(N+jstar(i)-1,N)
!        ind1 = 1+modulo(N+jstar(i),N)
!        dof(1) = output_bsl(ind)
!        dof(2) = output_bsl(ind1)
!        dof(3) = deriv(1,ind)
!        dof(4) = deriv(2,ind)
!
!        !we have jstar(i) = i
!        do ii=-1,2
!          !xi(ii) = jstar(i+ii)+alpha(i+ii)-i
!          xi(ii) = jstar(i+ii)+alpha(i+ii)-jstar(i)
!          !print *,'xi(ii)=',ii,xi(ii)
!        enddo
!        !stop
!        
!        do ii=-1,2
!          a = real(ii,f64)
!          b = xi(ii)
!          fxi(ii) = contribution_simpson_hermite(a,b,dof)
!        enddo    
!        flux(i) = (7._f64/12._f64)*(fxi(0)+fxi(1))
!        flux(i) = flux(i)-(1._f64/12._f64)*(fxi(-1)+fxi(2))
!        !if jstar(i)>i ->+sum(j=i+1..jstar(i))
!        do ii=i+1,jstar(i)
!          ind1 = 1+modulo(ii+N-1,N)
!          flux(i) = flux(i)+output_bsl(ind1)
!        enddo  
!
!      else
!        !SLL_ERROR('update_solution_csl_periodic','temporary we do not want this case') 
!        if(jstar(i)<i-1)then
!        !if(xstarj.le.(xj-delta))then
!          print *,'xstarj<=xj-delta'
!          print *,'xstarj=',xstarj
!          print *,'xj-delta=',xj-delta
!          print *,'i=',i
!          SLL_ERROR('update_solution_csl_periodic','dt may be too big')
!        endif
!
!        ind = 1+modulo(N+jstar(i)-1,N)
!        ind1 = 1+modulo(N+jstar(i),N)
!        if(ind1.ne.i)then
!          print *,'#ind1 should be equal to i'
!          print *,'#ind1=',ind1
!          print *,'#i=',i
!          print *,'#xstarj=',xstarj          
!          print *,'#jstar(i)=',jstar(i)
!          print *,'#alpha(i)=',alpha(i)
!          print *,'#xstarj2=',eta_min+(jstar(i)+alpha(i)-1)*delta
!          SLL_ERROR('update_solution_csl_periodic',"bad value of ind")
!        endif
!        
!        dof(1) = output_bsl(ind)
!        dof(2) = output_bsl(ind1)
!        dof(3) = deriv(1,ind)
!        dof(4) = deriv(2,ind)
!        
!        
!        !we have jstar(i) = i-1
!        do ii=-1,2
!          !xi(ii) = jstar(i+ii)+alpha(i+ii)-i+1
!          xi(ii) = jstar(i+ii)+alpha(i+ii)-jstar(i)
!          !i-jstar(i-ii)-alpha(i-ii)
!          !print *,'xi(ii)=',ii,xi(ii)
!        enddo
!        !stop
!        do ii=-1,2
!          a = real(ii,f64)
!          b = xi(ii)
!          fxi(ii) = contribution_simpson_hermite(a,b,dof)
!        enddo    
!        flux(i) = (7._f64/12._f64)*(fxi(0)+fxi(1))
!        flux(i) = flux(i)-(1._f64/12._f64)*(fxi(-1)+fxi(2))
!        !flux(i) = flux(i)-output_bsl(ind1)
!        do ii=jstar(i)+1,i
!          ind1 = 1+modulo(ii+N-1,N)
!          flux(i) = flux(i)-output_bsl(ind1)
!        enddo  
!        !ind1 = 1+modulo(N+jstar(i),N)
!        !if jstar(i)=i-2 -> -output_bsl(i)-output_bsl(i-1)
!        !if jstar(i)=i-1 -> -output_bsl(i)
!        !if jstar(i)<i ->-sum(j=jstar(i)+1..i)
!        !if jstar(i)=i -> 0
!        !if jstar(i)>i ->+sum(j=i+1..jstar(i))
!        !if jstar(i)=i+1 ->+output_bsl(i+1)
!        !if jstar(i)=i+2 ->+output_bsl(i+1)+output_bsl(i+2)        
!                
!      endif  
!    enddo
    



    
    

    
    
    
    




    
    do i=1,Npts

      ind = 1+modulo(N+jstar(i)-1,N)
      dof(1) = output_bsl(ind)
      dof(2) = output_bsl(ind+1)
      dof(3) = deriv(1,ind)
      dof(4) = deriv(2,ind)

      do ii=-1,2
        xi(ii) = jstar(i+ii)+alpha(i+ii)-jstar(i)
      enddo
      do ii=-1,2
        a = real(ii,f64)
        b = xi(ii)
        fxi(ii) = contribution_simpson_hermite(a,b,dof)
      enddo    
      flux(i) = (7._f64/12._f64)*(fxi(0)+fxi(1))
      flux(i) = flux(i)-(1._f64/12._f64)*(fxi(-1)+fxi(2))

    !and A = sum(Jf(tn,xk),k=j+1..jstar), if jstar>j
    !and A = -sum(Jf(tn,xk),k=jstar..j-1), if jstar<j
      do ii=i+1,jstar(i)
        ind = 1+modulo(N+jstar(ii)-1,N)
        flux(i) = flux(i)+output_bsl(ind)
      enddo
      do ii=jstar(i),i-1
        ind = 1+modulo(N+jstar(ii)-1,N)
        flux(i) = flux(i)-output_bsl(ind)
      enddo
      
      !print *,i,fxi
    enddo
    
    flux(0) = flux(N)
    do i=1,N
      output(i) = input(i)+(flux(i)-flux(i-1)) !/delta
    enddo
    
    output(N+1) = output(1)
    
    err = maxval(abs(output-output_bsl)) 
    
    if(err>1.e-14)then
      print *,'#diff with bsl=',err
      
      do i=1,Npts
        print *,i,output(i)
      enddo
      stop
    endif
    
    
        
    !output = output_bsl
    
    !stop
    
  end subroutine update_solution_csl_periodic
  

 function evaluate_hermite_1d(x,dof) result(res)
    sll_real64, intent(in)::x
    sll_real64, dimension(:), intent(in) :: dof
    sll_real64 :: res
    sll_real64 :: w(4)

    w(1)=(2._f64*x+1)*(1._f64-x)*(1._f64-x)
    w(2)=x*x*(3._f64-2._f64*x)
    w(3)=x*(1._f64-x)*(1._f64-x)
    w(4)=x*x*(x-1._f64)
    res=dof(1)*w(1)+dof(2)*w(2)+dof(3)*w(3)+dof(4)*w(4)

  end function evaluate_hermite_1d

  function contribution_simpson_hermite( &
    a, &
    b, &
    dof) &
    result(res)
    sll_real64, intent(in) :: a
    sll_real64, intent(in) :: b
    sll_real64, dimension(:), intent(in) :: dof
    !sll_real64 :: eta
    sll_real64 :: res
    sll_real64 :: nodes(3,2)
    sll_int32 :: j
    
    nodes(1,1) = a
    nodes(3,1) = b
    nodes(2,1) = 0.5_f64*(a+b)
    do j=1,3
      nodes(j,2) = evaluate_hermite_1d(nodes(j,1),dof)
    enddo
    res = nodes(1,2)+4._f64*nodes(2,2)+nodes(3,2) 
    res = res*(b-a)/6._f64
    
  end function contribution_simpson_hermite


  subroutine compute_csl_ww(ww,r,s)
    sll_int32, intent(in) :: r
    sll_int32, intent(in) :: s
    sll_real64, dimension(r:s-1), intent(out) :: ww
    sll_real64, dimension(:), allocatable :: w
    sll_real64 :: tmp
    sll_int32 :: i
    sll_int32 :: j
    sll_int32 :: ierr

    SLL_ALLOCATE(w(r:s),ierr)

    ww = 0._f64
    do i=r,-1
      tmp=1._f64
      do j=r,i-1
        tmp=tmp*real(i-j,f64)
      enddo
      do j=i+1,s
        tmp=tmp*real(i-j,f64)
      enddo
      tmp=1._f64/tmp
      do j=r,i-1
        tmp=tmp*real(-j,f64)
      enddo
      do j=i+1,-1
        tmp=tmp*real(-j,f64)
      enddo
      do j=1,s
        tmp=tmp*real(-j,f64)
      enddo
      w(i)=tmp      
    enddo

    do i=1,s
      tmp=1._f64
      do j=r,i-1
        tmp=tmp*real(i-j,f64)
      enddo
      do j=i+1,s
        tmp=tmp*real(i-j,f64)
      enddo
      tmp=1._f64/tmp
      do j=r,-1
        tmp=tmp*real(-j,f64)
      enddo
      do j=1,i-1
        tmp=tmp*real(-j,f64)
      enddo
      do j=i+1,s
        tmp=tmp*real(-j,f64)
      enddo
      w(i)=tmp      
    enddo

    tmp=0._f64
    do i=r,-1
      tmp=tmp+w(i)
    enddo
    do i=1,s
      tmp=tmp+w(i)
    enddo
    w(0)=-tmp

    tmp=0._f64
    do i=r,-1
      tmp=tmp+w(i)
      ww(i)=-tmp
    enddo
    tmp=0._f64
    do i=s,1,-1
      tmp=tmp+w(i)
      ww(i-1)=tmp
    enddo
    
    !print *,'r,s=',r,s
    !print *,'w=',w

    !print *,'ww=',ww
    !SLL_DEALLOCATE_ARRAY(w,ierr)
    
  end subroutine compute_csl_ww


  function contribution_gauss_lagrange(a,b,xval,fval,r,s,xw,ng) result(res)
    sll_real64, intent(in) :: a
    sll_real64, intent(in) :: b
    sll_int32, intent(in) :: r
    sll_int32, intent(in) :: s
    sll_real64, dimension(r:s), intent(in) :: xval
    sll_real64, dimension(r:s), intent(in) :: fval
    sll_real64, dimension(:,:), intent(in) ::xw
    sll_int32, intent(in) :: ng
    sll_real64 :: res
    sll_int32 :: i
    sll_int32 :: d
    sll_real64 :: x
    !sll_real64 :: nodes(3)
    
    !print *,'a=',a,b
    !print *,'xval=',xval(r:s)
    !print *,'r=',r,s
    d = s-r
    res = 0._f64
    !do i=r,s
    !  print *,lagrange_interpolate(xval(i),d,xval(r:s),fval(r:s))
    !enddo
		!stop
	!x = a	
	!nodes(1) = lagrange_interpolate(x,d,xval(r:s),fval(r:s))	
	!x = 0.5_f64*(a+b)	
	!nodes(2) = lagrange_interpolate(x,d,xval(r:s),fval(r:s))	
	!x = b	
	!nodes(3) = lagrange_interpolate(x,d,xval(r:s),fval(r:s))	

    !res = nodes(1)+4._f64*nodes(2)+nodes(3) 
    !res = res*(b-a)/6._f64
		
	!print *,'res simpson=',res	

    res = 0._f64
    do i=1,ng
      x = a+xw(1,i)*(b-a)
      !print *,'x=',i,x
      res = res+(b-a)*xw(2,i)*lagrange_interpolate(x,d,xval(r:s),fval(r:s))
    enddo

	!print *,'res gauss=',res	


    
    
    
    
  end function contribution_gauss_lagrange  
  

end module sll_m_advection_1d_CSL_periodic
