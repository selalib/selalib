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


module sll_module_advection_1d_CSL_periodic
#include "sll_working_precision.h"
#include "sll_memory.h"
#include "sll_assert.h"
#include "sll_errors.h"
use sll_boundary_condition_descriptors
use sll_module_advection_1d_base
use sll_module_characteristics_1d_base
use sll_module_interpolators_1d_base
use sll_hermite_interpolation_1d_module

implicit none

  type,extends(sll_advection_1d_base) :: CSL_periodic_1d_advector
  
    class(sll_interpolator_1d_base), pointer  :: interp
    class(sll_characteristics_1d_base), pointer  :: charac
    sll_real64, dimension(:), pointer :: eta_coords
    sll_real64, dimension(:), pointer :: charac_feet
    sll_real64, dimension(:), pointer :: charac_feet_inside
    sll_real64, dimension(:), pointer :: csl_mat_init
    sll_real64, dimension(:), pointer :: csl_mat
    sll_real64, dimension(:), pointer :: fft_buf
    sll_int32 :: Npts
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
    hermite_degree) &  
    result(adv)      
    type(CSL_periodic_1d_advector), pointer :: adv
    class(sll_interpolator_1d_base), pointer :: interp
    class(sll_characteristics_1d_base), pointer  :: charac
    sll_int32, intent(in) :: Npts
    sll_real64, intent(in), optional :: eta_min
    sll_real64, intent(in), optional :: eta_max
    sll_real64, dimension(:), pointer, optional :: eta_coords
    sll_int32, intent(in), optional :: hermite_degree
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
      hermite_degree)    
    
  end function  new_CSL_periodic_1d_advector


  subroutine initialize_CSL_periodic_1d_advector(&
    adv, &
    interp, &
    charac, &
    Npts, &
    eta_min, &
    eta_max, &
    eta_coords, &
    hermite_degree)    
    class(CSL_periodic_1d_advector), intent(inout) :: adv
    class(sll_interpolator_1d_base), pointer :: interp
    class(sll_characteristics_1d_base), pointer  :: charac
    sll_int32, intent(in) :: Npts
    sll_real64, intent(in), optional :: eta_min
    sll_real64, intent(in), optional :: eta_max
    sll_real64, dimension(:), pointer, optional :: eta_coords
    sll_int32, optional :: hermite_degree
    sll_int32 :: ierr
    sll_int32 :: i
    sll_real64 :: delta_eta
    sll_real64 :: hermite_inversibility     
    
    hermite_inversibility = 1._f64
    
    adv%Npts = Npts
    adv%interp => interp
    adv%charac => charac
    SLL_ALLOCATE(adv%eta_coords(Npts),ierr)

    SLL_ALLOCATE(adv%charac_feet(Npts),ierr)
    SLL_ALLOCATE(adv%charac_feet_inside(Npts),ierr)
    SLL_ALLOCATE(adv%csl_mat(Npts-1),ierr)
    SLL_ALLOCATE(adv%csl_mat_init(Npts-1),ierr)
    SLL_ALLOCATE(adv%fft_buf(2*(Npts-1)+15),ierr)
    
    call dffti(Npts-1,adv%fft_buf)
    
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
    
    call compute_csl_mat( &
      interp, &
      adv%eta_coords(1), &
      adv%eta_coords(Npts), &
      Npts, &
      adv%csl_mat)
    
    
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
      do i=1,Npts-1
        adv%csl_mat_init(i) = adv%csl_mat(i)
      enddo
      
      
      
      call dfftf(Npts-1,adv%csl_mat,adv%fft_buf)
      hermite_inversibility = abs(adv%csl_mat(1))
      !print *,0,adv%csl_mat(1)
      do i=1,(Npts-1)/2
        !print *,i-1,adv%csl_mat(2*(i-1)),adv%csl_mat(2*(i-1)+1), &
        hermite_inversibility = min(hermite_inversibility, &
          sqrt(adv%csl_mat(2*i)**2+adv%csl_mat(2*i+1)**2))
      enddo
      !if(mod(Npts-1,2)==0)then
      !  print *,(Npts-1)/2,adv%csl_mat(Npts-1)
      !endif
      
      print *,'#hermite_inversibility=',hermite_inversibility
      if(hermite_inversibility<1.-12)then
        print *,'#csl_mat not invertible',hermite_inversibility
        stop
      endif
      call compute_inverse_csl_hermite_mat(Npts-1,adv%csl_mat)



      call check_solve_csl_mat( &
        adv%csl_mat, &
        adv%csl_mat_init, &
        Npts-1, &
        adv%fft_buf)
      
      
      
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

    call compute_csl_integral( &
      adv%Npts, &
      adv%interp, &
      adv%eta_coords, &
      adv%charac_feet, &
      output(1:adv%Npts-1))
!    if(maxval(input)>1.e-1)then
!      print *,'before'
!      do i=1,adv%Npts-1
!        print *,i,input(i),output(i),0.5_f64*(input(i)+input(i+1))
!      enddo
!!      stop
!    endif             
    
    call mult_circulant_mat( &
      adv%Npts-1, &
      adv%csl_mat, &
      output(1:adv%Npts-1), &
      adv%fft_buf)
    
    output(adv%Npts) = output(1)
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
    sll_real64, dimension(:), intent(in) :: feet
    sll_int32, intent(in) :: Npts
    sll_real64, intent(in), optional :: epsilon
    sll_real64 :: length
    sll_real64 :: gap_min
    sll_real64 :: gap_max
    sll_int32 :: i
    sll_real64 :: eps
    
    if(size(origin)<Npts)then
      SLL_ERROR("check_charac_feet","bad size for origin")
    endif
    if(size(feet)<Npts)then
      SLL_ERROR("check_charac_feet","bad size for feet")
    endif
    length = origin(Npts)-origin(1)
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
    sll_real64 :: tmp
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
    sll_int32 :: j
    sll_real64, dimension(:), allocatable :: f
    sll_int32 :: ierr
    sll_real64 :: a    
    sll_real64 :: b
    sll_real64 :: delta
    sll_real64 :: df0    
    sll_real64 :: df1
    sll_real64 :: w_deriv_left(4)    
    sll_real64 :: w_deriv_right(4)    
    sll_real64 :: N
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
      eta = a+(real(i-1)/3._f64)*(b-a)
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
  

end module sll_module_advection_1d_CSL_periodic
