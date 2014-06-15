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

! for the moment mimic of sll_periodic_interpolator_1d.F90

module sll_module_advection_1d_non_uniform_cubic_splines
#include "sll_working_precision.h"
#include "sll_memory.h"
#include "sll_assert.h"
use sll_boundary_condition_descriptors
use sll_module_advection_1d_base
use cubic_non_uniform_splines

implicit none

  type,extends(sll_advection_1d_base) :: non_uniform_cubic_splines_1d_advector

     sll_int32                            :: num_cells
     sll_real64                           :: xmin
     sll_real64                           :: xmax
     sll_real64, dimension(:), pointer :: node_positions
     sll_real64, dimension(:), pointer :: buf
     sll_int32, dimension(:), pointer :: ibuf
     sll_real64, dimension(:), pointer :: node_pos
     sll_real64, dimension(:), pointer :: coeffs
     sll_real64, dimension(:), pointer :: Xstar
  contains
    procedure, pass(adv) :: initialize => &
      initialize_non_uniform_cubic_splines_1d_advector
    procedure, pass(adv) :: advect_1d_constant => &
      non_uniform_cubic_splines_advect_1d_constant
    procedure, pass(adv) :: advect_1d => &
      non_uniform_cubic_splines_advect_1d
  
  end type non_uniform_cubic_splines_1d_advector
   
contains
  

  function new_non_uniform_cubic_splines_1d_advector(&
    num_cells, &
    xmin, &
    xmax, &
    order, &
    node_positions &
    ) &
    result(adv)      
    type(non_uniform_cubic_splines_1d_advector), pointer :: adv
    sll_int32,  intent(in)               :: num_cells
    sll_real64, intent(in)               :: xmin
    sll_real64, intent(in)               :: xmax
    sll_int32,  intent(in)               :: order
    sll_real64, dimension(:), intent(in), optional :: node_positions 
    sll_int32 :: ierr

    SLL_ALLOCATE(adv,ierr)
    call initialize_non_uniform_cubic_splines_1d_advector(&
      adv, &
      num_cells, &
      xmin, &
      xmax, &
      order, &
      node_positions)
    
  end function new_non_uniform_cubic_splines_1d_advector

  
  subroutine initialize_non_uniform_cubic_splines_1d_advector(&
      adv, &
      num_cells, &
      xmin, &
      xmax, &
      order, &
      node_positions)
      
    class(non_uniform_cubic_splines_1d_advector) :: adv
    sll_int32,  intent(in)               :: num_cells
    sll_real64, intent(in)               :: xmin
    sll_real64, intent(in)               :: xmax
    sll_int32,  intent(in)               :: order
    sll_real64, dimension(:), intent(in), optional :: node_positions 
    sll_int32 :: ierr
    sll_int32 :: i
    sll_real64 :: dx

!    call initialize_periodic_interp( &
!      adv%per_interp, &
!      num_cells, &
!      type, &
!      order)

    adv%num_cells = num_cells 
    adv%xmin = xmin 
    adv%xmax = xmax
    
    dx = (xmax-xmin)/real(num_cells,f64)
    
    if(order.ne.4)then
      print *,'#Warning order=4 is enforced'
      print *,'#in initialize_non_uniform_cubic_splines_1d_advector'
    endif
    
    if(present(node_positions))then
      if(size(node_positions,1)<num_cells+1)then
        print *,'#size problem for node_positions'
        print *,'#in subroutine initialize_non_uniform_cubic_splines_1d_advector'
        stop
      endif
      SLL_ALLOCATE(adv%node_positions(num_cells+1),ierr)
      adv%node_positions(1:num_cells+1) = &
        (node_positions(1:num_cells+1)-xmin)/(xmax-xmin)      
    else
      SLL_ALLOCATE(adv%node_positions(num_cells+1),ierr)
      do i=1,num_cells+1
        !adv%node_positions(i) = xmin + real(i-1,f64)*dx
        adv%node_positions(i) = real(i-1,f64)/real(num_cells,f64)
      enddo  
    endif

    SLL_ALLOCATE(adv%buf(10*num_cells),ierr)
    SLL_ALLOCATE(adv%ibuf(num_cells),ierr)
    SLL_ALLOCATE(adv%node_pos(-2:num_cells+2),ierr)
    SLL_ALLOCATE(adv%coeffs(-2:num_cells+2),ierr)
    SLL_ALLOCATE(adv%Xstar(1:num_cells+1),ierr)

 
  
  end subroutine initialize_non_uniform_cubic_splines_1d_advector   








  subroutine non_uniform_cubic_splines_advect_1d_constant(&
    adv, &
    A, &
    dt, &
    input, &
    output)
    class(non_uniform_cubic_splines_1d_advector) :: adv
    sll_real64, intent(in) :: A
    sll_real64, intent(in) :: dt 
    sll_real64, dimension(:), intent(in) :: input
    sll_real64, dimension(:), intent(out) :: output      
    sll_real64 :: shift
    sll_real64 :: xmin
    sll_real64 :: xmax
    sll_int32  :: num_cells
    sll_real64 :: alpha
      
    num_cells = adv%num_cells
    xmin = adv%xmin
    xmax = adv%xmax
    shift = A*dt/(xmax-xmin)*real(num_cells,f64)
    
    alpha = A*dt/(xmax-xmin)
    
    output(1:num_cells+1) = input(1:num_cells+1)

    call constant_advection_spl_non_unif_per( &
      output, &
      alpha, &
      adv%node_positions, &
      num_cells, &
      adv%buf, &
      adv%Xstar, &
      adv%node_pos, &
      adv%coeffs, &
      adv%ibuf)
    !print *,'#not implemented for the moment'
    !print *,'#non_uniform_cubic_splines_advect_1d_constant'
    !stop
      
!    call periodic_interp( &
!      adv%per_interp, &
!      output, &
!      input, &
!      shift)
!    ! complete by periodicity
!    output(num_cells+1) = output(1)
      
  end subroutine non_uniform_cubic_splines_advect_1d_constant

  subroutine non_uniform_cubic_splines_advect_1d(&
    adv, &
    A, &
    dt, &
    input, &
    output)
    class(non_uniform_cubic_splines_1d_advector) :: adv
    sll_real64, dimension(:), intent(in) :: A
    sll_real64, intent(in) :: dt 
    sll_real64, dimension(:), intent(in) :: input
    sll_real64, dimension(:), intent(out) :: output      
    
    print *,'#non_uniform_cubic_splines_advect_1d'
    print *,'#not implemented for the moment'
    print *,maxval(A)
    print *,dt
    print *,maxval(input)
    output= 0._f64
    print *,adv%num_cells
    stop
           
  end subroutine non_uniform_cubic_splines_advect_1d




  subroutine constant_advection_spl_non_unif_per(f, &
    alpha, &
    node_positions, &
    N, &
    buf, &
    Xstar, &
    node_pos, &
    coeffs, &
    ibuf)
    !alpha and node_positions are normalized to [0,1]
    !use numeric_constants
    !use cubic_non_uniform_splines
    implicit none
    
    sll_real64, dimension(:), intent(inout) :: f
    sll_real64, dimension(:), intent(in) :: node_positions
    !type(cubic_nonunif_spline_1D), pointer :: spl_per
    sll_int32,intent(in):: N
    sll_real64,intent(in)::alpha
    sll_real64 :: dx
    sll_int32  :: i
    !sll_real64 :: M,tmp,tmp2
    !temporary allocations
    sll_real64,dimension(:),pointer :: buf,Xstar,node_pos,coeffs
    sll_int32,dimension(:),pointer :: ibuf 
    
    
    dx = 1._f64/real(N,f64)
    
    !allocate(buf(10*N))
    !allocate(ibuf(N))
    !allocate(node_pos(-2:N+2),coeffs(-2:N+2))
    !allocate(Xstar(1:N+1))
    !print *,loc(buf)
    
    
    node_pos(0:N)=node_positions(1:N+1)
    
    
    !do i=1,N+1
    !  print *,i,node_positions(i)
    !enddo

    !do i=1,N+1
    !  print *,i,Xstar(i),f(i)
    !enddo
    
    !print *,dx
    
    do i=1,N+1    
      Xstar(i) = node_positions(i)-alpha
    enddo
    
    
    do i=1,N+1
      do while (Xstar(i).gt.1._f64)
        Xstar(i) = Xstar(i)-1._f64
      end do
      do while (Xstar(i).lt.0._f64)
        Xstar(i) = Xstar(i)+1._f64
      end do    
    enddo

   

    !from f compute the mean
!    do i=1,N
!      f(i)=f(i)*(node_positions(i+1)-node_positions(i))/dx
!    enddo
    
    
    !we compute the splines coefficients by solving the LU decomposition
!    M=0._f64
!    do i=1,N
!      M=M+f(i)
!    enddo
!    !M=M/real(N,rk)
!    do i=1,N
!      f(i)=f(i)-M*(node_positions(i+1)-node_positions(i))!/dx
!    enddo    
!    tmp=f(1)
!    f(1)=0._f64
!    do i=2,N
!      tmp2=f(i)
!      f(i)=f(i-1)+tmp
!      tmp=tmp2
!    enddo
    
    call setup_spline_nonunif_1D_periodic_aux( node_pos, N, buf, ibuf)
    call compute_spline_nonunif_1D_periodic_aux2( f, N, buf, ibuf, coeffs )
    call interpolate_array_value_nonunif_aux( Xstar, f, N, node_pos, coeffs,N)
    
    
!    tmp=f(1)
!    do i=1,N-1
!      f(i)=f(i+1)-f(i)+M*(node_positions(i+1)-node_positions(i))
!    enddo
!    f(N)=tmp-f(N)+M*(node_positions(1)+1._f64-node_positions(N))


    !from mean compute f
!    do i=1,N
!      f(i)=f(i)*dx/(node_positions(i+1)-node_positions(i))
!    enddo

    f(N+1) = f(1)
    
    !deallocate(buf)
    !deallocate(ibuf)
    !deallocate(node_pos,coeffs)
    !deallocate(Xstar)
    
  end subroutine constant_advection_spl_non_unif_per




end module sll_module_advection_1d_non_uniform_cubic_splines