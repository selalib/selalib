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

!we take an initial function that is constant for displacement
!  X1'(t) = A1_0
!  X2'(t) = A2_0
!we work on periodic [0,1]^2 => we assume that A1_0 and A2_0 are integers
!we then consider a displacement
!  X1'(t) = A1
!  X2'(t) = A2
!where A1 and A2 are not so far from A1_0 and A2_0
!we use a lot of points in x1 and few points in x2

program aligned_translation_2d
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#include "sll_memory.h"
#include "sll_working_precision.h"

  use sll_m_advection_1d_base, only: &
    sll_advection_1d_base

  use sll_m_advection_1d_periodic, only: &
    new_periodic_1d_advector

  use sll_m_advection_2d_oblic, only: &
    new_oblic_2d_advector, &
    oblic_2d_advector, &
    oblic_advect_2d_constant

  use sll_m_boundary_condition_descriptors, only: &
    sll_periodic

  use sll_m_constants, only: &
    sll_pi

  use sll_m_cubic_spline_interpolator_2d, only: &
    new_cubic_spline_interpolator_2d

  use sll_m_fcisl_toroidal, only: &
    compute_modulo_vect2d_inplace, &
    interpolate2d_toroidal

  use sll_m_hdf5_io_serial, only: &
    sll_hdf5_file_close, &
    sll_hdf5_file_create, &
    sll_hdf5_write_array

  use sll_m_interpolators_2d_base, only: &
    sll_c_interpolator_2d

  use sll_m_periodic_interp, only: &
    lagrange, &
    spline

  use sll_m_timer, only: &
    sll_set_time_mark, &
    sll_time_elapsed_since, &
    sll_time_mark

  use sll_m_utilities, only: &
    int2string

  use sll_m_xdmf, only: &
    sll_xdmf_close, &
    sll_xdmf_open, &
    sll_xdmf_write_array

  implicit none
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  type(oblic_2d_advector), pointer :: adv  
  class(sll_advection_1d_base), pointer :: adv_x1
  class(sll_advection_1d_base), pointer :: adv_x2
  sll_int32 :: i1
  sll_int32 :: i2
  sll_int32 :: Nc_x1
  sll_int32 :: Nc_x2
  sll_real64 :: A1
  sll_real64 :: A2
  sll_real64 :: A1_0
  sll_real64 :: A2_0
  sll_int32 :: k_mode
  sll_real64, dimension(:,:), allocatable :: f
  sll_real64, dimension(:,:), allocatable :: f_init
  sll_real64, dimension(:,:), allocatable :: f_exact
  !sll_real64, dimension(:) :: buf_x1
  !sll_real64, dimension(:) :: buf_x2
  sll_int32 :: ierr
  sll_real64 :: x1
  sll_real64 :: x2
  sll_real64 :: x1_min
  sll_real64 :: x1_max
  sll_real64 :: x2_min
  sll_real64 :: x2_max
  sll_real64 :: delta_x1
  sll_real64 :: delta_x2
  sll_real64 :: dt
  sll_int32 :: nb_step
  sll_int32 :: num_dt
  sll_int32 :: step
  sll_int32 :: istep
  sll_real64 :: err
!  sll_real64 :: alpha
!  sll_int32 :: i0
  sll_real64 :: dt_loc
  sll_real64, dimension(:,:), allocatable :: buf
  sll_real64, dimension(:,:), allocatable :: f_new
  sll_int32 :: d
  sll_int32 :: r
  sll_int32 :: s
  sll_real64, dimension(:), allocatable :: xx
  sll_real64, dimension(:), allocatable :: x1_array
  sll_real64, dimension(:), allocatable :: x2_array
!  sll_int32 :: ell
!  sll_int32 :: i2_loc
  character(len=256) :: filename
  sll_int32 :: IO_stat
  sll_int32, parameter  :: input_file = 99
  sll_int32 :: i
  sll_real64 :: err0
!  sll_real64 :: err1
  sll_real64 :: err2
  sll_real64 :: dt_max0
  sll_real64 :: dt_max2
  sll_int32 :: num_dt1
  sll_int32 :: verbose
  type(sll_time_mark) :: t0
  sll_real64 :: time0
  sll_real64 :: time2
  sll_real64 :: time3
  sll_real64 :: time4
  sll_real64, dimension(:,:), allocatable :: feet_x1
  sll_real64, dimension(:,:), allocatable :: feet_x2
  sll_real64 :: dt_max3
  sll_real64 :: dt_max4
  sll_real64 :: err3  
  sll_real64 :: err4  
  class(sll_c_interpolator_2d), pointer :: interp_classic
  sll_real64, dimension(:), allocatable :: params_aligned
  sll_int32 :: hermite_p
  sll_int32 :: lag_p
  sll_int32 :: lag_r
  sll_int32 :: lag_s
  sll_real64 :: iota  
  sll_real64 :: R0
  sll_real64 :: F0
  sll_real64 :: smallr
  sll_real64 :: psipr





   
  ! namelists for data input
  namelist /params/ &
    k_mode, &
    Nc_x1, &
    Nc_x2, &
    x1_min, &
    x1_max, &
    x2_min, &
    x2_max, &
    dt, &
    nb_step, &
    num_dt, &
    d, &
    A1_0, &
    A2_0, &
    A1, &
    A2, &
    verbose
  
  !initialization
  k_mode = 3
  Nc_x1 = 512
  Nc_x2 = 16
  dt = 0.1_f64
  nb_step = 10
  num_dt = 1  
  d = 5
  x1_min = 0._f64
  x1_max = 1._f64
  x2_min = 0._f64
  x2_max = 1._f64
  verbose = 0
  
  
  A1_0 = 3._f64
  A2_0 = 7._f64  ! we should assume A2>A1>=0
  
  A1 = 2.8357_f64
  A2 = 7.18459_f64

  call get_command_argument(1, filename)
  if (len_trim(filename) .ne. 0)then
    if(verbose==1)then
      print*,'#read namelist'
    endif
    open(unit = input_file, file=trim(filename)//'.nml',IOStat=IO_stat)
      if( IO_stat /= 0 ) then
        print *, '#aligned_translation_2d failed to open file ', trim(filename)//'.nml'
        stop
      end if
      read(input_file, params) 
    close(input_file)    
  else
    if(verbose==1)then
      print *,'#use default parameters'
    endif    
  endif
  if(verbose==1)then
    print *,'#k_mode=',k_mode
    print *,'#Nc_x1=',Nc_x1
    print *,'#Nc_x2=',Nc_x2
    print *,'#x1_min x1_max=',x1_min,x1_max
    print *,'#x2_min x2_max=',x2_min,x2_max
    print *,'#nb_step=',nb_step
    print *,'#dt=',dt
    print *,'#num_dt=',num_dt
    print *,'#d=',d
    print *,'#A1_0',A1_0
    print *,'#A2_0',A2_0
    print *,'#A1=',A1
    print *,'#A2=',A2
  endif
  delta_x1 = (x1_max-x1_min)/real(Nc_x1,f64)
  delta_x2 = (x2_max-x2_min)/real(Nc_x2,f64)  
  r = -(d-1)/2
  s = (d+1)/2

  iota = A1/A2

  hermite_p = 6
  lag_p = d
  lag_r = -lag_p/2
  lag_s = (lag_p+1)/2
  
  R0 = 10._f64
  smallr = 2._f64
  psipr = 4._f64
  F0=-psipr*R0/(smallr*iota)
 








  SLL_ALLOCATE(xx(r:s),ierr)
  SLL_ALLOCATE(buf(r:s,Nc_x1+1),ierr)  
  SLL_ALLOCATE(f(Nc_x1+1,Nc_x2+1),ierr)
  SLL_ALLOCATE(f_init(Nc_x1+1,Nc_x2+1),ierr)
  SLL_ALLOCATE(f_exact(Nc_x1+1,Nc_x2+1),ierr)
  SLL_ALLOCATE(f_new(Nc_x1+1,Nc_x2+1),ierr)
  SLL_ALLOCATE(x1_array(Nc_x1+1),ierr)
  SLL_ALLOCATE(x2_array(Nc_x2+1),ierr)
  SLL_ALLOCATE(feet_x1(Nc_x1+1,Nc_x2+1),ierr)
  SLL_ALLOCATE(feet_x2(Nc_x1+1,Nc_x2+1),ierr)
  
  do i=1,Nc_x1+1
    x1_array(i) = x1_min+real(i-1,f64)*delta_x1
  enddo
  do i=1,Nc_x2+1
    x2_array(i) = x2_min+real(i-1,f64)*delta_x2
  enddo

  do i1=r,s
    xx(i1) = real(i1,f64)
  enddo
  
  
  adv_x1 => new_periodic_1d_advector( &
    Nc_x1, &
    x1_min, &
    x1_max, &
!    LAGRANGE, & 
    SPLINE, & 
    4) 
!  adv_x2 => new_periodic_1d_advector( &
!    Nc_x2, &
!    x2_min, &
!    x2_max, &
!!    LAGRANGE, & 
!    SPLINE, & 
!    4) 

  adv_x2 => new_periodic_1d_advector( &
    Nc_x2, &
    x2_min, &
    x2_max, &
    LAGRANGE, & 
    d+1)
!    SPLINE, & 
!    4) 


  interp_classic => new_cubic_spline_interpolator_2d( &
    Nc_x1+1, &
    Nc_x2+1, &
    x1_min, &
    x1_max, &
    x2_min, &
    x2_max, &
    SLL_PERIODIC, &
    SLL_PERIODIC)


  SLL_ALLOCATE(params_aligned(11),ierr)
  params_aligned(1) = R0
  params_aligned(2) = psipr
  params_aligned(3) = F0 
  params_aligned(4) = smallr
  params_aligned(5) = real(hermite_p,f64) 
  params_aligned(6) = real(lag_r,f64)
  params_aligned(7) = real(lag_s,f64)
  params_aligned(8) = 0._f64
  params_aligned(9) = 2._f64*sll_pi
  params_aligned(10) = 0._f64
  params_aligned(11) = 2._f64*sll_pi





  
  if(verbose==1)then
    print *,'#error fexact-finit=',maxval(f_exact-f_init)
  endif  
err0 = 0._f64
dt_max0 = dt
if(num_dt>1)then
  num_dt1=num_dt+1
else
  num_dt1 = num_dt  
endif 
do istep = 1,num_dt1  
  dt_loc = real(istep,f64)/real(num_dt,f64)*dt
  if(istep==num_dt+1)then
    dt_loc = dt_max2
  endif
  do i2=1,Nc_x2+1
    x2 = x2_min+real(i2-1,f64)*delta_x2
    do i1=1,Nc_x1+1
      x1 = x1_min+real(i1-1,f64)*delta_x1
      x2 = x2_min+real(i2-1,f64)*delta_x2
      f_init(i1,i2) = sin(2._f64*sll_pi*real(k_mode,f64) &
        *(-A2_0*(x1-x1_min)/(x1_max-x1_min)+A1_0*(x2-x2_min)/(x2_max-x2_min)))
      x1 = x1 - A1*real(nb_step,f64)*dt_loc
      x2 = x2 - A2*real(nb_step,f64)*dt_loc
      f_exact(i1,i2) = sin(2._f64*sll_pi*real(k_mode,f64) &
        *(-A2_0*(x1-x1_min)/(x1_max-x1_min)+A1_0*(x2-x2_min)/(x2_max-x2_min)))
    enddo
  enddo

  call sll_set_time_mark(t0)

  !classical method with splitting  
  f = f_init    
  err = 0._f64
     
  do step = 1,nb_step    
    !advection in x1
    do i2=1,Nc_x2+1
      call adv_x1%advect_1d_constant(A1, dt_loc, f(1:Nc_x1+1,i2), f(1:Nc_x1+1,i2))
    enddo    
    !advection in x2
    do i1=1,Nc_x1+1
      call adv_x2%advect_1d_constant(A2, dt_loc, f(i1,1:Nc_x2+1), f(i1,1:Nc_x2+1))
    enddo          
  enddo
  
  time0 = sll_time_elapsed_since(t0)
  print*,'#time for classical method', time0
    
  err = maxval(abs(f-f_exact))
  if(err>err0)then
    dt_max0 = dt_loc
    err0 = err 
  endif  
enddo
  if(verbose==1)then  
    print *,'#err for classical method=',err0
  endif

#ifndef NOHDF5
      call plot_f_cartesian( &
        0, &
        f, &
        x1_array, &
        Nc_x1+1, &
        x2_array, &
        Nc_x2+1, &
        'fold', 0._f64 )        
!      call plot_f_cartesian( &
!        iplot, &
!        f_visu_light, &
!        sim%x1_array_light, &
!        np_x1_light, &
!        node_positions_x2_light, &
!        sim%num_dof_x2_light, &
!        'light_f', time_init )        
#endif


#ifndef NOHDF5
      call plot_f_cartesian( &
        0, &
        f-f_exact, &
        x1_array, &
        Nc_x1+1, &
        x2_array, &
        Nc_x2+1, &
        'errorf_classic', 0._f64 )        
!      call plot_f_cartesian( &
!        iplot, &
!        f_visu_light, &
!        sim%x1_array_light, &
!        np_x1_light, &
!        node_positions_x2_light, &
!        sim%num_dof_x2_light, &
!        'light_f', time_init )        
#endif


  
  
!!new method
!  f = f_init  
!  err = 0._f64  
!  alpha = A2*dt/delta_x2
!  i0 = floor(alpha)
!  alpha = alpha-i0  
!  print *,'#i0=',i0,alpha
!    
!  do step =1,nb_step
!    do i2=1,Nc_x2+1
!      !choose several dt_loc so that advection in x2 is exact
!      do ell=r,s
!        dt_loc = real(ell+i0,f64)*delta_x2/A2         
!        i2_loc = modulo(i2-ell-i0-1,Nc_x2)+1
!        call adv_x1%advect_1d_constant( &
!          A1, &
!          dt_loc, &
!          f(1:Nc_x1+1,i2_loc), &
!          buf(ell,1:Nc_x1+1))
!      enddo
!      ! interpolate between these values 
!      do i1=1,Nc_x1+1
!        f_new(i1,i2) = lagrange_interpolate(alpha, d, xx, buf(r:s,i1) )
!      enddo
!    enddo    
!    f = f_new
!  enddo
!  err = maxval(abs(f-f_exact))
!  err1=err
!  print *,'#err with new method=',err

  !new method using oblic advector

  adv => new_oblic_2d_advector( &
    Nc_x1, &
    adv_x1, &
    Nc_x2, &
    x2_min, &
    x2_max, &
    r, &
    s )

  
err2 = 0._f64
dt_max2 = dt 
do istep = 1,num_dt1  
  f = f_init  
  dt_loc = real(istep,f64)/real(num_dt,f64)*dt
  if(istep==num_dt+1)then
    dt_loc = dt_max2
  endif

  do i2=1,Nc_x2+1
    x2 = x2_min+real(i2-1,f64)*delta_x2
    do i1=1,Nc_x1+1
      x1 = x1_min+real(i1-1,f64)*delta_x1
      x2 = x2_min+real(i2-1,f64)*delta_x2
      f_init(i1,i2) = sin(2._f64*sll_pi*real(k_mode,f64) &
        *(-A2_0*(x1-x1_min)/(x1_max-x1_min)+A1_0*(x2-x2_min)/(x2_max-x2_min)))
      x1 = x1 - A1*real(nb_step,f64)*dt_loc
      x2 = x2 - A2*real(nb_step,f64)*dt_loc
      f_exact(i1,i2) = sin(2._f64*sll_pi*real(k_mode,f64) &
        *(-A2_0*(x1-x1_min)/(x1_max-x1_min)+A1_0*(x2-x2_min)/(x2_max-x2_min)))

    enddo
  enddo

  call sll_set_time_mark(t0)


  do step =1,nb_step
    call oblic_advect_2d_constant( &
      adv, &
      A1, &
      A2, &
      dt_loc, &
      f, &
      f_new)
    f = f_new      
  enddo

  time2 = sll_time_elapsed_since(t0)
  print*,'#time for new method', time2
  
    
  err = maxval(abs(f-f_exact))
  if(err>err2)then
    dt_max2 = dt_loc
    err2 = err
  endif
enddo
  if(verbose==1)then  
    print *,'#err with new method using oblic advector=',err2
  endif
  



err3 = 0._f64
dt_max3 = dt 
do istep = 1,num_dt1  
  f = f_init  
  dt_loc = real(istep,f64)/real(num_dt,f64)*dt
  if(istep==num_dt+1)then
    dt_loc = dt_max3
  endif

  do i2=1,Nc_x2+1
    x2 = x2_min+real(i2-1,f64)*delta_x2
    do i1=1,Nc_x1+1
      x1 = x1_min+real(i1-1,f64)*delta_x1
      x2 = x2_min+real(i2-1,f64)*delta_x2
      f_init(i1,i2) = sin(2._f64*sll_pi*real(k_mode,f64) &
        *(-A2_0*(x1-x1_min)/(x1_max-x1_min)+A1_0*(x2-x2_min)/(x2_max-x2_min)))
      x1 = x1 - A1*real(nb_step,f64)*dt_loc
      x2 = x2 - A2*real(nb_step,f64)*dt_loc
      f_exact(i1,i2) = sin(2._f64*sll_pi*real(k_mode,f64) &
        *(-A2_0*(x1-x1_min)/(x1_max-x1_min)+A1_0*(x2-x2_min)/(x2_max-x2_min)))
      x1 = x1_min+real(i1-1,f64)*delta_x1
      x2 = x2_min+real(i2-1,f64)*delta_x2
      x1 = x1 - A1*dt_loc
      x2 = x2 - A2*dt_loc
      feet_x1(i1,i2) = x1
      feet_x2(i1,i2) = x2
    enddo
  enddo
  call compute_modulo_vect2d_inplace(feet_x1,Nc_x1+1,Nc_x2+1,x1_max-x1_min)
  call compute_modulo_vect2d_inplace(feet_x2,Nc_x1+1,Nc_x2+1,x2_max-x2_min)

  
  call sll_set_time_mark(t0)


  do step =1,nb_step
    call interp_classic%interpolate_array( &
      Nc_x1+1, &
      Nc_x2+1, &
      f, &
      feet_x1, &
      feet_x2, f_new)
    f = f_new      
  enddo

  time3 = sll_time_elapsed_since(t0)
  print*,'#time for classical method using charac', time3
  
    
  err = maxval(abs(f-f_exact))
  if(err>err3)then
    dt_max3 = dt_loc
    err3 = err
  endif
enddo


err4 = 0._f64
dt_max4 = dt 
do istep = 1,num_dt1  
  f = f_init  
  dt_loc = real(istep,f64)/real(num_dt,f64)*dt
  if(istep==num_dt+1)then
    dt_loc = dt_max3
  endif

  do i2=1,Nc_x2+1
    x2 = x2_min+real(i2-1,f64)*delta_x2
    do i1=1,Nc_x1+1
      x1 = x1_min+real(i1-1,f64)*delta_x1
      x2 = x2_min+real(i2-1,f64)*delta_x2
      f_init(i1,i2) = sin(2._f64*sll_pi*real(k_mode,f64) &
        *(-A2_0*(x1-x1_min)/(x1_max-x1_min)+A1_0*(x2-x2_min)/(x2_max-x2_min)))
      x1 = x1 - A1*real(nb_step,f64)*dt_loc
      x2 = x2 - A2*real(nb_step,f64)*dt_loc
      f_exact(i1,i2) = sin(2._f64*sll_pi*real(k_mode,f64) &
        *(-A2_0*(x1-x1_min)/(x1_max-x1_min)+A1_0*(x2-x2_min)/(x2_max-x2_min)))
      x1 = x1_min+real(i1-1,f64)*delta_x1
      x2 = x2_min+real(i2-1,f64)*delta_x2
      x1 = x1 - A1*dt_loc
      x2 = x2 - A2*dt_loc
      feet_x1(i1,i2) = x1
      feet_x2(i1,i2) = x2
    enddo
  enddo
  call compute_modulo_vect2d_inplace(feet_x1,Nc_x1+1,Nc_x2+1,x1_max-x1_min)
  call compute_modulo_vect2d_inplace(feet_x2,Nc_x1+1,Nc_x2+1,x2_max-x2_min)

  
  call sll_set_time_mark(t0)


  do step =1,nb_step
    f_new = interpolate2d_toroidal( &
      Nc_x1+1, &
      Nc_x2+1, &
      f, &
      feet_x1, &
      feet_x2, &
      params_aligned)
    f = f_new      
  enddo

  time4 = sll_time_elapsed_since(t0)
  print*,'#time for new method using charac', time4
  
    
  err = maxval(abs(f-f_exact))
  if(err>err4)then
    dt_max4 = dt_loc
    err4 = err
  endif
enddo






  
  
  
  print *,Nc_x1,Nc_x2,d,dt,nb_step,k_mode,A1,A2,A1_0,A2_0,err0,err2,dt_max0,dt_max2
  
  print *,"#err1=",err0,time0
  print *,"#err2=",err2,time2
  print *,"#err3=",err3,time3
  print *,"#err4=",err4,time4
  
#ifndef NOHDF5
      call plot_f_cartesian( &
        0, &
        f, &
        x1_array, &
        Nc_x1+1, &
        x2_array, &
        Nc_x2+1, &
        'fnew', 0._f64 )        
!      call plot_f_cartesian( &
!        iplot, &
!        f_visu_light, &
!        sim%x1_array_light, &
!        np_x1_light, &
!        node_positions_x2_light, &
!        sim%num_dof_x2_light, &
!        'light_f', time_init )        
#endif

#ifndef NOHDF5
      call plot_f_cartesian( &
        0, &
        f-f_exact, &
        x1_array, &
        Nc_x1+1, &
        x2_array, &
        Nc_x2+1, &
        'errorf_new', 0._f64 )        
!      call plot_f_cartesian( &
!        iplot, &
!        f_visu_light, &
!        sim%x1_array_light, &
!        np_x1_light, &
!        node_positions_x2_light, &
!        sim%num_dof_x2_light, &
!        'light_f', time_init )        
#endif

 
  
  
  


contains

#ifndef NOHDF5
!*********************
!*********************

  !---------------------------------------------------
  ! Save the mesh structure
  !---------------------------------------------------
  subroutine plot_f_cartesian( &
    iplot, &
    f, &
    node_positions_x1, &
    nnodes_x1, &
    node_positions_x2, &
    nnodes_x2, &
    array_name, time)    
    !mesh_2d)

    sll_int32 :: file_id
    sll_int32 :: error
    sll_real64, dimension(:), intent(in) :: node_positions_x1
    sll_real64, dimension(:), intent(in) :: node_positions_x2    
     character(len=*), intent(in) :: array_name !< field name
    sll_real64, dimension(:,:), allocatable :: x1
    sll_real64, dimension(:,:), allocatable :: x2
    sll_int32, intent(in) :: nnodes_x1
    sll_int32, intent(in) :: nnodes_x2
    sll_int32 :: i, j
    sll_int32, intent(in) :: iplot
    character(len=4)      :: cplot
    sll_real64, dimension(:,:), intent(in) :: f
    sll_real64 :: time
    
    if (iplot == 1) then

      SLL_ALLOCATE(x1(nnodes_x1,nnodes_x2), error)
      SLL_ALLOCATE(x2(nnodes_x1,nnodes_x2), error)
      do j = 1,nnodes_x2
        do i = 1,nnodes_x1
          x1(i,j) = node_positions_x1(i) !x1_min+real(i-1,f32)*dx1
          x2(i,j) = node_positions_x2(j) !x2_min+real(j-1,f32)*dx2
        end do
      end do
      call sll_hdf5_file_create("cartesian_mesh-x1.h5",file_id,error)
      call sll_hdf5_write_array(file_id,x1,"/x1",error)
      call sll_hdf5_file_close(file_id, error)
      call sll_hdf5_file_create("cartesian_mesh-x2.h5",file_id,error)
      call sll_hdf5_write_array(file_id,x2,"/x2",error)
      call sll_hdf5_file_close(file_id, error)
      deallocate(x1)
      deallocate(x2)

    end if

    call int2string(iplot,cplot)
    call sll_xdmf_open(trim(array_name)//cplot//".xmf","cartesian_mesh", &
      nnodes_x1,nnodes_x2,file_id,error)
    write(file_id,"(a,f8.3,a)") "<Time Value='",time,"'/>"
    call sll_xdmf_write_array(trim(array_name)//cplot,f,"values", &
      error,file_id,"Node")
    call sll_xdmf_close(file_id,error)
  end subroutine plot_f_cartesian

#endif

end program
