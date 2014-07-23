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
#include "sll_working_precision.h"
#include "sll_assert.h"
#include "sll_memory.h"
use sll_module_advection_1d_base
use sll_module_advection_1d_periodic
use lagrange_interpolation
use sll_module_advection_2d_oblic


implicit none
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
  sll_int32 :: step
  sll_real64 :: err
  sll_real64 :: alpha
  sll_int32 :: i0
  sll_real64 :: dt_loc
  sll_real64, dimension(:,:), allocatable :: buf
  sll_real64, dimension(:,:), allocatable :: f_new
  sll_int32 :: d
  sll_int32 :: r
  sll_int32 :: s
  sll_real64, dimension(:), allocatable :: xx
  sll_real64, dimension(:), allocatable :: x1_array
  sll_real64, dimension(:), allocatable :: x2_array
  sll_int32 :: ell
  sll_int32 :: i2_loc
  character(len=256) :: filename
  sll_int32 :: IO_stat
  sll_int32, parameter  :: input_file = 99
  sll_int32 :: i
  
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
    d, &
    A1_0, &
    A2_0, &
    A1, &
    A2
  
  !initialization
  k_mode = 3
  Nc_x1 = 512
  Nc_x2 = 16
  dt = 0.1_f64
  nb_step = 10  
  d = 5
  
  A1_0 = 3._f64
  A2_0 = 7._f64  ! we should assume A2>A1>=0
  
  A1 = 2.8357_f64
  A2 = 7.18459_f64

  call get_command_argument(1, filename)
  if (len_trim(filename) .ne. 0)then
    print*,'#read namelist'
    open(unit = input_file, file=trim(filename)//'.nml',IOStat=IO_stat)
      if( IO_stat /= 0 ) then
        print *, '#aligned_translation_2d failed to open file ', trim(filename)//'.nml'
        stop
      end if
      read(input_file, params) 
    close(input_file)    
  else
    print *,'#use default parameters'  
  endif
  print *,'#k_mode=',k_mode
  print *,'#Nc_x1=',Nc_x1
  print *,'#Nc_x2=',Nc_x2
  print *,'#x1_min x1_max=',x1_min,x1_max
  print *,'#x2_min x2_max=',x2_min,x2_max
  print *,'#nb_step=',nb_step
  print *,'#dt=',dt
  print *,'#d=',d
  print *,'#A1_0',A1_0
  print *,'#A2_0',A2_0
  print *,'#A1=',A1
  print *,'#A2=',A2
  
  delta_x1 = (x1_max-x1_min)/real(Nc_x1,f64)
  delta_x2 = (x2_max-x2_min)/real(Nc_x2,f64)  
  r = -(d-1)/2
  s = (d+1)/2
  SLL_ALLOCATE(xx(r:s),ierr)
  SLL_ALLOCATE(buf(r:s,Nc_x1+1),ierr)  
  SLL_ALLOCATE(f(Nc_x1+1,Nc_x2+1),ierr)
  SLL_ALLOCATE(f_init(Nc_x1+1,Nc_x2+1),ierr)
  SLL_ALLOCATE(f_exact(Nc_x1+1,Nc_x2+1),ierr)
  SLL_ALLOCATE(f_new(Nc_x1+1,Nc_x2+1),ierr)
  SLL_ALLOCATE(x1_array(Nc_x1+1),ierr)
  SLL_ALLOCATE(x2_array(Nc_x2+1),ierr)
  
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
  adv_x2 => new_periodic_1d_advector( &
    Nc_x2, &
    x2_min, &
    x2_max, &
!    LAGRANGE, & 
    SPLINE, & 
    4) 
  
  do i2=1,Nc_x2+1
    x2 = x2_min+real(i2-1,f64)*delta_x2
    do i1=1,Nc_x1+1
      x1 = x1_min+real(i1-1,f64)*delta_x1
      x2 = x2_min+real(i2-1,f64)*delta_x2
      f_init(i1,i2) = sin(2._f64*sll_pi*real(k_mode,f64)*(-A2_0*x1+A1_0*x2))
      x1 = x1 - A1*real(nb_step,f64)*dt
      x2 = x2 - A2*real(nb_step,f64)*dt
      f_exact(i1,i2) = sin(2._f64*sll_pi*real(k_mode,f64)*(-A2_0*x1+A1_0*x2))
    enddo
  enddo
  
  !classical method with splitting  
  f = f_init    
  err = 0._f64
     
  do step = 1,nb_step    
    !advection in x1
    do i2=1,Nc_x2+1
      call adv_x1%advect_1d_constant(A1, dt, f(1:Nc_x1+1,i2), f(1:Nc_x1+1,i2))
    enddo    
    !advection in x2
    do i1=1,Nc_x1+1
      call adv_x2%advect_1d_constant(A2, dt, f(i1,1:Nc_x2+1), f(i1,1:Nc_x2+1))
    enddo          
  enddo  
  err = maxval(abs(f-f_exact))  
  print *,'#err for classical method=',err


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

  
  
  !new method
  f = f_init  
  err = 0._f64  
  alpha = A2*dt/delta_x2
  i0 = floor(alpha)
  alpha = alpha-i0  
  print *,'#i0=',i0,alpha
    
  do step =1,nb_step
    do i2=1,Nc_x2+1
      !choose several dt_loc so that advection in x2 is exact
      do ell=r,s
        dt_loc = real(ell+i0,f64)*delta_x2/A2         
        i2_loc = modulo(i2-ell-i0-1,Nc_x2)+1
        call adv_x1%advect_1d_constant( &
          A1, &
          dt_loc, &
          f(1:Nc_x1+1,i2_loc), &
          buf(ell,1:Nc_x1+1))
      enddo
      ! interpolate between these values 
      do i1=1,Nc_x1+1
        f_new(i1,i2) = lagrange_interpolate(alpha, d, xx, buf(r:s,i1) )
      enddo
    enddo    
    f = f_new
  enddo
  err = maxval(abs(f-f_exact))
  print *,'#err with new method=',err

  !new method using oblic advector
  f = f_init  

  adv => new_oblic_2d_advector( &
    Nc_x1, &
    adv_x1, &
    Nc_x2, &
    x2_min, &
    x2_max, &
    r, &
    s )
  
  err = 0._f64  
  do step =1,nb_step
    call oblic_advect_2d_constant( &
      adv, &
      A1, &
      A2, &
      dt, &
      f, &
      f_new)
    f = f_new      
  enddo  
  err = maxval(abs(f-f_exact))
  print *,'#err with new method using oblic advector=',err
  
  

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
    use sll_xdmf
    use sll_hdf5_io_serial
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
