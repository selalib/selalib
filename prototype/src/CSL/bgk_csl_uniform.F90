program bgk_csl
#include "sll_working_precision.h"
#include "sll_mesh_types.h"
#include "sll_memory.h"
  use sll_constants
  use distribution_function
  use sll_diagnostics
  use sll_csl
  use sll_splines
  implicit none
  external compute_translate_nodes_periodic,compute_non_unif_integral
  !external poisson1dpertrap
  sll_real64 :: x2_min,x2_max,x1_min,x1_max,x1,x2,delta_x1,delta_x2
  sll_real64 :: mu,xi,L,H
  sll_int    :: i,j,N_phi,err,N_x1,N_x2,i1,i2,N,nb_step
  LOGICAL :: ex
  sll_real64,dimension(:), pointer :: phi,node_positions_x1,node_positions_x2
  sll_real64,dimension(:), pointer :: new_node_positions,buf_1d,rho,E,E_store
  sll_real64,dimension(:,:), pointer :: f,f_store
  sll_real64 :: phi_val,delta_x1_phi,xx,dt,alpha,val
  sll_int :: ii,step,test_case
  type(cubic_nonunif_spline_1D), pointer :: spl_per_x1, spl_per_x2
  
  

  N_x1 = 128!128
  N_x2 = 128!256
  dt = 0.1_f64
  nb_step = 600
  
  test_case = 4
  
  
  N = max(N_x1,N_x2)
  
  print *,'#max(N_x1,N_x2)=',N
  
  
  x2_max = 10._f64
  x2_min = -x2_max
  
  SLL_ALLOCATE(f(N_x1+1,N_x2+1),err)
  SLL_ALLOCATE(f_store(N_x1+1,N_x2+1),err)
  SLL_ALLOCATE(node_positions_x1(N_x1+1),err)
  SLL_ALLOCATE(node_positions_x2(N_x2+1),err)
  SLL_ALLOCATE(new_node_positions(N+1),err)
  SLL_ALLOCATE(buf_1d(N+1),err)
  SLL_ALLOCATE(rho(N_x1+1),err)
  SLL_ALLOCATE(E(N_x1+1),err)
  SLL_ALLOCATE(E_store(N_x1+1),err)
    
  
  spl_per_x1 =>  new_cubic_nonunif_spline_1D( N_x1, PERIODIC_SPLINE)
  spl_per_x2 =>  new_cubic_nonunif_spline_1D( N_x2, PERIODIC_SPLINE)
  
  
  !physical parameters
  mu=0.92_f64
  xi=0.90_f64
  L=14.71_f64
  
  inquire(file='half_phi.dat', exist=ex) 
  if(.not.(ex))then
    print *,'file half_phi.dat does not exist'
    stop
  endif  
  open(unit=900,file='half_phi.dat')  
    read(900,*) N_phi,L
    N_phi = 2*N_phi
    SLL_ALLOCATE(phi(N_phi+1),err)
    do j=1,N_phi/2+1
      read(900,*) i,x1,x2
      phi(i)=x1
    enddo
    do j=N_phi/2+2,N_phi+1
      phi(j)=phi(N_phi+2-j)
    enddo
  close(900)
  
  !L = 4._f64*sll_pi
  
  x1_min = 0._f64
  x1_max = L

  if(test_case==4.or.test_case==5)then
    L = 4._f64*sll_pi
    x1_min = 0._f64
    x1_max = L
    x2_min = -6._f64
    x2_max = -x2_min
  endif



  delta_x1_phi = (x1_max-x1_min)/real(N_phi,f64)
  
  

  open(unit=900,file='phi0.dat')  
    do i1=1,N_phi+1
      x1 = x1_min+real(i1-1,f64)*delta_x1_phi
      write(900,*) x1,phi(i1)
    enddo
  close(900)
  
  

  
  print *,'#N_phi=',N_phi
  
  delta_x1 = (x1_max-x1_min)/real(N_x1,f64)
  delta_x2 = (x2_max-x2_min)/real(N_x2,f64)

  do i1=1,N_x1+1
    x1 = x1_min+real(i1-1,f64)*delta_x1
    node_positions_x1(i1) = x1
  enddo
  do i2=1,N_x2+1
    x2 = x2_min+real(i2-1,f64)*delta_x2
    node_positions_x2(i2) = x2
  enddo
  
  
  
  do i1=1,N_x1+1
    do i2=1,N_x2+1
      x1 = x1_min+real(i1-1,f64)*delta_x1
      x2 = x2_min+real(i2-1,f64)*delta_x2
      phi_val = 0._f64
      xx = (x1-x1_min)/(x1_max-x1_min)
      if(xx<=0._f64)then
        xx = 0._f64
      endif
      if(xx>=1._f64)then
        xx = 1._f64-1e-15_f64
      endif
      xx = xx*real(N_phi,f64)
      ii = floor(xx)
      xx = xx-real(ii,f64)      
      phi_val = (1._f64-xx)*phi(ii+1)+xx*phi(ii+2)
      !phi_val = 0._f64
      H = 0.5_f64*x2*x2 + phi_val

      if(test_case==1)then
        val = mu/(sqrt(2._f64*sll_pi))*(2._f64-2._f64*xi)/(3._f64-2._f64*xi)*(1._f64+H/(1._f64-xi))*exp(-H)
      endif
      if(test_case==2)then
        val = mu/(sqrt(2._f64*sll_pi))*(2._f64-2._f64*xi)/(3._f64-2._f64*xi)*(1._f64+H/(1._f64-xi))*exp(-H)
        val = val*(1._f64+0.1_f64*cos(2._f64*sll_pi/L*x1))
      endif
      if(test_case==3)then
        val = exp(-0.5_f64*40._f64*((x1-.5_f64)**2+(x2-.5_f64)**2))
      endif
      if(test_case==4)then
        !linear landau damping
        val = 1._f64/(sqrt(2._f64*sll_pi))*exp(-0.5_f64*x2*x2)
        val = val*(1._f64+0.001_f64*cos(2._f64*sll_pi/L*x1))
      endif
      if(test_case==5)then
        !non linear landau damping
        val = 1._f64/(sqrt(2._f64*sll_pi))*exp(-0.5_f64*x2*x2)
        val = val*(1._f64+0.5_f64*cos(2._f64*sll_pi/L*x1))
      endif
      f(i1,i2) = val
    enddo
  enddo
  
  !print *,f(N_x1+1,:)
  
  !stop

  !splitting x1 dt/2
!  do i2=1,N_x2+1
!    buf_1d(1:N_x1+1) = f(1:N_x1+1,i2)
!    call compute_spline_nonunif( buf_1d, spl_per_x1, node_positions_x1)
!    alpha = node_positions_x2(i2)*0.5_f64*dt
!    new_node_positions(1:N_x1+1) = node_positions_x1(1:N_x1+1)
!    call compute_translate_nodes_periodic(-alpha,N_x1,node_positions_x1(1:N_x1+1),new_node_positions(1:N_x1+1))
!    call interpolate_array_value_nonunif( new_node_positions, buf_1d, N_x1+1, spl_per_x1)
!    f(1:N_x1+1,i2) = buf_1d(1:N_x1+1)
!  enddo
!  
!  !compute_rho
!  do i1=1,N_x1+1
!    buf_1d(1:N_x2+1) = f(i1,1:N_x2+1)
!    call compute_non_unif_integral(node_positions_x2(1:N_x2+1),buf_1d(1:N_x2+1),N_x2+1,val)
!    rho(i1)=val-1._f64
!  enddo
!  E=rho
!  call poisson1dpertrap(E,x1_max-x1_min,N_x1)
!
!  open(unit=900,file='field_1.dat')  
!    do i1=1,N_x1+1
!      x1 = x1_min+real(i1-1,f64)*delta_x1
!      write(900,*) x1,E(i1),rho(i1)
!    enddo
!  close(900)
  
  
  !time iteration
  
  do step=1,nb_step
    
    
    f_store = f
    !compute E
    do i1=1,N_x1+1
      buf_1d(1:N_x2+1) = f(i1,1:N_x2+1)
      call compute_non_unif_integral(node_positions_x2(1:N_x2+1),buf_1d(1:N_x2+1),N_x2+1,val)
      rho(i1)=val-1._f64
    enddo
    E=rho
    call poisson1dpertrap(E,x1_max-x1_min,N_x1)
    
    E_store = E
    

    !diagnostic
    val=0._f64
    do i1=1,N_x1
      val = val+E(i1)*E(i1)
    enddo
    val = val/real(N_x1,f64)
    print *,(real(step,f64)-1._f64)*dt,val



    !splitting x1 dt/2
    do i2=1,N_x2+1
      buf_1d(1:N_x1+1) = f(1:N_x1+1,i2)
      call compute_spline_nonunif( buf_1d, spl_per_x1, node_positions_x1)
      alpha = node_positions_x2(i2)*dt*0.5_f64
      new_node_positions(1:N_x1+1) = node_positions_x1(1:N_x1+1)
      call compute_translate_nodes_periodic(-alpha,N_x1,node_positions_x1(1:N_x1+1),new_node_positions(1:N_x1+1))
      call interpolate_array_value_nonunif( new_node_positions, buf_1d, N_x1+1, spl_per_x1)
      f(1:N_x1+1,i2) = buf_1d(1:N_x1+1)
    enddo
    !compute E
    !do i1=1,N_x1+1
    !  buf_1d(1:N_x2+1) = f(i1,1:N_x2+1)
    !  call compute_non_unif_integral(node_positions_x2(1:N_x2+1),buf_1d(1:N_x2+1),N_x2+1,val)
    !  rho(i1)=val-1._f64
    !enddo
    !E=rho
    !call poisson1dpertrap(E,x1_max-x1_min,N_x1)





    !diagnostic
    !val=0._f64
    !do i1=1,N_x1
    !  val = val+E(i1)*E(i1)
    !enddo
    !val = val/real(N_x1,f64)
    !print *,(real(step,f64)-0.5_f64)*dt,val







    
    !splitting x2 dt
    do i1=1,N_x1+1
      buf_1d(1:N_x2+1) = f(i1,1:N_x2+1)
      call compute_spline_nonunif( buf_1d, spl_per_x2, node_positions_x2)
      alpha = E(i1)*dt
      new_node_positions(1:N_x2+1) = node_positions_x2(1:N_x2+1)
      call compute_translate_nodes_periodic(-alpha,N_x2,node_positions_x2(1:N_x2+1),new_node_positions(1:N_x2+1))
      call interpolate_array_value_nonunif( new_node_positions, buf_1d, N_x2+1, spl_per_x2)
      f(i1,1:N_x2+1) = buf_1d(1:N_x2+1)
    enddo
    

    !splitting x1 dt/2
    do i2=1,N_x2+1
      buf_1d(1:N_x1+1) = f(1:N_x1+1,i2)
      call compute_spline_nonunif( buf_1d, spl_per_x1, node_positions_x1)
      alpha = node_positions_x2(i2)*dt*0.5_f64
      new_node_positions(1:N_x1+1) = node_positions_x1(1:N_x1+1)
      call compute_translate_nodes_periodic(-alpha,N_x1,node_positions_x1(1:N_x1+1),new_node_positions(1:N_x1+1))
      call interpolate_array_value_nonunif( new_node_positions, buf_1d, N_x1+1, spl_per_x1)
      f(1:N_x1+1,i2) = buf_1d(1:N_x1+1)
    enddo

    
    !compute E
    do i1=1,N_x1+1
      buf_1d(1:N_x2+1) = f(i1,1:N_x2+1)
      call compute_non_unif_integral(node_positions_x2(1:N_x2+1),buf_1d(1:N_x2+1),N_x2+1,val)
      rho(i1)=val-1._f64
    enddo
    E=rho
    call poisson1dpertrap(E,x1_max-x1_min,N_x1)

    !print *,N_x1, modulo(1,N_x1),modulo(N_x1,N_x1),modulo(0,N_x1),modulo(-1,N_x1)
    
    !stop

    !E = 0.5_f64*(E+E_store)
    do i=1,N_x1
      E(i) = 0.5_f64*(E_store(i)+1._f64*dt*E_store(i)*(E_store(1+modulo(i,N_x1))-E_store(1+modulo(i+N_x1-2,N_x1))/(2._f64*delta_x1))+E(i))
    enddo
    E(N_x1+1)=E(1)

    f = f_store

    !splitting x1 dt/2
    do i2=1,N_x2+1
      buf_1d(1:N_x1+1) = f(1:N_x1+1,i2)
      call compute_spline_nonunif( buf_1d, spl_per_x1, node_positions_x1)
      alpha = node_positions_x2(i2)*dt*0.5_f64
      new_node_positions(1:N_x1+1) = node_positions_x1(1:N_x1+1)
      call compute_translate_nodes_periodic(-alpha,N_x1,node_positions_x1(1:N_x1+1),new_node_positions(1:N_x1+1))
      call interpolate_array_value_nonunif( new_node_positions, buf_1d, N_x1+1, spl_per_x1)
      f(1:N_x1+1,i2) = buf_1d(1:N_x1+1)
    enddo
    !compute E
    !do i1=1,N_x1+1
    !  buf_1d(1:N_x2+1) = f(i1,1:N_x2+1)
    !  call compute_non_unif_integral(node_positions_x2(1:N_x2+1),buf_1d(1:N_x2+1),N_x2+1,val)
    !  rho(i1)=val-1._f64
    !enddo
    !E=rho
    !call poisson1dpertrap(E,x1_max-x1_min,N_x1)





    !diagnostic
    !val=0._f64
    !do i1=1,N_x1
    !  val = val+E(i1)*E(i1)
    !enddo
    !val = val/real(N_x1,f64)
    !print *,(real(step,f64)-0.5_f64)*dt,val







    
    !splitting x2 dt
    do i1=1,N_x1+1
      buf_1d(1:N_x2+1) = f(i1,1:N_x2+1)
      call compute_spline_nonunif( buf_1d, spl_per_x2, node_positions_x2)
      alpha = E(i1)*dt
      new_node_positions(1:N_x2+1) = node_positions_x2(1:N_x2+1)
      call compute_translate_nodes_periodic(-alpha,N_x2,node_positions_x2(1:N_x2+1),new_node_positions(1:N_x2+1))
      call interpolate_array_value_nonunif( new_node_positions, buf_1d, N_x2+1, spl_per_x2)
      f(i1,1:N_x2+1) = buf_1d(1:N_x2+1)
    enddo
    

    !splitting x1 dt/2
    do i2=1,N_x2+1
      buf_1d(1:N_x1+1) = f(1:N_x1+1,i2)
      call compute_spline_nonunif( buf_1d, spl_per_x1, node_positions_x1)
      alpha = node_positions_x2(i2)*dt*0.5_f64
      new_node_positions(1:N_x1+1) = node_positions_x1(1:N_x1+1)
      call compute_translate_nodes_periodic(-alpha,N_x1,node_positions_x1(1:N_x1+1),new_node_positions(1:N_x1+1))
      call interpolate_array_value_nonunif( new_node_positions, buf_1d, N_x1+1, spl_per_x1)
      f(1:N_x1+1,i2) = buf_1d(1:N_x1+1)
    enddo

    
    
    
    
    
    
  enddo

  open(unit=900,file='field_final.dat')  
    do i1=1,N_x1+1
      x1 = x1_min+real(i1-1,f64)*delta_x1
      write(900,*) x1,E(i1),rho(i1)
    enddo
  close(900)


  
end program


subroutine compute_translate_nodes_periodic(alpha,N_cells,old_node_positions,new_node_positions)
  ! compute displaced nodes in the case of a translation
  ! the nodes are put in [x_min,x_max] by periodicity
  use sll_constants
  implicit none

  sll_int,intent(in) :: N_cells
  sll_real64,dimension(1:N_cells+1) :: old_node_positions
  sll_real64,dimension(1:N_cells+1) :: new_node_positions
  sll_int :: i
  sll_real64 :: alpha,x_min,x_max,xx  
  x_min = old_node_positions(1)
  x_max = old_node_positions(N_cells+1)
  
  do i=1,N_cells+1
    xx = (old_node_positions(i)-alpha-x_min)/(x_max-x_min)
    do while(xx>=1._f64)
      xx = xx-1._f64
    enddo
    do while(xx<0._f64)
        xx = xx+1._f64
    enddo
    if(xx>1._f64)then
      print *,'Problem of localization',i,xx
      stop
    endif
    new_node_positions(i) = x_min+xx*(x_max-x_min)
  enddo
  
  
end subroutine compute_translate_nodes_periodic


subroutine compute_non_unif_integral(x_points,f_points,N_points,val)
  use sll_constants
  implicit none
  sll_real64,intent(out) :: val
  sll_int,intent(in) :: N_points
  sll_real64,dimension(1:N_points) :: x_points,f_points
  sll_int :: i
  sll_real64 :: tmp,x1,x2,fval1,fval2
  val = 0._f64
  if(N_points<=1)then
    print *,'bad value of N_points=',N_points
    stop
  endif
  do i=1,N_points-1
    x1 = x_points(i)
    x2 = x_points(i+1)
    if(x2<x1)then
      print *,i,'bad integration points x1=',x1,'x2=',x2
      stop
    endif
    fval1 = f_points(i)
    fval2 = f_points(i+1)
    tmp = 0.5_f64*(fval1+fval2)*(x2-x1)
    val=val+tmp
  enddo
  
  
end  subroutine compute_non_unif_integral


subroutine poisson1dpertrap(E,L,N)
  use sll_constants
  implicit none
  sll_int,intent(in)::N
  sll_real64,dimension(N+1),intent(inout)::E
  sll_real64,intent(in)::L
  sll_int::i
  sll_real64::eold,enew,dx2,tmp
  !ensures at first that Ein is of mean zero
  tmp=0._f64
  do i=1,N
    tmp=tmp+E(i)
  enddo
  tmp=-tmp/real(N,f64)
  do i=1,N
    E(i)=E(i)+tmp
  enddo
  dx2=0.5_f64*L/real(N,f64)
  eold=E(1)
  E(1)=0._f64
  tmp=0._f64
  do i=1,N-1
    enew=E(i+1)
    E(i+1)=E(i)+dx2*(eold+enew)
    tmp=tmp+E(i+1)
    eold=enew
  enddo
  tmp=-tmp/real(N,f64)
  do i=1,N
    E(i)=E(i)+tmp
  enddo
  E(N+1)=E(1)
end subroutine poisson1dpertrap

