program contrib_rho_tester

#include "sll_working_precision.h"
#include "sll_assert.h"
#include "sll_memory.h"


!function
!f(H)=mu/(sqrt(2._f64*sll_pi))*(2._f64-2._f64*xi)/(3._f64-2._f64*xi)*(1._f64+H/(1._f64-xi))*exp(-H)
!mu=0.92_f64
!xi=0.90_f64
!L=14.71_f64
!Hmin 
!Hmax
!integration
!1._f64-mu*(3._f64-2._f64*xi+2*eta_1)/(3._f64-2._f64*xi)*exp(-eta_1)

  use cubic_nonuniform_splines
  use sll_constants
  use contrib_rho_module

  !use init1d_module
  
  implicit none
  
  sll_int ::N_int1,N_int2,i,j,N_cells_eta1,N_cells_eta2,N_size,N_store,s,jj!,Nr,Nth
  
  sll_int :: N_int1_quarter,N_int2_quarter
  
  sll_real64,dimension(:),allocatable::store_contrib_rho
  integer,dimension(:),allocatable::store_index_contrib_rho
  integer,dimension(:,:),pointer::size_contrib_rho
  sll_real64,dimension(:,:,:),pointer::integration_points
  
  
  sll_real64::dx1,dx2,x1_min,x1_max,x2_min,x2_max,x1,x2,eta1_min,eta1_max,eta2_min,eta2_max
  
  sll_real64,dimension(:),pointer::x2_min_tab,x2_max_tab
  
  sll_int :: test_case
  sll_real64 :: tmp,tmp2
  
  sll_real64,dimension(:), pointer :: node_positions_eta1,node_positions_eta2,f_eta1,f_eta2
  sll_real64,dimension(:), pointer :: node_positions_eta1_non_unif
  sll_real64,dimension(:,:), pointer :: f,coef
  type(cubic_nonunif_spline_1D), pointer :: spl_eta1, spl_eta2
  sll_int :: err,f_case,N_eta1_non_unif,ii,N_phi,N_integ
  sll_real64,dimension(:), pointer :: rho, rho_exact,rho_exact_grid
  sll_real64,dimension(:), pointer :: phi
  sll_real64,dimension(:,:), pointer :: integration_points_non_unif(:,:)
  sll_real64 :: eta_1,eta_2
  sll_real64 :: h_min,h_max,f_val,delta_f
  sll_real64,dimension(:), pointer :: v_positions(:)
  sll_real64 :: L,mu,xi
  
  sll_real64 :: M,phi_val,v_min,v_max,delta_v_before,delta_v_after,delta_h
  
  sll_real64 :: f_min_after,f_max_after,f_min_before,f_max_before
  
  sll_int :: N_int2_before,N_int2_after,v_positions_case
  

  test_case=-1
  v_positions_case = 2
  
  mu=0.92_f64
  xi=0.90_f64
  L=14.71_f64
  
  M = 1._f64
  
  f_case = 5
  N_int1 = 500
  
  N_int2_before = 30
  N_int2_after = 30
  
  
  N_int2 = 2*N_int2_before+2*N_int2_after
  
  v_max = 10._f64
  v_min = -v_max
  
  N_cells_eta1 = 24
  N_cells_eta2 = 32
  
  N_store = N_cells_eta2 * N_int1 *5
  
  
  
  x1_min=0._f64
  x1_max=2._f64
  
  x2_min=-1._f64
  x2_max=3._f64
  
  x1_min=0._f64
  x1_max=1._f64
  
  x2_min=0._f64
  x2_max=1._f64
  
  
  dx1=(x1_max-x1_min)/real(N_int1,f64)
  dx2=(x2_max-x2_min)/real(N_int2,f64)
  
  SLL_ALLOCATE(integration_points(2,N_int2+1,N_int1+1),err)
  SLL_ALLOCATE(size_contrib_rho(-1:N_cells_eta2+1,N_int1+1),err)
  SLL_ALLOCATE(store_index_contrib_rho(N_store),err)
  SLL_ALLOCATE(store_contrib_rho(N_store),err)
  
  
  SLL_ALLOCATE(x2_min_tab(1:N_int1+1),err)
  SLL_ALLOCATE(x2_max_tab(1:N_int1+1),err)
  
  x2_min_tab = x2_min
  x2_max_tab = x2_max

  spl_eta1 =>  new_cubic_nonunif_spline_1D( N_cells_eta1, HERMITE_SPLINE)
  spl_eta2 =>  new_cubic_nonunif_spline_1D( N_cells_eta2, PERIODIC_SPLINE)

  SLL_ALLOCATE(node_positions_eta1(N_cells_eta1+1), err)
  SLL_ALLOCATE(node_positions_eta2(N_cells_eta2+1), err)
  SLL_ALLOCATE(f_eta1(N_cells_eta1+1), err)
  SLL_ALLOCATE(f_eta2(N_cells_eta2+1), err)
  SLL_ALLOCATE(f(N_cells_eta1+1,N_cells_eta2+1), err)
  SLL_ALLOCATE(coef(-1:N_cells_eta1+1,-1:N_cells_eta2+1), err)
  SLL_ALLOCATE(rho(N_int1+1), err)
  SLL_ALLOCATE(rho_exact(N_int1+1), err)
  SLL_ALLOCATE(rho_exact_grid(2*N_int1+1), err)
  SLL_ALLOCATE(integration_points_non_unif(2,N_int2+1), err)

  do i=1,N_cells_eta1+1
    node_positions_eta1(i)=real(i-1,f64)/real(N_cells_eta1,f64)
  enddo

  do i=1,N_cells_eta2+1
    node_positions_eta2(i)=real(i-1,f64)/real(N_cells_eta2,f64)
  enddo

  
  ! compute initial function and its spline coefficients
  
  f_eta1 = 0._f64
  f_eta2 = 0._f64
  call compute_spline_nonunif( f_eta1, spl_eta1, node_positions_eta1)
  call compute_spline_nonunif( f_eta2, spl_eta2, node_positions_eta2)


  if(test_case==0)then
    
    open(unit=900,file='integration_points.dat')
    read(900,*) x1,x2
    x1_min = x1
    x1_max = x2
    do s=1,N_int1+1
      read(900,*) i,x1
      x2_max_tab(i) =x1
      x2_min_tab(i) =-x1
    enddo
    do s=1,(N_int1+1)*(N_int2+1)    
      read(900,* ) i,j,x1,x2
      integration_points(1,j,i) = x1
      integration_points(2,j,i) = x2
      !print *,i,j,x1,x2
    enddo
    close(900)
  endif

  if(test_case==-1)then
    
    
    open(unit=900,file='phi.dat')
      read(900,*) N_phi,L
      SLL_ALLOCATE(phi(N_phi+1),err)
      do j=1,N_phi+1
        read(900,*) i,x1,x2
        phi(i)=x1
      enddo
    close(900)
    
    ! compute rho on uniform mesh
    

    jj= N_phi/N_int1
    if(jj==0)then
      print *,'bad compatibility between N_phi=',N_phi,' and N_int1=',N_int1
      print *,'N_phi/N_int1 should be a non zero integer'
      stop
    endif
    if(jj*N_int1/=N_phi)then
      print *,'bad compatibility between N_phi=',N_phi,' and N_int1=',N_int1
      print *,'N_phi/N_int1 should be a non zero integer'
      stop
    endif
    
    
    dx1 = (L/2._f64)/real(N_int1,f64)
    dx2 = (v_max-v_min)/real(N_int2,f64)
    do i=1,N_int1+1
      phi_val=phi(1+(i-1)*jj)
      x1=real(i-1)*dx1
      rho(i)=0._f64
      do j=1,N_int2+1
        x2=v_min+real(j-1,f64)*dx2
        eta_1 = phi_val+0.5_f64*x2*x2
        tmp = mu/(sqrt(2._f64*sll_pi))*(2._f64-2._f64*xi)/(3._f64-2._f64*xi)*(1._f64+eta_1/(1._f64-xi))*exp(-M*eta_1)
        rho(i)=rho(i)+tmp
      enddo
      rho(i) = rho(i)*dx2
      eta_1=phi_val
      rho_exact(i)=mu/(3._f64-2._f64*xi)*exp(-M*eta_1)*(2._f64*(1-xi+eta_1)/sqrt(M)+1._f64/(M*sqrt(M)))
      rho_exact_grid(i)=rho_exact(i)
    enddo
    
    tmp = 0._f64
    open(unit=900,file='rho_exact.dat')    
      do i=1,N_int1+1
        x1=real(i-1)*dx1
        write(900,*) x1,rho(i),rho_exact(i),rho_exact(i)-rho(i)
        if(abs(rho_exact(i)-rho(i))>tmp)then
          tmp = abs(rho_exact(i)-rho(i))
        endif
      enddo
    close(900)
    
    print *,'# error of rho on uniform grid',N_int2,tmp
    
    !now, compute on non uniform mesh given by h values
    delta_v_before = sqrt(4._f64*phi(N_phi+1))/real(N_int2_before,f64)
    delta_v_after = (v_max-sqrt(4._f64*phi(N_phi+1)))/real(N_int2_after,f64)
    
    print *,'#delta_v=',delta_v_before,delta_v_after
    
    
    !plot of fequil.dat
    SLL_ALLOCATE(v_positions(N_int2+1),err)
    
    
    h_max = xi
    h_min = 0._f64
    delta_h = (h_max-h_min)/real(N_int2_before,f64)
    do j=1,N_int2_before+1
      eta_1 = h_min+real(j-1,f64)*delta_h
      v_positions(j) = mu/(sqrt(2._f64*sll_pi))*(2._f64-2._f64*xi)/(3._f64-2._f64*xi)*(1._f64+eta_1/(1._f64-xi))*exp(-eta_1)
    enddo
    
    x1 = v_positions(2)-v_positions(1)
    x2 = v_positions(2)-v_positions(1)
     do j=1,N_int2_before
      if(v_positions(j+1)-v_positions(j)<x1)then
        x1 = v_positions(j+1)-v_positions(j)
      endif
      if(v_positions(j+1)-v_positions(j)>x2)then
        x2 = v_positions(j+1)-v_positions(j)
      endif
    enddo
      
   print *,'# check for monotonicity of fequil before',x1,x2
    
    
    
    f_min_before = v_positions(1)
    f_max_before = v_positions(N_int2_before+1)
    open(unit=900,file='fequil_before.dat')    
      do j=1,N_int2_before+1
        eta_1 = h_min+real(j-1,f64)*delta_h
        write(900,*) eta_1,v_positions(j)
      enddo
    close(900)




    h_max = 0.5_f64*v_max*v_max
    h_min = xi!phi(N_phi+1)
    delta_h = (h_max-h_min)/real(N_int2_after,f64)
    do j=1,N_int2_after+1
      eta_1 = h_min+real(j-1,f64)*delta_h
      v_positions(N_int2_before+j) = mu/(sqrt(2._f64*sll_pi))*(2._f64-2._f64*xi)/(3._f64-2._f64*xi)*(1._f64+eta_1/(1._f64-xi))*exp(-eta_1)
    enddo

    x1 = v_positions(N_int2_before+2)-v_positions(N_int2_before+1)
    x2 = x1
     do j=N_int2_before+1,N_int2_before+N_int2_after
      if(v_positions(j+1)-v_positions(j)<x1)then
        x1 = v_positions(j+1)-v_positions(j)
      endif
      if(v_positions(j+1)-v_positions(j)>x2)then
        x2 = v_positions(j+1)-v_positions(j)
      endif
    enddo
   
   print *,'# check for monotonicity of fequil after',x1,x2

    f_max_after = v_positions(N_int2_before+1)
    f_min_after = v_positions(N_int2_before+N_int2_after+1)
   
   open(unit=900,file='fequil_after.dat')    
     do j=1,N_int2_after+1
        eta_1 = h_min+real(j-1,f64)*delta_h
        write(900,*) eta_1,v_positions(N_int2_before+j)
     enddo
   close(900)

  print *,'# fmin/fmax',f_min_before,f_max_before,f_max_after,f_min_after,f_max_before-f_max_after
  
  print *,'#phi',phi(N_phi+1),2._f64* f_max_before
  
  eta_1 = phi(N_phi+1)
  tmp= mu/(sqrt(2._f64*sll_pi))*(2._f64-2._f64*xi)/(3._f64-2._f64*xi)*(1._f64+eta_1/(1._f64-xi))*exp(-M*eta_1)
  
  !print *,'#f(phi(L/2))=',tmp
  
    open(unit=900,file='fequil.dat')    
      h_max = xi
      h_min = 0._f64
      delta_h = (h_max-h_min)/real(N_int2_before,f64)      
      do j=1,N_int2_before+1
        eta_1 = h_min+real(j-1,f64)*delta_h
        write(900,*) eta_1,v_positions(j)
      enddo
      h_max = 0.5_f64*v_max*v_max
      h_min = xi!phi(N_phi+1)
      delta_h = (h_max-h_min)/real(N_int2_after,f64)
      do j=2,N_int2_after+1
        eta_1 = h_min+real(j-1,f64)*delta_h
        write(900,*) eta_1,v_positions(N_int2_before+j)
      enddo
    close(900)

    !compute of level lines of f for v_positions
    
    !stop
    v_positions = 0._f64

    !N_int2_before = 2*N_int2_before
    !N_int2_after = 2*N_int2_after

    
    delta_f = (f_max_before-f_min_before)/real(N_int2_before) 
    eta_1 = 0._f64
    do i=1,N_int2_before+1
      f_val = f_min_before + real(i-1,f64)*delta_f
       j=0
       tmp=1._f64
       do while ((abs(tmp).gt.1.e-15_f64).and.(j<=100))
         if(eta_1<xi)then
           tmp  = (mu/sqrt(2._f64*sll_pi))*(2._f64-2._f64*xi)/(3._f64-2._f64*xi)*(1._f64+eta_1/(1._f64-xi))*exp(-eta_1)-f_val
           tmp2 = (mu/sqrt(2._f64*sll_pi))*2._f64* (eta_1-xi )*exp(-eta_1)/(-3._f64 + 2._f64* xi)
           !print*,c,f,fp, En(i) 
           !alpha=1.d0/(1.d0-fp)
           eta_1=eta_1-tmp/tmp2
         else         
           !eta_1=0
           print *,'problem of eta_1 before',eta_1
           stop
         endif
        j=j+1
      enddo
      if(j>=100)then
        print *,'problem of convergence for Newton',eta_1,tmp,j,tmp2
      endif
      v_positions(i) = eta_1
      !print *,i,eta_1
    enddo
    
    
    
    
    delta_f = (f_max_after-f_min_after)/real(N_int2_after) 
    !eta_1 = 0._f64
    eta_1 = 1.5*xi
    do i=2,N_int2_after+1
      f_val = f_max_after - real(i-1,f64)*delta_f
       j=0
       tmp=1._f64
       do while ((abs(tmp).gt.1.e-15_f64).and.(j<=100))
         if(eta_1>=xi-1.e-15)then
           tmp  = (mu/sqrt(2._f64*sll_pi))*(2._f64-2._f64*xi)/(3._f64-2._f64*xi)*(1._f64+eta_1/(1._f64-xi))*exp(-eta_1)-f_val
           tmp2 = (mu/sqrt(2._f64*sll_pi))*2._f64* (eta_1-xi )*exp(-eta_1)/(-3._f64 + 2._f64* xi)
           !print*,c,f,fp, En(i) 
           !alpha=1.d0/(1.d0-fp)
           eta_1=eta_1-tmp/tmp2
         else
           !eta_1=0
           print *,'problem of eta_1 after',eta_1
           stop
         endif
        j=j+1
      enddo
      if(j>=100)then
        print *,'problem of convergence for Newton',eta_1,tmp,j,tmp2
      endif
      v_positions(i+N_int2_before) = eta_1
      !print *,i,eta_1
    enddo



    open(unit=900,file='fequil_inverse.dat')    
    do i=1,N_int2_before+N_int2_after+1
      write(900,*) i,v_positions(i)
    enddo
    close(900)

    dx1 = 0.5_f64*L/real(N_int1,f64)
    jj= N_phi/N_int1
    open(unit=900,file='h_lines.dat')    
    do j=1,N_int2_before+N_int2_after+1
      !eta_1=0.5_f64*v_positions(j)**2
      eta_1=v_positions(j)
      do i=1,N_int1+1
        phi_val=phi(1+(i-1)*jj)
        x1=real(i-1)*dx1
        if(eta_1>=phi_val)then
          x2=sqrt(2._f64*(eta_1-phi_val))
          write(900,*) x1,x2
        endif  
      enddo
    enddo
    close(900)

     
    
    !N_int2_before = N_int2_before/2
    !N_int2_after = N_int2_after/2
    
    
    !stop
      tmp=sqrt(2._f64*v_positions(1))
      do j=2,N_int2_before+N_int2_after+1
        v_positions(N_int2/2+j)=sqrt(2._f64*v_positions(j))      
      enddo
      v_positions(N_int2/2+1)=tmp
      do j=1,N_int2/2
        v_positions(N_int2/2+1-j) = -v_positions(N_int2/2+1+j)
      enddo
    
    !do j=1,N_int2+1
    !  print *,j,v_positions(j)
    !enddo
    
    !print *,N_int2_before+N_int2_after+1+N_int2/2
    !stop
    if(v_positions_case==1)then
      do j=1,N_int2_before+1
        v_positions(N_int2/2+j)=(real(j,f64)-0.5_f64)*delta_v_before      
      enddo
      do j=1,N_int2_after
        v_positions(N_int2/2+N_int2_before+1+j)=v_positions(N_int2/2+N_int2_before+1)+real(j,f64)*delta_v_after
      enddo
      do j=1,N_int2/2
        v_positions(N_int2/2+1-j) = -v_positions(N_int2/2+1+j)
      enddo
    endif
    
    
    !do i=1,N_int2+1
    !  print *,i,v_positions(i)
    !enddo
    !print *,sqrt(2._f64*phi(N_phi+1))
    !stop
    do i=1,N_int1+1
      phi_val=phi(1+(i-1)*jj)
      x1=real(i-1)*dx1
      N_integ=0
      do j=1,N_int2+1
        x2=v_positions(j)
        !x2=v_min+real(j-1,f64)*dx2
        eta_1 = 0.5_f64*x2*x2
        if(eta_1>=phi_val)then
          N_integ = N_integ+1
          integration_points(1,N_integ,i) = eta_1
          integration_points(2,N_integ,i) = x2
          if(x2>=0) then
            x2 = sqrt(2._f64*(eta_1-phi_val))
          else
            x2 = -sqrt(2._f64*(eta_1-phi_val))
          endif
          integration_points_non_unif(1,N_integ) = x2
          integration_points_non_unif(2,N_integ) = &
            &mu/(sqrt(2._f64*sll_pi))*(2._f64-2._f64*xi)/(3._f64-2._f64*xi)*(1._f64+eta_1/(1._f64-xi))*exp(-M*eta_1)
          
          !print *,i,j,integration_points_non_unif(1,N_integ),integration_points_non_unif(2,N_integ)
        else
          x2=-1._f64
        endif
        !print *,i,j,eta_1,x2,phi_val
        
      enddo
      !rho(i)= compute_non_unif_integral(integration_points_non_unif,N_integ)
      !rho(i)= compute_non_unif_integral_spline_old(integration_points_non_unif,N_integ,10000)
      !rho(i)= compute_non_unif_integral_spline(integration_points_non_unif,N_integ)
      rho(i)= compute_non_unif_integral_gaussian(integration_points_non_unif,N_integ)
      
    enddo
    tmp=0._f64
    open(unit=900,file='rho_non_unif.dat')    
      do i=1,N_int1+1
        x1=real(i-1)*dx1
        write(900,*) x1,rho(i),rho_exact(i),rho_exact(i)-rho(i)
        if(abs(rho_exact(i)-rho(i))>tmp)then
          tmp = abs(rho_exact(i)-rho(i))
        endif
      enddo
    close(900)

    print *,'# error of rho on non uniform grid',N_int2,tmp,abs(rho_exact(1)-rho(1))
    
    stop
    
    !now compute f on h grid
    !jj= N_phi/N_cells_eta1
    !if(jj==0)then
    !  print *,'bad compatibility between N_phi=',N_phi,' and N_cells_eta1=',N_cells_eta1
    !  print *,'N_phi/N_cells_eta1 should be a non zero integer'
    !  stop
    !endif
    !if(jj*N_cells_eta1/=N_phi)then
    !  print *,'bad compatibility between N_phi=',N_phi,' and N_cells_eta1=',N_cells_eta1
    !  print *,'N_phi/N_cells_eta1 should be a non zero integer'
    !  stop
    !endif
    
    
    open(unit=900,file='integration_points_store.dat')
      read(900,*) N_int1_quarter,N_int2_quarter,h_min,h_max,N_eta1_non_unif
    close(900)
    
    !tmp=dx2
    
    tmp=sqrt(2._f64*h_max)/real(N_cells_eta1,f64)
    do i=1,N_cells_eta1+1
      x2=real(i-1,f64)*tmp
      eta_1=0.5_f64*x2*x2
      do j=1,N_cells_eta2+1
        f(i,j)=mu/(sqrt(2._f64*sll_pi))*(2._f64-2._f64*xi)/(3._f64-2._f64*xi)*(1._f64+eta_1/(1._f64-xi))*exp(-M*eta_1)
      enddo
    enddo
    
   
    
    
    !print *,'N_phi=',N_phi,'L=',L
    
    !do i=1,N_phi+1
    !  print *,i,phi(i)
    !enddo
    
    open(unit=900,file='integration_points_store.dat')
    read(900,*) N_int1_quarter,N_int2_quarter,h_min,h_max,N_eta1_non_unif
    
    SLL_DEALLOCATE_ARRAY(integration_points,err)
    SLL_DEALLOCATE_ARRAY(x2_max_tab,err)
    SLL_DEALLOCATE_ARRAY(x2_min_tab,err)
    SLL_DEALLOCATE_ARRAY(size_contrib_rho,err)
    SLL_DEALLOCATE_ARRAY(rho, err)
    SLL_DEALLOCATE_ARRAY(rho_exact, err)
    
    N_int1 = 2*N_int1_quarter
    N_int2 = 2*N_int2_quarter
    
    SLL_ALLOCATE(integration_points(2,N_int2+1,N_int1+1),err)
    SLL_ALLOCATE(x2_max_tab(N_int1+1),err)
    SLL_ALLOCATE(x2_min_tab(N_int1+1),err)
    SLL_ALLOCATE(size_contrib_rho(-1:N_cells_eta2+1,N_int1+1),err)
    SLL_ALLOCATE(rho(N_int1+1), err)
    SLL_ALLOCATE(rho_exact(N_int1+1), err)
  
 
    !read(900,*) x1,x2
    !x1_min = x1
    !x1_max = x2
    x1_min = 0._f64
    x1_max= L
    do s=1,N_int1_quarter+1
      read(900,*) i,x1
      !print *,i,x1
      x2_max_tab(i) =x1
      x2_min_tab(i) =-x1
    enddo
    
    print *,'#x2_max:',x2_max_tab(1),x2_max_tab(N_int1_quarter+1)
    
    if(N_eta1_non_unif>=0)then
      SLL_ALLOCATE(node_positions_eta1_non_unif(N_eta1_non_unif+1),err)    
    endif
    
    do s=1,N_eta1_non_unif+1
      read(900,*) i, node_positions_eta1_non_unif(i)
    enddo
    
    !stop
    
    do s=1,(N_int1_quarter+1)*(N_int2_quarter+1)    
      read(900,* ) i,j,x1,x2
      !print *,i,j,x1,x2
      
      integration_points(1,j+N_int2_quarter,i) = x1
      integration_points(2,j+N_int2_quarter,i) = x2
      !print *,i,j,x1,x2
    enddo    
                    
    close(900)
    
    !print *,N_int1_quarter,N_int2_quarter
    
    
    
    do i = 1,N_int1/2
      integration_points(1,N_int2+1,i)=integration_points(1,N_int2,i)
      integration_points(2,N_int2+1,i)=integration_points(2,N_int2,i)
    enddo
    
    
    do j = 1,N_int2/2+1
      integration_points(1,j,N_int1/2+1)=integration_points(1,j,N_int1/2)
      integration_points(2,j,N_int1/2+1)=integration_points(2,j,N_int1/2)
    enddo
    
    integration_points(1,N_int2+1,N_int1/2+1)=integration_points(1,N_int2+1,N_int1/2)
    integration_points(2,N_int2+1,N_int1/2+1)=integration_points(2,N_int2+1,N_int1/2)
    
    
    !stop
    

    do s=1,N_int1_quarter
      x2_max_tab(N_int1-s+2) =x2_max_tab(s)
      x2_min_tab(N_int1-s+2) =x2_min_tab(s)
      !print *,s,N_int1-s+1,x2_max_tab(s),x2_max_tab(N_int1-s+1)
    enddo
    
    
    ! A  ! C !  vmax 
    ! B  ! D !   0
    !0   L/2 L  vmin
    
    ! A is filled 
    
    
    ! fill B
    do i=1,N_int1_quarter+1      
      do j=1,N_int2_quarter
        integration_points(1,j,i)=integration_points(1,N_int2+2-j,i)
        integration_points(2,j,i)=-integration_points(2,N_int2+2-j,i)+2._f64*sll_pi
      enddo
    enddo
    
    ! fill C

    do i=1,N_int1_quarter      
      do j=1,N_int2_quarter+1
        integration_points(1,j+N_int2_quarter,i+N_int1_quarter+1)=integration_points(1,j+N_int2_quarter,-i+N_int1_quarter+1)
        integration_points(2,j+N_int2_quarter,i+N_int1_quarter+1)=-integration_points(2,j+N_int2_quarter,-i+N_int1_quarter+1)+sll_pi
      enddo
    enddo

    ! fill D

    do i=1,N_int1_quarter      
      do j=1,N_int2_quarter
        integration_points(1,j,i+N_int1_quarter+1)=integration_points(1,N_int2+2-j,i+N_int1_quarter+1)
        integration_points(2,j,i+N_int1_quarter+1)=-integration_points(2,N_int2+2-j,i+N_int1_quarter+1)+2._f64*sll_pi
      enddo
    enddo

    
    !eta1_min = integration_points(1,1,1)
    !eta1_max = eta1_min
    !eta2_min = integration_points(2,1,1)
    !eta2_max = eta2_min
    ! normalization of integration points
    

    integration_points(1,:,:) = (integration_points(1,:,:)-h_min)/(h_max-h_min)
    integration_points(2,:,:) = (integration_points(2,:,:))/(2.0_f64*sll_pi)
    
    eta1_min = h_min
    eta1_max = h_max
    eta2_min = 0._f64
    eta2_max = 2.0_f64*sll_pi
    
    
    do i=1,N_int1+1
      do j=1,N_int2+1
        x1=integration_points(1,j,i)
        x2=integration_points(2,j,i)
        !if(x1<eta1_min)then
        !  eta1_min = x1
        !endif
        !if(x2<eta2_min)then
        !  eta2_min = x2
        !endif
        !if(x1>eta1_max)then
        !  eta1_max = x1
        !endif
        !if(x2>eta2_max)then
        !  eta2_max = x2
        !endif
        !print *,i,j,x1,x2
      enddo
    enddo
    
    
    
    
    print *,'#',eta1_min,eta1_max,eta2_min,eta2_max
    
    x2_max_tab(N_int1_quarter+1)=x2_max_tab(N_int1_quarter)
    x2_min_tab(N_int1_quarter+1)=x2_min_tab(N_int1_quarter)
    !do i=1,N_int1+1
    !  print *,i,x2_max_tab(i),x2_min_tab(i)
    !enddo
    
    !print *,h_min,h_max
    
  endif

  
  if(test_case==1)then
    do i=1,N_int1+1
      do j=1,N_int2+1
        integration_points(1,j,i) = x1_min+real(i-1,f64)*dx1
        integration_points(2,j,i) = x2_min+real(j-1,f64)*dx2
      enddo
    enddo
    eta1_min = x1_min
    eta1_max = x1_max
    eta2_min = x2_min
    eta2_max = x2_max
    
    integration_points(1,:,:) = (integration_points(1,:,:)-eta1_min)/(eta1_max-eta1_min)
    integration_points(2,:,:) = (integration_points(2,:,:)-eta2_min)/(eta2_max-eta2_min)

    open(unit=900,file='integration_points.dat')
    write(900,*) x1_min,x1_max
    do i=1,N_int1+1
      write(900,*) i,x2_max_tab(i)
    enddo
    
    do i=1,N_int1+1
      do j=1,N_int2+1
        write(900,* ) i,j,integration_points(1,j,i),integration_points(2,j,i)
      enddo  
    enddo
    close(900)

    
  endif
  
  if(test_case==2)then
    call random_number(integration_points(1:2,1:N_int2+1,1:N_int1+1))
  endif
  
  if(test_case==3)then
    x1_min=-5._f64
    x1_max=5._f64
  
    x2_min=-3._f64
    x2_max=3._f64

    x2_min_tab = x2_min
    x2_max_tab = x2_max

  
    dx1=(x1_max-x1_min)/real(N_int1,f64)
    dx2=(x2_max-x2_min)/real(N_int2,f64)
    !polar coordinates
    do i=1,N_int1+1
      do j=1,N_int2+1
        x1 = x1_min+real(i-1,f64)*dx1
        x2 = x2_min+real(j-1,f64)*dx2
        integration_points(1,j,i) = sqrt(x1*x1+x2*x2)
        if(abs(x1)>1.e-12)then
          integration_points(2,j,i) = datan(x2/x1)
        else
          integration_points(2,j,i) = 0._f64
        endif  
      enddo
    enddo
    
    eta1_min = integration_points(1,1,1) 
    eta1_max = eta1_min
    eta2_min = integration_points(2,1,1) 
    eta2_max = eta2_min
    do i=1,N_int1+1
      do j=1,N_int2+1
        x1=integration_points(1,j,i)
        if(x1>eta1_max)then
          eta1_max=x1
        endif
        if(x1<eta1_min)then
          eta1_min=x1
        endif
        x2=integration_points(2,j,i)
        if(x2>eta2_max)then
          eta2_max=x2
        endif
        if(x2<eta2_min)then
          eta2_min=x2
        endif        
      enddo
    enddo
    
    !print *,'#',eta1_min,eta1_max,eta2_min,eta2_max  
    
    integration_points(1,:,:) = (integration_points(1,:,:)-eta1_min)/(eta1_max-eta1_min)
    integration_points(2,:,:) = (integration_points(2,:,:)-eta2_min)/(eta2_max-eta2_min)
    
  endif



  do j=1,N_cells_eta2+1
    do i=1,N_cells_eta1+1
      eta_1 = eta1_min+(eta1_max-eta1_min)*node_positions_eta1(i)
      if(N_eta1_non_unif>=0)then      
        eta_1 = real(N_eta1_non_unif,f64)*node_positions_eta1(i)
        ii=floor(eta_1)
        eta_1 = eta_1-real(ii,f64)
        eta_1 = (1._f64-eta_1)*node_positions_eta1_non_unif(ii+1)+eta_1*node_positions_eta1_non_unif(ii+2)
      endif  
      eta_2 = eta2_min+(eta2_max-eta2_min)*node_positions_eta2(j)
      if(f_case==1) then      
        f(i,j) = 1._f64        
      endif
      if(f_case==2) then
        f(i,j) = eta_1
      endif
      if(f_case==3) then
        f(i,j) = eta_1*cos(eta_2)
      endif
      if(f_case==4) then
        f(i,j) = exp(-10*eta_1**2)
      endif
      if(f_case==5) then
        f(i,j)=mu/(sqrt(2._f64*sll_pi))*(2._f64-2._f64*xi)/(3._f64-2._f64*xi)*(1._f64+eta_1/(1._f64-xi))*exp(-M*eta_1)       
      endif
    enddo
  enddo  
  
  
  do j=1,N_cells_eta2+1
     f_eta1 = f(1:N_cells_eta1+1,j)     
     call compute_spline_nonunif( f_eta1, spl_eta1)
     coef(-1:N_cells_eta1+1,j) = spl_eta1%coeffs(-1:N_cells_eta1+1)
  enddo

  print *,'#',eta1_min,eta1_max,h_min,h_max,exp(-h_max*M)
  
  !stop

  do i=-1,N_cells_eta1+1
     f_eta2 = coef(i,1:N_cells_eta2+1)
     call compute_spline_nonunif( f_eta2, spl_eta2)
     coef(i,-1:N_cells_eta2+1) = spl_eta2%coeffs(-1:N_cells_eta2+1)
  enddo

  !do j=1,N_int2+1
  !  print *,j,integration_points(1,j,1),integration_points(2,j,1)
  !enddo
  
  !stop
  !integration_points(1,N_int2+1,1) =1._f64-1e-13
  !integration_points(1,1,1) =1._f64-1e-13
  
  do i=1,N_int1+1
    rho_exact_grid(i)=0._f64
    do j=1,N_int2+1
      eta_1 = integration_points(1,j,i)
      eta_2 = integration_points(2,j,i)
      !if(N_eta1_non_unif<0)then
      !  eta_1 = eta1_min+(eta1_max-eta1_min)*eta_1
      !else
      !  eta_1 = real(N_eta1_non_unif,f64)*eta_1
      !  ii=floor(eta_1)
      !  eta_1 = eta_1-real(ii,f64)
      !  eta_1 = (1._f64-eta_1)*node_positions_eta1_non_unif(ii+1)+eta_1*node_positions_eta1_non_unif(ii+2)      
      !endif
      eta_1 = eta1_min+(eta1_max-eta1_min)*eta_1      
      eta_2 = eta2_min+(eta2_max-eta2_min)*eta_2
      if(f_case==1) then      
        tmp = 1._f64        
      endif
      if(f_case==2) then
        tmp = eta_1
      endif
      if(f_case==3) then
        tmp = eta_1*cos(eta_2)
      endif
      if(f_case==4) then
        tmp = exp(-10*eta_1**2)
      endif
      if(f_case==5) then
        tmp = mu/(sqrt(2._f64*sll_pi))*(2._f64-2._f64*xi)/(3._f64-2._f64*xi)*(1._f64+eta_1/(1._f64-xi))*exp(-eta_1)
      endif      
      !x1 = x1_min+real(i-1,f64)*dx1
      !x2 = x2_min+real(j-1,f64)*dx2
      !tmp = exp(-10._f64*eta_1**2)
      !tmp = exp(-10._f64*(x1**2+x2**2))
      rho_exact_grid(i) = rho_exact_grid(i)+tmp
      !print *,i,j,tmp,eta_1
    enddo
    !print *,''
    !rho_exact_grid(i) = rho_exact_grid(i)*(x2_max_tab(i)-x2_min_tab(i))/real(N_int2+1,f64)
    rho_exact_grid(i) = rho_exact_grid(i)*(x2_max_tab(i)-x2_min_tab(i))/real(N_int2,f64)
  enddo
  
  
  do i=1,N_int1+1
    eta_1 = integration_points(1,N_int2/2+1,i)
    !if(N_eta1_non_unif<0)then
    !  eta_1 = eta1_min+(eta1_max-eta1_min)*eta_1
    !else
    !  eta_1 = real(N_eta1_non_unif,f64)*eta_1
    !  ii=floor(eta_1)
    !  eta_1 = eta_1-real(ii,f64)
    !  eta_1 = (1._f64-eta_1)*node_positions_eta1_non_unif(ii+1)+eta_1*node_positions_eta1_non_unif(ii+2)      
    !endif
    eta_1 = eta1_min+(eta1_max-eta1_min)*eta_1
    if(f_case==5)then
      !rho_exact(i)=1._f64-mu*(3._f64-2._f64*xi+2*eta_1)/(3._f64-2._f64*xi)*exp(-eta_1)
      !rho_exact(i)=mu*(3._f64-2._f64*xi+2*eta_1)/(3._f64-2._f64*xi)*exp(-eta_1)      
      rho_exact(i)=mu/(3._f64-2._f64*xi)*exp(-M*eta_1)*(2._f64*(1-xi+eta_1)/sqrt(M)+1._f64/(M*sqrt(M)))      
    endif  
  enddo
  
  !print *,node_positions_eta1_non_unif
  
  !print *,eta1_min,eta1_max
  
  
  !we redefine the integration_points on [0,1], bu using the composition of the two transforms:
  ! uniform mesh-> non uniform mesh -> mapping
  if(test_case==-1)then
    integration_points(1,:,:) = eta1_min+(eta1_max-eta1_min)*integration_points(1,:,:)
    do i=1,N_int1+1
      do j=1,N_int2+1
        x1=integration_points(1,j,i)
        ii=1
        do while(node_positions_eta1_non_unif(ii)<=x1)
          ii = ii+1
        enddo
        ii=ii-1
        x1=(x1-node_positions_eta1_non_unif(ii))/(node_positions_eta1_non_unif(ii+1)-node_positions_eta1_non_unif(ii))
        integration_points(1,j,i) = (real(ii,f64)+x1)/real(N_eta1_non_unif)
      enddo
    enddo
  endif
  
  !stop
  
  !compute splines
  
  size_contrib_rho =0
  
  
  
  !N_size =  
  N_size =  compute_contrib_rho(integration_points,N_int1,N_int2,size_contrib_rho,&
     & store_index_contrib_rho,store_contrib_rho,N_store,N_cells_eta1,N_cells_eta2)  
  
  print *,'#',N_size,N_store

  !stop
  
  if(N_size>N_store)then
    deallocate(store_index_contrib_rho)
    deallocate(store_contrib_rho)
    N_store=N_size
    allocate(store_index_contrib_rho(N_store))
    allocate(store_contrib_rho(N_store))    
    N_size =  compute_contrib_rho(integration_points,N_int1,N_int2,size_contrib_rho,&
     & store_index_contrib_rho,store_contrib_rho,N_store,N_cells_eta1,N_cells_eta2)  
  endif
  
  if(N_size>N_store)then
    print *,'N_size=',N_size,'N_store',N_store
    stop
  endif
  
  print *,'#',N_size,N_store
  
  !stop

  s=0
  do i=1,N_int1+1
    tmp=0._f64
    do j=-1,N_cells_eta2+1
      do jj=1,size_contrib_rho(j,i)
        s=s+1
        tmp = tmp+coef(store_index_contrib_rho(s),j)*store_contrib_rho(s)
      enddo
    enddo
    !rho(i) = tmp*(x2_max_tab(i)-x2_min_tab(i))/real(N_int2+1,f64)
    rho(i) = tmp*(x2_max_tab(i)-x2_min_tab(i))/real(N_int2,f64)
  enddo  

  !change for rho(L/2)
  rho(N_int1/2+1)=0.5_f64*(rho(N_int1/2)+rho(N_int1/2+2))
  

  
  s=0
  do i=1,N_int1+1
    tmp=0._f64
    do j=-1,N_cells_eta2+1
      !print *,i,j,size_contrib_rho(j,i)
      do jj=1,size_contrib_rho(j,i)
        s=s+1
        !print *,jj,store_index_contrib_rho(s),store_contrib_rho(s)
        tmp = tmp+coef(store_index_contrib_rho(s),j)*store_contrib_rho(s)
      enddo
    enddo
    !rho(i) = tmp*(x2_max_tab(i)-x2_min_tab(i))/real(N_int2+1,f64)
    print *,x1_min+real(i-1,f64)*(x1_max-x1_min)/real(N_int1,f64),rho(i),rho_exact(i),rho_exact_grid(i),&
    &2._f64*x2_max_tab(i),&
    &tmp/real(N_int2+1,f64)!tmp*dx2,tmp*(x2_max-x2_min)/real(N_int2+1,f64)   
  enddo  

    tmp=0._f64
    open(unit=900,file='rho_mapped_mesh.dat')    
      do i=1,N_int1+1
        x1=real(i-1)*dx1
        write(900,*) x1,rho(i),rho_exact(i),rho_exact(i)-rho(i)
        if(abs(rho_exact(i)-rho(i))>tmp)then
          tmp = abs(rho_exact(i)-rho(i))
        endif
      enddo
    close(900)

    print *,'# error of rho on mapped mesh',N_int1,N_int2,N_cells_eta1,N_cells_eta2,tmp,abs(rho_exact(1)-rho(1))

  
  !print *,'rho0'
  
  !do j=1,N_int2/2+1
  !  x2 = real(j-1,f64)*(x2_max_tab(1)-x2_min_tab(1))/real(N_int2,f64)    
  !  print *,j,integration_points(1,j+N_int2/2,1),integration_points(2,j+N_int2/2,1),(x2**2/2._f64)/h_max
  !enddo


end program