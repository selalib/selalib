
!******************************************************************************
!
! Selalib      
! Module: unit_test_nonuniform.F90
!
!> @brief 
!> cubic_nonuniform_splines unit test
!   
!> @authors                    
!> Michel Mehrenberger (mehrenbe@math.unistra.fr)
!                                  
!******************************************************************************


program nonuniform_spline_tester
#include "sll_working_precision.h"
#include "sll_assert.h"
#include "sll_memory.h"

  use cubic_non_uniform_splines
  use sll_constants
  !use sort_module
  implicit none
  
  
  
  logical                                :: test_passed
  sll_int32 :: err
  sll_int32 :: N,i,N_new,j
  sll_real64,dimension(:), pointer :: node_positions,f_per,f_hrmt,f
  sll_real64,dimension(:), pointer :: new_node_positions,f_new
  sll_real64,dimension(:), pointer :: fine_node_positions,f_fine
  sll_real64,dimension(:,:,:), pointer :: f_deriv
  type(cubic_nonunif_spline_1D), pointer :: spl_per, spl_hrmt, spl 
  sll_real64 :: x,xmin,xmax,tmp,linf_err(4),linf(4)
  sll_real64 :: xmin_val,xmax_val,slope_left,slope_right,fmin_val,fmax_val,local_xval(4)
  sll_real64 :: p(4),pp(4),w(4),fp(4)
  sll_real64 :: node_uniformity_min,node_uniformity_max,unif_val_min,unif_val_max
  sll_int32 :: test,nb_test,ival,N_test1,N_test2,index_max_err(6),unif_case,bdr_case
  sll_real64 :: max_err(6)


  test_passed = .true.
  



  
  xmin_val = -10._f64
  xmax_val = 10._f64
  
  fmin_val = -5._f64
  fmax_val = 5._f64
  
  unif_val_max = 2._f64
  unif_val_min= 1._f64/unif_val_max !has to be always the case
  
  ! results heavily depend on the value unif_val_max
  ! if unif_val_max=1, we return to the uniform case
  ! 
  
    
  unif_case = 1 ! only this case for the moment
  
  bdr_case = 2 !1=periodic 2=hermite
  
  nb_test = 10000
  
  N = 64
  
  N_new = 15
  
  max_err=0.0_f64
  
  !definition of the mesh
  SLL_ALLOCATE(node_positions(N+1), err)
  SLL_ALLOCATE(f_per(N+1), err)
  SLL_ALLOCATE(f_hrmt(N+1), err)
  
  SLL_ALLOCATE(fine_node_positions(4*N), err)
  SLL_ALLOCATE(f_fine(4*N), err)
  
  SLL_ALLOCATE(f_deriv(3,2,N), err)
  
  SLL_ALLOCATE(new_node_positions(N_new), err)
  SLL_ALLOCATE(f_new(N_new), err)
  
  spl_per =>  new_cubic_nonunif_spline_1D( N, SLL_PERIODIC)

  spl_hrmt =>  new_cubic_nonunif_spline_1D( N, SLL_HERMITE)
  
  
  do bdr_case=1,2
  
  if(bdr_case==1)then 
    spl => spl_per
    f   => f_per 
  endif  
  if(bdr_case==2)then 
    spl => spl_hrmt
    f   => f_hrmt
  endif  
  
  do test=1,nb_test
    
    call random_number(xmin)
    xmin = xmin_val+xmin*(xmax_val-xmin_val)
    call random_number(xmax)
    xmax = xmin_val+xmax*(xmax_val-xmin_val)  
    if(xmax<xmin)then
      tmp = xmin
      xmin = xmax
      xmax = tmp
    endif
    !if(unif_case==2)then  
    !  call random_number(node_positions(2:N))
    !  node_positions(1) = 0._f64
    !  node_positions(N+1) = 1._f64
    
      !do i=2,N
      !  node_positions(i)=real(i-1,f64)/real(N,f64)
      !enddo
      !xmin=0._f64
      !xmax=1._f64
     ! node_positions = xmin + (xmax-xmin)*node_positions
    
      !call sort(node_positions,N+1)    ! not for the moment in the repository
    !endif
    if(unif_case==1)then
      node_positions(1) = 0._f64
      call random_number(tmp)
      node_positions(2) = tmp
      do i=3,N+1
        call random_number(tmp)
        if(mod(i,2)==1)then
          tmp = unif_val_min + tmp*(1._f64-unif_val_min)
        else
          tmp = 1./(unif_val_min + tmp*(1._f64-unif_val_min))!unif_val_max + tmp*(1._f64-unif_val_max)        
        endif
        !tmp = unif_val_min + tmp*(unif_val_max-unif_val_min)
        !print *,i,tmp
        node_positions(i) = node_positions(i-1) + tmp*(node_positions(i-1)-node_positions(i-2))
      enddo
      !print *,node_positions
      !stop
      node_positions = node_positions/node_positions(N+1)
      node_positions(1) = 0._f64
      node_positions(N+1) = 1._f64
      node_positions = xmin + (xmax-xmin)*node_positions
    endif
    
    node_uniformity_min = 1._f64
    node_uniformity_max = 1._f64
    
    do i=2,N
      tmp = (node_positions(i+1)-node_positions(i))/(node_positions(i)-node_positions(i-1))
      if(tmp<node_uniformity_min)then
        node_uniformity_min = tmp
      endif
      if(tmp>node_uniformity_max)then
        node_uniformity_max = tmp
      endif
    enddo
    

    
    call random_number(new_node_positions(1:N_new))
    new_node_positions = xmin + (xmax-xmin)*new_node_positions
    
    f = 0.0_f64
    
    if(test<=N) then
      f(test) = 1._f64
    else
      call random_number(f)
    endif
    
    if(test<=N) then
      f(test) = 1._f64
    endif
    N_test1=N
    if(test==N_test1+1)then 
      f =1._f64
    endif
    if(test==N_test1+2)then
      do i=1,N+1
        f(i)=sin(2._f64*sll_pi/(xmax-xmin)*node_positions(i))
      enddo
    endif
    
    N_test2=N_test1+2
    if(test>N_test2)then
      call random_number(f)
    endif

    call random_number(slope_left) 
    call random_number(slope_right) 

    !call compute_spline_nonunif( f, spl, node_positions)
    if (bdr_case==1)then
      call compute_spline_nonunif( f_per, spl_per, node_positions)
    endif 
    if(bdr_case==2)then
      call compute_spline_nonunif( f_hrmt, spl_hrmt, node_positions,slope_left,slope_right)
      !call compute_spline_nonunif( f_hrmt, spl_hrmt, sl=slope_left,sr=slope_right)
    endif 
    
    
    call random_number(tmp)
    ival = floor(tmp*real(N,f64))+1
    
    ival = 1
    !print *,ival+1,N
    do i=1,N
      
      !call random_number(local_xval)      
      !if(i==ival)then
        local_xval(1) = 0._f64
        local_xval(2) = 1._f64/3._f64
        local_xval(3) = 2._f64/3._f64
        local_xval(4) = 1._f64-1e-12!(2._f64*local_xval(4)-1.0_f64)*1e-12        
      !endif
      x = node_positions(i)
      tmp = node_positions(i+1)-x
      fine_node_positions(4*(i-1)+1) = x + tmp*local_xval(1)
      fine_node_positions(4*(i-1)+2) = x + tmp*local_xval(2)      
      fine_node_positions(4*(i-1)+3) = x + tmp*local_xval(3)      
      fine_node_positions(4*(i-1)+4) = x + tmp*local_xval(4)      
    enddo
    
    
    call interpolate_array_value_nonunif( fine_node_positions, f_fine,4*N, spl)
    
    
    !print *,'ival=',ival
    
    !print *,f_fine_per
    
    !stop
    
    do i=1,N
      x = node_positions(i)
      tmp = node_positions(i+1)-x
      p(1)=(fine_node_positions(4*(i-1)+1)-x)/tmp
      p(2)=(fine_node_positions(4*(i-1)+2)-x)/tmp
      p(3)=(fine_node_positions(4*(i-1)+3)-x)/tmp
      p(4)=(fine_node_positions(4*(i-1)+4)-x)/tmp
      
      pp(1:4)=1._f64-p(1:4)
      
      fp(1:4)=f_fine((4*(i-1)+1):(4*(i-1)+4))
      if(fp(1).ne.fp(1))then
        print *,'detection of NaN'
        stop
      endif
      if(fp(2).ne.fp(2))then
        print *,'detection of NaN'
        stop
      endif
      if(fp(3).ne.fp(3))then
        print *,'detection of NaN'
        stop
      endif
      if(fp(4).ne.fp(4))then
        print *,'detection of NaN'
        stop
      endif
      
      !compute of f(x_i)
      w(1) = (p(2)*p(3)*p(4))/((p(2)-p(1))*(p(3)-p(1))*(p(4)-p(1)))
      w(2) = (p(1)*p(3)*p(4))/((p(1)-p(2))*(p(3)-p(2))*(p(4)-p(2)))
      w(3) = (p(1)*p(2)*p(4))/((p(1)-p(3))*(p(2)-p(3))*(p(4)-p(3)))
      w(4) = (p(1)*p(2)*p(3))/((p(1)-p(4))*(p(2)-p(4))*(p(3)-p(4)))
      f_deriv(1,1,i) = w(1)*fp(1)+w(2)*fp(2)+w(3)*fp(3)+w(4)*fp(4)     
      
      !compute f(x_{i+1})
      w(1) = -(pp(2)*pp(3)*pp(4))/((p(2)-p(1))*(p(3)-p(1))*(p(4)-p(1)))
      w(2) = -(pp(1)*pp(3)*pp(4))/((p(1)-p(2))*(p(3)-p(2))*(p(4)-p(2)))
      w(3) = -(pp(1)*pp(2)*pp(4))/((p(1)-p(3))*(p(2)-p(3))*(p(4)-p(3)))
      w(4) = -(pp(1)*pp(2)*pp(3))/((p(1)-p(4))*(p(2)-p(4))*(p(3)-p(4)))
      f_deriv(1,2,i) = w(1)*fp(1)+w(2)*fp(2)+w(3)*fp(3)+w(4)*fp(4)     
      
      !compute f'(x_{i})
      w(1) = -(p(2)*p(3)+p(2)*p(4)+p(3)*p(4))/((p(2)-p(1))*(p(3)-p(1))*(p(4)-p(1)))
      w(2) = -(p(1)*p(3)+p(1)*p(4)+p(3)*p(4))/((p(1)-p(2))*(p(3)-p(2))*(p(4)-p(2)))
      w(3) = -(p(1)*p(2)+p(1)*p(4)+p(2)*p(4))/((p(1)-p(3))*(p(2)-p(3))*(p(4)-p(3)))
      w(4) = -(p(1)*p(2)+p(1)*p(3)+p(2)*p(3))/((p(1)-p(4))*(p(2)-p(4))*(p(3)-p(4)))

      f_deriv(2,1,i) = (w(1)*fp(1)+w(2)*fp(2)+w(3)*fp(3)+w(4)*fp(4))/tmp     

      !compute f'(x_{i+1})
      w(1) = -(pp(2)*pp(3)+pp(2)*pp(4)+pp(3)*pp(4))/((p(2)-p(1))*(p(3)-p(1))*(p(4)-p(1)))
      w(2) = -(pp(1)*pp(3)+pp(1)*pp(4)+pp(3)*pp(4))/((p(1)-p(2))*(p(3)-p(2))*(p(4)-p(2)))
      w(3) = -(pp(1)*pp(2)+pp(1)*pp(4)+pp(2)*pp(4))/((p(1)-p(3))*(p(2)-p(3))*(p(4)-p(3)))
      w(4) = -(pp(1)*pp(2)+pp(1)*pp(3)+pp(2)*pp(3))/((p(1)-p(4))*(p(2)-p(4))*(p(3)-p(4)))

      f_deriv(2,2,i) = (w(1)*fp(1)+w(2)*fp(2)+w(3)*fp(3)+w(4)*fp(4))/tmp   

      !compute f''(x_{i})
      w(1) = 2._f64*(p(2)+p(3)+p(4))/((p(2)-p(1))*(p(3)-p(1))*(p(4)-p(1)))
      w(2) = 2._f64*(p(1)+p(3)+p(4))/((p(1)-p(2))*(p(3)-p(2))*(p(4)-p(2)))
      w(3) = 2._f64*(p(1)+p(2)+p(4))/((p(1)-p(3))*(p(2)-p(3))*(p(4)-p(3)))
      w(4) = 2._f64*(p(1)+p(2)+p(3))/((p(1)-p(4))*(p(2)-p(4))*(p(3)-p(4)))
      f_deriv(3,1,i) = (w(1)*fp(1)+w(2)*fp(2)+w(3)*fp(3)+w(4)*fp(4))/(tmp*tmp)     

      !compute f''(x_{i+1})
      w(1) = -2._f64*(pp(2)+pp(3)+pp(4))/((p(2)-p(1))*(p(3)-p(1))*(p(4)-p(1)))
      w(2) = -2._f64*(pp(1)+pp(3)+pp(4))/((p(1)-p(2))*(p(3)-p(2))*(p(4)-p(2)))
      w(3) = -2._f64*(pp(1)+pp(2)+pp(4))/((p(1)-p(3))*(p(2)-p(3))*(p(4)-p(3)))
      w(4) = -2._f64*(pp(1)+pp(2)+pp(3))/((p(1)-p(4))*(p(2)-p(4))*(p(3)-p(4)))
      f_deriv(3,2,i) = (w(1)*fp(1)+w(2)*fp(2)+w(3)*fp(3)+w(4)*fp(4))/(tmp*tmp)     
      
      !print *,i,tmp
      
    enddo
    
    !call compute_spline_nonunif( f_hrmt, spl_hrmt, sl=f_deriv(2,1,1), sr=f_deriv(2,2,N))
    !call compute_spline_nonunif( f_per, spl_per)
    
    !print *,f_deriv(2,1,1),f_deriv(2,2,N)
    
    !do i=-1,N+1
    !  print *,i,spl_per%coeffs(i),spl_hrmt%coeffs(i)
    !  print *,i,spl_per%node_positions(i),spl_hrmt%node_positions(i)
    !enddo
    
    
    !check for interpolation of f
    linf_err(1)=0._f64
    linf(1)=0._f64
    do i =1,N
      tmp = abs(f_deriv(1,1,i)-f(i))        
      if(tmp>linf_err(1))then
        linf_err(1) = tmp
      endif
      tmp = abs(f(i))        
      if(tmp>linf(1))then
        linf(1) = tmp
      endif      
    enddo
    if(bdr_case==2) then
      tmp = abs(f_deriv(1,2,N)-f(N+1))
      if(tmp>linf_err(1))then
        linf_err(1) = tmp
      endif
    endif
    
    
        
     !check for continuity of f
    linf_err(2)=0._f64
    linf(2)=0._f64
    do i =1,N
      if(i==1)then
        if(bdr_case==1)then
          tmp = abs(f_deriv(1,1,i)-f_deriv(1,2,N))
        else
          tmp = 0._f64
        endif  
      else
        tmp = abs(f_deriv(1,1,i)-f_deriv(1,2,i-1))        
      endif
      if(tmp>linf_err(2))then
        linf_err(2) = tmp
      endif
      tmp = abs(f_deriv(1,1,i))        
      if(tmp>linf(2))then
        linf(2) = tmp
      endif      
      !if(i>1)then
      !  print *,i,f_deriv(1,1,i),f_deriv(1,2,i-1),f_per(i)
      !endif  
    enddo

    !check for continuity of f'
    linf_err(3)=0._f64
    linf(3)=0._f64
    do i =1,N
      if(i==1)then
        if(bdr_case==1)then
          tmp = abs(f_deriv(2,1,i)-f_deriv(2,2,N))
        else
          tmp = 0._f64
        endif
      else
        tmp = abs(f_deriv(2,1,i)-f_deriv(2,2,i-1))        
      endif
      if(tmp>linf_err(3))then
        linf_err(3) = tmp
      endif
      tmp=abs(f_deriv(2,1,i))
      if(tmp>linf(3))then
        linf(3) = tmp
      endif      
      !if(i>1)then
      !  print *,i,f_deriv(2,1,i),f_deriv(2,2,i-1)
      !endif  
      
    enddo
    if(bdr_case==2)then
      tmp = abs(f_deriv(2,1,1)-slope_left)
      if(tmp>linf_err(3))then
        linf_err(3) = tmp
      endif
      tmp = abs(f_deriv(2,2,N)-slope_right)
      if(tmp>linf_err(3))then
        linf_err(3) = tmp
      endif
    endif
    

    !check for continuity of f''
    linf_err(4)=0._f64
    linf(4)=0._f64
    do i =1,N
      if(i==1)then
        if(bdr_case==1)then
          tmp = abs(f_deriv(3,1,i)-f_deriv(3,2,N))
        else
          tmp=0._f64
        endif
      else
        tmp = abs(f_deriv(3,1,i)-f_deriv(3,2,i-1))        
      endif
      if(tmp>linf_err(4))then
        linf_err(4) = tmp
      endif
      tmp=abs(f_deriv(3,1,i))
      if(tmp>linf(4))then
        linf(4) = tmp
      endif
      !if(i>1)then
      !  print *,i,f_deriv(3,1,i),f_deriv(3,2,i-1)
      !endif  
    enddo
    
    !print *,test,min(linf_err(1)/linf(1),linf(1)),min(linf_err(2)/linf(2),linf(2)),min(linf_err(3)/linf(3),linf(3)),min(linf_err(4)/linf(4),linf(4)),1._f64/node_uniformity_min,node_uniformity_max
    
    do j=1,4
      tmp=min(linf_err(j)/linf(j),linf(j))
      if(tmp.ge.max_err(j))then        
        max_err(j)=tmp
        index_max_err(j)=test
      endif
    enddo
      tmp=1._f64/node_uniformity_min
      if(tmp.ge.max_err(5))then
        max_err(5)=tmp
        index_max_err(5)=test
      endif
      tmp=node_uniformity_max
      if(tmp.ge.max_err(6))then
        max_err(6)=tmp
        index_max_err(6)=test
      endif
    
    call interpolate_array_value_nonunif( new_node_positions, f_new,N_new, spl)
    
  enddo
  
  print *,'#boundary_case',bdr_case
  print *,'#error: node,C0,C1,C2', max_err(1),max_err(2),max_err(3),max_err(4)
  print *,'#non uniformity',max_err(5),max_err(6)
  print *,'#index',index_max_err(1),index_max_err(2),index_max_err(3),index_max_err(4),index_max_err(5),index_max_err(6)
  
  if(max_err(1)>1.e-13)then
    print *,'#problem with node interpolation', max_err(1)
    print *,'#bdr_case=',bdr_case
    test_passed = .false.
  endif

  if(max_err(2)>1.e-12)then
    print *,'#problem with continuity', max_err(2)  
    print *,'#bdr_case=',bdr_case
    test_passed = .false.
  endif

  if(max_err(3)>1.e-10)then
    print *,'#problem with continuity of derivative', max_err(3)  
    print *,'#bdr_case=',bdr_case
    test_passed = .false.
  endif

  if(max_err(4)>1.e-10)then
    print *,'#problem with continuity of second derivative', max_err(4)  
    print *,'#bdr_case=',bdr_case
    test_passed = .false.
  endif
  
  if(abs(max_err(5)/unif_val_max-1._f64)>1.e-2)then
    print *,'#problem with random non uniformity',unif_val_max,max_err(5)
    print *,'#bdr_case=',bdr_case
    test_passed = .false.
  endif

  if(abs(max_err(6)/unif_val_max-1._f64)>1.e-2)then
    print *,'#problem with random non uniformity',unif_val_max,max_err(6)
    print *,'#bdr_case=',bdr_case
    test_passed = .false.
  endif
  
  
enddo  
  
   if (test_passed .eqv. .true.) then
     print *, '#'
     print *, '# cubic_non_uniform_splines unit test: PASSED'
     print *, '# '
  else
     print *, ' '
     print *, 'cubic_non_uniform_splines unit test: FAILED'
     print *, ' '
  endif
 
  
  
end program nonuniform_spline_tester
