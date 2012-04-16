program radial_1d_SL
#include "sll_working_precision.h"
#include "sll_assert.h"
#include "sll_memory.h"

  use sll_splines
  use cubic_nonuniform_splines
  use numeric_constants
  implicit none

  type(sll_spline_2D), pointer :: spl_bsl,spl_bsl_nc,spl_fsl,spl_fsl_nc
  sll_int32  :: N,Neta1,Neta2,mesh_case,test_case,step,nb_step,visu_step,field_case
  sll_int32  :: i,j,bc1_type,bc2_type,err
  sll_real64 :: eta1,delta_eta1,eta1_min,eta1_max,eta2,delta_eta2,eta2_min,eta2_max
  sll_real64 :: x1,x2,x1_min,x2_min,x1_max,x2_max,dt
  sll_real64 :: alpha_mesh
  sll_real64 :: val,val_bsl,val_bsl_nc,val_fsl,val_fsl_nc,tmp1,tmp2
  sll_real64 :: a1,a2,eta1t,eta2t,k1eta1,k2eta1,k3eta1,k4eta1,k1eta2,k2eta2,k3eta2,k4eta2
  sll_real64,dimension(:,:), pointer :: f,fh_bsl,fh_bsl_nc,fh_fsl,fh_fsl_nc
  sll_real64,dimension(:,:), pointer :: x1_array,x2_array,x1c_array,x2c_array,jac_array,eta1feet,eta2feet
  sll_real64,dimension(:,:), pointer :: diag
	
  ! for the python script polar-ex1e.py
  !namelist /param/ N
  N = 64
  !read(*,NML=param);open(unit=900,file="gyrof.param");write(900,NML=param);close(900)
 
  ! ---- physical parameters ----
  Neta1 = N
  Neta2 = N
    
  ! ---- mesh ----
  mesh_case = 1
  ! domain : disc of radius eta1_max with a hole of radius eta1_min
  ! BC     : hermite-periodic
  if (mesh_case==1) then
    eta1_min = 0.2_f64
    eta1_max = 0.8_f64
    eta2_min = 0._f64
    eta2_max = 2._f64*sll_pi
    
    bc1_type = PERIODIC_SPLINE
    bc2_type = PERIODIC_SPLINE
  end if
  ! domain : Collela curvilinear mesh 
  ! BC     : periodic-periodic
  if(mesh_case==2)then
    eta1_min = 0._f64
    eta1_max = 1._f64
    eta2_min = 0._f64
    eta2_max = 2._f64
    
    bc1_type = PERIODIC_SPLINE
    bc2_type = PERIODIC_SPLINE
  end if
  
  delta_eta1 = (eta1_max-eta1_min)/real(Neta1,f64)
  delta_eta2 = (eta2_max-eta2_min)/real(Neta2,f64)  
  
  ! visualization parameters
  visu_step = 1
  
  ! ---- distribution function ----
  ! 1 : gaussian in r
  ! 2 : gaussian in r
  ! 3 : cos eta2
  ! 4 : truncated gaussian in eta2
  ! 5 : f=1
  ! 6 : gaussian in r and eta2
  ! 7 : dirac 
  test_case = 5
  
  ! ---- field ---- 
  ! 1 : divergence free complex1 field
  ! 2 : "rotation" (constant coef advection in eta2)
  ! 3 : translation of vector (a1,a2)
  field_case = 2
  a1=0.001_f64
  a2=0.002_f64
  
  ! ---- time-scheme parameters ----
  dt = 0._f64/10._f64
  nb_step = 1
  
	
  print *,'# N=',N
  print *,'# T=',real(nb_step,f64)*dt
  print *,'# test_case=',test_case
  print *,'# field_case=',field_case
  
	
  ! allocations
  SLL_ALLOCATE(f(Neta1+1,Neta2+1), err)
  SLL_ALLOCATE(fh_bsl(Neta1+1,Neta2+1), err)
  SLL_ALLOCATE(fh_bsl_nc(Neta1+1,Neta2+1), err)
  SLL_ALLOCATE(fh_fsl(Neta1+1,Neta2+1), err)
  SLL_ALLOCATE(fh_fsl_nc(Neta1+1,Neta2+1), err)
  SLL_ALLOCATE(x1_array(Neta1+1,Neta2+1), err)
  SLL_ALLOCATE(x2_array(Neta1+1,Neta2+1), err)
  SLL_ALLOCATE(x1c_array(Neta1+1,Neta2+1), err)
  SLL_ALLOCATE(x2c_array(Neta1+1,Neta2+1), err)
  SLL_ALLOCATE(eta1feet(Neta1+1,Neta2+1), err)
  SLL_ALLOCATE(eta2feet(Neta1+1,Neta2+1), err)
  SLL_ALLOCATE(jac_array(Neta1+1,Neta2+1), err)
  
  SLL_ALLOCATE(diag(10,0:nb_step), err)
	
  ! creation spline 
  spl_bsl => new_spline_2D(Neta1+1, Neta2+1, &
    eta1_min, eta1_max, &
    0._f64, 2._f64*sll_pi, &
    bc1_type, bc2_type)!,&
    !const_slope_x1_min = 0._f64,const_slope_x1_max = 0._f64)
  spl_bsl_nc => new_spline_2D(Neta1+1, Neta2+1, &
    eta1_min, eta1_max, &
    0._f64, 2._f64*sll_pi, &
    bc1_type, bc2_type)!,&
    !const_slope_x1_min = 0._f64,const_slope_x1_max = 0._f64)
  spl_fsl => new_spline_2D(Neta1+1, Neta2+1, &
    eta1_min, eta1_max, &
    0._f64, 2._f64*sll_pi, &
    bc1_type, bc2_type)!,&
    !const_slope_x1_min = 0._f64,const_slope_x1_max = 0._f64)
  spl_fsl_nc => new_spline_2D(Neta1+1, Neta2+1, &
    eta1_min, eta1_max, &
    0._f64, 2._f64*sll_pi, &
    bc1_type, bc2_type)!,&
    !const_slope_x1_min = 0._f64,const_slope_x1_max = 0._f64)
  
  ! Analytic distribution function
  ! and data for the curvilinear mesh
  do i=1,Neta1+1
    eta1 = eta1_min + real(i-1,f64)*delta_eta1
    do j=1,Neta2+1
      eta2 = eta2_min + real(j-1,f64)*delta_eta2
      
      ! polar mesh
      if (mesh_case==1) then
        x1_array(i,j) = eta1*cos(eta2)
        x2_array(i,j) = eta1*sin(eta2)
        jac_array(i,j) = eta1
      end if
      
      ! Collela mesh
      if (mesh_case==2)then
        alpha_mesh = 1.e-1_f64 !0.1_f64
        
        eta1t = eta1 + real(i-1,f64)*delta_eta1/2._f64
        eta2t = eta2 + real(j-1,f64)*delta_eta2/2._f64
       
        x1_min = eta1_min + alpha_mesh * sin(2*sll_pi*eta1_min) * sin(2*sll_pi*eta2_min)
        x1_max = eta1_max + alpha_mesh * sin(2*sll_pi*eta1_max) * sin(2*sll_pi*eta2_max)
        
        x2_min = eta2_min + alpha_mesh * sin(2*sll_pi*eta1_min) * sin(2*sll_pi*eta2_min)
        x2_max = eta2_max + alpha_mesh * sin(2*sll_pi*eta1_max) * sin(2*sll_pi*eta2_max)
       
        x1_array(i,j) = eta1 + alpha_mesh * sin(2*sll_pi*eta1) * sin(2*sll_pi*eta2)
        x2_array(i,j) = eta2 + alpha_mesh * sin(2*sll_pi*eta1) * sin(2*sll_pi*eta2)
             
        jac_array(i,j) = (1.0_f64 + alpha_mesh *2._f64 *sll_pi * cos (2*sll_pi*eta1t) * sin (2*sll_pi*eta2t)) * &
            (1.0_f64 + alpha_mesh *2._f64 * sll_pi * sin (2*sll_pi*eta1t) * cos (2*sll_pi*eta2t)) - &
            alpha_mesh *2._f64 *sll_pi * sin (2*sll_pi*eta1t) * cos (2*sll_pi*eta2t) * &
            alpha_mesh *2._f64 * sll_pi * cos (2*sll_pi*eta1t) * sin (2*sll_pi*eta2t)

        x1_array(i,j) = x1_min + x1_array(i,j)*(x1_max-x1_min)
        x2_array(i,j) = x2_min + x2_array(i,j)*(x2_max-x2_min)
             
        jac_array(i,j) = jac_array(i,j)*(x1_max-x1_min)*(x2_max-x2_min)
      end if

      ! test-function
      if (test_case==1) then
        f(i,j) = exp(-3000._f64*(eta1-1.5_f64*eta1_min)**2)
      endif
      if (test_case==2) then
        f(i,j) = exp(-100._f64*(eta1-0.5_f64*(eta1_max+eta1_min))**2)
      endif
      if (test_case==3) then
        f(i,j) = cos(eta2)
      endif
      if (test_case==4) then
        if(eta1>=0.4.and.eta1<=0.6)then
          f(i,j) = exp(-30._f64*(eta2-0.5*sll_pi)**2)
        else
          f(i,j) = 0._f64
        endif  
      endif
      if (test_case==5) then
        f(i,j) = 1._f64
      endif
      if (test_case==6) then
        f(i,j) = exp(-100._f64*(eta1-0.5_f64*(eta1_max+eta1_min))**2)*exp(-30._f64*(eta2-0.5*sll_pi)**2)
      endif
      if (test_case==7) then
        f(i,j) = 0._f64
        if ((i==17).and.(j==12)) then
          f(i,j) = real(32,f64)
        end if
      endif
    enddo  
  enddo

  ! intialization of the numerical distribution function
  fh_bsl = f
  fh_bsl_nc = f
  fh_fsl = f
  fh_fsl_nc = f
  
  ! diagnostic with the weighted (temp1) ad the "classical" (temp2) masses 
  tmp1=0._f64
  tmp2=0._f64
  do i=1,Neta1+1
    do j=1,Neta2+1
      tmp1 = tmp1 + f(i,j)*jac_array(i,j)
      tmp2 = tmp2 + f(i,j)
    enddo
  enddo
  tmp1 = tmp1*delta_eta1*delta_eta2
  tmp2 = tmp2*delta_eta1*delta_eta2
  
  ! initialization of the diagnostics
  ! with the weighted mass 
  diag(1,0) = tmp1  ! analytic solution
  diag(3,0) = tmp1  ! BSL
  diag(5,0) = tmp1  ! BSL NC
  diag(7,0) = tmp1  ! FSL
  diag(9,0) = tmp1  ! FSL NC
  ! with the classical mass
  diag(2,0) = tmp2  ! analytic solution
  diag(4,0) = tmp2  ! BSL
  diag(6,0) = tmp2  ! BSL NC
  diag(8,0) = tmp2  ! FSL
  diag(10,0) = tmp2 ! FSL NC
  
  ! checking that the function is close to 0 at the r boundary 
  print *,'#f0 boundary:',f(1,1),f(Neta1+1,1)
  
  ! saving f
  open(unit=900,file='f0.dat')  
  do i=1,Neta1+1
    eta1 = eta1_min + real(i-1,f64)*delta_eta1
    do j=1,Neta2+1
      eta2 = eta2_min + real(j-1,f64)*delta_eta2
      write(900,*) x1_array(i,j),x2_array(i,j),eta1,eta2,f(i,j)
    enddo
    write(900,*) ' '
  end do
  close(900)

  ! Evolution in time - Analytic solution
  open(unit=900,file='carac.dat')  
  open(unit=800,file='fn.dat')
  do i=1,Neta1+1
    eta1 = eta1_min + real(i-1,f64)*delta_eta1
    do j=1,Neta2+1
      eta2 = eta2_min + real(j-1,f64)*delta_eta2
      write(800,'(f0.3,a,f0.3)',advance='no') eta1,' ',eta2
      write(800,'(a,f0.3,a,f0.3)',advance='no') ' ',x1_array(i,j),' ',x2_array(i,j)
      do step=1,nb_step
          
        !RK4
        if(field_case==1)then
          eta1t = eta1
          eta2t = eta2
          k1eta1 = eta1t*(eta1t-eta1_min)*(eta1_max-eta1t)*cos(eta2t)
          k1eta2 = -(eta1t*(eta1_max+eta1_min-2._f64*eta1t)+2._f64*(eta1t-eta1_min)*(eta1_max-eta1t))*sin(eta2t)
          eta1t = eta1 + 0.5_f64*dt*k1eta1
          eta2t = eta2 +  0.5_f64*dt*k1eta2
          k2eta1 = eta1t*(eta1t-eta1_min)*(eta1_max-eta1t)*cos(eta2t)
          k2eta2 = -(eta1t*(eta1_max+eta1_min-2._f64*eta1t)+2._f64*(eta1t-eta1_min)*(eta1_max-eta1t))*sin(eta2t)
          eta1t = eta1 + 0.5_f64*dt*k2eta1
          eta2t = eta2 +  0.5_f64*dt*k2eta2
          k3eta1 = eta1t*(eta1t-eta1_min)*(eta1_max-eta1t)*cos(eta2t)
          k3eta2 = -(eta1t*(eta1_max+eta1_min-2._f64*eta1t)+2._f64*(eta1t-eta1_min)*(eta1_max-eta1t))*sin(eta2t)
          eta1t = eta1 + dt*k3eta1
          eta2t = eta2 +  dt*k3eta2
          k4eta1 = eta1t*(eta1t-eta1_min)*(eta1_max-eta1t)*cos(eta2t)
          k4eta2 = -(eta1t*(eta1_max+eta1_min-2._f64*eta1t)+2._f64*(eta1t-eta1_min)*(eta1_max-eta1t))*sin(eta2t)
          eta1 = eta1 + dt/6._f64*(k1eta1+2._f64*k2eta1+2._f64*k3eta1+k4eta1)
          eta2 = eta2 + dt/6._f64*(k1eta2+2._f64*k2eta2+2._f64*k3eta2+k4eta2)
        endif
        ! rotation
        if(field_case==2)then
          eta2 = eta2 + dt
        endif
        ! translation
        if(field_case==3)then
          x1 = x1_array(i,j) + a1*dt
          x2 = x2_array(i,j) + a2*dt
          eta1 = sqrt(x1**2+x2**2)
          if(x2>=0)then
            eta2 = acos(x1/eta1)
          else
            eta2 = 2._f64*sll_pi-acos(x1/eta1)
          endif
        endif
          
        ! correction of eta2
        if(eta2>2._f64*sll_pi)then
          eta2= eta2-2._f64*sll_pi
        endif
        if(eta2<0._f64)then
          eta2= eta2+2._f64*sll_pi
        endif
          					
        if (test_case==1) then
          val = exp(-3000._f64*(eta1-1.5_f64*eta1_min)**2)
        endif
        if (test_case==2) then
          val = exp(-100._f64*(eta1-0.5_f64*(eta1_max+eta1_min))**2)
        endif
        if (test_case==3) then
          val = cos(eta2)
        endif
        if (test_case==4) then
          if(eta1>=0.4.and.eta1<=0.6)then
            val = exp(-30._f64*(eta2-0.5*sll_pi)**2)
          else
            val = 0._f64
          endif
        endif
        if (test_case==5) then
          val = 1._f64
        endif
        if (test_case==6) then
          val = exp(-100._f64*(eta1-0.5_f64*(eta1_max+eta1_min))**2)*exp(-30._f64*(eta2-0.5*sll_pi)**2)
        endif
        if (test_case==7) then
          val = 0._f64
          if ((i==17).and.(j==12)) then
            val = real(32,f64)
          end if
        endif
			
        diag(1,step) = diag(1,step) + val*delta_eta1*delta_eta2*jac_array(i,j)
        diag(2,step) = diag(2,step) + val*delta_eta1*delta_eta2
					
        if(mod(step,visu_step)==0)then
          write(900,*) x1_array(i,j),x2_array(i,j)
          write(800,'(a,f0.3)',advance='no') ' ',val
        endif
        f(i,j) = val
      enddo
      write(800,*) ' '
    enddo
    write(900,*) ' '
    write(800,*) ' '
  enddo  
  close(900)
  close(800)
  
  open(unit=900,file='diag0.dat')
  do step=0,nb_step
    write(900,*) real(step,f64)*dt, diag(1,step), diag(1,step)/diag(1,0)-1._f64, &
                  diag(2,step), diag(2,step)/diag(2,0)-1._f64
  enddo
  close(900)
  

	! Evolution in time - Numerical solution	
  do step=1,nb_step
  
    do i=1,Neta1+1
      do j=1,Neta2+1
        fh_bsl_nc(i,j) = fh_bsl_nc(i,j)*jac_array(i,j)
        fh_fsl(i,j)    = fh_fsl(i,j)*jac_array(i,j)
      enddo
      fh_bsl(i,Neta2+1)    = fh_bsl(i,1)
      fh_bsl_nc(i,Neta2+1) = fh_bsl_nc(i,1)
      fh_fsl(i,Neta2+1)    = fh_fsl(i,1)
      fh_fsl_nc(i,Neta2+1) = fh_fsl_nc(i,1)
    enddo
      
    call compute_spline_2D(fh_bsl,spl_bsl)
    call compute_spline_2D(fh_bsl_nc,spl_bsl_nc)
    call compute_spline_2D(fh_fsl,spl_fsl)
    call compute_spline_2D(fh_fsl_nc,spl_fsl_nc)
    
    do i=1,Neta1+1
      do j=1,Neta2+1
      				    
        ! ------------- BSL part -----------------
        
        eta1 = eta1_min + real(i-1,f64)*delta_eta1
        eta2 = eta2_min + real(j-1,f64)*delta_eta2
        
        ! RK4
        if(field_case==1)then
          eta1t = eta1
          eta2t = eta2
          k1eta1 = eta1t*(eta1t-eta1_min)*(eta1_max-eta1t)*cos(eta2t)
          k1eta2 = -(eta1t*(eta1_max+eta1_min-2._f64*eta1t)+2._f64*(eta1t-eta1_min)*(eta1_max-eta1t))*sin(eta2t)
          eta1t = eta1 + 0.5_f64*dt*k1eta1
          eta2t = eta2 +  0.5_f64*dt*k1eta2
          k2eta1 = eta1t*(eta1t-eta1_min)*(eta1_max-eta1t)*cos(eta2t)
          k2eta2 = -(eta1t*(eta1_max+eta1_min-2._f64*eta1t)+2._f64*(eta1t-eta1_min)*(eta1_max-eta1t))*sin(eta2t)
          eta1t = eta1 + 0.5_f64*dt*k2eta1
          eta2t = eta2 +  0.5_f64*dt*k2eta2
          k3eta1 = eta1t*(eta1t-eta1_min)*(eta1_max-eta1t)*cos(eta2t)
          k3eta2 = -(eta1t*(eta1_max+eta1_min-2._f64*eta1t)+2._f64*(eta1t-eta1_min)*(eta1_max-eta1t))*sin(eta2t)
          eta1t = eta1 + dt*k3eta1
          eta2t = eta2 +  dt*k3eta2
          k4eta1 = eta1t*(eta1t-eta1_min)*(eta1_max-eta1t)*cos(eta2t)
          k4eta2 = -(eta1t*(eta1_max+eta1_min-2._f64*eta1t)+2._f64*(eta1t-eta1_min)*(eta1_max-eta1t))*sin(eta2t)
          eta1 = eta1 + dt/6._f64*(k1eta1+2._f64*k2eta1+2._f64*k3eta1+k4eta1)
          eta2 = eta2 + dt/6._f64*(k1eta2+2._f64*k2eta2+2._f64*k3eta2+k4eta2)
        endif
        if(field_case==2)then
          eta2 = eta2 + dt
        endif
        if(field_case==3)then
          x1 = x1_array(i,j) + a1*dt
          x2 = x2_array(i,j) + a2*dt
          eta1 = sqrt(x1**2+x2**2)
          if(x2>=0)then
            eta2 = acos(x1/eta1)
          else
            eta2 = 2._f64*sll_pi-acos(x1/eta1)
          endif
        endif
        
        ! correction of eta2
        if(eta2>2._f64*sll_pi)then
          eta2= eta2-2._f64*sll_pi
        endif
        if(eta2<0._f64)then
          eta2= eta2+2._f64*sll_pi
        endif
        
        fh_bsl(i,j)    = interpolate_value_2D(eta1,eta2,spl_bsl)
        fh_bsl_nc(i,j) = interpolate_value_2D(eta1,eta2,spl_bsl_nc)/jac_array(i,j)
    
        diag(3,step) = diag(3,step) + fh_bsl(i,j) * delta_eta1*delta_eta2*jac_array(i,j)
        diag(4,step) = diag(4,step) + fh_bsl(i,j) * delta_eta1*delta_eta2
        
        diag(5,step) = diag(5,step) + fh_bsl_nc(i,j) * delta_eta1*delta_eta2*jac_array(i,j)
        diag(6,step) = diag(6,step) + fh_bsl_nc(i,j) * delta_eta1*delta_eta2
        
        ! ------------- FSL part -----------------
        
        dt = -dt
        
        eta1 = eta1_min + real(i-1,f64)*delta_eta1
        eta2 = eta2_min + real(j-1,f64)*delta_eta2
        
        ! RK4
        if(field_case==1)then
          eta1t = eta1
          eta2t = eta2
          k1eta1 = eta1t*(eta1t-eta1_min)*(eta1_max-eta1t)*cos(eta2t)
          k1eta2 = -(eta1t*(eta1_max+eta1_min-2._f64*eta1t)+2._f64*(eta1t-eta1_min)*(eta1_max-eta1t))*sin(eta2t)
          eta1t = eta1 + 0.5_f64*dt*k1eta1
          eta2t = eta2 +  0.5_f64*dt*k1eta2
          k2eta1 = eta1t*(eta1t-eta1_min)*(eta1_max-eta1t)*cos(eta2t)
          k2eta2 = -(eta1t*(eta1_max+eta1_min-2._f64*eta1t)+2._f64*(eta1t-eta1_min)*(eta1_max-eta1t))*sin(eta2t)
          eta1t = eta1 + 0.5_f64*dt*k2eta1
          eta2t = eta2 +  0.5_f64*dt*k2eta2
          k3eta1 = eta1t*(eta1t-eta1_min)*(eta1_max-eta1t)*cos(eta2t)
          k3eta2 = -(eta1t*(eta1_max+eta1_min-2._f64*eta1t)+2._f64*(eta1t-eta1_min)*(eta1_max-eta1t))*sin(eta2t)
          eta1t = eta1 + dt*k3eta1
          eta2t = eta2 +  dt*k3eta2
          k4eta1 = eta1t*(eta1t-eta1_min)*(eta1_max-eta1t)*cos(eta2t)
          k4eta2 = -(eta1t*(eta1_max+eta1_min-2._f64*eta1t)+2._f64*(eta1t-eta1_min)*(eta1_max-eta1t))*sin(eta2t)
          eta1 = eta1 + dt/6._f64*(k1eta1+2._f64*k2eta1+2._f64*k3eta1+k4eta1)
          eta2 = eta2 + dt/6._f64*(k1eta2+2._f64*k2eta2+2._f64*k3eta2+k4eta2)
        endif
        if(field_case==2)then
          eta2 = eta2 + dt
        endif
        if(field_case==3)then
          x1 = x1_array(i,j) + a1*dt
          x2 = x2_array(i,j) + a2*dt
          eta1 = sqrt(x1**2+x2**2)
          if(x2>=0)then
            eta2 = acos(x1/eta1)
          else
            eta2 = 2._f64*sll_pi-acos(x1/eta1)
          endif
        endif
        
        ! correction of eta2
        if(eta2>2._f64*sll_pi)then
          eta2= eta2-2._f64*sll_pi
        endif
        if(eta2<0._f64)then
          eta2= eta2+2._f64*sll_pi
        endif
        
        ! conditions for the mass preservation
        if ((i==1).or.(i==Neta1+1)) then
          eta1 = eta1_min + real(i-1,f64)*delta_eta1
        else
          eta1 = min(max(eta1, eta1_min + delta_eta1), eta1_min + real(Neta1,f64)*delta_eta1)
        end if
        
        eta1feet(i,j) = eta1
        eta2feet(i,j) = eta2
				
        dt = -dt
        
      enddo
    enddo
    
    call deposit_value_2D(eta1feet,eta2feet,spl_fsl,fh_fsl)
    call deposit_value_2D(eta1feet,eta2feet,spl_fsl_nc,fh_fsl_nc)
    
    do i=1,Neta1+1
      do j=1,Neta2+1
        fh_fsl(i,j) = fh_fsl(i,j)/jac_array(i,j)
        
        if (abs(fh_fsl_nc(i,j)-1._f64).ge.1E-12) then
          print*,i,j,fh_fsl_nc(i,j)
        end if
      enddo
    enddo
    
    do i=1,Neta1+1
      do j=1,Neta2+1
        diag(7,step)  = diag(7,step)  + fh_fsl(i,j)*delta_eta1*delta_eta2*jac_array(i,j)
        diag(8,step)  = diag(8,step)  + fh_fsl(i,j)*delta_eta1*delta_eta2
        diag(9,step)  = diag(9,step)  + fh_fsl_nc(i,j)*delta_eta1*delta_eta2*jac_array(i,j)
        diag(10,step) = diag(10,step) + fh_fsl_nc(i,j)*delta_eta1*delta_eta2
      enddo
    enddo
    
  enddo
  
  ! diagnostics
  open(unit=900,file='diag.dat')
  do step=0,nb_step
    write(900,*) real(step,f64)*dt, &
      diag(3,step),diag(3,step)/diag(1,0)-1._f64,diag(4,step),diag(4,step)/diag(1,0)-1._f64, &
      diag(5,step),diag(5,step)/diag(1,0)-1._f64,diag(6,step),diag(6,step)/diag(1,0)-1._f64, &
      diag(7,step),diag(7,step)/diag(1,0)-1._f64,diag(8,step),diag(8,step)/diag(1,0)-1._f64, &
      diag(9,step),diag(9,step)/diag(1,0)-1._f64,diag(10,step),diag(10,step)/diag(1,0)-1._f64
  end do
  close(900)

  ! L^\infty
  val = 0._f64
  val_bsl    = 0._f64
  val_bsl_nc = 0._f64
  val_fsl    = 0._f64
  val_fsl_nc = 0._f64
  do i=1,Neta1+1
    do j=1,Neta2+1
      val_bsl    = max(val_bsl,   abs(f(i,j)-fh_bsl(i,j)))
      val_bsl_nc = max(val_bsl_nc,abs(f(i,j)-fh_bsl_nc(i,j)))
      val_fsl    = max(val_fsl,   abs(f(i,j)-fh_fsl(i,j)))
      val_fsl_nc = max(val_fsl_nc,abs(f(i,j)-fh_fsl_nc(i,j)))
      val        = max(val,       abs(fh_fsl(i,j)-fh_bsl(i,j)))
    enddo
  enddo  
  print*,N,val_bsl,val_bsl_nc,val_fsl,val_fsl_nc,val

  open(unit=900,file='fh.dat')  
  do i=1,Neta1+1
    eta1 = eta1_min+real(i-1,f64)*delta_eta1
    do j=1,Neta2+1
      eta2 = eta2_min + real(j-1,f64)*delta_eta2
      write(900,*) eta1,eta2,x1_array(i,j),x2_array(i,j),f(i,j),fh_bsl(i,j),fh_bsl_nc(i,j),fh_fsl(i,j),fh_fsl_nc(i,j)
    enddo
    write(900,*) ' '
  enddo
  close(900)
  
end program