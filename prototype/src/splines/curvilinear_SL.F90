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
  sll_real64 :: a1,a2,eta1t,eta2t,eta1c,eta2c,k1eta1,k2eta1,k3eta1,k4eta1,k1eta2,k2eta2,k3eta2,k4eta2
  sll_real64,dimension(:,:), pointer :: f,fh_bsl,fh_bsl_nc,fh_fsl,fh_fsl_nc
  sll_real64,dimension(:,:), pointer :: x1_array,x2_array,x1c_array,x2c_array,jac_array,eta1feet,eta2feet,eta1tot,eta2tot
  sll_real64,dimension(:,:), pointer :: diag
	 
  ! for the python script polar-ex1e.py
  namelist /param/ N
  N = 32
  read(*,NML=param);open(unit=900,file="gyrof.param");write(900,NML=param);close(900)
  
  ! ---- * mesh type * ----
  ! 1 : cartesian
  ! 2 : polar
  ! 3 : polar like (polar with r^2)
  ! 4 : Collela  
  mesh_case = 1
  if (mesh_case > 4) then
    print*,'Non existing case'
    print*,'mesh_case = 1 (cartesian), 2 (polar), 3 (polar-like) or 4 (Collela)'
    STOP
  end if
  
  ! mesh type : cartesian
  ! domain    : square [eta1_min eta1_max] x [eta2_min eta2_max]
  ! BC        : periodic-periodic
  if (mesh_case==1) then
    eta1_min = 0._f64
    eta1_max = 2._f64*sll_pi
    eta2_min = 0._f64
    eta2_max = 2._f64*sll_pi
    
    bc1_type = HERMITE_SPLINE
    bc2_type = HERMITE_SPLINE
  end if
  ! mesh type : polar or polar-like
  ! domain    : disc of radius eta1_max with a hole of radius eta1_min
  ! BC        : hermite-periodic
  if ((mesh_case==2).or.(mesh_case==3)) then
    eta1_min = 0.2_f64
    eta1_max = 0.8_f64
    eta2_min = 0._f64
    eta2_max = 2._f64*sll_pi

    bc1_type = HERMITE_SPLINE
    bc2_type = PERIODIC_SPLINE
  end if
  ! mesh type : Collela
  ! domain    : square [eta1_min eta1_max] x [eta2_min eta2_max]
  ! BC        : periodic-periodic
  if (mesh_case==4) then
    eta1_min = 0._f64
    eta1_max = 1._f64
    eta2_min = 0._f64
    eta2_max = 2._f64
    
    bc1_type = PERIODIC_SPLINE
    bc2_type = PERIODIC_SPLINE
  end if
    
  Neta1 = N
  Neta2 = N
  
  delta_eta1 = (eta1_max-eta1_min)/real(Neta1,f64)
  delta_eta2 = (eta2_max-eta2_min)/real(Neta2,f64)  
  
  ! ---- * distribution function * ----
  ! 1 : f=1
  ! 2 : cos eta2
  ! 3 : centered gaussian in eta1 and eta2 (mesh_case = 1)
  ! 4 : centered gaussian in eta1 and eta2 (mesh_case = 2,3) 
  ! 5 : gaussian in eta1
  ! 6 : centered dirac 
  test_case = 3
  if (test_case > 6) then
    print*,'Non existing case'
    print*,'test_case = 1 (f=1), 2 (cos(eta2)), 3 (centered gaussian in eta1 and eta2 for cartesian mesh),'
    print*,'4 (centered gaussian in eta1 and eta2 for polar/polar-like mesh), 5 (gaussian in eta1) or 6 (centered dirac)'
    STOP
  end if
  
  ! ---- * advecton field * ---- 
  ! 1 : translation of vector (a1,a2)
  ! 2 : rotation
  ! 3 : divergence free complex field (mesh_case = 1)
  ! 4 : divergence free complex field (mesh_case = 2,3)
  ! 5 : divergence free non-symmetric field (mesh_case = 1)
  field_case = 5
  if (field_case > 5) then
    print*,'Non existing case'
    print*,'field_case = 1 (transaltion), 2 (rotation),'
    print*,'3 (complex divergence free field for cartesian mesh), 4 (complex divergence free field for polar/polar-like mesh)'
    print*,'or 5 (divergence free non-symmetric field)'
    STOP
  end if
  a1=0.001_f64
  a2=0.002_f64
  
  ! ---- * time-scheme parameters * ----
  dt = 0.001_f64*delta_eta1 ! 1._f64/1._f64
  nb_step = 1 !floor(1._f64/dt)
  
  ! visualization parameters
  visu_step = 1
  
  ! ---- * Beginning of the program * ----
	
  print *,'# N=',N
  print *,'# T=',real(nb_step,f64)*dt
  print *,'# mesh_case=',mesh_case
  print *,'# test_case=',test_case
  print *,'# field_case=',field_case
	
  ! allocations of the arrays
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
  SLL_ALLOCATE(eta1tot(Neta1+1,Neta2+1), err)
  SLL_ALLOCATE(eta2tot(Neta1+1,Neta2+1), err)
  SLL_ALLOCATE(jac_array(Neta1+1,Neta2+1), err)
  SLL_ALLOCATE(diag(10,0:nb_step), err)
	
  ! creation of the splines
  spl_bsl => new_spline_2D(Neta1+1, Neta2+1, &
    eta1_min, eta1_max, &
    0._f64, 2._f64*sll_pi, &
    bc1_type, bc2_type)
  spl_bsl_nc => new_spline_2D(Neta1+1, Neta2+1, &
    eta1_min, eta1_max, &
    0._f64, 2._f64*sll_pi, &
    bc1_type, bc2_type)
  spl_fsl => new_spline_2D(Neta1+1, Neta2+1, &
    eta1_min, eta1_max, &
    0._f64, 2._f64*sll_pi, &
    bc1_type, bc2_type)
  spl_fsl_nc => new_spline_2D(Neta1+1, Neta2+1, &
    eta1_min, eta1_max, &
    0._f64, 2._f64*sll_pi, &
    bc1_type, bc2_type)
  
  ! ---- Initializations ----
  
  ! Analytic distribution function and data for the mesh
  open(unit=900,file='f0.dat')
  do i=1,Neta1+1
    eta1 = eta1_min + real(i-1,f64)*delta_eta1
    do j=1,Neta2+1
      eta2 = eta2_min + real(j-1,f64)*delta_eta2
      
      ! cartesian mesh
      if (mesh_case==1) then
        x1_array(i,j) = eta1
        x2_array(i,j) = eta2
        jac_array(i,j) = 1._f64
      end if
      
      ! polar mesh
      if (mesh_case==2) then
        x1_array(i,j) = eta1*cos(eta2)
        x2_array(i,j) = eta1*sin(eta2)
        jac_array(i,j) = eta1
      end if
      
      ! polar-like mesh
      if (mesh_case==3) then
        x1_array(i,j) = eta1*eta1*cos(eta2)
        x2_array(i,j) = eta1*eta1*sin(eta2)
        jac_array(i,j) = 2._f64*eta1*eta1*eta1
      end if
      
      ! Collela mesh
      if (mesh_case==4) then
        alpha_mesh = 1.e-1_f64
        
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
      
      ! total displacement
      eta1tot(i,j) = eta1
      eta2tot(i,j) = eta2

      ! test-function
      if (test_case==1) then
        f(i,j) = 1._f64
      endif
      if (test_case==2) then
        f(i,j) = cos(eta2)
      endif
      if (test_case==3) then
        eta1c = 0.5_f64*(eta1_max+eta1_min)
        eta2c = 0.5_f64*(eta2_max+eta2_min)
        f(i,j) = exp(-2._f64*(eta1-eta1c)**2)*exp(-2._f64*(eta2-eta2c)**2)
      endif
      if (test_case==4) then
        eta1c = 0.5_f64*(eta1_max+eta1_min)
        eta2c = 0.5_f64*(eta2_max+eta2_min)
        f(i,j) = exp(-100._f64*(eta1-eta1c)**2)*exp(-30._f64*(eta2-eta2c)**2)
      endif
      if (test_case==5) then
        f(i,j) = exp(-100._f64*(eta1-0.5_f64*(eta1_max+eta1_min))**2)
      endif
      if (test_case==6) then
        f(i,j) = 0._f64
        if ((i==(Neta1+1)/2).and.(j==(Neta2+1)/2)) then
          f(i,j) = 1._f64
        end if
      endif
      write(900,*) x1_array(i,j),x2_array(i,j),eta1,eta2,f(i,j)
    enddo
    write(900,*) ' ' 
  enddo
  close(900)

  ! Distribution functions for the four methods
  fh_bsl    = f
  fh_bsl_nc = f
  fh_fsl    = f
  fh_fsl_nc = f
  
  ! Diagnostic with the weighted mass (temp1) and the "classical" mass (temp2) 
  diag = 0._f64
  
  tmp1 = 0._f64
  tmp2 = 0._f64
  do i=1,Neta1+bc1_type
    do j=1,Neta2+bc2_type
      tmp1 = tmp1 + f(i,j)*delta_eta1*delta_eta2*jac_array(i,j)
      tmp2 = tmp2 + f(i,j)*delta_eta1*delta_eta2
    enddo
  enddo
  
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
  
  ! ---- Evolution in time ----

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
      				    
        ! ------------ Analytic part -----------------
        
        eta1 = eta1tot(i,j)
        eta2 = eta2tot(i,j)
        
        if (field_case==1) then          
          if (mesh_case==1) then
            eta1 = eta1 + a1*dt 
            eta2 = eta2 + a2*dt
          endif
          if (mesh_case==2) then
            x1 = x1_array(i,j) + a1*dt
            x2 = x2_array(i,j) + a2*dt
            eta1 = sqrt(x1**2+x2**2)
            if (x2>=0) then
              eta2 = acos(x1/eta1)
            else
              eta2 = 2._f64*sll_pi-acos(x1/eta1)
            endif
          endif
          if (mesh_case==3) then
            x1 = x1_array(i,j) + a1*dt
            x2 = x2_array(i,j) + a2*dt
            eta1 = (x1**2+x2**2)**(0.25_f64)
            if (x2>=0) then
              eta2 = acos(x1/(eta1*eta1))
            else
              eta2 = 2._f64*sll_pi-acos(x1/(eta1*eta1))
            endif
          end if
        end if
        if (field_case==2) then
          if (mesh_case==1) then            
            eta1c = 0.5_f64*(eta1_max+eta1_min)
            eta2c = 0.5_f64*(eta2_max+eta2_min)
        
            eta1t = cos(dt)*(eta1-eta1c)  + sin(dt)*(eta2-eta2c) + eta1c
            eta2t = -sin(dt)*(eta1-eta1c) + cos(dt)*(eta2-eta2c) + eta2c
 
            eta1 = eta1t
            eta2 = eta2t
          end if
          if ((mesh_case==2).or.(mesh_case==3)) then
            eta2 = eta2 + dt
          end if
        endif
        if (field_case==3) then
          eta1t = eta1
          eta2t = eta2
          k1eta1 = -(3*eta1t*eta1t+6*eta1t*eta2t+3*eta2*eta2t+1)
          k1eta2 = 3*eta1t*eta1t+6*eta1t*eta2t+3*eta2*eta2t+1
          eta1t = eta1 + 0.5_f64*dt*k1eta1
          eta2t = eta2 + 0.5_f64*dt*k1eta2
          k2eta1 = -(3*eta1t*eta1t+6*eta1t*eta2t+3*eta2*eta2t+1)
          k2eta2 = 3*eta1t*eta1t+6*eta1t*eta2t+3*eta2*eta2t+1
          eta1t = eta1 + 0.5_f64*dt*k2eta1
          eta2t = eta2 + 0.5_f64*dt*k2eta2
          k3eta1 = -(3*eta1t*eta1t+6*eta1t*eta2t+3*eta2*eta2t+1)
          k3eta2 = 3*eta1t*eta1t+6*eta1t*eta2t+3*eta2*eta2t+1
          eta1t = eta1 + dt*k3eta1
          eta2t = eta2 + dt*k3eta2
          k4eta1 = -(3*eta1t*eta1t+6*eta1t*eta2t+3*eta2*eta2t+1)
          k4eta2 = 3*eta1t*eta1t+6*eta1t*eta2t+3*eta2*eta2t+1
          eta1 = eta1 + dt/6._f64*(k1eta1+2._f64*k2eta1+2._f64*k3eta1+k4eta1)
          eta2 = eta2 + dt/6._f64*(k1eta2+2._f64*k2eta2+2._f64*k3eta2+k4eta2)
        endif
        ! RK4
        if (field_case==4) then
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
        if (field_case==5) then
          eta1t = eta1*cos(2*sqrt(a1*a2)*step*dt) - eta2*sqrt(a2/a1)*sin(2*sqrt(a1*a2)*step*dt)
          eta2t = eta1*sqrt(a1/a2)*sin(2*sqrt(a1*a2)*step*dt) - eta2*sin(2*sqrt(a1*a2)*step*dt)
          eta1 = eta1t
          eta2 = eta2t
        end if
                
        ! Corrections in case of periodic BC
        if (bc1_type==PERIODIC_SPLINE) then
          if (eta1>eta1_max) then
            eta1 = eta1-(eta1_max-eta1_min)
          endif
          if (eta1<eta1_min) then
            eta1 = eta1+(eta1_max-eta1_min)
          endif
        end if
        if (bc2_type==PERIODIC_SPLINE) then
          if (eta2>eta2_max) then
            eta2 = eta2-(eta2_max-eta2_min)
          endif
          if (eta2<eta2_min) then
            eta2 = eta2+(eta2_max-eta2_min)
          endif
        end if
          					
        if (test_case==1) then
          f(i,j) = 1._f64
        endif
        if (test_case==2) then
          f(i,j) = cos(eta2)
        endif
        if (test_case==4) then
          eta1c = 0.5_f64*(eta1_max+eta1_min)
          eta2c = 0.5_f64*(eta2_max+eta2_min)
          f(i,j) = exp(-100._f64*(eta1-eta1c)**2)*exp(-30._f64*(eta2-eta2c)**2)
        endif
        if (test_case==3) then
          eta1c = 0.5_f64*(eta1_max+eta1_min)
          eta2c = 0.5_f64*(eta2_max+eta2_min)
          f(i,j) = exp(-2._f64*(eta1-eta1c)**2)*exp(-2._f64*(eta2-eta2c)**2)
        endif
        if (test_case==5) then
          f(i,j) = exp(-100._f64*(eta1-0.5_f64*(eta1_max+eta1_min))**2)
        endif
        if (test_case==6) then
          f(i,j) = 0._f64
          if ((i==(Neta1+1)/2).and.(j==(Neta2+1)/2)) then
            f(i,j) = 1._f64
          end if
        endif
        
        eta1tot(i,j) = eta1
        eta2tot(i,j) = eta2
        
        ! ------------ BSL part -----------------
        
        eta1 = eta1_min + real(i-1,f64)*delta_eta1
        eta2 = eta2_min + real(j-1,f64)*delta_eta2
        
        if (field_case==1) then          
          if (mesh_case==1) then
            eta1 = eta1 + a1*dt 
            eta2 = eta2 + a2*dt
          endif
          if (mesh_case==2) then
            x1 = x1_array(i,j) + a1*dt
            x2 = x2_array(i,j) + a2*dt
            eta1 = sqrt(x1**2+x2**2)
            if (x2>=0) then
              eta2 = acos(x1/eta1)
            else
              eta2 = 2._f64*sll_pi-acos(x1/eta1)
            endif
          endif
          if (mesh_case==3) then
            x1 = x1_array(i,j) + a1*dt
            x2 = x2_array(i,j) + a2*dt
            eta1 = (x1**2+x2**2)**(0.25_f64)
            if (x2>=0) then
              eta2 = acos(x1/(eta1*eta1))
            else
              eta2 = 2._f64*sll_pi-acos(x1/(eta1*eta1))
            endif
          end if
        end if
        if (field_case==2) then
          if (mesh_case==1) then            
            eta1c = 0.5_f64*(eta1_max+eta1_min)
            eta2c = 0.5_f64*(eta2_max+eta2_min)
        
            eta1t = cos(dt)*(eta1-eta1c)  + sin(dt)*(eta2-eta2c) + eta1c
            eta2t = -sin(dt)*(eta1-eta1c) + cos(dt)*(eta2-eta2c) + eta2c
 
            eta1 = eta1t
            eta2 = eta2t
          end if
          if ((mesh_case==2).or.(mesh_case==3)) then
            eta2 = eta2 + dt
          end if
        endif
        if (field_case==3) then
          eta1t = eta1
          eta2t = eta2
          k1eta1 = -(3*eta1t*eta1t+6*eta1t*eta2t+3*eta2*eta2t+1)
          k1eta2 = 3*eta1t*eta1t+6*eta1t*eta2t+3*eta2*eta2t+1
          eta1t = eta1 + 0.5_f64*dt*k1eta1
          eta2t = eta2 + 0.5_f64*dt*k1eta2
          k2eta1 = -(3*eta1t*eta1t+6*eta1t*eta2t+3*eta2*eta2t+1)
          k2eta2 = 3*eta1t*eta1t+6*eta1t*eta2t+3*eta2*eta2t+1
          eta1t = eta1 + 0.5_f64*dt*k2eta1
          eta2t = eta2 + 0.5_f64*dt*k2eta2
          k3eta1 = -(3*eta1t*eta1t+6*eta1t*eta2t+3*eta2*eta2t+1)
          k3eta2 = 3*eta1t*eta1t+6*eta1t*eta2t+3*eta2*eta2t+1
          eta1t = eta1 + dt*k3eta1
          eta2t = eta2 + dt*k3eta2
          k4eta1 = -(3*eta1t*eta1t+6*eta1t*eta2t+3*eta2*eta2t+1)
          k4eta2 = 3*eta1t*eta1t+6*eta1t*eta2t+3*eta2*eta2t+1
          eta1 = eta1 + dt/6._f64*(k1eta1+2._f64*k2eta1+2._f64*k3eta1+k4eta1)
          eta2 = eta2 + dt/6._f64*(k1eta2+2._f64*k2eta2+2._f64*k3eta2+k4eta2)
        endif
        ! RK4
        if (field_case==4) then
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
        if (field_case==5) then
          eta1t = eta1*cos(2*sqrt(a1*a2)*step*dt) - eta2*sqrt(a2/a1)*sin(2*sqrt(a1*a2)*step*dt)
          eta2t = eta1*sqrt(a1/a2)*sin(2*sqrt(a1*a2)*step*dt) - eta2*sin(2*sqrt(a1*a2)*step*dt)
          eta1 = eta1t
          eta2 = eta2t
        end if
                
        ! Corrections in case of periodic BC
        if (bc1_type==PERIODIC_SPLINE) then
          if (eta1>eta1_max) then
            eta1 = eta1-(eta1_max-eta1_min)
          endif
          if (eta1<eta1_min) then
            eta1 = eta1+(eta1_max-eta1_min)
          endif
        end if
        if (bc2_type==PERIODIC_SPLINE) then
          if (eta2>eta2_max) then
            eta2 = eta2-(eta2_max-eta2_min)
          endif
          if (eta2<eta2_min) then
            eta2 = eta2+(eta2_max-eta2_min)
          endif
        end if
        
        fh_bsl(i,j)    = interpolate_value_2D(eta1,eta2,spl_bsl)
        fh_bsl_nc(i,j) = interpolate_value_2D(eta1,eta2,spl_bsl_nc)/jac_array(i,j)
        
        ! ------------ FSL part -----------------
        
        dt = -dt
        
        eta1 = eta1_min + real(i-1,f64)*delta_eta1
        eta2 = eta2_min + real(j-1,f64)*delta_eta2
        
        if (field_case==1) then          
          if (mesh_case==1) then
            eta1 = eta1 + a1*dt 
            eta2 = eta2 + a2*dt
          endif
          if (mesh_case==2) then
            x1 = x1_array(i,j) + a1*dt
            x2 = x2_array(i,j) + a2*dt
            eta1 = sqrt(x1**2+x2**2)
            if (x2>=0) then
              eta2 = acos(x1/eta1)
            else
              eta2 = 2._f64*sll_pi-acos(x1/eta1)
            endif
          endif
          if (mesh_case==3) then
            x1 = x1_array(i,j) + a1*dt
            x2 = x2_array(i,j) + a2*dt
            eta1 = (x1**2+x2**2)**(0.25_f64)
            if (x2>=0) then
              eta2 = acos(x1/(eta1*eta1))
            else
              eta2 = 2._f64*sll_pi-acos(x1/(eta1*eta1))
            endif
          end if
        end if
        if (field_case==2) then
          if (mesh_case==1) then            
            eta1c = 0.5_f64*(eta1_max+eta1_min)
            eta2c = 0.5_f64*(eta2_max+eta2_min)
        
            eta1t = cos(dt)*(eta1-eta1c)  + sin(dt)*(eta2-eta2c) + eta1c
            eta2t = -sin(dt)*(eta1-eta1c) + cos(dt)*(eta2-eta2c) + eta2c
 
            eta1 = eta1t
            eta2 = eta2t
          end if
          if ((mesh_case==2).or.(mesh_case==3)) then
            eta2 = eta2 + dt
          end if
        endif
        if (field_case==3) then
          eta1t = eta1
          eta2t = eta2
          k1eta1 = -(3*eta1t*eta1t+6*eta1t*eta2t+3*eta2*eta2t+1)
          k1eta2 = 3*eta1t*eta1t+6*eta1t*eta2t+3*eta2*eta2t+1
          eta1t = eta1 + 0.5_f64*dt*k1eta1
          eta2t = eta2 +  0.5_f64*dt*k1eta2
          k2eta1 = -(3*eta1t*eta1t+6*eta1t*eta2t+3*eta2*eta2t+1)
          k2eta2 = 3*eta1t*eta1t+6*eta1t*eta2t+3*eta2*eta2t+1
          eta1t = eta1 + 0.5_f64*dt*k2eta1
          eta2t = eta2 +  0.5_f64*dt*k2eta2
          k3eta1 = -(3*eta1t*eta1t+6*eta1t*eta2t+3*eta2*eta2t+1)
          k3eta2 = 3*eta1t*eta1t+6*eta1t*eta2t+3*eta2*eta2t+1
          eta1t = eta1 + dt*k3eta1
          eta2t = eta2 +  dt*k3eta2
          k4eta1 = -(3*eta1t*eta1t+6*eta1t*eta2t+3*eta2*eta2t+1)
          k4eta2 = 3*eta1t*eta1t+6*eta1t*eta2t+3*eta2*eta2t+1
          eta1 = eta1 + dt/6._f64*(k1eta1+2._f64*k2eta1+2._f64*k3eta1+k4eta1)
          eta2 = eta2 + dt/6._f64*(k1eta2+2._f64*k2eta2+2._f64*k3eta2+k4eta2)
        endif
        ! RK4
        if (field_case==4) then
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
        if (field_case==5) then
          eta1t = eta1*cos(2*sqrt(a1*a2)*step*dt) - eta2*sqrt(a2/a1)*sin(2*sqrt(a1*a2)*step*dt)
          eta2t = eta1*sqrt(a1/a2)*sin(2*sqrt(a1*a2)*step*dt) - eta2*sin(2*sqrt(a1*a2)*step*dt)
          eta1 = eta1t
          eta2 = eta2t
        end if
                
        ! Corrections in case of periodic BC
        if (bc1_type==PERIODIC_SPLINE) then
          if (eta1>eta1_max) then
            eta1 = eta1-(eta1_max-eta1_min)
          endif
          if (eta1<eta1_min) then
            eta1 = eta1+(eta1_max-eta1_min)
          endif
        end if
        if (bc2_type==PERIODIC_SPLINE) then
          if (eta2>eta2_max) then
            eta2 = eta2-(eta2_max-eta2_min)
          endif
          if (eta2<eta2_min) then
            eta2 = eta2+(eta2_max-eta2_min)
          endif
        end if
                                
        ! conditions for the mass preservation
        !if ((i==1).or.(i==Neta1+1)) then
        !  eta1 = eta1_min + real(i-1,f64)*delta_eta1
        !else
        !  eta1 = min(max(eta1, eta1_min + delta_eta1), eta1_min + real(Neta1,f64)*delta_eta1)
        !end if
        
        eta1feet(i,j) = eta1
        eta2feet(i,j) = eta2
				
        dt = -dt
        
        ! Actual coordinates 
        if (mesh_case==2) then 
          x1 = eta1
          x2 = eta2
        end if
        if (mesh_case==2) then 
          x1 = eta1*cos(eta2)
          x2 = eta1*sin(eta2)
        end if
        if (mesh_case==3) then 
          x1 = eta1*eta1*cos(eta2)
          x2 = eta1*eta1*sin(eta2)
        end if
        if (mesh_case==4) then 
          print*,'Maillage de Collela non encore disponible'
          STOP
        end if
        
      enddo
    enddo
    
    print*,'avant',fh_fsl(1,1),fh_fsl(17,1)
    
    call deposit_value_2D(eta1feet,eta2feet,spl_fsl,fh_fsl)
    call deposit_value_2D(eta1feet,eta2feet,spl_fsl_nc,fh_fsl_nc)
    
    print*,'après',fh_fsl(1,1),fh_fsl(17,1)
    
    ! some adding operations
    
    ! BC for BSL 
    if (bc1_type==PERIODIC_SPLINE) then
      fh_bsl(Neta1+1,:)    = fh_bsl(1,:)
      fh_bsl_nc(Neta1+1,:) = fh_bsl_nc(1,:)
    end if
    if (bc2_type==PERIODIC_SPLINE) then
      fh_bsl(:,Neta2+1)    = fh_bsl(:,1)
      fh_bsl_nc(:,Neta2+1) = fh_bsl_nc(:,1)
    end if
    
    do i=1,Neta1+1
      do j=1,Neta2+1
        fh_fsl(i,j) = fh_fsl(i,j)/jac_array(i,j)
      enddo
    enddo
    
    ! ------------ Diagnostics -----------------
    
    do i=1,Neta1+bc1_type
      do j=1,Neta2+bc2_type
        diag(1,step)  = diag(1,step)  + f(i,j)*delta_eta1*delta_eta2*jac_array(i,j)
        diag(2,step)  = diag(2,step)  + f(i,j)*delta_eta1*delta_eta2
        diag(3,step)  = diag(3,step)  + fh_bsl(i,j)*delta_eta1*delta_eta2*jac_array(i,j)
        diag(4,step)  = diag(4,step)  + fh_bsl(i,j)*delta_eta1*delta_eta2
        diag(5,step)  = diag(5,step)  + fh_bsl_nc(i,j)*delta_eta1*delta_eta2*jac_array(i,j)
        diag(6,step)  = diag(6,step)  + fh_bsl_nc(i,j)*delta_eta1*delta_eta2
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
      diag(3,step),diag(3,step)/diag(1,step)-1._f64,diag(4,step),diag(4,step)/diag(2,step)-1._f64, &
      diag(5,step),diag(5,step)/diag(1,step)-1._f64,diag(6,step),diag(6,step)/diag(2,step)-1._f64, &
      diag(7,step),diag(7,step)/diag(1,step)-1._f64,diag(8,step),diag(8,step)/diag(2,step)-1._f64, &
      diag(9,step),diag(9,step)/diag(1,step)-1._f64,diag(10,step),diag(10,step)/diag(2,step)-1._f64
  end do
  close(900)
  open(unit=950,file='carac.dat')
  do i=1,Neta1+1
    do j=1,Neta2+1
      write(950,*) eta1tot(i,j),eta2tot(i,j)
    end do
    end do
  close(950)

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

! conservation de la masse
!plot 'diag.dat' u 1:3 w lp title 'BSL', 'diag.dat' u 1:7 w lp title 'BSL NC', 'diag.dat' u 1:11 w lp title 'FSL', 'diag.dat' u 1:15 w lp title 'FSL NC'
! convergence en espace
!plot 'Conv_cart_nonsym_f3.dat' u 1:2 w lp title 'BSL', 'Conv_cart_nonsym_f3.dat' u 1:3 w lp title 'BSL NC', 'Conv_cart_nonsym_f3.dat' u 1:4 w lp title 'FSL', 'Conv_cart_nonsym_f3.dat' u 1:5 w lp title 'FSL NC'