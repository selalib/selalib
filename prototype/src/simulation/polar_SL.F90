program radial_1d_SL
#include "sll_working_precision.h"
#include "sll_assert.h"
#include "sll_memory.h"

  use sll_cubic_splines
  use cubic_non_uniform_splines
  use sll_constants
  implicit none
  
  type(sll_cubic_spline_2D), pointer :: spl_bsl,spl_bsl_nc,spl_fsl,spl_fsl_nc
  sll_int32  :: N,Nr,Ntheta,i,j,err,test_case,step,nb_step,visu_step,field_case
  sll_real64 :: rmin,rmax,r,dr,x,y,dt,theta,dtheta,val,val_bsl,val_bsl_nc,val_fsl,val_fsl_nc,tmp1,tmp2
  sll_real64 :: a1,a2,rr,thetath,k1r,k2r,k3r,k4r,k1theta,k2theta,k3theta,k4theta
  sll_real64,dimension(:,:), pointer :: f,fh_bsl,fh_bsl_nc,fh_fsl,fh_fsl_nc
  sll_real64,dimension(:,:), pointer :: rfeet,thetafeet
  sll_real64,dimension(:,:), pointer :: diag

  ! for the python script polar-exe.py
  !namelist /param/ N
  N = 256
  !read(*,NML=param);open(unit=900,file="gyrof.param");write(900,NML=param);close(900)
 
  ! physical parameters
  Nr = N
  Ntheta = N

  ! domain : disc of radius rmax with a hole of radius rmin
  rmin = 0.2_f64
  rmax = 0.8_f64
  
  ! time-scheme parameters
  dt = 1._f64/10._f64
  nb_step = 10
  
  ! visualization parameters
  visu_step = 1
  
  ! distribution function
  ! 1 : gaussian in r
  ! 2 : gaussian in r
  ! 3 : cos theta
  ! 4 : truncated gaussian in theta
  ! 5 : f=1
  ! 6 : gaussian in r and theta
  ! 7 : dirac 
  test_case = 6
  ! field
  ! 1 : divergence free complex field
  ! 2 : "rotation" (constant coef advection in theta)
  ! 3 : translation of vector (a1,a2)
  field_case = 2
  a1=0.001_f64
  a2=0.002_f64
	
  print *,'#N=',N
  print *,'#T=',real(nb_step,f64)*dt
  
  dr = (rmax-rmin)/real(Nr,f64)
  dtheta = 2._f64*sll_pi/real(Ntheta,f64)
	
  ! allocations
  SLL_ALLOCATE(f(Nr+1,Ntheta+1), err)
  SLL_ALLOCATE(fh_bsl(Nr+1,Ntheta+1), err)
  SLL_ALLOCATE(fh_bsl_nc(Nr+1,Ntheta+1), err)
  SLL_ALLOCATE(fh_fsl(Nr+1,Ntheta+1), err)
  SLL_ALLOCATE(fh_fsl_nc(Nr+1,Ntheta+1), err)
  SLL_ALLOCATE(rfeet(Nr+1,Ntheta+1), err)
  SLL_ALLOCATE(thetafeet(Nr+1,Ntheta+1), err)
  SLL_ALLOCATE(diag(10,0:nb_step), err)
	
  ! creation spline 
  spl_bsl => new_cubic_spline_2D(Nr+1, Ntheta+1, &
    rmin, rmax, &
    0._f64, 2._f64*sll_pi, &
    SLL_HERMITE, SLL_PERIODIC,&
    const_slope_x1_min = 0._f64,const_slope_x1_max = 0._f64)
  spl_bsl_nc => new_cubic_spline_2D(Nr+1, Ntheta+1, &
    rmin, rmax, &
    0._f64, 2._f64*sll_pi, &
    SLL_HERMITE, SLL_PERIODIC,&
    const_slope_x1_min = 0._f64,const_slope_x1_max = 0._f64)
  spl_fsl => new_cubic_spline_2D(Nr+1, Ntheta+1, &
    rmin, rmax, &
    0._f64, 2._f64*sll_pi, &
    SLL_HERMITE, SLL_PERIODIC,&
    const_slope_x1_min = 0._f64,const_slope_x1_max = 0._f64)
  spl_fsl_nc => new_cubic_spline_2D(Nr+1, Ntheta+1, &
    rmin, rmax, &
    0._f64, 2._f64*sll_pi, &
    SLL_HERMITE, SLL_PERIODIC,&
    const_slope_x1_min = 0._f64,const_slope_x1_max = 0._f64)

 
  
  ! Analytic distribution function
  do i=1,Nr+1
    r = rmin+real(i-1,f64)*dr
    do j=1,Ntheta+1
      theta = real(j-1,f64)*dtheta
      x = r*cos(theta)
      y = r*sin(theta)

      ! test-function
      if (test_case==1) then
        f(i,j) = exp(-3000._f64*(r-1.5_f64*rmin)**2)
      endif
      if (test_case==2) then
        f(i,j) = exp(-100._f64*(r-0.5_f64*(rmax+rmin))**2)
      endif
      if (test_case==3) then
        f(i,j) = cos(theta)
      endif
      if (test_case==4) then
        if(r>=0.4.and.r<=0.6)then
          f(i,j) = exp(-30._f64*(theta-0.5*sll_pi)**2)
        else
          f(i,j) = 0._f64
        endif  
      endif
      if (test_case==5) then
        f(i,j) = 1._f64
      endif
      if (test_case==6) then
        f(i,j) = exp(-100._f64*(r-0.5_f64*(rmax+rmin))**2)*exp(-30._f64*(theta-0.5*sll_pi)**2)
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
  do i=1,Nr
    r = rmin+real(i-1,f64)*dr
    do j=1,Ntheta
      tmp1 = tmp1 + f(i,j)*r
      tmp2 = tmp2 + f(i,j)
    enddo
  enddo
  tmp1 = tmp1*dr*dtheta
  tmp2 = tmp2*dr*dtheta
  
  ! initialization of the diagnostics
  ! with the weighted mass 
  diag(1,0) = tmp1  ! analytic solution
  diag(3,0) = tmp1  ! BSL
  diag(5,0) = tmp1  ! BSL NC
  diag(7,0) = tmp1  ! FSL
  diag(9,0) = tmp1  ! FSL NC
  ! with the classical mass
  diag(2,0) = tmp2  ! analytic solution
  diag(4,0) = tmp2  ! BSL NC
  diag(6,0) = tmp2  ! BSL
  diag(8,0) = tmp2  ! FSL
  diag(10,0) = tmp2 ! FSL NC
  
  ! checking that the function is close to 0 at the r boundary 
  print *,'#f0 boundary:',f(1,1),f(Nr+1,1)
  
  ! saving f
  open(unit=900,file='f0.dat')  
  do i=1,Nr+1
    r = rmin+real(i-1,f64)*dr
    do j=1,Ntheta+1
      theta = real(j-1,f64)*dtheta
      x = r*cos(theta)
      y = r*sin(theta)
      write(900,*) x,y,r,theta,f(i,j)
    enddo
    write(900,*) ' '
  end do
  close(900)

  ! Evolution in time - Analytic solution
  open(unit=900,file='carac.dat')  
  open(unit=800,file='fn.dat')
  do i=1,Nr+1
    do j=1,Ntheta+1
      r=rmin+real(i-1,f64)*dr
      theta=real(j-1,f64)*dtheta
      x = r*cos(theta)
      y = r*sin(theta)
      write(800,'(f0.3,a,f0.3)',advance='no') r,' ',theta
      write(800,'(a,f0.3,a,f0.3)',advance='no') ' ',x,' ',y
      do step=1,nb_step
          
        !RK4
        if(field_case==1)then
          rr = r
          thetath = theta
          k1r = rr*(rr-rmin)*(rmax-rr)*cos(thetath)
          k1theta = -(rr*(rmax+rmin-2._f64*rr)+2._f64*(rr-rmin)*(rmax-rr))*sin(thetath)
          rr = r + 0.5_f64*dt*k1r
          thetath = theta +  0.5_f64*dt*k1theta
          k2r = rr*(rr-rmin)*(rmax-rr)*cos(thetath)
          k2theta = -(rr*(rmax+rmin-2._f64*rr)+2._f64*(rr-rmin)*(rmax-rr))*sin(thetath)
          rr = r + 0.5_f64*dt*k2r
          thetath = theta +  0.5_f64*dt*k2theta
          k3r = rr*(rr-rmin)*(rmax-rr)*cos(thetath)
          k3theta = -(rr*(rmax+rmin-2._f64*rr)+2._f64*(rr-rmin)*(rmax-rr))*sin(thetath)
          rr = r + dt*k3r
          thetath = theta +  dt*k3theta
          k4r = rr*(rr-rmin)*(rmax-rr)*cos(thetath)
          k4theta = -(rr*(rmax+rmin-2._f64*rr)+2._f64*(rr-rmin)*(rmax-rr))*sin(thetath)
          r = r+dt/6._f64*(k1r+2._f64*k2r+2._f64*k3r+k4r)
          theta = theta+dt/6._f64*(k1theta+2._f64*k2theta+2._f64*k3theta+k4theta)
        endif
        ! rotation
        if(field_case==2)then
          theta = theta+dt
        endif
        ! translation
        if(field_case==3)then
          x = x + a1*dt
          y = y + a2*dt
          r = sqrt(x**2+y**2)
          if(y>=0)then
            theta = acos(x/r)
          else
            theta = 2._f64*sll_pi-acos(x/r)
          endif
        endif
          
        ! correction of theta
        if(theta>2._f64*sll_pi)then
          theta= theta-2._f64*sll_pi
        endif
        if(theta<0._f64)then
          theta= theta+2._f64*sll_pi
        endif
          
        x = r*cos(theta)
        y = r*sin(theta)
          					
        if (test_case==1) then
          val = exp(-3000._f64*(r-1.5_f64*rmin)**2)
        endif
        if (test_case==2) then
          val = exp(-100._f64*(r-0.5_f64*(rmax+rmin))**2)
        endif
        if (test_case==3) then
          val = cos(theta)
        endif
        if (test_case==4) then
          if(r>=0.4.and.r<=0.6)then
            val = exp(-30._f64*(theta-0.5*sll_pi)**2)
          else
            val = 0._f64
          endif
        endif
        if (test_case==5) then
          val = 1._f64
        endif
        if (test_case==6) then
          val = exp(-100._f64*(r-0.5_f64*(rmax+rmin))**2)*exp(-30._f64*(theta-0.5*sll_pi)**2)
        endif
        if (test_case==7) then
          val = 0._f64
          if ((i==17).and.(j==12)) then
            val = real(32,f64)
          end if
        endif
			
        diag(1,step) = diag(1,step) + val*(rmin+real(i-1,f64)*dr)*dr*dtheta
        diag(2,step) = diag(2,step) + val*dr*dtheta
					
        if(mod(step,visu_step)==0)then
          write(900,*) x,y
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
  
    do i=1,Nr+1
      r=rmin+real(i-1,f64)*dr
      do j=1,Ntheta+1
        fh_bsl_nc(i,j) = fh_bsl_nc(i,j)*r
        fh_fsl(i,j) = fh_fsl(i,j)*r
      enddo
      fh_bsl(i,Ntheta+1)=fh_bsl(i,1)
      fh_bsl_nc(i,Ntheta+1)=fh_bsl_nc(i,1)
      fh_fsl(i,Ntheta+1)=fh_fsl(i,1)
      fh_fsl_nc(i,Ntheta+1)=fh_fsl_nc(i,1)
    enddo
      
    call compute_cubic_spline_2D(fh_bsl,spl_bsl)
    call compute_cubic_spline_2D(fh_bsl_nc,spl_bsl_nc)
    call compute_cubic_spline_2D(fh_fsl,spl_fsl)
    call compute_cubic_spline_2D(fh_fsl_nc,spl_fsl_nc)
    
    do i=1,Nr+1
      do j=1,Ntheta+1
      				    
        ! ------------- BSL part -----------------
        
        r=rmin+real(i-1,f64)*dr
        theta=real(j-1,f64)*dtheta
        x= r*cos(theta)
        y = r*sin(theta) 
        
        ! RK4
        if(field_case==1)then
          rr = r
          thetath = theta
          k1r = rr*(rr-rmin)*(rmax-rr)*cos(thetath)
          k1theta = -(rr*(rmax+rmin-2._f64*rr)+2._f64*(rr-rmin)*(rmax-rr))*sin(thetath)
          rr = r + 0.5_f64*dt*k1r
          thetath = theta +  0.5_f64*dt*k1theta
          k2r = rr*(rr-rmin)*(rmax-rr)*cos(thetath)
          k2theta = -(rr*(rmax+rmin-2._f64*rr)+2._f64*(rr-rmin)*(rmax-rr))*sin(thetath)
          rr = r + 0.5_f64*dt*k2r
          thetath = theta +  0.5_f64*dt*k2theta
          k3r = rr*(rr-rmin)*(rmax-rr)*cos(thetath)
          k3theta = -(rr*(rmax+rmin-2._f64*rr)+2._f64*(rr-rmin)*(rmax-rr))*sin(thetath)
          rr = r + dt*k3r
          thetath = theta +  dt*k3theta
          k4r = rr*(rr-rmin)*(rmax-rr)*cos(thetath)
          k4theta = -(rr*(rmax+rmin-2._f64*rr)+2._f64*(rr-rmin)*(rmax-rr))*sin(thetath)
          r = r+dt/6._f64*(k1r+2._f64*k2r+2._f64*k3r+k4r)
          theta = theta+dt/6._f64*(k1theta+2._f64*k2theta+2._f64*k3theta+k4theta)
        endif
        if(field_case==2)then
          theta = theta+dt
        endif
        if(field_case==3)then
          x = x + a1*dt
          y = y + a2*dt
          r = sqrt(x**2+y**2)
          if(y>=0)then
            theta = acos(x/r)
          else
            theta = 2._f64*sll_pi-acos(x/r)
          endif
        endif
        
        ! correction of theta
        if(theta>2._f64*sll_pi)then
          theta= theta-2._f64*sll_pi
        endif
        if(theta<0._f64)then
          theta= theta+2._f64*sll_pi
        endif
        
        fh_bsl(i,j)    = interpolate_value_2D(r,theta,spl_bsl)
        fh_bsl_nc(i,j) = interpolate_value_2D(r,theta,spl_bsl_nc)/(rmin+real(i-1,f64)*dr)
    
        diag(3,step) = diag(3,step) + fh_bsl(i,j) * (rmin+real(i-1,f64)*dr)*dr*dtheta
        diag(4,step) = diag(4,step) + fh_bsl(i,j) * dr*dtheta
        
        diag(5,step) = diag(5,step) + fh_bsl_nc(i,j) * (rmin+real(i-1,f64)*dr)*dr*dtheta
        diag(6,step) = diag(6,step) + fh_bsl_nc(i,j) * dr*dtheta
        
        ! ------------- FSL part -----------------
    
        r=rmin+real(i-1,f64)*dr
        theta=real(j-1,f64)*dtheta
        x= r*cos(theta)
        y = r*sin(theta) 
        
        dt = -dt
        
        ! RK4
        if(field_case==1)then
          rr = r
          thetath = theta
          k1r = rr*(rr-rmin)*(rmax-rr)*cos(thetath)
          k1theta = -(rr*(rmax+rmin-2._f64*rr)+2._f64*(rr-rmin)*(rmax-rr))*sin(thetath)
          rr = r + 0.5_f64*dt*k1r
          thetath = theta +  0.5_f64*dt*k1theta
          k2r = rr*(rr-rmin)*(rmax-rr)*cos(thetath)
          k2theta = -(rr*(rmax+rmin-2._f64*rr)+2._f64*(rr-rmin)*(rmax-rr))*sin(thetath)
          rr = r + 0.5_f64*dt*k2r
          thetath = theta +  0.5_f64*dt*k2theta
          k3r = rr*(rr-rmin)*(rmax-rr)*cos(thetath)
          k3theta = -(rr*(rmax+rmin-2._f64*rr)+2._f64*(rr-rmin)*(rmax-rr))*sin(thetath)
          rr = r + dt*k3r
          thetath = theta +  dt*k3theta
          k4r = rr*(rr-rmin)*(rmax-rr)*cos(thetath)
          k4theta = -(rr*(rmax+rmin-2._f64*rr)+2._f64*(rr-rmin)*(rmax-rr))*sin(thetath)
          r = r+dt/6._f64*(k1r+2._f64*k2r+2._f64*k3r+k4r)
          theta = theta+dt/6._f64*(k1theta+2._f64*k2theta+2._f64*k3theta+k4theta)
        endif
        if(field_case==2)then
          theta = theta+dt
        endif
        if(field_case==3)then
          x = x + a1*dt
          y = y + a2*dt
          r = sqrt(x**2+y**2)
          if(y>=0)then
            theta = acos(x/r)
          else
            theta = 2._f64*sll_pi-acos(x/r)
          endif
        endif
        
        ! correction of theta
        if(theta>2._f64*sll_pi)then
          theta= theta-2._f64*sll_pi
        endif
        if(theta<0._f64)then
          theta= theta+2._f64*sll_pi
        endif
        
        rfeet(i,j) = r
        thetafeet(i,j) = theta
				
        dt = -dt
        
      enddo
    enddo
    
    call deposit_value_2D(rfeet,thetafeet,spl_fsl,fh_fsl)
    call deposit_value_2D(rfeet,thetafeet,spl_fsl_nc,fh_fsl_nc)
    
    do i=1,Nr+1
      r=rmin+real(i-1,f64)*dr
      do j=1,Ntheta+1
        fh_fsl(i,j) = fh_fsl(i,j)/r
      enddo
    enddo
    
    do i=1,Nr+1
      r = rmin+real(i-1,f64)*dr
      do j=1,Ntheta+1
        theta = real(j-1,f64)*dtheta
        diag(7,step) = diag(7,step)+ fh_fsl(i,j)*(rmin+real(i-1,f64)*dr)*dr*dtheta
        diag(8,step) = diag(8,step)+ fh_fsl(i,j)*dr*dtheta
        diag(9,step) = diag(9,step)+ fh_fsl_nc(i,j)*(rmin+real(i-1,f64)*dr)*dr*dtheta
        diag(10,step) = diag(10,step)+ fh_fsl_nc(i,j)*dr*dtheta
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
  val_bsl = 0._f64
  val_bsl_nc = 0._f64
  val_fsl = 0._f64
  val_fsl_nc = 0._f64
  do i=1,Nr+1
    do j=1,Ntheta+1
      val_bsl = max(val_bsl,abs(f(i,j)-fh_bsl(i,j)))
      val_bsl_nc = max(val_bsl_nc,abs(f(i,j)-fh_bsl_nc(i,j)))
      val_fsl = max(val_fsl,abs(f(i,j)-fh_fsl(i,j)))
      val_fsl_nc = max(val_fsl_nc,abs(f(i,j)-fh_fsl_nc(i,j)))
      val = max(val,abs(fh_fsl(i,j)-fh_bsl(i,j)))
    enddo
  enddo  
  print*,N,val_bsl,val_bsl_nc,val_fsl,val_fsl_nc,val

  open(unit=900,file='fh.dat')  
  do i=1,Nr!+1
    r = rmin+real(i-1,f64)*dr
    do j=1,Ntheta!+1
      theta = real(j-1,f64)*dtheta
      x = r*cos(theta)
      y = r*sin(theta)
      write(900,*) r,theta,x,y,f(i,j),fh_bsl(i,j),fh_bsl_nc(i,j),fh_fsl(i,j),fh_fsl_nc(i,j)
    enddo
    write(900,*) ' '
  enddo
  close(900)
	
  

end program
